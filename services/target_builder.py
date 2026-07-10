"""Agentically assemble a new resistance target from just a protein name/structure.

This is the "bring your own target" core: a user uploads a receptor (or names a
protein), and instead of hand-curating a :class:`~targets.Target`, a Claude
tool-use agent **researches** the system over PubMed and returns a structured
spec — the drug class + inhibitor panel, the known resistance/variant mutations
(with the wildtype→mutant residues), the binding site, a suggested reference PDB,
and whether a public genotype–phenotype dataset exists (which decides whether the
target can be *validated* or only *triaged*).

The structural half (chains, pocket center) comes from the uploaded PDB
(``services/pdb_intake.py``); this module supplies the biology. Together they
build a Target the rest of the pipeline already knows how to run.

Same model + pattern as the literature agent (Haiku 4.5, forced-tool structured
output, PubMed grounding), so it degrades gracefully and never invents a PMID it
did not retrieve.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import config
from services.literature_agent import pubmed_search

_SEARCH_TOOL = {
    "name": "pubmed_search",
    "description": ("Search PubMed for drug-resistance / inhibitor literature on "
                    "this protein and return abstracts. Use focused queries "
                    "(protein + 'resistance mutations', + 'inhibitor', + a drug "
                    "name). Call more than once to cover drugs and mutations."),
    "input_schema": {
        "type": "object",
        "properties": {
            "query": {"type": "string"},
            "max_results": {"type": "integer", "minimum": 1, "maximum": 6},
        },
        "required": ["query"],
    },
}

_SUBMIT_TOOL = {
    "name": "submit_target",
    "description": ("Submit the assembled resistance-target spec. Only include "
                    "mutations and drugs you are confident are real and "
                    "documented; cite PMIDs you actually retrieved."),
    "input_schema": {
        "type": "object",
        "properties": {
            "target_id": {"type": "string", "description": "Short UPPER_SNAKE id, e.g. FLU_NA."},
            "label": {"type": "string", "description": "Human label, e.g. 'Influenza A neuraminidase'."},
            "organism": {"type": "string"},
            "drug_class": {"type": "string", "description": "e.g. 'neuraminidase inhibitor'."},
            "binding_site": {"type": "string", "description": "e.g. 'catalytic sialic-acid pocket'."},
            "drugs": {
                "type": "array",
                "items": {"type": "object", "properties": {
                    "abbrev": {"type": "string"}, "name": {"type": "string"}},
                    "required": ["abbrev", "name"]},
                "description": "Approved/known inhibitors of this target.",
            },
            "resistance_mutations": {
                "type": "array",
                "items": {"type": "object", "properties": {
                    "mutation": {"type": "string", "description": "e.g. H275Y (wt+pos+mut)."},
                    "note": {"type": "string"}},
                    "required": ["mutation"]},
                "description": "Documented resistance/variant mutations.",
            },
            "reference_pdb": {"type": "string", "description": "Suggested PDB id with an inhibitor bound."},
            "dataset": {
                "type": "object",
                "properties": {
                    "exists": {"type": "boolean"},
                    "name": {"type": "string"},
                    "notes": {"type": "string"},
                },
                "required": ["exists"],
                "description": "Public genotype-phenotype resistance dataset (decides validate vs triage-only).",
            },
            "confidence": {"type": "string", "enum": ["high", "medium", "low"]},
            "citation_pmids": {"type": "array", "items": {"type": "string"}},
        },
        "required": ["target_id", "label", "drug_class", "binding_site", "drugs",
                     "resistance_mutations", "dataset", "confidence"],
    },
}


@dataclass
class TargetSpec:
    protein: str
    spec: dict = field(default_factory=dict)   # the submit_target payload
    citations: list = field(default_factory=list)
    n_searches: int = 0
    error: str = ""

    @property
    def validatable(self) -> bool:
        return bool(self.spec.get("dataset", {}).get("exists"))


def _system_prompt() -> str:
    return (
        "You are a structural-pharmacology research assistant assembling a "
        "drug-resistance target definition for a docking pipeline. Given a "
        "protein, research its approved inhibitors and its documented resistance "
        "mutations using the pubmed_search tool, then call submit_target. Mutations "
        "must be written as wildtype+position+mutant (e.g. H275Y). Ground specific "
        "claims in retrieved literature; if you are unsure a mutation or dataset is "
        "real, omit it or lower the confidence rather than inventing it."
    )


def assemble_target(protein: str, hint: str = "", model: str | None = None,
                    max_iters: int = 7, max_searches: int = 5) -> TargetSpec:
    """Run the agent to assemble a resistance-target spec for ``protein``.

    ``hint`` is optional free text (e.g. a drug class or the PDB header title from
    an uploaded structure). Returns a :class:`TargetSpec`; ``spec`` is the
    submit_target payload, empty with ``error`` set on failure.
    """
    import anthropic

    model = model or config.CLAUDE_MODEL
    res = TargetSpec(protein=protein)
    client = anthropic.Anthropic()

    user = (f"Protein / target: {protein}\n"
            + (f"Context hint: {hint}\n" if hint else "")
            + "Assemble its resistance-target spec (inhibitors + resistance "
              "mutations + reference PDB + whether a public resistance dataset "
              "exists). Research first, then submit_target.")
    messages = [{"role": "user", "content": user}]
    seen: dict[str, dict] = {}
    submission = None

    for i in range(max_iters):
        force = res.n_searches >= max_searches or i >= max_iters - 1
        extra = {"tool_choice": {"type": "tool", "name": "submit_target"}} if force else {}
        resp = client.messages.create(
            model=model, max_tokens=1600, system=_system_prompt(),
            tools=[_SEARCH_TOOL, _SUBMIT_TOOL], messages=messages, **extra,
        )
        messages.append({"role": "assistant", "content": resp.content})
        tool_results, stop = [], False
        for b in resp.content:
            if b.type != "tool_use":
                continue
            if b.name == "submit_target":
                submission = b.input
                tool_results.append({"type": "tool_result", "tool_use_id": b.id, "content": "recorded"})
                stop = True
            elif b.name == "pubmed_search":
                res.n_searches += 1
                hits = pubmed_search(b.input["query"], int(b.input.get("max_results", 4)))
                for h in hits:
                    if h["pmid"]:
                        seen[h["pmid"]] = h
                tool_results.append({"type": "tool_result", "tool_use_id": b.id,
                                     "content": json.dumps([{k: h[k] for k in
                                     ("pmid", "title", "year", "journal", "abstract")} for h in hits])
                                     if hits else "No results; broaden the query."})
        if tool_results:
            messages.append({"role": "user", "content": tool_results})
        if stop or resp.stop_reason == "end_turn":
            break

    if submission is None:
        res.error = "agent did not submit a target spec"
        return res
    res.spec = submission
    cited = [p for p in (submission.get("citation_pmids") or []) if p in seen]
    res.citations = [{k: seen[p][k] for k in ("pmid", "title", "year", "journal")} for p in cited]
    return res
