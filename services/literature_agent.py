"""An agentic, literature-grounded builder for resistance-mechanism ground truth.

The faithfulness eval (``scripts/07``) grades Claude's generated explanations
against a curated reference of *known* structural mechanisms
(``mechanism_ground_truth.json``). That reference was hand-written and covered
only 15 protease mutations — so RT had nothing to score against, and PI coverage
was partial.

This module builds that reference **agentically**. For each (drug, mutation)
pair, Claude is given a ``pubmed_search`` tool (real NCBI E-utilities, abstracts
included) and the structural context, and runs a multi-step research loop:
search → read abstracts → refine → submit a concise mechanism cited to the
papers it actually read. A second Claude pass then **verifies** that the
submitted mechanism is supported by those abstracts and that the cited PMIDs are
real and relevant — unsupported entries are flagged ``verified=false`` and kept
out of the scored reference. This is deliberately a *different* task and prompt
from explanation generation (extract established mechanism from literature vs.
hypothesise from structure), and every entry carries its provenance + citations,
so it augments rather than launders the hand-curated seed.

Uses the Anthropic SDK tool-use loop on the same model as the rest of the
pipeline (Haiku 4.5 — the larger models refuse HIV-resistance content).
"""

from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime, timezone

import requests

import config
from services.explanation import (
    _DRUG_CLASS,
    _enzyme_noun,
    build_structural_context,
)
from targets import Target

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


# --- PubMed tool backend (search + abstracts) --------------------------------

def pubmed_search(query: str, max_results: int = 4, timeout: int = 20) -> list[dict]:
    """Search PubMed and return records **with abstracts** (best-effort).

    Backs the agent's ``pubmed_search`` tool. Returns dicts with keys
    ``pmid, title, year, journal, abstract``. Empty list on any failure so the
    agent degrades gracefully rather than crashing mid-loop.
    """
    try:
        r = requests.get(f"{_EUTILS}/esearch.fcgi", timeout=timeout, params={
            "db": "pubmed", "term": query, "retmax": max_results,
            "retmode": "json", "sort": "relevance"})
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        fetch = requests.get(f"{_EUTILS}/efetch.fcgi", timeout=timeout, params={
            "db": "pubmed", "id": ",".join(ids), "retmode": "xml"})
        root = ET.fromstring(fetch.text)
        out = []
        for art in root.findall(".//PubmedArticle"):
            pmid = art.findtext(".//PMID", default="")
            title = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
            # AbstractText may be split into labelled sections; concatenate.
            abst = " ".join("".join(node.itertext()) for node in art.findall(".//AbstractText"))
            journal = art.findtext(".//Journal/ISOAbbreviation") or art.findtext(".//Journal/Title") or ""
            year = art.findtext(".//JournalIssue/PubDate/Year") or art.findtext(".//PubDate/MedlineDate") or ""
            out.append({
                "pmid": pmid,
                "title": title.rstrip("."),
                "year": (year or "")[:4],
                "journal": journal,
                "abstract": abst[:1800],  # cap tokens; abstracts are enough signal
            })
        return out
    except Exception:  # noqa: BLE001 - tool is best-effort by design
        return []


# --- Agent tool definitions --------------------------------------------------

_PUBMED_TOOL = {
    "name": "pubmed_search",
    "description": (
        "Search PubMed for peer-reviewed literature and return matching articles "
        "with their abstracts. Use focused queries combining the enzyme, the "
        "mutation, the drug, and terms like 'resistance mechanism' or 'crystal "
        "structure'. Call it more than once with refined queries if the first "
        "results are thin."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "query": {"type": "string", "description": "PubMed query string."},
            "max_results": {"type": "integer", "description": "1-6 articles.",
                            "minimum": 1, "maximum": 6},
        },
        "required": ["query"],
    },
}

_SUBMIT_TOOL = {
    "name": "submit_mechanism",
    "description": (
        "Submit the final, literature-grounded structural mechanism for this "
        "mutation. Only cite PMIDs you actually retrieved and read via "
        "pubmed_search. Keep the mechanism to 2-3 sentences, specific about the "
        "molecular interaction (van der Waals, H-bond, sterics, electrostatics)."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "mechanism": {"type": "string",
                          "description": "2-3 sentence structural mechanism of resistance."},
            "key_interaction": {"type": "string",
                                "description": "The single most important interaction changed."},
            "affects_drugs": {"type": "array", "items": {"type": "string"},
                              "description": "Drug abbreviations this mutation affects."},
            "confidence": {"type": "string", "enum": ["high", "medium", "low"]},
            "citation_pmids": {"type": "array", "items": {"type": "string"},
                               "description": "PMIDs supporting the mechanism."},
        },
        "required": ["mechanism", "key_interaction", "affects_drugs", "confidence", "citation_pmids"],
    },
}


@dataclass
class MechanismResult:
    mutation: str
    mechanism: str = ""
    key_interaction: str = ""
    affects_drugs: list = field(default_factory=list)
    confidence: str = "low"
    citations: list = field(default_factory=list)   # [{pmid,title,year,journal}]
    verified: bool = False
    verify_note: str = ""
    n_searches: int = 0
    error: str = ""


def _system_prompt(t: Target) -> str:
    return (
        f"You are a structural virologist compiling a citation-grounded reference "
        f"of {_enzyme_noun(t)} drug-resistance mechanisms for a benchmark. For the "
        f"given mutation you must ground every claim in primary literature you "
        f"retrieve with the pubmed_search tool — do NOT rely on memory for "
        f"specific findings, and never cite a PMID you did not retrieve. Search "
        f"(refining queries as needed), read the abstracts, then call "
        f"submit_mechanism with a concise, mechanistically specific answer. If the "
        f"literature is thin, say so via a 'low' confidence rather than inventing "
        f"detail."
    )


def research_mechanism(
    drug: str,
    mutation: str,
    target: Target | None = None,
    model: str | None = None,
    max_iters: int = 8,
    max_searches: int = 5,
    verify: bool = True,
) -> MechanismResult:
    """Run the agentic research loop for one mutation and return a MechanismResult.

    ``drug`` seeds the query framing but the agent decides the full drug list it
    affects. The loop runs Claude with the pubmed_search + submit_mechanism tools,
    collecting every abstract it reads. Once it has done ``max_searches`` searches
    (or is on the last iteration) submission is *forced* via ``tool_choice`` so a
    thorough searcher can't run past the budget without answering.
    """
    import anthropic

    t = target if target is not None else config.ACTIVE_TARGET
    model = model or config.CLAUDE_MODEL
    res = MechanismResult(mutation=mutation)

    try:
        ctx = build_structural_context(mutation, target=t)
    except Exception:  # noqa: BLE001
        # Structure not built yet (e.g. RT before docking): fall back to a
        # minimal context from the mutation string so the agent can still
        # research the mechanism from literature.
        from services.explanation import AA_PROPERTIES, _parse_mutation, _region_for
        wt, pos, mut = _parse_mutation(mutation)
        ctx = {
            "position": pos, "wt_aa": wt, "mut_aa": mut,
            "wt_residue": AA_PROPERTIES.get(wt, {}).get("name", wt),
            "mut_residue": AA_PROPERTIES.get(mut, {}).get("name", mut),
            "region": _region_for(pos, t),
            "note": "structure not built; distances unavailable",
        }

    client = anthropic.Anthropic()
    drug_full = t.drugs.get(drug, drug)
    drug_class = _DRUG_CLASS.get(t.name, "inhibitor")
    user = (
        f"Enzyme: {_enzyme_noun(t)}\n"
        f"Mutation: {mutation}\n"
        f"Example affected drug: {drug_full} ({drug}), a {drug_class}\n"
        f"Structural context (computed from the wildtype structure):\n"
        f"{json.dumps(ctx, indent=2)}\n\n"
        f"Research and submit the structural resistance mechanism for {mutation}."
    )
    messages = [{"role": "user", "content": user}]
    seen_abstracts: dict[str, dict] = {}   # pmid -> record (for verification)
    submission = None

    for i in range(max_iters):
        force_submit = res.n_searches >= max_searches or i >= max_iters - 1
        extra = ({"tool_choice": {"type": "tool", "name": "submit_mechanism"}}
                 if force_submit else {})
        resp = client.messages.create(
            model=model, max_tokens=1400, system=_system_prompt(t),
            tools=[_PUBMED_TOOL, _SUBMIT_TOOL], messages=messages, **extra,
        )
        messages.append({"role": "assistant", "content": resp.content})
        tool_results = []
        stop = False
        for block in resp.content:
            if block.type != "tool_use":
                continue
            if block.name == "submit_mechanism":
                submission = block.input
                tool_results.append({"type": "tool_result", "tool_use_id": block.id,
                                     "content": "recorded"})
                stop = True
            elif block.name == "pubmed_search":
                res.n_searches += 1
                hits = pubmed_search(block.input["query"],
                                     int(block.input.get("max_results", 4)))
                for h in hits:
                    if h["pmid"]:
                        seen_abstracts[h["pmid"]] = h
                payload = [{"pmid": h["pmid"], "title": h["title"], "year": h["year"],
                            "journal": h["journal"], "abstract": h["abstract"]} for h in hits]
                tool_results.append({"type": "tool_result", "tool_use_id": block.id,
                                     "content": json.dumps(payload) if payload
                                     else "No results; try a broader query."})
        if tool_results:
            messages.append({"role": "user", "content": tool_results})
        if stop or resp.stop_reason == "end_turn":
            break

    if submission is None:
        res.error = "agent did not submit a mechanism"
        return res

    res.mechanism = (submission.get("mechanism") or "").strip()
    res.key_interaction = (submission.get("key_interaction") or "").strip()
    res.affects_drugs = list(submission.get("affects_drugs") or [])
    res.confidence = submission.get("confidence", "low")
    # Keep only citations the agent actually retrieved (drops hallucinated PMIDs).
    cited = [p for p in (submission.get("citation_pmids") or []) if p in seen_abstracts]
    res.citations = [{k: seen_abstracts[p][k] for k in ("pmid", "title", "year", "journal")}
                     for p in cited]

    if verify and res.mechanism:
        res.verified, res.verify_note = _verify_grounding(
            res, seen_abstracts, t, model, client)
    return res


def _verify_grounding(res: MechanismResult, seen: dict, t: Target,
                      model: str, client) -> tuple[bool, str]:
    """Second-pass check: is the mechanism supported by the retrieved abstracts?

    Guards against a confidently-worded but unsupported mechanism or an
    irrelevant citation slipping into the scored reference.
    """
    # Structured output via a forced tool call (this SDK has no output_config).
    verify_tool = {
        "name": "report_verification",
        "description": "Report whether the mechanism is grounded in the abstracts.",
        "input_schema": {
            "type": "object",
            "properties": {
                "supported": {"type": "boolean"},
                "citations_relevant": {"type": "boolean"},
                "note": {"type": "string"},
            },
            "required": ["supported", "citations_relevant", "note"],
        },
    }
    abstracts = "\n\n".join(
        f"PMID {p} ({r['journal']} {r['year']}): {r['title']}. {r['abstract']}"
        for p, r in seen.items()
    ) or "(no abstracts were retrieved)"
    prompt = (
        f"Mutation: {res.mutation} ({_enzyme_noun(t)}).\n\n"
        f"Proposed mechanism:\n{res.mechanism}\n\n"
        f"Abstracts the author read:\n{abstracts}\n\n"
        f"Is the proposed structural mechanism supported by these abstracts (not "
        f"contradicted, not invented beyond what a domain expert would accept)? "
        f"Are the retrieved citations relevant to this mutation's resistance?"
    )
    try:
        resp = client.messages.create(
            model=model, max_tokens=400,
            system=("You verify that a proposed drug-resistance mechanism is "
                    "grounded in the provided abstracts and general structural "
                    "biology. Be strict about invented specifics."),
            messages=[{"role": "user", "content": prompt}],
            tools=[verify_tool],
            tool_choice={"type": "tool", "name": "report_verification"},
        )
        v = next((b.input for b in resp.content
                  if b.type == "tool_use" and b.name == "report_verification"), None)
        if not v:
            return False, "no verification returned"
        return bool(v["supported"] and v["citations_relevant"]), v.get("note", "")
    except Exception as exc:  # noqa: BLE001
        return False, f"verification failed: {exc}"


def to_ground_truth_entry(res: MechanismResult) -> dict:
    """Convert a MechanismResult to a mechanism_ground_truth.json value.

    Keeps the fields the faithfulness eval reads (``mechanism``,
    ``affects_drugs``, ``source``) and adds provenance so agent-built entries are
    distinguishable from hand-curated ones.
    """
    src = "literature agent (Claude + PubMed)"
    if res.citations:
        src += "; " + ", ".join(f"PMID {c['pmid']}" for c in res.citations[:4])
    return {
        "mechanism": res.mechanism,
        "affects_drugs": res.affects_drugs,
        "source": src,
        "key_interaction": res.key_interaction,
        "confidence": res.confidence,
        "verified": res.verified,
        "provenance": "agent",
        "citations": res.citations,
        "built_at": datetime.now(timezone.utc).isoformat(),
    }
