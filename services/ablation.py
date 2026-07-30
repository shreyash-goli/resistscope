"""Faithfulness attribution ablation: does the *structural pipeline* make the
explanations faithful, or does Claude's parametric knowledge of (famous)
resistance mutations?

The rigorous benchmark (``services/benchmark.py``) shows the docking ΔΔG is a
chance-level ranker (ROC-AUC ≈ 0.51) that loses to a zero-cost prevalence
baseline. That raises a sharp question for the interpretability finding: when
Claude explains *why* V82A resists darunavir at 72% faithfulness, is that the
structural context we feed it, or biology it already memorized?

This module runs the *same* explanation model, system prompt, and (downstream)
Claude judge across three conditions that differ **only** in the context handed
to the explanation step:

- ``full``      : real structural context + real ΔΔG (the production pipeline).
- ``minimal``   : drug + mutation identity only — no structural facts, no ΔΔG.
- ``corrupted`` : production-shaped context, but the *pipeline-derived* fields
                  (ligand distance, direct-contact flag, subpocket, ΔΔG/ΔG) are
                  deranged across pairs, while the mutation's intrinsic chemistry
                  (which the mutation string itself fixes — volume/charge/
                  hydrophobicity change) is left correct so the context never
                  self-contradicts.

Contrasts: ``minimal`` vs ``full`` = the total value of stating structural
context; ``corrupted`` vs ``full`` = the value of *correct* geometry/energetics.
Flat faithfulness across conditions ⇒ the pipeline is not what makes the
explanations faithful (Claude's prior is), which predicts the pipeline's value
lies on obscure / novel mutations (RT, BYO targets), not textbook DRMs.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import config
from services.explanation import _DRUG_CLASS, _resolve, _system_prompt
from targets import Target

CONDITIONS = ("full", "minimal", "corrupted")

# Pipeline-derived (structure/docking) context fields the corruption deranges.
# Intrinsic-chemistry fields (position, wt/mut aa & residue, volume/charge/
# hydrophobicity change) are NOT here: the mutation string fixes them, so
# scrambling them would make the context contradict the mutation itself.
_CORRUPT_FIELDS = (
    "distance_from_ligand_centroid_angstrom",
    "min_distance_to_ligand_angstrom",
    "contacts_ligand_directly",
    "region",
    "delta_delta_g",
    "delta_g",
)


def _ddg_str(ddg: float | None) -> str:
    return (
        f"{ddg:+.2f} kcal/mol (positive = weaker binding = resistance)"
        if ddg is not None
        else "not available"
    )


def build_user_prompt(condition: str, drug: str, drug_full: str, drug_class: str,
                      mutation: str, ddg: float | None, context: dict) -> str:
    """User prompt for one condition.

    ``full`` / ``corrupted`` reproduce the production explanation prompt
    (``services.explanation.generate_explanation``, cite=False) verbatim so the
    control is the real pipeline; ``minimal`` drops the ΔΔG line and the
    structural-context JSON entirely.
    """
    head = f"Drug: {drug_full} ({drug}), a {drug_class}\nMutation: {mutation}\n"
    tail = (f"\nExplain the structural mechanism by which {mutation} likely "
            f"reduces {drug_full} binding.")
    if condition == "minimal":
        return head + tail
    return (
        head
        + f"Predicted delta-delta-G: {_ddg_str(ddg)}\n"
        + f"Structural context:\n{json.dumps(context, indent=2)}\n"
        + tail
    )


def derange_contexts(reals: list[tuple[dict, float | None]], seed: int = 0
                     ) -> list[tuple[dict, float | None]]:
    """Corrupt each (context, ΔΔG) by donating another pair's pipeline fields.

    Uses a derangement (no pair keeps its own donor), so every corrupted context
    has realistic-but-wrong geometry/energetics for its mutation. The ΔΔG shown
    in the prompt follows the donor too. Intrinsic-chemistry fields are untouched.
    """
    n = len(reals)
    if n < 2:
        return [(dict(c), d) for c, d in reals]
    rng = np.random.default_rng(seed)
    perm = np.arange(n)
    while True:
        rng.shuffle(perm)
        if not np.any(perm == np.arange(n)):
            break
    out = []
    for i, (ctx, _ddg) in enumerate(reals):
        donor_ctx, donor_ddg = reals[perm[i]]
        c = dict(ctx)
        for f in _CORRUPT_FIELDS:
            if f in donor_ctx:
                c[f] = donor_ctx[f]
            else:
                c.pop(f, None)  # donor lacked it → don't keep our own
        out.append((c, donor_ddg))
    return out


def generate_for_condition(drug: str, mutation: str, ddg: float | None,
                           context: dict, condition: str, cache_dir: Path,
                           model: str | None = None,
                           target: Target | None = None) -> str:
    """Generate (or return cached) an explanation for one pair under one condition.

    Cached per condition dir as ``{drug}_{mutation}.json`` so re-runs are free
    and each condition's dir is directly consumable by
    ``services.explanation.evaluate_faithfulness``.
    """
    t = _resolve(target)
    cache_dir = Path(cache_dir)
    cache_path = cache_dir / f"{drug}_{mutation}.json"
    if cache_path.exists():
        return json.loads(cache_path.read_text())["explanation"]

    import anthropic

    model = model or config.CLAUDE_MODEL
    drug_full = t.drugs.get(drug, drug)
    drug_class = _DRUG_CLASS.get(t.name, "inhibitor")
    prompt = build_user_prompt(condition, drug, drug_full, drug_class,
                               mutation, ddg, context)

    client = anthropic.Anthropic()
    resp = client.messages.create(
        model=model,
        max_tokens=config.CLAUDE_MAX_TOKENS,
        system=_system_prompt(t),
        messages=[{"role": "user", "content": prompt}],
    )
    expl = next((b.text for b in resp.content if b.type == "text"), "").strip()

    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps({
        "drug": drug,
        "mutation": mutation,
        "condition": condition,
        "delta_delta_g": (float(ddg) if ddg is not None else None),
        "structural_context": (context if condition != "minimal" else None),
        "explanation": expl,
        "model": model,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }, indent=2))
    return expl
