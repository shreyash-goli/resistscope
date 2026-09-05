"""LLM-free structural-consistency check on generated explanations (paper Appendix).

Motivation
----------
The faithfulness evaluation is a closed loop over LLMs: Claude generates, Claude
judges, and RT ground truth is agent-built. Cross-model judge agreement shows
*consistency*, not *correctness*. This script adds an out-of-family check that
uses no language model at all.

Method
------
One claim in every explanation is mechanically verifiable against the receptor:
**does the mutated residue touch the ligand?** We
  1. read the ground-truth geometry from the pipeline's own structural context
     (``min_distance_to_ligand_angstrom``; contact iff < 4 Å), and
  2. classify what the *explanation text* asserts about contact, using a
     deterministic regex lexicon with negation resolved first (so "does not
     contact the inhibitor" is a NON-contact claim, not a contact one).
Explanations making no explicit claim are counted as ABSTAIN and excluded from
accuracy, not scored as wrong.

Three analyses fall out:
  A. Contact-claim accuracy per ablation condition. `minimal` has no distance
     field and `corrupted` has a flipped one, so if the pipeline's geometry is
     what drives the explanations, accuracy should track full > minimal >
     corrupted -- reproducing the ablation *without* an LLM judge.
  B. Judge score vs. structural correctness. If the LLM judge is tracking
     something real, structurally-wrong explanations should score lower.
  C. Overall accuracy and abstention rate.

Outputs
-------
data/validation/structural_consistency.csv    (per explanation)
data/validation/structural_consistency.json   (summary)

Usage: python scripts/16_structural_consistency.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
CONTACT_THRESHOLD_A = 4.0

DIRS = {
    "protease": ROOT / "data/ablation",
    "RT": ROOT / "data/ablation_HIV1_RT",
}
SCORES = {
    "protease": ROOT / "data/validation/faithfulness_ablation.parquet",
    "RT": ROOT / "data/rt/validation/faithfulness_ablation.parquet",
}

# --- deterministic claim lexicon -------------------------------------------
# Negation patterns are tested FIRST; each asserts the residue does NOT touch
# the ligand. Ordering matters because "does not contact" contains "contact".
NON_CONTACT = [
    r"(?:does|do|is|are|was|were)\s+not\s+(?:directly\s+)?(?:contact|touch|interact)",
    r"without\s+(?:directly\s+)?(?:contact|touching|interacting)",
    r"no\s+direct\s+(?:contact|interaction|van\s+der\s+waals)",
    r"not\s+in\s+direct\s+contact",
    r"rather\s+than\s+direct\s+(?:contact|steric)",
    r"not\s+directly\s+contact",
    r"(?:distal|remote)\s+(?:to|from)\s+the\s+(?:ligand|inhibitor|pocket|binding)",
    r"indirect(?:ly)?\s+(?:allosteric|mechanism|effect|modulat)",
    r"allosteric(?:ally)?\s+(?:mechanism|modulat|effect|rather)",
]
CONTACT = [
    r"direct\s+(?:van\s+der\s+waals\s+)?contact",
    r"directly\s+contact(?:s|ing)?\b",
    r"makes?\s+(?:critical\s+|extensive\s+|favorable\s+)?van\s+der\s+waals\s+contacts?\s+with",
    r"in\s+direct\s+contact\s+with",
    r"direct\s+steric\s+(?:clash|contact)",
    r"packs?\s+(?:directly\s+)?against\s+the\s+(?:ligand|inhibitor)",
]

_NON = [re.compile(p) for p in NON_CONTACT]
_CON = [re.compile(p) for p in CONTACT]


def classify_claim(text: str) -> str | None:
    """Return 'contact', 'no_contact', or None (abstain). Negation wins."""
    t = " ".join(text.lower().split())
    if any(p.search(t) for p in _NON):
        return "no_contact"
    if any(p.search(t) for p in _CON):
        return "contact"
    return None


def load_truth(target: str) -> dict[tuple[str, str], float]:
    """True min-distance per (drug, mutation), read from the `full` condition."""
    truth = {}
    for f in sorted((DIRS[target] / "full").glob("*.json")):
        d = json.loads(f.read_text())
        ctx = d.get("structural_context") or {}
        dist = ctx.get("min_distance_to_ligand_angstrom")
        if dist is not None:
            truth[(d["drug"], d["mutation"])] = float(dist)
    return truth


def build_rows(target: str) -> pd.DataFrame:
    truth = load_truth(target)
    rows = []
    for cond in ("full", "minimal", "corrupted"):
        for f in sorted((DIRS[target] / cond).glob("*.json")):
            d = json.loads(f.read_text())
            key = (d["drug"], d["mutation"])
            if key not in truth:
                continue
            dist = truth[key]
            actual = "contact" if dist < CONTACT_THRESHOLD_A else "no_contact"
            claim = classify_claim(d["explanation"])
            rows.append({
                "target": target, "condition": cond,
                "drug": d["drug"], "mutation": d["mutation"],
                "true_min_dist": dist, "actual": actual, "claim": claim,
                "correct": (None if claim is None else claim == actual),
            })
    return pd.DataFrame(rows)


def attach_judge_scores(df: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for target, path in SCORES.items():
        if not path.exists():
            continue
        s = pd.read_parquet(path)[["drug", "mutation", "condition", "score"]]
        s["target"] = target
        frames.append(s)
    if not frames:
        return df
    sc = pd.concat(frames, ignore_index=True)
    return df.merge(sc, on=["target", "drug", "mutation", "condition"], how="left")


def main() -> None:
    df = pd.concat([build_rows(t) for t in DIRS], ignore_index=True)
    df = attach_judge_scores(df)
    df.to_csv(ROOT / "data/validation/structural_consistency.csv", index=False)

    scored = df[df["claim"].notna()]
    summary = {
        "n_explanations": int(len(df)),
        "n_explicit_claims": int(len(scored)),
        "abstention_rate": float(1 - len(scored) / len(df)),
        "overall_accuracy": float(scored["correct"].mean()),
        "contact_threshold_A": CONTACT_THRESHOLD_A,
        "by_condition": {}, "by_target_condition": {},
    }

    print(f"n={len(df)} explanations | explicit contact claim in {len(scored)} "
          f"({len(scored)/len(df):.0%}) | overall accuracy {scored['correct'].mean():.3f}\n")

    print("=== A. contact-claim accuracy by condition (no LLM involved) ===")
    for cond in ("full", "minimal", "corrupted"):
        g = scored[scored.condition == cond]
        acc = g["correct"].mean()
        summary["by_condition"][cond] = {
            "n": int(len(g)), "accuracy": float(acc),
            "n_correct": int(g["correct"].sum()),
        }
        print(f"  {cond:10s} n={len(g):3d}  accuracy={acc:.3f}  ({int(g['correct'].sum())}/{len(g)})")

    # full vs corrupted, and full vs minimal: 2x2 Fisher exact
    for a, b in (("full", "minimal"), ("full", "corrupted")):
        ga, gb = scored[scored.condition == a], scored[scored.condition == b]
        tab = [[int(ga.correct.sum()), int((~ga.correct.astype(bool)).sum())],
               [int(gb.correct.sum()), int((~gb.correct.astype(bool)).sum())]]
        odds, p = fisher_exact(tab)
        summary[f"fisher_{a}_vs_{b}"] = {"p": float(p), "odds_ratio": float(odds),
                                         "table": tab}
        print(f"  {a} vs {b}: Fisher p={p:.4g} (OR={odds:.2f})")

    print("\n=== per target x condition ===")
    for (t, c), g in scored.groupby(["target", "condition"]):
        summary["by_target_condition"][f"{t}/{c}"] = {
            "n": int(len(g)), "accuracy": float(g["correct"].mean())}
        print(f"  {t:9s} {c:10s} n={len(g):3d}  accuracy={g['correct'].mean():.3f}")

    if "score" in scored.columns and scored["score"].notna().any():
        print("\n=== B. does the LLM judge track structural correctness? ===")
        ok = scored[scored.correct == True]["score"].dropna()
        bad = scored[scored.correct == False]["score"].dropna()
        if len(ok) and len(bad):
            u, p = mannwhitneyu(ok, bad, alternative="greater")
            summary["judge_vs_structure"] = {
                "mean_score_structurally_correct": float(ok.mean()),
                "mean_score_structurally_wrong": float(bad.mean()),
                "n_correct": int(len(ok)), "n_wrong": int(len(bad)),
                "mannwhitney_p": float(p),
            }
            print(f"  judge score | structurally CORRECT: {ok.mean():.2f} (n={len(ok)})")
            print(f"  judge score | structurally WRONG:   {bad.mean():.2f} (n={len(bad)})")
            print(f"  Mann-Whitney one-sided p = {p:.4g}")

    (ROOT / "data/validation/structural_consistency.json").write_text(
        json.dumps(summary, indent=2))
    print("\nwrote structural_consistency.{csv,json}")


if __name__ == "__main__":
    main()
