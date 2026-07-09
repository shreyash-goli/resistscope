"""Scripted end-to-end CLI demo for ResistScope.

Docks a compound (default: darunavir) against wildtype HIV-1 protease and the top
clinically-relevant resistance mutations, prints a scored table, the robustness
score, a Claude mechanistic explanation (with PubMed citations) for the worst
mutation, and the validation summary.

Usage::

    python demo.py                       # darunavir, top 10 mutations (~5-10 min)
    python demo.py --smiles "CC(C)..."   # a custom compound
    python demo.py --mutations 5         # fewer mutations = faster
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd  # noqa: E402

import config  # noqa: E402
from services.docking import (  # noqa: E402
    dock_against_panel, dock_single, get_benchmark_smiles, smiles_to_pdbqt,
)
from services.explanation import build_structural_context, generate_explanation  # noqa: E402
from services.mutation_panel import load_panel  # noqa: E402
from services.scoring import compute_robustness_scores  # noqa: E402

BANNER = r"""
================================================================
  ResistScope - HIV-1 Protease Resistance Triage
================================================================
"""


def _severity(ddg):
    if ddg is None:
        return "?"
    if ddg >= config.DDG_DANGER_THRESHOLD:
        return "HIGH"
    if ddg >= config.DDG_WARNING_THRESHOLD:
        return "MED"
    return "low"


def clinically_relevant_panel(n: int) -> pd.DataFrame:
    """Top-N most prevalent major DRMs across all drugs (deduplicated)."""
    frames = [load_panel(p.stem, config.PANELS_DIR)
              for p in sorted(config.PANELS_DIR.glob("*.parquet"))]
    allp = pd.concat(frames, ignore_index=True)
    agg = allp[allp["is_primary"]].groupby("mutation").agg(
        position=("position", "first"),
        wildtype_aa=("wildtype_aa", "first"),
        mutant_aa=("mutant_aa", "first"),
        mean_log_fold_resistance=("mean_log_fold_resistance", "mean"),
        n_isolates=("n_isolates", "max"),
        is_primary=("is_primary", "any"),
    ).reset_index()
    # Only mutations we have a receptor structure for.
    agg = agg[agg["mutation"].apply(lambda m: (config.MUTANTS_DIR / f"{m}.pdbqt").exists())]
    return agg.sort_values("n_isolates", ascending=False).head(n).reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", default=None, help="Compound SMILES (default: darunavir).")
    parser.add_argument("--mutations", type=int, default=10, help="How many top DRMs to dock.")
    args = parser.parse_args()

    print(BANNER)
    smiles = args.smiles or get_benchmark_smiles()["DRV"]
    # Use the real drug id for the darunavir default (nicer explanation + citations).
    drug_id = "DRV" if args.smiles is None else "QUERY"
    label = "darunavir" if args.smiles is None else "custom compound"
    print(f"Compound : {label}")
    print(f"SMILES   : {smiles}\n")

    # (a) Wildtype dock
    print("[1/5] Docking against wildtype HIV-1 protease ...")
    ligand = smiles_to_pdbqt(smiles)
    wt = dock_single(ligand, config.STRUCTURES_DIR / "wildtype.pdbqt")
    wt_dg = wt["delta_g"]
    print(f"      Wildtype binding affinity: {wt_dg:.2f} kcal/mol\n")

    # (b) Dock against the top clinically-relevant mutations
    panel = clinically_relevant_panel(args.mutations)
    print(f"[2/5] Docking against {len(panel)} top clinical resistance mutations "
          f"(this is the slow part) ...")
    results = dock_against_panel(smiles, drug_id, panel, progress=False)

    # (c) Scored table
    muts = results[results["mutation"] != "WT"].dropna(subset=["delta_delta_g"]).copy()
    muts = muts.sort_values("delta_delta_g", ascending=False)
    print("\n[3/5] Predicted binding change per mutation "
          "(ddG > 0 = weaker binding = resistance):\n")
    print(f"      {'mutation':<10}{'dG':>8}{'ddG':>8}  severity")
    print(f"      {'-'*34}")
    for _, r in muts.iterrows():
        ddg = float(r["delta_delta_g"])
        print(f"      {r['mutation']:<10}{r['delta_g']:>8.2f}{ddg:>+8.2f}  {_severity(ddg)}")

    # (d) Robustness score
    scores = compute_robustness_scores(results, panel)
    print(f"\n[4/5] Robustness score: {scores['robustness_0_100']:.0f} / 100 "
          f"(100 = binding unaffected across the panel)")
    print(f"      worst-case ddG = {scores['worst_case_ddg']:+.2f} kcal/mol, "
          f"prevalence-weighted = {scores['prevalence_weighted_ddg']:+.2f}")

    # (e) Explain the worst mutation (with PubMed citations)
    worst = muts.iloc[0]
    print(f"\n[5/5] Mechanistic explanation for the worst mutation "
          f"({worst['mutation']}, ddG {float(worst['delta_delta_g']):+.2f}):\n")
    try:
        ctx = build_structural_context(
            worst["mutation"], {"delta_delta_g": float(worst["delta_delta_g"]),
                                "delta_g": float(worst["delta_g"])})
        expl = generate_explanation(drug_id, worst["mutation"],
                                    float(worst["delta_delta_g"]), ctx, cite=True)
        print("      " + expl.replace("\n", "\n      "))
        cache = config.EXPLANATIONS_DIR / f"{drug_id}_{worst['mutation']}.json"
        if cache.exists():
            import json
            cites = json.loads(cache.read_text()).get("citations", [])
            if cites:
                print("\n      Literature (PubMed):")
                for c in cites:
                    print(f"        - {c['journal']} {c['year']} (PMID {c['pmid']}): {c['title']}")
    except Exception as exc:  # noqa: BLE001
        print(f"      (explanation unavailable: {exc})")
        print("      Set ANTHROPIC_API_KEY to enable Claude explanations.")

    # (f) Validation summary
    corr_p = config.VALIDATION_DIR / "validation_correlations.parquet"
    faith_p = config.VALIDATION_DIR / "faithfulness_scores.parquet"
    if corr_p.exists():
        corr = pd.read_parquet(corr_p)
        ov = corr[(corr.drug == "OVERALL") & (corr.scoring_method == "docking_ddg")]
        print("\n" + "=" * 64)
        print("Validation (precomputed benchmark):")
        if len(ov):
            print(f"  per-mutation Spearman rho (pooled) = {ov.spearman_rho.iloc[0]:+.3f}")
        if faith_p.exists():
            f = pd.read_parquet(faith_p)
            print(f"  explanation faithfulness = {100*(f.score==2).mean():.0f}% correct "
                  f"(mean {f.score.mean():.2f}/2, n={len(f)})")
        print("=" * 64)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
