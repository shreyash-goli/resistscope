"""05: Correlation analysis of predicted vs. real fold-resistance.

Loads the docking benchmark, correlates predicted delta_delta_g (and baselines)
against measured Rhee fold-resistance per drug and pooled, saves the plot, and
prints the headline Spearman number.

Usage::

    python scripts/05_validate.py
    python scripts/05_validate.py --docking data/docking_results/benchmark_docking.parquet
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

import config  # noqa: E402
from services.validation import (  # noqa: E402
    build_merged,
    compute_enrichment,
    per_drug_correlation,
    plot_validation,
    run_full_validation,
)

pd.set_option("display.width", 120)
pd.set_option("display.max_rows", 200)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--docking", type=Path,
                        default=config.DOCKING_DIR / "benchmark_docking.parquet")
    parser.add_argument("--out-plot", type=Path,
                        default=config.VALIDATION_DIR / "validation_plot.png")
    args = parser.parse_args()

    if not args.docking.exists():
        print(f"ERROR: docking results not found at {args.docking}\n"
              f"Run scripts/04_gpu_batch.py (or 04_dock_benchmark.py) first.")
        return 1

    print(f"Validating {args.docking}\n")
    results = run_full_validation(docking_parquet=args.docking)

    # Per-drug table for the docking method.
    dd = results[results.scoring_method == "docking_ddg"].copy()
    print("=== Docking ddG vs measured fold-resistance (Spearman) ===")
    show = dd[["drug", "spearman_rho", "spearman_pvalue", "pearson_r", "n_mutations"]]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # Method comparison (pooled): does docking beat the data/knowledge baselines?
    print("\n=== Method comparison (OVERALL, pooled) — docking vs baselines ===")
    ov = results[results.drug == "OVERALL"][
        ["scoring_method", "spearman_rho", "spearman_pvalue", "pearson_r", "n_mutations"]]
    print(ov.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    merged = build_merged(docking_parquet=args.docking)

    # Per-drug signal, DRMs only (where the confounded target is least noisy).
    print("\n=== Per-drug Spearman on primary DRMs ===")
    pdc = per_drug_correlation(merged)
    prim = pdc[pdc.scope == "primary"][["drug", "spearman_rho", "spearman_pvalue", "n_mutations"]]
    print(prim.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # Enrichment: are the top-ddG predictions enriched for known major DRMs?
    print("\n=== Enrichment: top-N ddG recovers known major DRMs ===")
    enr = compute_enrichment(merged)
    print(enr.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    plot_path = plot_validation(results, output_path=args.out_plot)
    print(f"\nPlot saved: {plot_path}")

    headline = results[(results.drug == "OVERALL") &
                       (results.scoring_method == "docking_ddg")].iloc[0]
    drv = prim[prim.drug == "DRV"]
    best_enr = enr.loc[enr["enrichment"].idxmax()]
    print("\n" + "=" * 64)
    print("HEADLINE")
    print(f"  Per-mutation Spearman (pooled) : {headline['spearman_rho']:+.3f} "
          f"(p={headline['spearman_pvalue']:.2e}, n={int(headline['n_mutations'])})")
    if len(drv):
        print(f"  Darunavir primary-DRM Spearman : {drv['spearman_rho'].iloc[0]:+.3f} "
              f"(p={drv['spearman_pvalue'].iloc[0]:.2f}, n={int(drv['n_mutations'].iloc[0])})")
    print(f"  DRM enrichment @ top-{int(best_enr['top_n'])} ddG      : "
          f"{best_enr['enrichment']:.2f}x (vs {best_enr['base_rate']:.0%} base rate)")
    print("=" * 64)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
