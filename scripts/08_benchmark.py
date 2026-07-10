"""08: Rigorous benchmarking of the docking-ΔΔG predictor.

Beyond the plain correlation (script 05), this puts a significance test or
confidence interval on every claim: permutation p-values + bootstrap CIs on the
top-N DRM enrichment, threshold-free ROC/PR-AUC for DRM recovery, and an honest
de-confounding analysis (does a cleaner single-mutation target rescue the
magnitude correlation? — for this dataset, no).

Writes ``benchmark_metrics.json`` (+ parquets) to the target's validation dir,
consumed by the API ``/benchmark`` endpoint and the Validation tab.

Usage::

    python scripts/08_benchmark.py                 # HIV-1 protease (default)
    python scripts/08_benchmark.py --target rt     # reverse transcriptase
    python scripts/08_benchmark.py --fast          # fewer permutations (dev)
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

import config  # noqa: E402
from services.benchmark import run_benchmark  # noqa: E402

pd.set_option("display.width", 120)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", default="HIV1_PR",
                    help="Docking target: HIV1_PR / pr (default) or HIV1_RT / rt.")
    ap.add_argument("--docking", type=Path, default=None,
                    help="Docking parquet (default: target's benchmark_docking.parquet).")
    ap.add_argument("--fast", action="store_true",
                    help="Fewer permutations/bootstraps for a quick dev run.")
    args = ap.parse_args()

    t = config.set_active_target(args.target)
    print(f"Target: {t.name} ({t.label})")

    docking = args.docking or (config.DOCKING_DIR / "benchmark_docking.parquet")
    if not Path(docking).exists():
        print(f"ERROR: docking results not found at {docking}\n"
              f"Run scripts/04_gpu_batch.py (or 04_dock_benchmark.py) first.")
        return 1

    n_perm, n_boot = (1000, 1000) if args.fast else (5000, 5000)
    summary = run_benchmark(docking_parquet=Path(docking),
                            n_perm=n_perm, n_boot=n_boot)

    enr = pd.DataFrame(summary["enrichment"])
    rank = pd.DataFrame(summary["ranking"])
    decon = pd.DataFrame(summary["deconfounding"])

    print("\n=== Top-N DRM enrichment (permutation p, bootstrap 95% CI) ===")
    print(enr[["top_n", "precision", "enrichment", "ci_low", "ci_high", "perm_p"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    print("\n=== DRM-recovery ranking (ROC-AUC / PR-AUC; 0.5 = chance) ===")
    print(rank[["drug", "roc_auc", "pr_auc", "base_rate", "n"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    print("\n=== De-confounding: does a cleaner target rescue the correlation? ===")
    print(decon[["subset", "n_isolates_used", "n_pairs", "spearman_rho", "spearman_pvalue"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    h = summary["headline"]
    print("\n" + "=" * 66)
    print("HEADLINE (honest)")
    print(f"  Best DRM enrichment : {h['best_enrichment']:.2f}× at top-"
          f"{h['best_enrichment_top_n']} "
          f"(95% CI {h['best_enrichment_ci'][0]:.2f}–{h['best_enrichment_ci'][1]:.2f}, "
          f"perm p={h['best_enrichment_perm_p']:.1e})")
    print(f"  Pooled DRM ROC-AUC  : {h['pooled_roc_auc']:.3f} "
          f"(95% CI {h['pooled_roc_auc_ci'][0]:.3f}–{h['pooled_roc_auc_ci'][1]:.3f}) "
          f"→ ~chance as a global ranker")
    print(f"  De-confounding rescues magnitude correlation? "
          f"{'YES' if h['deconfounding_rescues_correlation'] else 'NO'}")
    print("  Read: docking ΔΔG is a coarse DRM-triage flag, not a quantitative")
    print("  resistance predictor.")
    print("=" * 66)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
