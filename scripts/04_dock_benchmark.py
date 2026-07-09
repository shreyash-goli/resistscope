"""04: Dock the benchmark PIs against their mutation panels.

For each protease inhibitor with a panel, docks the drug against wildtype + the
mutations in its panel and records delta_g / delta_delta_g. Results are written
to ``data/docking_results/benchmark_docking.parquet``.

Because AutoDock Vina is CPU-only and each dock takes ~15-25 s, the full panel
(~1,600 docks) is a multi-hour job. Use ``--subset`` to scope it:

    python scripts/04_dock_benchmark.py --subset primary       # ~180 docks, ~1 h
    python scripts/04_dock_benchmark.py --subset prevalence --min-isolates 20
    python scripts/04_dock_benchmark.py --subset all           # full panel, ~9 h

Other options: ``--drugs DRV NFV`` (subset of drugs), ``--replicates 3``
(average seeded runs for noise-robust ddG), ``--exhaustiveness``, ``--workers``,
``--cpu-per-worker``, ``--out``.
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

import config  # noqa: E402
from services.docking import dock_against_panel, get_benchmark_smiles  # noqa: E402
from services.mutation_panel import load_panel  # noqa: E402


def _select_subset(panel: pd.DataFrame, subset: str, min_isolates: int) -> pd.DataFrame:
    """Filter a panel to the requested docking subset."""
    if subset == "all":
        return panel
    if subset == "primary":
        return panel[panel["is_primary"]].reset_index(drop=True)
    if subset == "prevalence":
        return panel[panel["n_isolates"] >= min_isolates].reset_index(drop=True)
    raise ValueError(f"unknown subset: {subset}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subset", choices=["primary", "prevalence", "all"],
                        default="primary", help="Which mutations to dock.")
    parser.add_argument("--min-isolates", type=int, default=20,
                        help="Prevalence threshold for --subset prevalence.")
    parser.add_argument("--drugs", nargs="+", default=None,
                        help="Restrict to these drug abbreviations (default: all with panels).")
    parser.add_argument("--replicates", type=int, default=1,
                        help="Seeded docks per pair, averaged (noise robustness).")
    parser.add_argument("--exhaustiveness", type=int, default=config.VINA_EXHAUSTIVENESS)
    parser.add_argument("--backend", choices=["vina", "unidock"], default="vina",
                        help="Docking engine: CPU Vina (default) or GPU Uni-Dock.")
    parser.add_argument("--unidock-bin", default="unidock",
                        help="Uni-Dock executable name/path (--backend unidock).")
    parser.add_argument("--search-mode", choices=["fast", "balance", "detail"],
                        default=None,
                        help="Uni-Dock effort preset (overrides --exhaustiveness).")
    parser.add_argument("--workers", type=int, default=1,
                        help="Concurrent docks (Vina only; default 1 = sequential, "
                             "each Vina using all cores; robust for long runs).")
    parser.add_argument("--cpu-per-worker", type=int, default=None,
                        help="Vina threads per worker when --workers > 1.")
    parser.add_argument("--out", type=Path,
                        default=config.DOCKING_DIR / "benchmark_docking.parquet")
    args = parser.parse_args()

    smiles = get_benchmark_smiles()

    # Drugs with both a SMILES and a panel on disk.
    drugs = args.drugs or [
        d for d in config.PI_DRUGS
        if (config.PANELS_DIR / f"{d}.parquet").exists() and d in smiles
    ]

    print(f"Subset: {args.subset}"
          + (f" (n_isolates>={args.min_isolates})" if args.subset == "prevalence" else "")
          + f" | replicates: {args.replicates} | exhaustiveness: {args.exhaustiveness}")
    print(f"Drugs: {drugs}\n")

    all_results = []
    t_start = time.time()
    for drug in drugs:
        panel = load_panel(drug, config.PANELS_DIR)
        sub = _select_subset(panel, args.subset, args.min_isolates)
        print(f"=== {drug} ({config.PI_DRUGS.get(drug, drug)}): "
              f"{len(sub)} mutations + WT ===")
        t0 = time.time()
        res = dock_against_panel(
            smiles[drug], drug, sub,
            exhaustiveness=args.exhaustiveness,
            replicates=args.replicates,
            n_workers=args.workers,
            cpu_per_worker=args.cpu_per_worker,
            backend=args.backend,
            unidock_bin=args.unidock_bin,
            search_mode=args.search_mode,
        )
        all_results.append(res)
        dt = time.time() - t0
        n_fail = int(res["delta_g"].isna().sum())
        print(f"    done in {dt/60:.1f} min | {n_fail} failed dock(s)\n")

    results = pd.concat(all_results, ignore_index=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    results.to_parquet(args.out, index=False)

    # Summary.
    print("=" * 64)
    print(f"{'drug':<6}{'docked':>10}{'failed':>10}{'mean_ddG':>12}{'max_ddG':>12}")
    print("-" * 64)
    for drug in drugs:
        r = results[results["drug"] == drug]
        muts = r[r["mutation"] != "WT"]
        docked = int(muts["delta_g"].notna().sum())
        failed = int(muts["delta_g"].isna().sum())
        mean_ddg = muts["delta_delta_g"].mean()
        max_ddg = muts["delta_delta_g"].max()
        print(f"{drug:<6}{docked:>10}{failed:>10}{mean_ddg:>12.3f}{max_ddg:>12.3f}")
    print("=" * 64)
    print(f"Total wall time: {(time.time()-t_start)/60:.1f} min")
    print(f"Saved: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
