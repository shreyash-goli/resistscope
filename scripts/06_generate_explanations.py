"""06: Generate per-mutation mechanistic explanations via Claude.

For every (drug, mutation) pair in the docking results with |delta_delta_g|
above a threshold, build the structural context and generate a cached Claude
explanation. Re-runs are cheap: cached pairs are skipped without an API call.

Requires ANTHROPIC_API_KEY (or an ``ant auth login`` profile).

Usage::

    python scripts/06_generate_explanations.py
    python scripts/06_generate_explanations.py --min-ddg 0.5 --limit 20
    python scripts/06_generate_explanations.py --drugs DRV --primary-only
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

import config  # noqa: E402
from services.explanation import build_structural_context, generate_explanation  # noqa: E402
from services.mutation_panel import load_panel  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--docking", type=Path,
                        default=config.DOCKING_DIR / "benchmark_docking.parquet")
    parser.add_argument("--min-ddg", type=float, default=config.DDG_WARNING_THRESHOLD,
                        help="Only explain pairs with |delta_delta_g| >= this (kcal/mol).")
    parser.add_argument("--drugs", nargs="+", default=None,
                        help="Restrict to these drug abbreviations.")
    parser.add_argument("--primary-only", action="store_true",
                        help="Only explain primary (major) DRMs.")
    parser.add_argument("--ground-truth", action="store_true",
                        help="Explain exactly the (drug, mutation) pairs in "
                             "mechanism_ground_truth.json (ignores --min-ddg) so "
                             "the faithfulness eval (scripts/07) has full coverage.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Cap the number of pairs (for a cheap test run).")
    parser.add_argument("--sleep", type=float, default=0.5,
                        help="Seconds to sleep between API calls (rate limiting).")
    parser.add_argument("--cite", action="store_true",
                        help="Ground each explanation in real PubMed literature "
                             "(NCBI E-utilities) and store citations.")
    args = parser.parse_args()

    if not args.docking.exists():
        print(f"ERROR: docking results not found at {args.docking}. Run scripts/04 first.")
        return 1

    docking = pd.read_parquet(args.docking)
    docking = docking[docking["mutation"] != "WT"]
    docking = docking[docking["delta_delta_g"].notna()]
    if args.drugs:
        docking = docking[docking["drug"].isin(args.drugs)]

    # Optionally restrict to primary DRMs (join is_primary from panels).
    if args.primary_only:
        primary = set()
        for drug in docking["drug"].unique():
            panel = load_panel(drug, config.PANELS_DIR)
            for m in panel[panel["is_primary"]]["mutation"]:
                primary.add((drug, m))
        docking = docking[
            docking.apply(lambda r: (r["drug"], r["mutation"]) in primary, axis=1)
        ]

    if args.ground_truth:
        import json
        gt = json.loads((config.DATA_DIR / "mechanism_ground_truth.json").read_text())
        gt_pairs = {(d, mut) for mut, v in gt.items()
                    for d in (v.get("affects_drugs") or list(config.PI_DRUGS))}
        pairs = docking[
            docking.apply(lambda r: (r["drug"], r["mutation"]) in gt_pairs, axis=1)
        ].copy()
    else:
        pairs = docking[docking["delta_delta_g"].abs() >= args.min_ddg].copy()
    pairs = pairs.sort_values("delta_delta_g", ascending=False)
    if args.limit:
        pairs = pairs.head(args.limit)

    print(f"Explaining {len(pairs)} (drug, mutation) pairs "
          f"(|ddG| >= {args.min_ddg}, model={config.CLAUDE_MODEL})\n")

    n_new = n_cached = n_failed = 0
    for _, row in pairs.iterrows():
        drug, mut, ddg = row["drug"], row["mutation"], float(row["delta_delta_g"])
        cache_path = config.EXPLANATIONS_DIR / f"{drug}_{mut}.json"
        cached = cache_path.exists()
        try:
            context = build_structural_context(
                mut, {"delta_delta_g": ddg, "delta_g": row.get("delta_g")}
            )
            expl = generate_explanation(drug, mut, ddg, context, cite=args.cite)
            if cached:
                n_cached += 1
            else:
                n_new += 1
                print(f"  [{drug} {mut}] ddG={ddg:+.2f}: {expl[:100]}...")
                time.sleep(args.sleep)
        except Exception as exc:  # noqa: BLE001 - keep the batch alive
            msg = str(exc)
            if "authentication" in msg.lower() or "api_key" in msg.lower():
                print("\nERROR: no Anthropic credentials found. Set ANTHROPIC_API_KEY "
                      "(or run `ant auth login`) and re-run — cached pairs will be "
                      "skipped so no work is repeated.")
                return 1
            n_failed += 1
            print(f"  FAILED {drug} {mut}: {type(exc).__name__}: {exc}")

    print(f"\nDone. Generated {n_new} new, {n_cached} from cache, {n_failed} failed.")
    print(f"Explanations cached in {config.EXPLANATIONS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
