"""03: Generate all mutant PDB + PDBQT structures.

Prepares the wildtype receptor (clean 3OXC -> add H + fix Asp25 protonation ->
PDBQT), then generates a point-mutant receptor (PDB + PDBQT) for every unique
mutation across all per-drug panels. Existing PDBQTs are skipped so the build is
resumable.

Usage::

    python scripts/03_build_mutant_cache.py              # full build (~286 mutants)
    python scripts/03_build_mutant_cache.py --limit 5    # quick smoke test
    python scripts/03_build_mutant_cache.py --force-wt   # rebuild wildtype too
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402
from services.structure_prep import (  # noqa: E402
    build_mutant_cache,
    collect_unique_mutations,
    prepare_wildtype,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Only process the first N unique mutations (for testing).",
    )
    parser.add_argument(
        "--force-wt", action="store_true",
        help="Force rebuild of the wildtype receptor even if it exists.",
    )
    args = parser.parse_args()

    print("=== Wildtype preparation ===")
    prepare_wildtype(force=args.force_wt)

    n_total = len(collect_unique_mutations(config.PANELS_DIR))
    n_run = n_total if args.limit is None else min(args.limit, n_total)
    print("\n=== Mutant cache ===")
    print(f"Unique mutations across all panels: {n_total}")
    print(f"Processing: {n_run}  (est. ~5-8 s each => ~{n_run * 7 // 60} min uncached)\n")

    pdbqt_paths = build_mutant_cache(config.PANELS_DIR, limit=args.limit)

    n_have = sum(1 for p in pdbqt_paths if Path(p).exists())
    print(f"\nDone. {n_have}/{n_run} mutant PDBQTs present in {config.MUTANTS_DIR}")
    if n_have < n_run:
        print(f"  {n_run - n_have} mutation(s) failed — see FAILED lines above.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
