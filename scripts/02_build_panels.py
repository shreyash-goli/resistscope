"""02: Parse a Stanford HIVdb dataset into per-drug mutation panels.

Builds one parquet panel per drug in the target's panel (protease inhibitors or
NNRTIs) from the target's raw dataset and prints a summary (drug, #mutations,
#isolates).

Run from anywhere::

    python scripts/02_build_panels.py              # protease (default)
    python scripts/02_build_panels.py --target rt  # reverse transcriptase (NNRTIs)
"""

import argparse
import sys
from pathlib import Path

# Make the project root importable so ``import config`` / ``services`` work
# regardless of the current working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402
from services.mutation_panel import (  # noqa: E402
    _drug_column_for,
    build_all_panels,
    load_panel,
    parse_dataset,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target", default="HIV1_PR",
        help="Docking target: HIV1_PR / pr (default) or HIV1_RT / rt.",
    )
    args = parser.parse_args()
    t = config.set_active_target(args.target)

    raw_path = config.RAW_DIR / t.dataset_filename
    print(f"Target: {t.name} ({t.label})")
    print(f"Reading dataset: {raw_path}")
    print(f"Writing panels to: {t.panels_dir}\n")

    print("Building panels...")
    paths = build_all_panels(target=t)

    # Isolate counts per drug (isolates with a usable fold-resistance value).
    df = parse_dataset(raw_path, t)

    print("\nSummary")
    print("-" * 60)
    print(f"{'drug':<6}{'name':<16}{'mutations':>12}{'isolates':>12}")
    print("-" * 60)
    for drug, name in t.drugs.items():
        if drug not in paths:
            print(f"{drug:<6}{name:<16}{'—':>12}{'(no data)':>12}")
            continue
        panel = load_panel(drug, t.panels_dir)
        col = _drug_column_for(drug, t)
        n_isolates = int((df[col].notna() & (df[col] > 0)).sum())
        print(f"{drug:<6}{name:<16}{len(panel):>12}{n_isolates:>12}")
    print("-" * 60)
    print(f"\nWrote {len(paths)} panel(s) to {t.panels_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
