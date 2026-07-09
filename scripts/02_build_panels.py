"""02: Parse the Rhee dataset into per-drug mutation panels.

Builds one parquet panel per protease inhibitor in ``config.PI_DRUGS`` from
``data/raw/PI_DataSet.txt`` and prints a summary (drug, #mutations, #isolates).

Run from anywhere::

    python scripts/02_build_panels.py
"""

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
    parse_pi_dataset,
)


def main() -> int:
    print(f"Reading dataset: {config.RAW_DIR / 'PI_DataSet.txt'}")
    print(f"Writing panels to: {config.PANELS_DIR}\n")

    print("Building panels...")
    paths = build_all_panels(config.PANELS_DIR)

    # Isolate counts per drug (isolates with a usable fold-resistance value).
    df = parse_pi_dataset(config.RAW_DIR / "PI_DataSet.txt")

    print("\nSummary")
    print("-" * 60)
    print(f"{'drug':<6}{'name':<16}{'mutations':>12}{'isolates':>12}")
    print("-" * 60)
    for drug, name in config.PI_DRUGS.items():
        if drug not in paths:
            print(f"{drug:<6}{name:<16}{'—':>12}{'(no data)':>12}")
            continue
        panel = load_panel(drug, config.PANELS_DIR)
        col = _drug_column_for(drug)
        n_isolates = int((df[col].notna() & (df[col] > 0)).sum())
        print(f"{drug:<6}{name:<16}{len(panel):>12}{n_isolates:>12}")
    print("-" * 60)
    print(f"\nWrote {len(paths)} panel(s) to {config.PANELS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
