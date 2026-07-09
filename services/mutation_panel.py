"""Parse the Rhee PI dataset and build per-drug mutation panels.

The Stanford HIVdb PI genotype-phenotype dataset (``data/raw/PI_DataSet.txt``)
is a tab-separated file with one viral isolate per row and these columns:

- ``SeqID``: isolate identifier.
- Eight drug columns of **linear** fold-resistance: ``FPV, ATV, IDV, LPV,
  NFV, SQV, TPV, DRV``. Missing values are encoded as ``NA``.
- Position columns ``P1``..``P99``: the amino acid observed at each HIV-1
  protease residue. ``-`` means "same as the consensus/reference" (wildtype);
  a letter is the observed residue; a cell may contain several letters for a
  sequence mixture (e.g. ``IV``). Non-amino-acid codes (``.``, ``X``, ``#``
  insertion, ``~`` deletion, ``*`` stop) are ignored.
- ``CompMutList``: a pre-formatted mutation list, kept only as a cross-check.

Two things differ from a naive reading of the dataset and are handled here:

1. Mutations are derived from the ``P1``..``P99`` columns using the HIV-1
   protease reference sequence (NOT from binary ``10F``-style columns, which
   this file does not have).
2. Fold-resistance is stored linearly, so it is ``log10``-transformed to
   produce ``mean_log_fold_resistance``.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import config

# HIV-1 protease consensus/reference sequence, residues 1-99 (1-indexed).
# Used to determine the wildtype amino acid at each position: a value "A" in
# column "P82" with reference V at position 82 means the mutation V82A.
REFERENCE_SEQUENCE = (
    "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQ"
    "ILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
)
assert len(REFERENCE_SEQUENCE) == 99, "HIV-1 protease reference must be 99 residues"

# The 20 standard amino acids; anything else in a position cell is not a point
# substitution we can name (mixtures are split into individual letters first).
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# Drug fold-resistance columns present in the raw dataset.
DATASET_DRUG_COLUMNS = ["FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"]

# Column-name prefix for the per-position genotype columns.
POSITION_COLUMNS = [f"P{i}" for i in range(1, 100)]

# Alias map from a config drug abbreviation to the dataset column name, for the
# cases where they differ. Most match directly; RTV (ritonavir) is not tested
# in the Rhee 2006 dataset and has no column.
DRUG_COLUMN_ALIASES: dict[str, str] = {
    # config_abbrev: dataset_column
}


def _reference_aa(position: int) -> str:
    """Return the wildtype (reference) amino acid at a 1-indexed position."""
    return REFERENCE_SEQUENCE[position - 1]


def _drug_column_for(drug: str) -> str | None:
    """Map a config drug abbreviation to its dataset column, or None if absent."""
    if drug in DATASET_DRUG_COLUMNS:
        return drug
    return DRUG_COLUMN_ALIASES.get(drug)


def _mutations_in_cell(position: int, cell: str) -> list[str]:
    """Return the mutation strings implied by one position cell.

    ``cell`` is the raw value of a ``P{position}`` column. ``-`` (wildtype),
    empty, or non-amino-acid codes yield no mutations. A mixture like ``IV``
    yields one mutation per distinct residue that differs from wildtype.
    """
    if not isinstance(cell, str):
        return []
    cell = cell.strip().upper()
    if cell in ("", "-", "."):
        return []
    wt = _reference_aa(position)
    muts: list[str] = []
    seen: set[str] = set()
    for aa in cell:
        if aa not in STANDARD_AMINO_ACIDS:
            continue  # skip X, insertion/deletion markers, stops, etc.
        if aa == wt:
            continue  # matches wildtype -> not a mutation
        if aa in seen:
            continue
        seen.add(aa)
        muts.append(f"{wt}{position}{aa}")
    return muts


def parse_pi_dataset(filepath: Path) -> pd.DataFrame:
    """Parse the raw PI TSV into a clean per-isolate DataFrame.

    Returns one row per isolate with columns:

    - ``isolate_id`` (str)
    - one numeric column per dataset drug (``FPV``..``DRV``), holding the raw
      **linear** fold-resistance with ``NA`` parsed as ``NaN``
    - ``mutations``: a ``list[str]`` of mutation strings (e.g. ``["D30N",
      "M46I"]``) derived from the ``P1``..``P99`` genotype columns

    The ``log10`` transform is applied later, in :func:`build_mutation_panel`.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(
            f"PI dataset not found at {filepath}. "
            f"Run scripts/01_download_data.py first."
        )

    # Read everything as string so the position columns keep '-'/letters intact
    # and 'NA'/'' become NaN.
    df = pd.read_csv(
        filepath,
        sep="\t",
        dtype=str,
        na_values=["NA", ""],
        keep_default_na=True,
    )

    present_drugs = [c for c in DATASET_DRUG_COLUMNS if c in df.columns]
    present_positions = [c for c in POSITION_COLUMNS if c in df.columns]

    out = pd.DataFrame()
    out["isolate_id"] = df.get("SeqID", pd.Series(range(len(df)))).astype(str)

    # Coerce drug columns to numeric linear fold-resistance.
    for drug in present_drugs:
        out[drug] = pd.to_numeric(df[drug], errors="coerce")

    # Derive per-isolate mutation lists from the position columns.
    def _row_mutations(row: pd.Series) -> list[str]:
        muts: list[str] = []
        for col in present_positions:
            pos = int(col[1:])
            muts.extend(_mutations_in_cell(pos, row[col]))
        return muts

    out["mutations"] = df[present_positions].apply(_row_mutations, axis=1)

    return out


def build_mutation_panel(
    df: pd.DataFrame,
    drug: str,
    min_isolates: int = 3,
) -> pd.DataFrame:
    """Compute per-mutation summary statistics for one drug.

    ``df`` is the output of :func:`parse_pi_dataset`. For every mutation, this
    averages the ``log10`` fold-resistance over the isolates that (a) carry the
    mutation and (b) have a non-missing, positive fold-resistance for ``drug``.

    Returns a DataFrame with columns:
    ``mutation`` (e.g. ``"V82A"``), ``position`` (int), ``wildtype_aa`` (str),
    ``mutant_aa`` (str), ``mean_log_fold_resistance`` (float),
    ``n_isolates`` (int), ``is_primary`` (bool, True for major PI DRMs per
    ``config.PRIMARY_PI_MUTATIONS``) — keeping only mutations seen in
    >= ``min_isolates`` qualifying isolates, sorted by position then mutant
    residue.

    Returns an empty (correctly-typed) DataFrame if ``drug`` has no column in
    the dataset (e.g. RTV in the Rhee 2006 data).
    """
    empty = pd.DataFrame(
        {
            "mutation": pd.Series(dtype="object"),
            "position": pd.Series(dtype="int64"),
            "wildtype_aa": pd.Series(dtype="object"),
            "mutant_aa": pd.Series(dtype="object"),
            "mean_log_fold_resistance": pd.Series(dtype="float64"),
            "n_isolates": pd.Series(dtype="int64"),
            "is_primary": pd.Series(dtype="bool"),
        }
    )

    col = _drug_column_for(drug)
    if col is None or col not in df.columns:
        return empty

    # Keep only isolates with a usable (positive) fold-resistance for this drug.
    sub = df[["mutations", col]].copy()
    sub = sub[sub[col].notna() & (sub[col] > 0)]
    if sub.empty:
        return empty
    sub["log_fold"] = np.log10(sub[col].astype(float))

    # Explode to one (isolate, mutation) pair per row, then aggregate.
    exploded = sub.explode("mutations").dropna(subset=["mutations"])
    if exploded.empty:
        return empty

    grouped = exploded.groupby("mutations")["log_fold"].agg(
        mean_log_fold_resistance="mean",
        n_isolates="size",
    )
    grouped = grouped[grouped["n_isolates"] >= min_isolates].reset_index()
    grouped = grouped.rename(columns={"mutations": "mutation"})

    if grouped.empty:
        return empty

    # Split the mutation string into components (wt / position / mutant).
    grouped["wildtype_aa"] = grouped["mutation"].str[0]
    grouped["mutant_aa"] = grouped["mutation"].str[-1]
    grouped["position"] = grouped["mutation"].str[1:-1].astype(int)
    grouped["is_primary"] = grouped["mutation"].isin(config.PRIMARY_PI_MUTATIONS)

    grouped = grouped[
        [
            "mutation",
            "position",
            "wildtype_aa",
            "mutant_aa",
            "mean_log_fold_resistance",
            "n_isolates",
            "is_primary",
        ]
    ]
    grouped["n_isolates"] = grouped["n_isolates"].astype("int64")
    grouped = grouped.sort_values(["position", "mutant_aa"]).reset_index(drop=True)
    return grouped


def build_all_panels(output_dir: Path = config.PANELS_DIR) -> dict[str, Path]:
    """Build and save a mutation panel for every drug in ``config.PI_DRUGS``.

    Parses the raw dataset once, builds one panel per drug, and writes each to
    ``{output_dir}/{drug}.parquet``. Drugs with no column in the dataset (e.g.
    RTV) are skipped with a warning and omitted from the returned mapping.

    Returns ``{drug: parquet_path}`` for the panels actually written.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_path = config.RAW_DIR / "PI_DataSet.txt"
    df = parse_pi_dataset(raw_path)

    paths: dict[str, Path] = {}
    for drug in config.PI_DRUGS:
        col = _drug_column_for(drug)
        if col is None or col not in df.columns:
            print(
                f"  [skip] {drug} ({config.PI_DRUGS[drug]}): "
                f"no fold-resistance column in the dataset"
            )
            continue
        panel = build_mutation_panel(df, drug)
        out_path = output_dir / f"{drug}.parquet"
        panel.to_parquet(out_path, index=False)
        n_isolates = int((df[col].notna() & (df[col] > 0)).sum())
        print(
            f"  [ok]   {drug} ({config.PI_DRUGS[drug]}): "
            f"{len(panel)} mutations from {n_isolates} isolates -> {out_path.name}"
        )
        paths[drug] = out_path

    return paths


def load_panel(drug: str, panels_dir: Path = config.PANELS_DIR) -> pd.DataFrame:
    """Load a precomputed per-drug mutation panel from parquet."""
    path = Path(panels_dir) / f"{drug}.parquet"
    if not path.exists():
        raise FileNotFoundError(
            f"No panel for {drug} at {path}. Run scripts/02_build_panels.py first."
        )
    return pd.read_parquet(path)
