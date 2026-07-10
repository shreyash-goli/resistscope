"""Parse a Stanford HIVdb genotype-phenotype dataset into per-drug panels.

The dataset (e.g. ``data/raw/PI_DataSet.txt`` for protease,
``data/rt/../NNRTI_DataSet.txt`` for reverse transcriptase) is a tab-separated
file with one viral isolate per row and these columns:

- ``SeqID``: isolate identifier.
- Drug columns of **linear** fold-resistance (protease: ``FPV, ATV, IDV, LPV,
  NFV, SQV, TPV, DRV``; NNRTI: ``NVP, EFV, ETR, RPV``). Missing values are ``NA``.
- Position columns ``P1``..``P{n}`` (n=99 for protease, 240 for RT): the amino
  acid observed at each residue. ``-`` means "same as the consensus/reference"
  (wildtype); a letter is the observed residue; a cell may contain several
  letters for a sequence mixture (e.g. ``IV``). Non-amino-acid codes (``.``,
  ``X``, ``#`` insertion, ``~`` deletion, ``*`` stop) are ignored.
- ``CompMutList``: a pre-formatted mutation list, kept only as a cross-check.

Everything target-specific — the reference sequence, the drug columns, the
number of position columns, and the "primary" mutation set — comes from the
active :class:`~targets.Target` (``config.ACTIVE_TARGET`` by default), so the
same code builds protease *or* RT panels. Two things differ from a naive
reading and are handled here:

1. Mutations are derived from the ``P1``..``P{n}`` columns using the target's
   reference sequence (NOT from binary ``10F``-style columns, which these files
   do not have).
2. Fold-resistance is stored linearly, so it is ``log10``-transformed to
   produce ``mean_log_fold_resistance``.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import config
from targets import Target

# The 20 standard amino acids; anything else in a position cell is not a point
# substitution we can name (mixtures are split into individual letters first).
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def _resolve(target: Target | None) -> Target:
    """Default an optional target to the active one."""
    return target if target is not None else config.ACTIVE_TARGET


def _reference_aa(position: int, target: Target) -> str:
    """Return the wildtype (reference) amino acid at a 1-indexed position."""
    return target.reference_seq[position - 1]


def _drug_column_for(drug: str, target: Target | None = None) -> str | None:
    """Map a drug abbreviation to its dataset column, or None if absent.

    A drug in the target's panel (e.g. RTV, DOR) that has no fold-resistance
    column in the raw dataset returns None and is skipped downstream.
    """
    t = _resolve(target)
    return drug if drug in t.drug_columns else None


def _mutations_in_cell(position: int, cell: str, target: Target) -> list[str]:
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
    wt = _reference_aa(position, target)
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


def parse_dataset(filepath: Path, target: Target | None = None) -> pd.DataFrame:
    """Parse a raw genotype-phenotype TSV into a clean per-isolate DataFrame.

    Returns one row per isolate with columns:

    - ``isolate_id`` (str)
    - one numeric column per dataset drug (from ``target.drug_columns``),
      holding the raw **linear** fold-resistance with ``NA`` parsed as ``NaN``
    - ``mutations``: a ``list[str]`` of mutation strings (e.g. ``["D30N",
      "M46I"]``) derived from the ``P1``..``P{n}`` genotype columns

    The ``log10`` transform is applied later, in :func:`build_mutation_panel`.
    """
    t = _resolve(target)
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(
            f"{t.label} dataset not found at {filepath}. "
            f"Run scripts/01_download_data.py first."
        )

    drug_columns = list(t.drug_columns)
    position_columns = [f"P{i}" for i in range(1, t.n_positions + 1)]

    # Read everything as string so the position columns keep '-'/letters intact
    # and 'NA'/'' become NaN.
    df = pd.read_csv(
        filepath,
        sep="\t",
        dtype=str,
        na_values=["NA", ""],
        keep_default_na=True,
    )

    present_drugs = [c for c in drug_columns if c in df.columns]
    present_positions = [c for c in position_columns if c in df.columns]

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
            muts.extend(_mutations_in_cell(pos, row[col], t))
        return muts

    out["mutations"] = df[present_positions].apply(_row_mutations, axis=1)

    return out


# Backward-compatible alias (protease-era name).
def parse_pi_dataset(filepath: Path, target: Target | None = None) -> pd.DataFrame:
    """Deprecated alias for :func:`parse_dataset`."""
    return parse_dataset(filepath, target)


def build_mutation_panel(
    df: pd.DataFrame,
    drug: str,
    min_isolates: int = 3,
    target: Target | None = None,
) -> pd.DataFrame:
    """Compute per-mutation summary statistics for one drug.

    ``df`` is the output of :func:`parse_dataset`. For every mutation, this
    averages the ``log10`` fold-resistance over the isolates that (a) carry the
    mutation and (b) have a non-missing, positive fold-resistance for ``drug``.

    Returns a DataFrame with columns:
    ``mutation`` (e.g. ``"V82A"``), ``position`` (int), ``wildtype_aa`` (str),
    ``mutant_aa`` (str), ``mean_log_fold_resistance`` (float),
    ``n_isolates`` (int), ``is_primary`` (bool, True for the target's major DRMs
    per ``target.primary_mutations``) — keeping only mutations seen in
    >= ``min_isolates`` qualifying isolates, sorted by position then mutant
    residue.

    Returns an empty (correctly-typed) DataFrame if ``drug`` has no column in
    the dataset (e.g. RTV / DOR).
    """
    t = _resolve(target)
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

    col = _drug_column_for(drug, t)
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
    grouped["is_primary"] = grouped["mutation"].isin(t.primary_mutations)

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


def build_all_panels(
    output_dir: Path | None = None,
    target: Target | None = None,
) -> dict[str, Path]:
    """Build and save a mutation panel for every drug in the target's panel.

    Parses the raw dataset once, builds one panel per drug, and writes each to
    ``{output_dir}/{drug}.parquet``. Drugs with no column in the dataset (e.g.
    RTV for protease, DOR for RT) are skipped with a warning and omitted from
    the returned mapping. ``output_dir`` defaults to the target's panels dir.

    Returns ``{drug: parquet_path}`` for the panels actually written.
    """
    t = _resolve(target)
    output_dir = Path(output_dir) if output_dir is not None else t.panels_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_path = config.RAW_DIR / t.dataset_filename
    df = parse_dataset(raw_path, t)

    paths: dict[str, Path] = {}
    for drug in t.drugs:
        col = _drug_column_for(drug, t)
        if col is None or col not in df.columns:
            print(
                f"  [skip] {drug} ({t.drugs[drug]}): "
                f"no fold-resistance column in the dataset"
            )
            continue
        panel = build_mutation_panel(df, drug, target=t)
        out_path = output_dir / f"{drug}.parquet"
        panel.to_parquet(out_path, index=False)
        n_isolates = int((df[col].notna() & (df[col] > 0)).sum())
        print(
            f"  [ok]   {drug} ({t.drugs[drug]}): "
            f"{len(panel)} mutations from {n_isolates} isolates -> {out_path.name}"
        )
        paths[drug] = out_path

    return paths


def load_panel(drug: str, panels_dir: Path | None = None) -> pd.DataFrame:
    """Load a precomputed per-drug mutation panel from parquet.

    ``panels_dir`` defaults to the active target's panels dir.
    """
    if panels_dir is None:
        panels_dir = config.PANELS_DIR
    path = Path(panels_dir) / f"{drug}.parquet"
    if not path.exists():
        raise FileNotFoundError(
            f"No panel for {drug} at {path}. Run scripts/02_build_panels.py first."
        )
    return pd.read_parquet(path)
