"""Delta-delta-G robustness scoring, baselines, and aggregation.

Two distinct jobs live here:

1. **Compound-level robustness score** (:func:`compute_robustness_scores`) — the
   tool's user-facing output: aggregate a compound's per-mutation ddG across the
   panel into a single 0-100 robustness number for triage.

2. **Baselines** for the validation ablation (:func:`compute_baseline_*`) — cheap
   reference predictions the docking method must beat to be worth anything.

Sign convention: ``delta_delta_g > 0`` means the mutant binds WORSE than
wildtype (predicted resistance). ``mean_log_fold_resistance`` is the measured
clinical readout (higher = more resistant).
"""

import numpy as np
import pandas as pd

import config


def _merge_docking_panel(docking_results: pd.DataFrame, panel: pd.DataFrame) -> pd.DataFrame:
    """Inner-join a drug's docking rows with its panel on mutation (excl. WT).

    Keeps only mutations with a usable (non-null) delta_delta_g.
    """
    muts = docking_results[docking_results["mutation"] != "WT"].copy()
    merged = muts.merge(
        panel[["mutation", "mean_log_fold_resistance", "n_isolates", "is_primary"]],
        on="mutation",
        how="inner",
    )
    return merged


def compute_robustness_scores(
    docking_results: pd.DataFrame,
    panel: pd.DataFrame,
) -> dict:
    """Aggregate one compound's docking results into robustness metrics.

    ``docking_results`` holds the rows for a single drug (WT + mutants).
    Returns a dict with simple/prevalence-weighted/worst-case ddG, a normalized
    0-100 robustness score (100 = binding unaffected across the panel), and
    scored/failed counts.
    """
    merged = _merge_docking_panel(docking_results, panel)
    scored = merged.dropna(subset=["delta_delta_g"])

    n_failed = int(merged["delta_delta_g"].isna().sum())
    n_scored = int(len(scored))

    if n_scored == 0:
        return {
            "simple_mean_ddg": None, "prevalence_weighted_ddg": None,
            "worst_case_ddg": None, "robustness_0_100": None,
            "n_mutations_scored": 0, "n_mutations_failed": n_failed,
        }

    ddg = scored["delta_delta_g"].to_numpy(dtype=float)
    weights = scored["n_isolates"].to_numpy(dtype=float)

    simple_mean = float(np.mean(ddg))
    weighted = float(np.average(ddg, weights=weights)) if weights.sum() > 0 else simple_mean
    worst = float(np.max(ddg))
    # Map 0 kcal/mol -> 100, 3 kcal/mol -> 0, clamped to [0, 100].
    robustness = max(0.0, min(100.0, 100.0 - weighted * 33.3))

    return {
        "simple_mean_ddg": simple_mean,
        "prevalence_weighted_ddg": weighted,
        "worst_case_ddg": worst,
        "robustness_0_100": robustness,
        "n_mutations_scored": n_scored,
        "n_mutations_failed": n_failed,
    }


def compute_baseline_mutation_count(panel: pd.DataFrame) -> float:
    """Baseline 1: number of high-resistance mutations known for the drug.

    Counts panel mutations with ``mean_log_fold_resistance > 1.0`` (>10x fold
    resistance). A docking-free, data-only reference.
    """
    return float((panel["mean_log_fold_resistance"] > 1.0).sum())


def compute_baseline_wt_only(docking_results: pd.DataFrame) -> float:
    """Baseline 2: the wildtype binding affinity only (mutation-blind).

    Tests whether simply knowing how well a compound binds wildtype is
    informative, ignoring all mutation data. Returns the WT delta_g.
    """
    wt = docking_results[docking_results["mutation"] == "WT"]
    if wt.empty or wt["delta_g"].isna().all():
        return float("nan")
    return float(wt["delta_g"].iloc[0])
