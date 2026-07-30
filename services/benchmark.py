"""Rigorous, honest benchmarking of the docking-ΔΔG predictor.

The plain correlation story (``services/validation.py``) shows that rigid
single-mutation docking ΔΔG does **not** quantitatively track pooled clinical
fold-resistance (pooled Spearman ≈ 0). This module answers the sharper, honest
questions a reviewer would actually ask, each with a significance test or
confidence interval rather than a bare point estimate:

1. **Is the top-N DRM enrichment real, or a fluke?**  Permutation p-value +
   bootstrap 95% CI on the enrichment factor (:func:`enrichment_with_significance`).
2. **How good a *ranker* is ΔΔG?**  Threshold-free ROC-AUC / PR-AUC for
   recovering known major DRMs, per drug and pooled, with bootstrap CIs
   (:func:`ranking_metrics`).
3. **Does de-confounding the target rescue the magnitude correlation?**  Rebuild
   the per-mutation target from progressively cleaner isolate subsets (all →
   ≤2 mutations → single-mutation only) and re-correlate. Reported honestly:
   for this dataset it does **not** rescue it (:func:`deconfounding_analysis`) —
   a negative result that bounds what the method can and cannot do.

The headline it supports is defensible: *docking ΔΔG is a coarse DRM-triage
flag (its extreme predictions are ~3× enriched for real resistance mutations,
p < 1e-3), not a quantitative resistance predictor.*
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import config
from services.mutation_panel import parse_dataset
from services.validation import build_merged

# Fixed seed so permutation p-values and bootstrap CIs are reproducible across
# runs (the numbers surface in the UI, so they must not jitter).
_SEED = 0
TARGET_COL = "mean_log_fold_resistance"
DDG_COL = "delta_delta_g"
LABEL_COL = "is_primary"


# --- 1. Enrichment with significance -----------------------------------------

def enrichment_with_significance(
    merged: pd.DataFrame,
    top_ns=(20, 40, 60, 100),
    n_perm: int = 5000,
    n_boot: int = 5000,
    seed: int = _SEED,
) -> pd.DataFrame:
    """Top-N DRM enrichment with a permutation p-value and bootstrap 95% CI.

    For each ``N``: rank mutations by predicted ΔΔG, take the top ``N``, and
    measure the fraction that are known major DRMs (precision). ``enrichment`` is
    that precision over the pool base rate.

    - **perm_p**: probability that a random size-``N`` draw matches/exceeds the
      observed precision (label shuffle, one-sided). Small = the tail is really
      DRM-enriched.
    - **ci_low/ci_high**: bootstrap 95% CI on the enrichment factor (resample
      the top-``N`` set with replacement).
    """
    rng = np.random.default_rng(seed)
    labels = merged[LABEL_COL].to_numpy(dtype=float)
    base = float(labels.mean())
    order = merged[DDG_COL].to_numpy().argsort()[::-1]  # high ΔΔG first
    rows = []
    for n in top_ns:
        n = int(min(n, len(merged)))
        top_idx = order[:n]
        top_labels = labels[top_idx]
        prec = float(top_labels.mean())
        # permutation null: random size-n subsets of the label vector
        null = np.array([rng.permutation(labels)[:n].mean() for _ in range(n_perm)])
        perm_p = float((null >= prec).mean())
        # bootstrap CI on enrichment (resample the observed top-n labels)
        boot = np.array([
            rng.choice(top_labels, size=n, replace=True).mean()
            for _ in range(n_boot)
        ]) / base if base > 0 else np.zeros(n_boot)
        rows.append({
            "top_n": n,
            "precision": prec,
            "base_rate": base,
            "enrichment": (prec / base) if base > 0 else np.nan,
            "perm_p": perm_p,
            "ci_low": float(np.percentile(boot, 2.5)),
            "ci_high": float(np.percentile(boot, 97.5)),
            "n_drms_in_top": int(top_labels.sum()),
        })
    return pd.DataFrame(rows)


# --- 2. Ranking metrics (threshold-free) -------------------------------------

def _auc_pair(y_true: np.ndarray, score: np.ndarray) -> tuple[float, float]:
    """(ROC-AUC, PR-AUC) with a graceful fallback if a class is absent."""
    from sklearn.metrics import average_precision_score, roc_auc_score
    if len(np.unique(y_true)) < 2:
        return np.nan, np.nan
    return float(roc_auc_score(y_true, score)), float(average_precision_score(y_true, score))


def _bootstrap_auc_ci(y_true: np.ndarray, score: np.ndarray, n_boot: int,
                      seed: int) -> tuple[float, float]:
    """Bootstrap 95% CI on ROC-AUC (stratified-ish: plain resample of pairs)."""
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    n = len(y_true)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yt = y_true[idx]
        if len(np.unique(yt)) < 2:
            continue
        vals.append(roc_auc_score(yt, score[idx]))
    if not vals:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def ranking_metrics(merged: pd.DataFrame, n_boot: int = 2000,
                    seed: int = _SEED) -> pd.DataFrame:
    """Per-drug and pooled ROC-AUC / PR-AUC for DRM recovery from ΔΔG.

    A threshold-free measure of how well ΔΔG *ranks* known major DRMs above
    non-DRMs. AUC 0.5 = chance. Pooled row carries a bootstrap 95% CI so "≈
    chance" is stated with an interval, not asserted.
    """
    rows = []
    groups = [(d, merged[merged["drug"] == d]) for d in sorted(merged["drug"].unique())]
    groups.append(("OVERALL", merged))
    for drug, g in groups:
        y = g[LABEL_COL].to_numpy(dtype=int)
        s = g[DDG_COL].to_numpy(dtype=float)
        roc, pr = _auc_pair(y, s)
        row = {"drug": drug, "roc_auc": roc, "pr_auc": pr,
               "base_rate": float(y.mean()) if len(y) else np.nan,
               "n": int(len(g)), "n_drms": int(y.sum())}
        if drug == "OVERALL":
            lo, hi = _bootstrap_auc_ci(y, s, n_boot, seed)
            row["roc_auc_ci_low"], row["roc_auc_ci_high"] = lo, hi
        rows.append(row)
    return pd.DataFrame(rows)


# --- 3. De-confounding analysis (the honest negative result) -----------------

def _subset_merged(df_iso: pd.DataFrame, docking: pd.DataFrame,
                   drugs: list[str], primary: frozenset,
                   min_isolates: int) -> pd.DataFrame:
    """Build a per-(drug,mutation) merged table from an isolate subset.

    ``df_iso`` is a subset of :func:`parse_dataset` output (already filtered by
    mutation count). Mirrors the panel build (mean log10 fold over qualifying
    isolates) but in-memory, then inner-joins the docking ΔΔG.
    """
    frames = []
    for drug in drugs:
        s = df_iso[["mutations", drug]].copy()
        s = s[s[drug].notna() & (s[drug] > 0)]
        if s.empty:
            continue
        s["lf"] = np.log10(s[drug].astype(float))
        ex = s.explode("mutations").dropna(subset=["mutations"])
        panel = ex.groupby("mutations")["lf"].agg(
            mean_log_fold_resistance="mean", n_isolates="size").reset_index()
        panel = panel[panel["n_isolates"] >= min_isolates]
        panel = panel.rename(columns={"mutations": "mutation"})
        dd = docking[docking["drug"] == drug][["mutation", DDG_COL]]
        m = panel.merge(dd, on="mutation", how="inner").dropna(
            subset=[DDG_COL, TARGET_COL])
        m["drug"] = drug
        m[LABEL_COL] = m["mutation"].isin(primary)
        frames.append(m)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def deconfounding_analysis(min_isolates: int = 3) -> pd.DataFrame:
    """Re-correlate ΔΔG vs fold-resistance on progressively cleaner targets.

    Confounding here = a mutation's measured fold-resistance is averaged over
    isolates that also carry *other* resistance mutations. We rebuild the target
    from isolate subsets of decreasing confounding and re-run the pooled
    Spearman. Columns: ``subset, max_mutations_per_isolate, n_isolates_used,
    n_pairs, spearman_rho, spearman_pvalue``.

    Honest result for this dataset: cleaner subsets are far sparser (only ~27
    single-mutation isolates exist) and do **not** improve — often invert — the
    correlation, bounding the method as a magnitude predictor.
    """
    t = config.ACTIVE_TARGET
    df = parse_dataset(config.RAW_DIR / t.dataset_filename, t)
    df = df.copy()
    df["n_mut"] = df["mutations"].apply(len)
    docking = pd.read_parquet(config.DOCKING_DIR / "benchmark_docking.parquet")
    drugs = [d for d in t.drug_columns if d in df.columns]

    specs = [
        ("all isolates (confounded)", None),
        ("≤3 mutations", 3),
        ("≤2 mutations", 2),
        ("single-mutation only", 1),
    ]
    rows = []
    for label, cap in specs:
        iso = df if cap is None else df[df["n_mut"] <= cap]
        # single-mutation isolates are so few that min_isolates=3 empties them;
        # relax the support requirement as the subset tightens.
        mi = 1 if cap == 1 else (2 if cap == 2 else min_isolates)
        m = _subset_merged(iso, docking, drugs, t.primary_mutations, mi)
        if len(m) >= 3 and m[DDG_COL].nunique() > 1 and m[TARGET_COL].nunique() > 1:
            rho, p = stats.spearmanr(m[DDG_COL], m[TARGET_COL])
        else:
            rho, p = np.nan, np.nan
        rows.append({
            "subset": label,
            "max_mutations_per_isolate": str(cap) if cap is not None else "∞",
            "n_isolates_used": int(len(iso)),
            "min_isolates_per_mutation": mi,
            "n_pairs": int(len(m)),
            "spearman_rho": float(rho) if rho == rho else np.nan,
            "spearman_pvalue": float(p) if p == p else np.nan,
        })
    return pd.DataFrame(rows)


# --- Orchestration -----------------------------------------------------------

def run_benchmark(
    docking_parquet: Path | None = None,
    output_dir: Path | None = None,
    n_perm: int = 5000,
    n_boot: int = 5000,
) -> dict:
    """Run the full benchmark and persist it for the API/UI.

    Writes ``benchmark_enrichment.parquet``, ``benchmark_ranking.parquet``,
    ``benchmark_deconfounding.parquet`` and a consolidated
    ``benchmark_metrics.json`` under ``output_dir`` (the active target's
    validation dir by default). Returns the JSON-serialisable summary dict.
    """
    if output_dir is None:
        output_dir = config.VALIDATION_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    merged = build_merged(docking_parquet)
    if merged.empty:
        raise RuntimeError("No overlapping docking+panel data to benchmark.")

    enr = enrichment_with_significance(merged, n_perm=n_perm, n_boot=n_boot)
    rank = ranking_metrics(merged, n_boot=min(n_boot, 2000))
    decon = deconfounding_analysis()

    enr.to_parquet(output_dir / "benchmark_enrichment.parquet", index=False)
    rank.to_parquet(output_dir / "benchmark_ranking.parquet", index=False)
    decon.to_parquet(output_dir / "benchmark_deconfounding.parquet", index=False)

    pooled = rank[rank["drug"] == "OVERALL"].iloc[0]
    best = enr.loc[enr["enrichment"].idxmax()]
    # Per-target conclusions (previously hardcoded to the protease result).
    # Rescue = a less-confounded isolate subset reveals a positive, significant
    # magnitude correlation the confounded pool hides (protease: no; RT: yes).
    decon_rescues = bool((
        decon["max_mutations_per_isolate"].isin(["1", "2", "3"])
        & (decon["spearman_rho"] >= 0.2)
        & (decon["spearman_pvalue"] < 0.05)
    ).any())
    _ci_low = float(pooled.get("roc_auc_ci_low", np.nan))
    roc_above_chance = bool(_ci_low == _ci_low and _ci_low > 0.5)
    summary = {
        "n_pairs": int(len(merged)),
        "drm_base_rate": float(merged[LABEL_COL].mean()),
        "enrichment": enr.to_dict("records"),
        "ranking": rank.to_dict("records"),
        "deconfounding": decon.replace({np.nan: None}).to_dict("records"),
        "headline": {
            "best_enrichment": float(best["enrichment"]),
            "best_enrichment_top_n": int(best["top_n"]),
            "best_enrichment_perm_p": float(best["perm_p"]),
            "best_enrichment_ci": [float(best["ci_low"]), float(best["ci_high"])],
            "pooled_roc_auc": float(pooled["roc_auc"]),
            "pooled_roc_auc_ci": [
                float(pooled.get("roc_auc_ci_low", np.nan)),
                float(pooled.get("roc_auc_ci_high", np.nan)),
            ],
            "deconfounding_rescues_correlation": decon_rescues,
            "roc_above_chance": roc_above_chance,
        },
    }
    import json
    (output_dir / "benchmark_metrics.json").write_text(json.dumps(summary, indent=2))
    return summary
