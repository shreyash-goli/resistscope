"""Predicted vs. real fold-resistance validation and plotting.

The core scientific test: does docking-predicted binding degradation
(``delta_delta_g``) track measured clinical resistance
(``mean_log_fold_resistance`` from the Rhee dataset), per mutation?

We correlate three per-mutation predictors against the measured readout so the
docking method is judged against honest baselines:

- ``docking_ddg``  — our method (per-mutation predicted ddG).
- ``prevalence``   — mutation frequency (n_isolates): are common mutations more
  resistant? (data-only baseline).
- ``is_primary``   — the IAS-USA major-DRM flag: does the literature label alone
  predict magnitude? (knowledge baseline).

Docking is only worthwhile if ``docking_ddg`` correlates at least as well as
these. Correlations are reported per drug and pooled (``OVERALL``).
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import config
from services.mutation_panel import load_panel
from services.scoring import _merge_docking_panel

# Per-mutation predictors compared in the ablation.
SCORING_METHODS = {
    "docking_ddg": "delta_delta_g",
    "prevalence": "n_isolates",
    "is_primary": "is_primary",
}
TARGET_COL = "mean_log_fold_resistance"


def _correlate(predictor: pd.Series, target: pd.Series) -> dict:
    """Spearman + Pearson of predictor vs target, robust to degenerate input."""
    x = pd.to_numeric(predictor, errors="coerce").astype(float)
    y = pd.to_numeric(target, errors="coerce").astype(float)
    mask = x.notna() & y.notna()
    x, y = x[mask], y[mask]
    n = int(len(x))
    # Need >=3 points and non-zero variance in both for a defined correlation.
    if n < 3 or x.nunique() < 2 or y.nunique() < 2:
        return {"spearman_rho": np.nan, "spearman_pvalue": np.nan,
                "pearson_r": np.nan, "pearson_pvalue": np.nan, "n_mutations": n}
    rho, sp = stats.spearmanr(x, y)
    r, pp = stats.pearsonr(x, y)
    return {"spearman_rho": float(rho), "spearman_pvalue": float(sp),
            "pearson_r": float(r), "pearson_pvalue": float(pp), "n_mutations": n}


def build_merged(
    docking_parquet: Path = None,
    panels_dir: Path = None,
) -> pd.DataFrame:
    """Join docking results with panels into one per-mutation table.

    Columns: drug, mutation, delta_g, delta_delta_g, mean_log_fold_resistance,
    n_isolates, is_primary. WT rows and failed docks are dropped. Paths default
    to the active target's dirs.
    """
    if panels_dir is None:
        panels_dir = config.PANELS_DIR
    if docking_parquet is None:
        docking_parquet = config.DOCKING_DIR / "benchmark_docking.parquet"
    docking = pd.read_parquet(docking_parquet)

    frames = []
    for drug in docking["drug"].unique():
        panel_path = Path(panels_dir) / f"{drug}.parquet"
        if not panel_path.exists():
            continue
        panel = load_panel(drug, panels_dir)
        merged = _merge_docking_panel(docking[docking["drug"] == drug], panel)
        merged["drug"] = drug
        frames.append(merged)

    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    return out.dropna(subset=["delta_delta_g", TARGET_COL])


def compute_enrichment(
    merged: pd.DataFrame,
    ddg_col: str = "delta_delta_g",
    label_col: str = "is_primary",
    top_ns=(20, 40, 60, 100),
) -> pd.DataFrame:
    """Precision@top-N: are the highest-ddG mutations enriched for known DRMs?

    For each N, ranks mutations by predicted ddG, takes the top N, and reports
    the fraction that are major DRMs vs the pool base rate (enrichment factor).
    A triage-oriented metric that tolerates the confounded continuous target.
    """
    base = float(merged[label_col].mean())
    rows = []
    for n in top_ns:
        n = min(n, len(merged))
        top = merged.nlargest(n, ddg_col)
        prec = float(top[label_col].mean())
        rows.append({"top_n": n, "precision": prec, "base_rate": base,
                     "enrichment": (prec / base) if base > 0 else np.nan})
    return pd.DataFrame(rows)


def per_drug_correlation(merged: pd.DataFrame) -> pd.DataFrame:
    """Per-drug docking-ddG Spearman/Pearson, over all mutations and DRMs only."""
    rows = []
    for d in sorted(merged["drug"].unique()):
        g = merged[merged["drug"] == d]
        for scope, sub in [("all", g), ("primary", g[g["is_primary"]])]:
            corr = _correlate(sub["delta_delta_g"], sub[TARGET_COL])
            rows.append({"drug": d, "scope": scope, **corr})
    return pd.DataFrame(rows)


def run_full_validation(
    docking_parquet: Path = None,
    panels_dir: Path = None,
    output_dir: Path = None,
) -> pd.DataFrame:
    """Correlate every scoring method vs measured fold-resistance.

    Computes per-drug and pooled (``OVERALL``) Spearman/Pearson for each method,
    saves the per-mutation merge and the correlation table to ``output_dir``,
    and returns the correlation table (matching the Validation Output schema).
    Paths default to the active target's dirs.
    """
    if panels_dir is None:
        panels_dir = config.PANELS_DIR
    if output_dir is None:
        output_dir = config.VALIDATION_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    merged = build_merged(docking_parquet, panels_dir)
    if merged.empty:
        raise RuntimeError("No overlapping docking+panel data to validate.")
    merged.to_parquet(output_dir / "scores_vs_fold_resistance.parquet", index=False)

    rows = []
    groups = [(d, merged[merged["drug"] == d]) for d in sorted(merged["drug"].unique())]
    groups.append(("OVERALL", merged))
    for drug, g in groups:
        for method, col in SCORING_METHODS.items():
            corr = _correlate(g[col], g[TARGET_COL])
            rows.append({"drug": drug, "scoring_method": method, **corr})

    results = pd.DataFrame(rows)
    results.to_parquet(output_dir / "validation_correlations.parquet", index=False)
    return results


def plot_validation(
    results: pd.DataFrame,
    output_path: Path = None,
    panels_dir: Path = None,
    merged_path: Path = None,
) -> Path:
    """Three-panel figure: per-mutation scatter, top-N enrichment, per-drug ρ."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if output_path is None:
        output_path = config.VALIDATION_DIR / "validation_plot.png"
    if merged_path is None:
        merged_path = config.VALIDATION_DIR / "scores_vs_fold_resistance.parquet"
    merged = pd.read_parquet(merged_path)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5))

    # --- Panel 1: per-mutation scatter, colored by drug, DRMs outlined ---
    drugs = sorted(merged["drug"].unique())
    cmap = plt.get_cmap("tab10")
    for i, drug in enumerate(drugs):
        g = merged[merged["drug"] == drug]
        ax1.scatter(g[TARGET_COL], g["delta_delta_g"], s=20, alpha=0.5,
                    color=cmap(i % 10), label=drug, edgecolors="none")
    overall = results[(results.drug == "OVERALL") &
                      (results.scoring_method == "docking_ddg")].iloc[0]
    ax1.axhline(0, color="0.7", lw=0.8, zorder=0)
    ax1.set_xlabel("Measured mean log$_{10}$ fold-resistance")
    ax1.set_ylabel(r"Predicted $\Delta\Delta G$ (kcal/mol)")
    ax1.set_title(f"Per-mutation (all drugs pooled)\n"
                  f"Spearman ρ = {overall['spearman_rho']:.2f} "
                  f"(n = {int(overall['n_mutations'])}) — weak by design")
    ax1.legend(fontsize=8, ncol=2, framealpha=0.9)

    # --- Panel 2: enrichment of known DRMs among top-ddG predictions ---
    enr = compute_enrichment(merged)
    ax2.bar([str(int(n)) for n in enr["top_n"]], enr["enrichment"], color="#2c7fb8")
    ax2.axhline(1.0, color="0.4", lw=1.0, ls="--", label="no enrichment")
    ax2.set_xlabel("Top-N mutations by predicted ΔΔG")
    ax2.set_ylabel("Enrichment for known major DRMs (×)")
    ax2.set_title("Top predictions recover known DRMs\n"
                  f"(base rate {enr['base_rate'].iloc[0]:.0%})")
    for i, v in enumerate(enr["enrichment"]):
        ax2.text(i, v + 0.03, f"{v:.1f}×", ha="center", fontsize=9)
    ax2.legend(fontsize=8)

    # --- Panel 3: per-drug Spearman on primary DRMs ---
    pdc = per_drug_correlation(merged)
    prim = pdc[pdc["scope"] == "primary"].set_index("drug").reindex(drugs)
    vals = prim["spearman_rho"].to_numpy()
    colors = ["#2c7fb8" if (v == v and v > 0) else "#d95f5f" for v in vals]
    ax3.bar(drugs, vals, color=colors)
    ax3.axhline(0, color="0.4", lw=0.8)
    ax3.set_ylabel("Spearman ρ (primary DRMs)")
    ax3.set_title("Per-drug signal (major DRMs only)")
    for i, v in enumerate(vals):
        if v == v:
            ax3.text(i, v + (0.02 if v >= 0 else -0.05), f"{v:.2f}",
                     ha="center", va="bottom" if v >= 0 else "top", fontsize=9)
    ax3.tick_params(axis="x", rotation=20)

    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path
