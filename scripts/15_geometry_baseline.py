"""Geometry-only predictive baselines + docking-dropout sensitivity (paper Appendix).

Two reviewer-driven analyses, both docking-free where possible.

1. **Geometry-only baselines.** The field-level ablation (scripts/13) says the
   *explanations* are carried by the mutation's geometry, not the docking energy.
   The obvious follow-up is whether geometry alone also carries the *ranking*.
   We score two predictors that need a receptor structure but no docking run:
     - ``-min_distance_to_ligand``  (closer to the ligand => more DRM-like)
     - ``contacts_ligand_directly`` (the binary contact flag)
   and compare them against docking ΔΔG and |Δ side-chain volume| on the same
   pairs, with bootstrap 95% CIs.

2. **Dropout sensitivity.** Pairs whose docking returned no valid pose are
   excluded from the benchmark. If those were enriched for DRMs, the reported
   AUC would be optimistic. We report the DRM rate among excluded vs kept pairs
   (Fisher exact), then bound the damage by re-scoring the pooled ROC-AUC with
   the excluded pairs imputed at the *worst case* (max ΔΔG => ranked top, all
   false positives) and at the median.

Outputs
-------
data/validation/geometry_baseline.csv     (both targets, all predictors)
data/validation/dropout_sensitivity.json  (both targets)

Usage: python scripts/15_geometry_baseline.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config  # noqa: E402
from services.explanation import build_structural_context  # noqa: E402

DDG, LABEL, TARGET = "delta_delta_g", "is_primary", "mean_log_fold_resistance"
SEED, N_BOOT = 0, 5000


def load_pairs(target_key: str) -> pd.DataFrame:
    """All attempted drug-mutation pairs for a target, ΔΔG NaN where docking failed."""
    if target_key == "protease":
        dock = pd.read_parquet(ROOT / "data/docking_results/benchmark_docking.parquet")
        panels = ROOT / "data/panels"
    else:
        dock = pd.read_parquet(ROOT / "data/rt/docking_results/benchmark_docking.parquet")
        panels = ROOT / "data/rt/panels"
    dock = dock[dock["mutation"] != "WT"]
    frames = []
    for p in sorted(panels.glob("*.parquet")):
        x = pd.read_parquet(p)
        x["drug"] = p.stem
        frames.append(x)
    pan = pd.concat(frames, ignore_index=True)
    cols = ["drug", "mutation", LABEL, TARGET]
    return dock.merge(pan[cols], on=["drug", "mutation"], how="inner")


def add_geometry(df: pd.DataFrame, target_name: str) -> pd.DataFrame:
    """Attach min-distance-to-ligand and the contact flag for every mutation."""
    t = config.get_target(target_name)
    cache: dict[str, tuple] = {}
    dists, contacts = [], []
    for mut in df["mutation"]:
        if mut not in cache:
            try:
                ctx = build_structural_context(mut, target=t)
                cache[mut] = (ctx["min_distance_to_ligand_angstrom"],
                              ctx["contacts_ligand_directly"])
            except Exception:
                cache[mut] = (np.nan, np.nan)
        d, c = cache[mut]
        dists.append(d)
        contacts.append(c)
    df = df.copy()
    df["min_dist"] = pd.to_numeric(pd.Series(dists), errors="coerce").to_numpy()
    df["contact_flag"] = pd.to_numeric(pd.Series(contacts), errors="coerce").to_numpy()
    df["absdvol"] = [_absdvol(m) for m in df["mutation"]]
    return df


def _absdvol(mut: str) -> float:
    from services.explanation import AA_PROPERTIES
    try:
        wt, mt = mut[0], mut[-1]
        return abs(AA_PROPERTIES[mt]["volume"] - AA_PROPERTIES[wt]["volume"])
    except Exception:
        return np.nan


def _roc_ci(y, s, rng):
    """ROC-AUC with a bootstrap 95% CI, NaN-safe."""
    m = ~(pd.isna(y) | pd.isna(s))
    y, s = np.asarray(y, float)[m], np.asarray(s, float)[m]
    if len(np.unique(y)) < 2:
        return np.nan, np.nan, np.nan, 0
    auc = roc_auc_score(y, s)
    boots = []
    n = len(y)
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        if len(np.unique(y[i])) < 2:
            continue
        boots.append(roc_auc_score(y[i], s[i]))
    lo, hi = np.percentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan)
    return auc, lo, hi, int(m.sum())


def geometry_table(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Predictor head-to-head on the scored pairs (docking ΔΔG available)."""
    rng = np.random.default_rng(SEED)
    d = df.dropna(subset=[DDG])          # same rows the benchmark uses
    preds = {
        "docking ΔΔG": d[DDG].to_numpy(),
        "-min distance to ligand": -d["min_dist"].to_numpy(),
        "contacts ligand (binary)": d["contact_flag"].to_numpy(),
        "|Δvolume|": d["absdvol"].to_numpy(),
    }
    rows = []
    for name, s in preds.items():
        auc, lo, hi, n = _roc_ci(d[LABEL].to_numpy(), s, rng)
        mm = ~(pd.isna(s) | pd.isna(d[TARGET].to_numpy()))
        rho, p = spearmanr(np.asarray(s, float)[mm],
                           d[TARGET].to_numpy()[mm]) if mm.sum() > 2 else (np.nan, np.nan)
        rows.append({"target": target, "predictor": name, "n": n,
                     "roc_auc": auc, "ci_low": lo, "ci_high": hi,
                     "magnitude_rho": float(rho), "rho_p": float(p)})
    return pd.DataFrame(rows)


def precision_at_10(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Per-drug precision@10 for docking ΔΔG vs. the distance baseline.

    ROC-AUC scores the whole ranking; the tool is read at the top of the list.
    The two can disagree: a predictor that orders the full list well may still
    fill its top-10 with near-but-negative candidates.
    """
    from scipy.stats import hypergeom
    d = df.dropna(subset=[DDG])
    rows = []
    for drug, g in d.groupby("drug", sort=True):
        n, n_drm = len(g), int(g[LABEL].sum())
        y = g[LABEL].to_numpy(dtype=int)
        for name, score in (("docking ΔΔG", g[DDG].to_numpy()),
                            ("-min distance", -g["min_dist"].to_numpy())):
            order = np.argsort(score)[::-1][:10]
            hits = int(y[order].sum())
            rows.append({
                "target": target, "drug": drug, "predictor": name,
                "hits_at_10": hits, "precision_at_10": hits / 10,
                "base_rate": n_drm / n,
                "hypergeom_p": float(hypergeom.sf(hits - 1, n, n_drm, 10)),
            })
    return pd.DataFrame(rows)


def dropout_sensitivity(df: pd.DataFrame, target: str) -> dict:
    """DRM enrichment of the excluded pairs, and worst-case AUC bounds."""
    excl, kept = df[df[DDG].isna()], df[df[DDG].notna()]
    out = {"target": target, "n_attempted": len(df), "n_excluded": len(excl),
           "excl_frac": len(excl) / len(df) if len(df) else 0.0,
           "drm_rate_excluded": float(excl[LABEL].mean()) if len(excl) else None,
           "drm_rate_kept": float(kept[LABEL].mean()),
           "n_drms_excluded": int(excl[LABEL].sum()) if len(excl) else 0}
    if len(excl):
        tab = [[int(excl[LABEL].sum()), int((~excl[LABEL].astype(bool)).sum())],
               [int(kept[LABEL].sum()), int((~kept[LABEL].astype(bool)).sum())]]
        odds, p = fisher_exact(tab)
        out["fisher_p"] = float(p)
        out["odds_ratio"] = float(odds)
    rng = np.random.default_rng(SEED)
    out["roc_observed"] = _roc_ci(kept[LABEL], kept[DDG], rng)[0]
    # worst case: every excluded pair scores above all kept pairs (top-ranked FPs)
    for label, fill in (("roc_worst_case", df[DDG].max() + 1.0),
                        ("roc_median_impute", df[DDG].median())):
        imp = df.copy()
        imp[DDG] = imp[DDG].fillna(fill)
        out[label] = _roc_ci(imp[LABEL], imp[DDG], rng)[0]
    return out


def main() -> None:
    tables, sens, patk = [], [], []
    for key, tname in (("protease", "HIV1_PR"), ("RT", "HIV1_RT")):
        df = load_pairs(key)
        df = add_geometry(df, tname)
        tables.append(geometry_table(df, key))
        patk.append(precision_at_10(df, key))
        sens.append(dropout_sensitivity(df, key))
        print(f"{key}: {len(df)} attempted, {df[DDG].isna().sum()} no-pose")

    tab = pd.concat(tables, ignore_index=True)
    tab.to_csv(ROOT / "data/validation/geometry_baseline.csv", index=False)
    pk = pd.concat(patk, ignore_index=True)
    pk.to_csv(ROOT / "data/validation/geometry_precision_at_10.csv", index=False)
    print("\n=== precision@10: docking vs distance ===")
    print(pk.pivot_table(index=["target", "drug"], columns="predictor",
                         values="hits_at_10").to_string())
    (ROOT / "data/validation/dropout_sensitivity.json").write_text(
        json.dumps(sens, indent=2))

    with pd.option_context("display.width", 200):
        print("\n=== geometry-only baselines ===")
        print(tab.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("\n=== dropout sensitivity ===")
    print(json.dumps(sens, indent=2))


if __name__ == "__main__":
    main()
