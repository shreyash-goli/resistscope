"""Per-drug precision@k for docking-ΔΔG DRM triage (paper Appendix).

Why this exists
---------------
``scripts/08_benchmark.py`` reports *pooled* top-N enrichment and ROC-AUC. Neither
matches how the tool is actually used: a chemist holds **one** compound, looks at
the top few flagged mutations for **that** compound, and asks "are these real?"

This script computes the decision-relevant metric — per-drug **precision@k**: rank
a single drug's own mutation panel by predicted ΔΔG (descending, matching
``services.benchmark.enrichment_with_significance``), take the top ``k``, and
report the fraction that are known major DRMs, against that drug's own base rate.

Significance is an exact one-sided hypergeometric test (probability that drawing
``k`` at random from the drug's panel yields at least as many DRMs). At k=5 a
bootstrap CI is uninformative, so we report the raw count and the exact p.

Also emits the per-drug ROC-AUC / magnitude-Spearman table used in the paper.

Outputs
-------
data/validation/precision_at_k.csv          (protease)
data/rt/validation/precision_at_k.csv       (RT)
data/validation/per_drug_table.csv          (both targets, paper table)

Usage: python scripts/14_precision_at_k.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]

DDG_COL = "delta_delta_g"
LABEL_COL = "is_primary"
TARGET_COL = "mean_log_fold_resistance"
KS = (5, 10, 20)


def load_protease() -> pd.DataFrame:
    """Protease pairs already carry ΔΔG, the DRM label and measured resistance."""
    df = pd.read_parquet(ROOT / "data/validation/scores_vs_fold_resistance.parquet")
    return df[df["mutation"] != "WT"].copy()


def load_rt() -> pd.DataFrame:
    """RT: join docked ΔΔG onto the per-drug panels that carry the DRM labels."""
    dock = pd.read_parquet(ROOT / "data/rt/docking_results/benchmark_docking.parquet")
    dock = dock[dock["mutation"] != "WT"]
    frames = []
    for path in sorted((ROOT / "data/rt/panels").glob("*.parquet")):
        panel = pd.read_parquet(path)
        panel["drug"] = path.stem
        frames.append(panel)
    panels = pd.concat(frames, ignore_index=True)
    merged = dock.merge(
        panels[["drug", "mutation", LABEL_COL, TARGET_COL, "n_isolates"]],
        on=["drug", "mutation"],
        how="inner",
    )
    return merged


def precision_at_k(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Per-drug precision@k with an exact hypergeometric one-sided p-value."""
    rows = []
    for drug, g in df.groupby("drug", sort=True):
        g = g.dropna(subset=[DDG_COL, LABEL_COL])
        n, n_drm = len(g), int(g[LABEL_COL].sum())
        if n == 0 or n_drm == 0:
            continue
        base = n_drm / n
        # high ΔΔG first — same convention as services.benchmark
        order = g[DDG_COL].to_numpy().argsort()[::-1]
        labels = g[LABEL_COL].to_numpy(dtype=int)
        for k in KS:
            kk = min(k, n)
            hits = int(labels[order[:kk]].sum())
            prec = hits / kk
            # P(X >= hits) drawing kk from n with n_drm successes
            p = float(hypergeom.sf(hits - 1, n, n_drm, kk))
            rows.append({
                "target": target, "drug": drug, "k": kk,
                "hits": hits, "precision": prec, "base_rate": base,
                "lift": prec / base if base > 0 else np.nan,
                "hypergeom_p": p, "n_panel": n, "n_drms": n_drm,
            })
    return pd.DataFrame(rows)


def pooled_precision_at_k(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Same metric on the pooled cross-drug ranking, for comparison."""
    g = df.dropna(subset=[DDG_COL, LABEL_COL])
    n, n_drm = len(g), int(g[LABEL_COL].sum())
    base = n_drm / n
    order = g[DDG_COL].to_numpy().argsort()[::-1]
    labels = g[LABEL_COL].to_numpy(dtype=int)
    rows = []
    for k in KS:
        hits = int(labels[order[:k]].sum())
        rows.append({
            "target": target, "drug": "POOLED", "k": k, "hits": hits,
            "precision": hits / k, "base_rate": base,
            "lift": (hits / k) / base, "n_panel": n, "n_drms": n_drm,
            "hypergeom_p": float(hypergeom.sf(hits - 1, n, n_drm, k)),
        })
    return pd.DataFrame(rows)


def per_drug_table(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Per-drug DRM-recovery ROC-AUC and magnitude Spearman ρ vs measured data."""
    rows = []
    for drug, g in df.groupby("drug", sort=True):
        g = g.dropna(subset=[DDG_COL, LABEL_COL])
        y, s = g[LABEL_COL].to_numpy(dtype=int), g[DDG_COL].to_numpy()
        roc = float(roc_auc_score(y, s)) if len(np.unique(y)) > 1 else np.nan
        gm = g.dropna(subset=[TARGET_COL])
        if len(gm) > 2:
            rho, pval = spearmanr(gm[DDG_COL], gm[TARGET_COL])
        else:
            rho, pval = np.nan, np.nan
        rows.append({
            "target": target, "drug": drug, "n": len(g),
            "n_drms": int(y.sum()), "base_rate": y.mean(),
            "roc_auc": roc, "spearman_rho": float(rho),
            "spearman_p": float(pval),
        })
    return pd.DataFrame(rows)


def latex_precision_table(pk: pd.DataFrame) -> str:
    """Compact LaTeX body: one row per drug, precision@5/10/20 as hits/k."""
    lines = []
    for target in ("protease", "RT"):
        sub = pk[(pk.target == target) & (pk.drug != "POOLED")]
        if sub.empty:
            continue
        lines.append(rf"\multicolumn{{5}}{{l}}{{\emph{{{target}}}}} \\")
        for drug in sorted(sub.drug.unique()):
            d = sub[sub.drug == drug].set_index("k")
            cells = []
            for k in KS:
                if k in d.index:
                    r = d.loc[k]
                    star = "*" if r.hypergeom_p < 0.05 else ""
                    cells.append(rf"{int(r.hits)}/{k}{star}")
                else:
                    cells.append("--")
            base = d.iloc[0]["base_rate"]
            lines.append(f"\\quad {drug} & " + " & ".join(cells) + f" & {base:.2f} \\\\")
    return "\n".join(lines)


def main() -> None:
    pr, rt = load_protease(), load_rt()
    print(f"protease pairs: {len(pr)}   RT pairs: {len(rt)}")

    pk = pd.concat([
        precision_at_k(pr, "protease"), pooled_precision_at_k(pr, "protease"),
        precision_at_k(rt, "RT"), pooled_precision_at_k(rt, "RT"),
    ], ignore_index=True)

    tbl = pd.concat([per_drug_table(pr, "protease"), per_drug_table(rt, "RT")],
                    ignore_index=True)

    (ROOT / "data/validation").mkdir(parents=True, exist_ok=True)
    (ROOT / "data/rt/validation").mkdir(parents=True, exist_ok=True)
    pk[pk.target == "protease"].to_csv(ROOT / "data/validation/precision_at_k.csv", index=False)
    pk[pk.target == "RT"].to_csv(ROOT / "data/rt/validation/precision_at_k.csv", index=False)
    tbl.to_csv(ROOT / "data/validation/per_drug_table.csv", index=False)

    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print("\n=== per-drug precision@k ===")
        print(pk.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print("\n=== per-drug ROC / rho ===")
        print(tbl.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    print("\n=== LaTeX (precision@k body) ===")
    print(latex_precision_table(pk))

    summary = {
        "ks": list(KS),
        "n_pairs": {"protease": len(pr), "RT": len(rt)},
        "median_lift_at_10": {
            t: float(pk[(pk.target == t) & (pk.k == 10) & (pk.drug != "POOLED")].lift.median())
            for t in ("protease", "RT")
        },
    }
    (ROOT / "data/validation/precision_at_k_summary.json").write_text(
        json.dumps(summary, indent=2))
    print("\nwrote precision_at_k.csv (both targets), per_drug_table.csv, summary json")


if __name__ == "__main__":
    main()
