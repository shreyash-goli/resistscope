"""Generate the paper's data figures (F2 predictive, F3 interpretability) from the
committed results. Colorblind-safe Okabe-Ito palette, publication style."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, ".")
import config
from services.explanation import AA_PROPERTIES
from services.validation import build_merged

FIG = Path("paper/figures"); FIG.mkdir(parents=True, exist_ok=True)
# Okabe-Ito colorblind-safe palette (fixed roles, never cycled)
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "purple": "#CC79A7",
      "vermillion": "#D55E00", "gray": "#9a9a9a"}
plt.rcParams.update({
    "savefig.dpi": 300, "font.size": 9, "font.family": "sans-serif",
    "axes.titlesize": 9.5, "axes.titleweight": "bold", "axes.labelsize": 8.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True,
    "axes.grid": True, "grid.color": "#ebebeb", "grid.linewidth": 0.6,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "legend.frameon": False,
})
pen = json.loads((config.DATA_DIR / "hivdb_penalty_scores.json").read_text())


def dvol(m):
    a, b = m[0], m[-1]
    if a in AA_PROPERTIES and b in AA_PROPERTIES:
        return abs(AA_PROPERTIES[b]["volume"] - AA_PROPERTIES[a]["volume"])
    return np.nan


def merged(tgt, gene):
    config.set_active_target(tgt)
    m = build_merged()
    m["hivdb"] = m.apply(lambda r: pen[gene].get(r.mutation, {}).get(r.drug, 0.0), axis=1)
    m["absdvol"] = m["mutation"].map(dvol)
    return m


def roc_ci(y, x, nb=1000):
    ok = pd.notna(x); y, x = y[ok], np.asarray(x)[ok]; n = len(y)
    rng = np.random.default_rng(0)
    boot = [roc_auc_score(y[i], x[i]) for i in (rng.integers(0, n, n) for _ in range(nb))
            if len(np.unique(y[i])) > 1]
    return roc_auc_score(y, x), np.percentile(boot, 2.5), np.percentile(boot, 97.5)


def fig2():
    preds = [("docking ΔΔG", "delta_delta_g"), ("prevalence", "n_isolates"),
             ("|Δvolume|", "absdvol"), ("HIVdb penalty", "hivdb")]
    D = {}
    for tgt, gene in [("HIV1_PR", "PR"), ("HIV1_RT", "RT")]:
        m = merged(tgt, gene); y = m["is_primary"].astype(int).to_numpy()
        D[gene] = {n: roc_ci(y, pd.to_numeric(m[c], errors="coerce").to_numpy()) for n, c in preds}

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9, 3.5))
    names = [p[0] for p in preds]; yp = np.arange(len(names))[::-1]; h = 0.38
    for gene, col, off in [("PR", OI["orange"], h / 2), ("RT", OI["blue"], -h / 2)]:
        pt = [D[gene][n][0] for n in names]
        hi = [D[gene][n][2] for n in names]
        err = [[D[gene][n][0] - D[gene][n][1] for n in names], [D[gene][n][2] - D[gene][n][0] for n in names]]
        axA.barh(yp + off, pt, height=h, color=col, label={"PR": "protease", "RT": "RT"}[gene],
                 xerr=err, error_kw=dict(lw=0.8, ecolor="#555", capsize=2))
        for y_, v, ci in zip(yp + off, pt, hi):
            axA.text(ci + 0.012, y_, f"{v:.2f}", va="center", ha="left", fontsize=6.5, color="#333")
    axA.axvline(0.5, color="#888", ls="--", lw=1)
    axA.text(0.5, len(names) - 0.42, "chance", fontsize=6.5, color="#888", ha="center")
    axA.set_yticks(yp); axA.set_yticklabels(names); axA.set_xlim(0.4, 1.08)
    axA.set_xlabel("DRM-recovery ROC-AUC (95% CI)")
    axA.set_title("a   Docking is chance on protease,\nabove-chance on RT; HIVdb wins both")
    axA.legend(loc="upper right")

    order = ["all isolates (confounded)", "≤3 mutations", "≤2 mutations", "single-mutation only"]
    for gene, path, col in [("PR", "data/validation/benchmark_deconfounding.parquet", OI["orange"]),
                            ("RT", "data/rt/validation/benchmark_deconfounding.parquet", OI["blue"])]:
        d = pd.read_parquet(path).set_index("subset").reindex(order)
        axB.plot(range(len(order)), d["spearman_rho"], "-o", color=col, lw=2, ms=6,
                 label={"PR": "protease", "RT": "RT"}[gene])
    axB.axhline(0, color="#888", lw=1)
    axB.set_xticks(range(len(order))); axB.set_xticklabels(["all", "≤3 mut", "≤2 mut", "single"])
    axB.set_ylabel("Spearman ρ  (ΔΔG vs. measured fold-resistance)")
    axB.set_xlabel("isolate subset  (less confounded →)")
    axB.set_title("b   De-confounding rescues RT (ρ→0.4),\nnot protease (ρ≤0)")
    axB.legend(loc="center left")
    fig.tight_layout(); fig.savefig(FIG / "fig2_benchmark.png", bbox_inches="tight"); plt.close(fig)
    print("wrote fig2_benchmark.png")


def fig3():
    abl = {"PR": pd.read_parquet("data/validation/faithfulness_ablation.parquet"),
           "RT": pd.read_parquet("data/rt/validation/faithfulness_ablation.parquet")}
    conds = ["full", "minimal", "corrupted"]
    ccol = {"full": OI["blue"], "minimal": OI["gray"], "corrupted": OI["vermillion"]}
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(12.5, 3.6))

    genes = ["PR", "RT"]; x = np.arange(len(genes)); w = 0.26
    for j, cond in enumerate(conds):
        means = [abl[g][abl[g].condition == cond].score.mean() for g in genes]
        axA.bar(x + (j - 1) * w, means, width=w, color=ccol[cond], label=cond)
        for xi, mv in zip(x + (j - 1) * w, means):
            axA.text(xi, mv + 0.03, f"{mv:.2f}", ha="center", fontsize=6.5, color="#333")
    axA.set_xticks(x); axA.set_xticklabels(["protease", "RT"]); axA.set_ylim(0, 2.1)
    axA.set_ylabel("mean faithfulness (0–2)")
    axA.set_title("a   Context helps, corruption hurts\n(both targets)"); axA.legend(loc="upper right")

    pooled = []
    for g in genes:
        wpv = abl[g].pivot_table(index=["drug", "mutation"], columns="condition", values="score")
        wpv["gain"] = wpv["full"] - wpv["minimal"]; pooled.append(wpv)
    P = pd.concat(pooled)
    rows = [(m, P[P.minimal == m].gain.mean(), (P[P.minimal == m].gain > 0).mean() * 100) for m in [0, 1, 2]]
    gains = [r[1] for r in rows]
    axB.bar([0, 1, 2], gains, width=0.62, color=[OI["green"] if v > 0 else OI["vermillion"] for v in gains])
    for xi, r in zip([0, 1, 2], rows):
        va = "bottom" if r[1] >= 0 else "top"
        axB.text(xi, r[1] + (0.04 if r[1] >= 0 else -0.03), f"{r[1]:+.2f}\n{r[2]:.0f}% fixed",
                 ha="center", va=va, fontsize=6.5, color="#333")
    axB.set_ylim(-0.42, 1.42)
    axB.axhline(0, color="#888", lw=1); axB.set_xticks([0, 1, 2])
    axB.set_xticklabels(["0\nprior wrong", "1\nprior vague", "2\nprior correct"])
    axB.set_ylabel("faithfulness gain from context (full − minimal)")
    axB.set_xlabel("model's prior strength (minimal-only score)")
    axB.set_title("b   Context helps most where\nthe prior is weakest")

    fa = pd.read_parquet("data/validation/field_ablation.parquet")
    w3 = fa.pivot_table(index=["drug", "mutation"], columns="condition", values="score")
    fields = [("ligand distance", "drop_distance"), ("subpocket", "drop_subpocket"),
              ("docking ΔΔG", "drop_ddg"), ("chemistry", "drop_chem")]
    dl = [(n, (w3["full"] - w3[c]).mean()) for n, c in fields]
    yp = np.arange(len(fields))[::-1]
    axC.barh(yp, [d for _, d in dl], height=0.6,
             color=[OI["purple"] if d > 0.05 else OI["gray"] for _, d in dl])
    for y_, (_, d) in zip(yp, dl):
        axC.text(d + 0.004, y_, f"{d:+.2f}", va="center", fontsize=6.5, color="#333")
    axC.axvline(0, color="#888", lw=1); axC.set_yticks(yp); axC.set_yticklabels([n for n, _ in dl])
    axC.set_xlabel("faithfulness lost when field dropped (full − drop)")
    axC.set_title("c   The lift is geometric\n(distance/subpocket), not energy")
    fig.tight_layout(); fig.savefig(FIG / "fig3_interpretability.png", bbox_inches="tight"); plt.close(fig)
    print("wrote fig3_interpretability.png")


fig2(); fig3()
