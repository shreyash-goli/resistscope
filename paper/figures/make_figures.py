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


def fig_toplist():
    """Per-drug precision@10, docking ΔΔG vs. the docking-free distance baseline.

    The complement to Table 1: the table is the global ranking, this is the top of
    the list, where the ordering reverses between targets. Reads only committed
    CSVs (scripts/15), so it needs no docking stack.
    """
    pk = pd.read_csv("data/validation/geometry_precision_at_10.csv")
    pk = pk.pivot(index=["target", "drug"], columns="predictor",
                  values=["precision_at_10", "hypergeom_p", "base_rate"])
    rows = []  # top-to-bottom: protease block, gap, RT block
    for tgt, label in [("protease", "protease  (PIs)"), ("RT", "RT  (NNRTIs)")]:
        block = pk.loc[tgt].sort_index()
        rows.append((label, None, None, None, None, None))
        for drug, r in block.iterrows():
            rows.append((drug, r[("precision_at_10", "docking ΔΔG")],
                         r[("precision_at_10", "-min distance")],
                         r[("hypergeom_p", "docking ΔΔG")],
                         r[("hypergeom_p", "-min distance")],
                         r[("base_rate", "docking ΔΔG")]))

    # Sized to print near 1:1 at \linewidth (5.5in) so the fonts survive.
    fig, ax = plt.subplots(figsize=(5.9, 3.6))
    yp = np.arange(len(rows))[::-1]
    for y_, (name, dock, dist, pd_, pdist, base) in zip(yp, rows):
        if dock is None:  # group header
            ax.text(-0.055, y_, name, ha="right", va="center", fontsize=8.5,
                    fontweight="bold", color="#333")
            continue
        ax.plot([dock, dist], [y_, y_], color="#cfcfcf", lw=2.2, zorder=1,
                solid_capstyle="round")
        ax.scatter([base], [y_], marker="|", s=90, color="#9a9a9a", lw=1.4, zorder=2)
        for v, p, col in [(dock, pd_, OI["blue"]), (dist, pdist, OI["green"])]:
            ax.scatter([v], [y_], s=62, color=col, zorder=3)
            if p < 0.05:
                ax.text(v, y_ + 0.26, "*", ha="center", va="center", fontsize=8, color=col)
        ax.text(-0.055, y_, name, ha="right", va="center", fontsize=8)

    ax.scatter([], [], s=62, color=OI["blue"], label="docking ΔΔG")
    ax.scatter([], [], s=62, color=OI["green"], label="−min distance  (no docking)")
    ax.scatter([], [], marker="|", s=90, color="#9a9a9a", lw=1.4, label="drug's base rate")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.30), ncol=3)

    ax.set_yticks([]); ax.set_ylim(-0.9, len(rows) - 0.1)
    ax.set_xlim(-0.06, 1.06); ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_xlabel("precision@10  (genuine clinical DRMs in the drug's top 10)")
    ax.set_title("The ordering reverses at the top of the list:\n"
                 "distance wins on protease, the energy wins in the NNRTI pocket")
    ax.spines["left"].set_visible(False)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig_toplist.png", bbox_inches="tight"); plt.close(fig)
    print("wrote fig_toplist.png")


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
    # Body figure at the printed size (\linewidth = 5.5in), so nothing is
    # downscaled and the default 9pt style holds. The former panel (a), mean
    # faithfulness by condition, is dropped: 4.3 states those deltas in text.
    abl = {"PR": pd.read_parquet("data/validation/faithfulness_ablation.parquet"),
           "RT": pd.read_parquet("data/rt/validation/faithfulness_ablation.parquet")}
    fig, (axB, axC) = plt.subplots(1, 2, figsize=(5.6, 2.5))

    genes = ["PR", "RT"]
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
        axB.text(xi, r[1] + (0.05 if r[1] >= 0 else -0.04), f"{r[1]:+.2f}\n{r[2]:.0f}% fixed",
                 ha="center", va=va, fontsize=6, color="#333")
    axB.set_ylim(-0.5, 1.62)
    axB.axhline(0, color="#888", lw=1); axB.set_xticks([0, 1, 2])
    axB.set_xticklabels(["0\nwrong", "1\nvague", "2\ncorrect"])
    axB.set_ylabel("gain (full − minimal)")
    axB.set_xlabel("model's prior strength")
    axB.set_title("a   Context helps most\nwhere the prior is weakest")

    fa = pd.read_parquet("data/validation/field_ablation.parquet")
    w3 = fa.pivot_table(index=["drug", "mutation"], columns="condition", values="score")
    fields = [("ligand distance", "drop_distance"), ("subpocket", "drop_subpocket"),
              ("docking ΔΔG", "drop_ddg"), ("chemistry", "drop_chem")]
    dl = [(n, (w3["full"] - w3[c]).mean()) for n, c in fields]
    yp = np.arange(len(fields))[::-1]
    axC.barh(yp, [d for _, d in dl], height=0.6,
             color=[OI["purple"] if d > 0.05 else OI["gray"] for _, d in dl])
    for y_, (_, d) in zip(yp, dl):
        axC.text(d + (0.004 if d >= 0 else -0.004), y_, f"{d:+.2f}", va="center",
                 ha="left" if d >= 0 else "right", fontsize=6.5, color="#333")
    axC.axvline(0, color="#888", lw=1); axC.set_yticks(yp)
    axC.set_yticklabels([n for n, _ in dl]); axC.set_xlim(-0.095, 0.305)
    axC.set_xlabel("faithfulness lost (full − drop)")
    axC.set_title("b   The lift is geometric,\nnot energetic")
    fig.tight_layout(w_pad=1.6)
    fig.savefig(FIG / "fig3_interpretability.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig3_interpretability.png")


fig2(); fig3()
