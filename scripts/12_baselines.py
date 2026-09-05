"""12: Baseline head-to-head (Table 1) — is docking-ΔΔG worth anything?

Compares four per-mutation predictors of resistance on both targets:
  - docking ΔΔG            (our method)
  - mutation prevalence    (n_isolates; docking-free, data-only)
  - |Δ side-chain volume|  (docking-free, structure-free physicochemical)
  - Stanford HIVdb penalty (expert rule-based system; fetched via Sierra GraphQL)

on two tasks: DRM recovery (ROC-AUC vs is_primary) and magnitude (Spearman vs
measured fold-resistance), with bootstrap 95% CIs, plus a leave-one-drug-out
jackknife of the docking finding.

HIVdb penalties are fetched once and cached to data/hivdb_penalty_scores.json.
Key gotcha: Sierra must be queried with POSITION-DISJOINT mutation groups (one AA
per position) or the per-mutation partialScores collapse.

Usage::

    python scripts/12_baselines.py            # fetch (or use cache) + print Table 1
    python scripts/12_baselines.py --refetch  # force re-fetch HIVdb penalties
"""

import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

import config  # noqa: E402
from services.explanation import AA_PROPERTIES  # noqa: E402
from services.validation import build_merged  # noqa: E402

API = "https://hivdb.stanford.edu/graphql"
GENES = {"HIV1_PR": "PR", "HIV1_RT": "RT"}


def _disjoint_groups(muts):
    """Partition mutations so each group has at most one AA per position."""
    groups: list[dict] = []
    for m in sorted(muts):
        pos = m[1:-1]
        for g in groups:
            if pos not in g:
                g[pos] = m
                break
        else:
            groups.append({pos: m})
    return [list(g.values()) for g in groups]


def _fetch_group(gene, group):
    q = ("{ mutationsAnalysis(mutations:[%s]){ drugResistance{ drugScores{ "
         "drug{displayAbbr} partialScores{ mutations{text} score } } } } }"
         % ",".join(f'"{gene}:{m}"' for m in group))
    req = urllib.request.Request(API, data=json.dumps({"query": q}).encode(),
                                 headers={"Content-Type": "application/json",
                                          "User-Agent": "Mozilla/5.0 (research)"})
    r = json.load(urllib.request.urlopen(req, timeout=120))
    out: dict = {}
    for g in r["data"]["mutationsAnalysis"]["drugResistance"]:
        for ds in g["drugScores"]:
            drug = ds["drug"]["displayAbbr"].replace("/r", "")
            for p in ds["partialScores"]:
                if len(p["mutations"]) == 1:
                    out.setdefault(p["mutations"][0]["text"], {})[drug] = float(p["score"])
    return out


def hivdb_penalties(gene, muts, cache_path, refetch=False):
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    if gene not in cache or refetch:
        pen: dict = {}
        for grp in _disjoint_groups(muts):
            for mut, sc in _fetch_group(gene, grp).items():
                pen.setdefault(mut, {}).update(sc)
            time.sleep(0.25)
        cache[gene] = pen
        cache_path.write_text(json.dumps(cache))
    return cache[gene]


def _absdvol(mut):
    wt, mt = mut[0], mut[-1]
    if wt in AA_PROPERTIES and mt in AA_PROPERTIES:
        return abs(AA_PROPERTIES[mt]["volume"] - AA_PROPERTIES[wt]["volume"])
    return np.nan


def _roc(y, x):
    m = pd.notna(x)
    return roc_auc_score(y[m], np.asarray(x)[m]) if len(np.unique(y[m])) > 1 else np.nan


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refetch", action="store_true")
    ap.add_argument("--n-boot", type=int, default=2000)
    args = ap.parse_args()
    cache_path = config.DATA_DIR / "hivdb_penalty_scores.json"
    rng = np.random.default_rng(0)
    preds = {"docking ΔΔG": "delta_delta_g", "prevalence": "n_isolates",
             "|Δvolume|": "absdvol", "HIVdb penalty": "hivdb"}

    rows = []
    print("======== TABLE 1: baseline head-to-head (bootstrap 95% CIs) ========")
    print(f"{'target':9s} {'predictor':14s} {'ROC-AUC [95% CI]':>24s} {'ρ [95% CI]':>24s}")
    for tgt, gene in GENES.items():
        config.set_active_target(tgt)
        m = build_merged()
        pen = hivdb_penalties(gene, sorted(m["mutation"].unique()), cache_path, args.refetch)
        m["hivdb"] = m.apply(lambda r: pen.get(r.mutation, {}).get(r.drug, 0.0), axis=1)
        m["absdvol"] = m["mutation"].map(_absdvol)
        y = m["is_primary"].astype(int).to_numpy()
        fold = m["mean_log_fold_resistance"].to_numpy()
        n = len(m)
        lbl = "protease" if gene == "PR" else "RT"
        for name, col in preds.items():
            x = pd.to_numeric(m[col], errors="coerce").to_numpy()
            roc_pt = _roc(y, x)
            xm = pd.notna(x)
            rho_pt = stats.spearmanr(x[xm], fold[xm])[0]
            rocs, rhos = [], []
            for _ in range(args.n_boot):
                idx = rng.integers(0, n, n)
                yb, xb, fb = y[idx], x[idx], fold[idx]
                mb = pd.notna(xb)
                if len(np.unique(yb[mb])) > 1:
                    rocs.append(roc_auc_score(yb[mb], xb[mb]))
                if pd.Series(xb[mb]).nunique() > 1 and pd.Series(fb[mb]).nunique() > 1:
                    rhos.append(stats.spearmanr(xb[mb], fb[mb])[0])
            rows.append({"target": lbl, "predictor": name,
                         "roc_auc": round(roc_pt, 3),
                         "roc_ci_low": round(np.percentile(rocs, 2.5), 3),
                         "roc_ci_high": round(np.percentile(rocs, 97.5), 3),
                         "magnitude_rho": round(rho_pt, 3),
                         "rho_ci_low": round(np.percentile(rhos, 2.5), 3),
                         "rho_ci_high": round(np.percentile(rhos, 97.5), 3)})
            rc = f"{roc_pt:.2f} [{np.percentile(rocs,2.5):.2f},{np.percentile(rocs,97.5):.2f}]"
            rh = f"{rho_pt:+.2f} [{np.percentile(rhos,2.5):+.2f},{np.percentile(rhos,97.5):+.2f}]"
            print(f"{lbl:9s} {name:14s} {rc:>24s} {rh:>24s}")

        # leave-one-drug-out jackknife of the docking DRM-recovery finding
        loo = [(_roc((m[m.drug != d]["is_primary"].astype(int).to_numpy()),
                     m[m.drug != d]["delta_delta_g"].to_numpy()))
               for d in sorted(m["drug"].unique())]
        print(f"{'':9s} {'(LODO docking ROC)':14s} range [{min(loo):.3f},{max(loo):.3f}] "
              f"over {len(loo)} drugs")
        print("-" * 74)

    out = config.VALIDATION_DIR / "baseline_table.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"saved {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
