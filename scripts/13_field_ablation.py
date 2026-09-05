"""13: Field-level attribution ablation — which structural-context fields carry
the faithfulness lift?

Drops one field group at a time from the full context (same model + same judge)
and re-scores faithfulness. Larger drop-vs-full gap = that field matters more.
Since the docking ΔΔG is a chance-level predictor, the hypothesis is that the
geometry (distance/subpocket) carries the lift, not the energy.

The 'full' condition reuses the scripts/11 ablation cache (identical generation
path), so only the 4 drop conditions are generated.

Usage::

    python scripts/13_field_ablation.py                # protease
    python scripts/13_field_ablation.py --target rt
"""

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import config  # noqa: E402
from services.ablation import generate_for_condition  # noqa: E402
from services.explanation import _faithfulness_pairs, evaluate_faithfulness  # noqa: E402

# (fields to strip from the context, whether to also null the ΔΔG prompt line)
FIELD_CONDS = {
    "full": ([], False),
    "drop_ddg": (["delta_delta_g", "delta_g"], True),
    "drop_distance": (["distance_from_ligand_centroid_angstrom",
                       "min_distance_to_ligand_angstrom", "contacts_ligand_directly"], False),
    "drop_subpocket": (["region"], False),
    "drop_chem": (["volume_change", "volume_change_A3", "charge_change",
                   "hydrophobicity_change"], False),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", default="HIV1_PR")
    ap.add_argument("--sleep", type=float, default=0.25)
    args = ap.parse_args()
    t = config.set_active_target(args.target)

    gt = json.loads(t.ground_truth_path.read_text())
    pairs = _faithfulness_pairs(t.explanations_dir, gt, target=t)
    print(f"Field ablation: {len(pairs)} pairs x {len(FIELD_CONDS)} conditions ({t.name})")

    def load(drug, mut):
        rec = json.loads((t.explanations_dir / f"{drug}_{mut}.json").read_text())
        return rec.get("structural_context") or {}, rec.get("delta_delta_g")

    abl_root = config.DATA_DIR / ("ablation" if t.name == "HIV1_PR" else f"ablation_{t.name}")
    base = config.DATA_DIR / ("field_ablation" if t.name == "HIV1_PR"
                              else f"field_ablation_{t.name}")

    frames = []
    for cond, (drop, null_ddg) in FIELD_CONDS.items():
        cond_dir = (abl_root / "full") if cond == "full" else (base / cond)
        n_new = 0
        for p in pairs:
            ctx, ddg = load(p["drug"], p["mutation"])
            if cond != "full":
                ctx = {k: v for k, v in ctx.items() if k not in drop}
                if null_ddg:
                    ddg = None
                existed = (cond_dir / f"{p['drug']}_{p['mutation']}.json").exists()
                generate_for_condition(p["drug"], p["mutation"], ddg, ctx, "full", cond_dir, target=t)
                if not existed:
                    n_new += 1
                    time.sleep(args.sleep)
        res = evaluate_faithfulness(explanations_dir=cond_dir,
                                    ground_truth_path=t.ground_truth_path, target=t)
        res["condition"] = cond
        frames.append(res)
        print(f"  [{cond:14s}] mean {res['score'].mean():.2f}  "
              f"correct {100*(res['score']==2).mean():.0f}%  (n={len(res)}, {n_new} new)")

    allr = pd.concat(frames, ignore_index=True)
    out = config.VALIDATION_DIR / "field_ablation.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    allr.to_parquet(out, index=False)

    w = allr.pivot_table(index=["drug", "mutation"], columns="condition", values="score")
    print("\n=== drop vs full (paired; larger Δ = field carries more of the lift) ===")
    for cond in FIELD_CONDS:
        if cond == "full" or cond not in w.columns:
            continue
        ww = w.dropna(subset=["full", cond])
        d = ww["full"] - ww[cond]
        pv = stats.wilcoxon(ww["full"], ww[cond]).pvalue if d.abs().sum() else float("nan")
        print(f"  full - {cond:14s}: Δ={d.mean():+.2f}  p={pv:.3f}  (n={len(ww)})")
    print(f"saved {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
