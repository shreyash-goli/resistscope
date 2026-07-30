"""11: Faithfulness attribution ablation — pipeline vs. parametric knowledge.

Regenerates the ground-truth explanation set under three context conditions
(full / minimal / corrupted; see ``services/ablation.py``) with the same model
and system prompt, then scores all three with the same Claude judge against the
same curated mechanisms. Answers: is the 72% explanation faithfulness driven by
our structural pipeline, or by Claude's prior knowledge of these (famous) DRMs?

Real per-pair context + ΔΔG are read from the existing cached explanation
records (``data/explanations/{drug}_{mutation}.json``), so no receptor structure
or docking stack is needed. Generations are cached per condition, so re-runs are
free and interruptible.

Usage::

    python scripts/11_faithfulness_ablation.py --limit 3    # cheap smoke test
    python scripts/11_faithfulness_ablation.py              # full n=46 x 3
    python scripts/11_faithfulness_ablation.py --target rt
"""

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import config  # noqa: E402
from services.ablation import CONDITIONS, derange_contexts, generate_for_condition  # noqa: E402
from services.explanation import _faithfulness_pairs, evaluate_faithfulness  # noqa: E402

pd.set_option("display.width", 140)


def _load_real_context(explanations_dir: Path, drug: str, mutation: str):
    """(structural_context, ΔΔG) from a cached production explanation record."""
    rec = json.loads((explanations_dir / f"{drug}_{mutation}.json").read_text())
    return rec.get("structural_context") or {}, rec.get("delta_delta_g")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", default="HIV1_PR",
                    help="HIV1_PR / pr (default) or HIV1_RT / rt.")
    ap.add_argument("--limit", type=int, default=None,
                    help="Cap pairs (smoke test).")
    ap.add_argument("--sleep", type=float, default=0.3,
                    help="Seconds between API calls.")
    ap.add_argument("--seed", type=int, default=0, help="Corruption derangement seed.")
    args = ap.parse_args()

    t = config.set_active_target(args.target)
    print(f"Target: {t.name} ({t.label})")

    if not t.ground_truth_path.exists():
        print(f"No ground truth at {t.ground_truth_path}. Run scripts/09 first.")
        return 1
    gt = json.loads(t.ground_truth_path.read_text())
    pairs = _faithfulness_pairs(t.explanations_dir, gt, target=t)
    if not pairs:
        print("No cached explanations overlap the ground truth. Run scripts/06 first.")
        return 1
    if args.limit:
        pairs = pairs[: args.limit]
    print(f"Ablating {len(pairs)} (drug, mutation) pairs x {len(CONDITIONS)} conditions "
          f"(model={config.CLAUDE_MODEL})")

    # Real context + ΔΔG per pair (from cached records); corrupted = derangement.
    reals = [_load_real_context(t.explanations_dir, p["drug"], p["mutation"]) for p in pairs]
    corrupted = derange_contexts(reals, seed=args.seed)
    ctx_by_cond = {"full": reals, "minimal": reals, "corrupted": corrupted}

    base = config.DATA_DIR / ("ablation" if t.name == "HIV1_PR" else f"ablation_{t.name}")
    for cond in CONDITIONS:
        cond_dir = base / cond
        n_new = 0
        for i, p in enumerate(pairs):
            ctx, ddg = ctx_by_cond[cond][i]
            existed = (cond_dir / f"{p['drug']}_{p['mutation']}.json").exists()
            try:
                generate_for_condition(p["drug"], p["mutation"], ddg, ctx, cond,
                                       cond_dir, target=t)
            except Exception as exc:  # noqa: BLE001
                msg = str(exc).lower()
                if "authentication" in msg or "api_key" in msg:
                    print("ERROR: no Anthropic credentials (set ANTHROPIC_API_KEY).")
                    return 1
                print(f"  FAILED [{cond}] {p['drug']} {p['mutation']}: {exc}")
                continue
            if not existed:
                n_new += 1
                time.sleep(args.sleep)
        print(f"  [{cond:9s}] {n_new} new generated -> {cond_dir}")

    # Judge every condition with the same judge against the same ground truth.
    frames = []
    for cond in CONDITIONS:
        res = evaluate_faithfulness(explanations_dir=base / cond,
                                    ground_truth_path=t.ground_truth_path, target=t)
        res["condition"] = cond
        frames.append(res)
    allres = pd.concat(frames, ignore_index=True)

    out = config.VALIDATION_DIR / "faithfulness_ablation.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    allres.to_parquet(out, index=False)

    # --- Per-condition summary ---
    print("\n" + "=" * 66)
    print("FAITHFULNESS BY CONDITION (0=contradicts, 1=vague, 2=correct)")
    rows = []
    for cond in CONDITIONS:
        g = allres[allres["condition"] == cond]
        rows.append({"condition": cond, "n": len(g), "mean": g["score"].mean(),
                     "pct_correct(2)": 100 * (g["score"] == 2).mean(),
                     "pct_ok(>=1)": 100 * (g["score"] >= 1).mean()})
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    # --- Paired analysis: same (drug,mutation) across conditions ---
    wide = allres.pivot_table(index=["drug", "mutation"], columns="condition",
                              values="score")
    if {"full", "minimal"}.issubset(wide.columns):
        paired = wide.dropna(subset=["full", "minimal"])
        d = paired["full"] - paired["minimal"]
        print("\n" + "-" * 66)
        print(f"full vs minimal (paired, n={len(paired)}): "
              f"mean Δ={d.mean():+.2f}  "
              f"full>minimal:{(d > 0).sum()}  tie:{(d == 0).sum()}  full<minimal:{(d < 0).sum()}")
        if d.abs().sum() > 0:
            w = stats.wilcoxon(paired["full"], paired["minimal"], zero_method="wilcox")
            print(f"    Wilcoxon signed-rank p = {w.pvalue:.3f}")
    if {"full", "corrupted"}.issubset(wide.columns):
        paired = wide.dropna(subset=["full", "corrupted"])
        d = paired["full"] - paired["corrupted"]
        print(f"full vs corrupted (paired, n={len(paired)}): "
              f"mean Δ={d.mean():+.2f}  "
              f"full>corrupt:{(d > 0).sum()}  tie:{(d == 0).sum()}  full<corrupt:{(d < 0).sum()}")
        if d.abs().sum() > 0:
            w = stats.wilcoxon(paired["full"], paired["corrupted"], zero_method="wilcox")
            print(f"    Wilcoxon signed-rank p = {w.pvalue:.3f}")
    print("=" * 66)
    print(f"Saved per-pair scores: {out}")
    print("Read: if full ≈ minimal ≈ corrupted, faithfulness is Claude's prior, "
          "not the pipeline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
