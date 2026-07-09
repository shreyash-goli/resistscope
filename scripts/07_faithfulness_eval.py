"""07: Claude-as-judge faithfulness evaluation of explanations.

Scores the cached Claude explanations (from scripts/06) against the curated
``data/mechanism_ground_truth.json`` mechanisms, using Claude as a judge on a
0-2 scale (0=contradicts, 1=consistent but vague, 2=correct primary mechanism).

Requires ANTHROPIC_API_KEY (or an ``ant auth login`` profile), and that
scripts/06 has generated explanations for at least some ground-truth mutations.

Usage::

    python scripts/07_faithfulness_eval.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

import config  # noqa: E402
from services.explanation import evaluate_faithfulness  # noqa: E402
from services.explanation import _faithfulness_pairs  # noqa: E402
import json  # noqa: E402

pd.set_option("display.width", 140)
pd.set_option("display.max_colwidth", 80)
pd.set_option("display.max_rows", 200)


def main() -> int:
    gt_path = config.DATA_DIR / "mechanism_ground_truth.json"
    ground_truth = json.loads(gt_path.read_text())
    pairs = _faithfulness_pairs(config.EXPLANATIONS_DIR, ground_truth)

    print(f"Ground-truth mechanisms: {len(ground_truth)}")
    print(f"Cached explanations matching ground truth: {len(pairs)}")
    if not pairs:
        print("\nNothing to evaluate. Generate explanations first:")
        print("  python scripts/06_generate_explanations.py --primary-only")
        return 1

    print(f"Judging {len(pairs)} explanation(s) with model={config.CLAUDE_MODEL} ...\n")
    try:
        results = evaluate_faithfulness()
    except Exception as exc:  # noqa: BLE001
        msg = str(exc)
        if "authentication" in msg.lower() or "api_key" in msg.lower():
            print("ERROR: no Anthropic credentials found. Set ANTHROPIC_API_KEY "
                  "(or run `ant auth login`) and re-run.")
            return 1
        raise

    # Per-explanation table.
    print(results.to_string(index=False))

    # Summary.
    n = len(results)
    mean_score = results["score"].mean()
    dist = results["score"].value_counts().reindex([0, 1, 2], fill_value=0)
    pct_correct = 100.0 * (results["score"] == 2).mean()
    pct_ok = 100.0 * (results["score"] >= 1).mean()

    out = config.VALIDATION_DIR / "faithfulness_scores.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    results.to_parquet(out, index=False)

    print("\n" + "=" * 60)
    print(f"Explanations judged : {n}")
    print(f"Score distribution  : 0={dist[0]}  1={dist[1]}  2={dist[2]}")
    print(f"Mean faithfulness   : {mean_score:.2f} / 2")
    print(f"Correct mechanism   : {pct_correct:.0f}% scored 2")
    print(f"Non-contradictory   : {pct_ok:.0f}% scored >= 1")
    print("=" * 60)
    print(f"Saved: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
