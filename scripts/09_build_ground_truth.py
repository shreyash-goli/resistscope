"""09: Agentically build the resistance-mechanism ground truth with Claude.

For each mutation, Claude runs a PubMed research loop (search → read abstracts →
submit) and a verification pass, producing a citation-grounded mechanism entry
(see services/literature_agent.py). Entries are merged into the target-scoped
``mechanism_ground_truth.json`` that the faithfulness eval (script 07) scores
against — this is what gives RT anything to score, and broadens PI coverage.

Hand-curated seed entries are never overwritten; agent entries are added
incrementally and cached (re-running skips mutations already present unless
--refresh). Unverified entries are written to a side log, not the scored file,
unless --include-unverified.

Usage::

    # protease — a couple of mutations to try it out
    python scripts/09_build_ground_truth.py --target pr --mutations I54V,V82T

    # reverse transcriptase — the canonical NNRTI DRMs (no structure needed)
    python scripts/09_build_ground_truth.py --target rt --primary --limit 10

    python scripts/09_build_ground_truth.py --target pr --primary   # all PI DRMs
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402
from services.literature_agent import research_mechanism, to_ground_truth_entry  # noqa: E402


def _load(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default="HIV1_PR", help="pr (default) or rt.")
    ap.add_argument("--mutations", default=None,
                    help="Comma-separated mutations, e.g. 'I54V,V82T'.")
    ap.add_argument("--primary", action="store_true",
                    help="Use the target's full primary-DRM set.")
    ap.add_argument("--limit", type=int, default=None, help="Cap how many to build.")
    ap.add_argument("--refresh", action="store_true",
                    help="Re-research mutations that already have an agent entry.")
    ap.add_argument("--include-unverified", action="store_true",
                    help="Write entries that failed the grounding check too.")
    args = ap.parse_args()

    t = config.set_active_target(args.target)
    print(f"Target: {t.name} ({t.label})")

    if args.mutations:
        muts = [m.strip() for m in args.mutations.split(",") if m.strip()]
    elif args.primary:
        muts = sorted(t.primary_mutations, key=lambda m: int(m[1:-1]))
    else:
        ap.error("give --mutations or --primary")
    if args.limit:
        muts = muts[: args.limit]

    gt_path = t.ground_truth_path
    log_path = t.validation_dir / "ground_truth_agent_log.json"
    gt = _load(gt_path)
    run_log = _load(log_path) if log_path.exists() else {}

    seed_drug = t.drug_columns[0] if t.drug_columns else next(iter(t.drugs))
    print(f"Seed drug for query framing: {seed_drug}; writing to {gt_path}\n")

    built, skipped, unverified = 0, 0, 0
    for mut in muts:
        existing = gt.get(mut)
        if existing and existing.get("provenance") != "agent" and not args.refresh:
            print(f"  [keep]  {mut}: hand-curated entry present — not touching.")
            skipped += 1
            continue
        if existing and existing.get("provenance") == "agent" and not args.refresh:
            print(f"  [cache] {mut}: agent entry present (use --refresh to redo).")
            skipped += 1
            continue

        print(f"  [research] {mut} …", end=" ", flush=True)
        res = research_mechanism(seed_drug, mut, target=t)
        run_log[mut] = to_ground_truth_entry(res) | {"verify_note": res.verify_note,
                                                      "n_searches": res.n_searches,
                                                      "error": res.error}
        if res.error:
            print(f"error: {res.error}")
            continue
        tag = "✓ verified" if res.verified else "⚠ unverified"
        print(f"{tag} · {res.n_searches} searches · {len(res.citations)} cites · {res.confidence}")

        if res.verified or args.include_unverified:
            gt[mut] = to_ground_truth_entry(res)
            built += 1
            if not res.verified:
                unverified += 1

    gt_path.parent.mkdir(parents=True, exist_ok=True)
    gt_path.write_text(json.dumps(gt, indent=2))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(run_log, indent=2))

    print(f"\nWrote {built} entries ({unverified} unverified) to {gt_path.name}; "
          f"{skipped} skipped. Full run log: {log_path}")
    print(f"Ground truth now covers {len(gt)} mutations for {t.name}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
