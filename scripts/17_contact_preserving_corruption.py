"""17: Contact-preserving corruption — does the pipeline add more than a boolean?

The `corrupted` arm of scripts/11 deranges every pipeline-derived field, which
includes flipping ``contacts_ligand_directly``. A skeptic can say that for a
mechanism-classification task, flipping that boolean hands the model the answer
key upside down: the drop in faithfulness would then show only that the model
reads its input, not that the *finer* geometry carries scientific value.

This script runs the cleaner control. We derange the pipeline fields **within
contact class** — contacting residues donate only to contacting residues and
non-contacting to non-contacting — so that:

  * ``contacts_ligand_directly`` keeps its TRUE value for every pair, and
  * the deranged ``min_distance_to_ligand`` stays consistent with it (both sides
    of the 4 Å threshold agree), avoiding a self-contradictory prompt.

The model therefore still knows whether the residue touches the ligand, but gets
the wrong distance, wrong subpocket and wrong energy. If faithfulness still
drops relative to `full`, the pipeline contributes beyond the binary contact
call. If it does not, the skeptic is right and we report that.

Usage::

    python scripts/17_contact_preserving_corruption.py --limit 3   # smoke test
    python scripts/17_contact_preserving_corruption.py             # full n=46
    python scripts/17_contact_preserving_corruption.py --target rt
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402
from services.ablation import _CORRUPT_FIELDS, generate_for_condition  # noqa: E402
from services.explanation import _faithfulness_pairs, evaluate_faithfulness  # noqa: E402

CONDITION = "corrupted_geom"
CONTACT_FIELD = "contacts_ligand_directly"
# Everything the standard corruption deranges, minus the contact boolean itself.
GEOM_FIELDS = tuple(f for f in _CORRUPT_FIELDS if f != CONTACT_FIELD)


def _load_real(explanations_dir: Path, drug: str, mutation: str):
    rec = json.loads((explanations_dir / f"{drug}_{mutation}.json").read_text())
    return rec.get("structural_context") or {}, rec.get("delta_delta_g")


def stratified_derange(reals: list[tuple[dict, float | None]],
                       muts: list[str], seed: int = 0):
    """Derange geometry/energy within contact class, at the level of *mutations*.

    Preserves ``contacts_ligand_directly`` exactly and keeps the donated distance
    on the same side of the contact threshold, so the corrupted context is wrong
    but not self-contradictory.

    The derangement is over distinct mutations rather than (drug, mutation)
    pairs: the same mutation recurs across drugs with identical geometry, so a
    pair-level derangement can donate a pair's own geometry back to it and leave
    the context effectively uncorrupted.
    """
    rng = np.random.default_rng(seed)
    out: list[tuple[dict, float | None]] = [(dict(c), d) for c, d in reals]

    groups: dict[bool, list[str]] = {}
    for (ctx, _), m in zip(reals, muts):
        groups.setdefault(bool(ctx.get(CONTACT_FIELD)), []).append(m)
    groups = {k: sorted(set(v)) for k, v in groups.items()}

    donor_of: dict[str, str] = {}
    for flag, uniq in groups.items():
        if len(uniq) < 2:
            print(f"  contact={flag}: only {len(uniq)} distinct mutation(s); unchanged")
            continue
        perm = np.arange(len(uniq))
        for _ in range(10_000):
            rng.shuffle(perm)
            if not np.any(perm == np.arange(len(uniq))):
                break
        for a, b in enumerate(perm):
            donor_of[uniq[a]] = uniq[int(b)]
        print(f"  contact={flag}: {len(uniq)} distinct mutations deranged within class")

    # one representative (context, ΔΔG) per donor mutation
    rep: dict[str, tuple[dict, float | None]] = {}
    for (ctx, ddg), m in zip(reals, muts):
        rep.setdefault(m, (ctx, ddg))

    for i, ((ctx, _ddg), m) in enumerate(zip(reals, muts)):
        d = donor_of.get(m)
        if d is None:
            continue
        donor_ctx, donor_ddg = rep[d]
        new = dict(ctx)
        for f in GEOM_FIELDS:
            if f in donor_ctx:
                new[f] = donor_ctx[f]
            else:
                new.pop(f, None)
        new[CONTACT_FIELD] = ctx.get(CONTACT_FIELD)      # keep the TRUE value
        out[i] = (new, donor_ddg)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", default="HIV1_PR")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--sleep", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    t = config.set_active_target(args.target)
    gt = json.loads(t.ground_truth_path.read_text())
    pairs = _faithfulness_pairs(t.explanations_dir, gt, target=t)
    if args.limit:
        pairs = pairs[: args.limit]
    print(f"Target {t.name}: {len(pairs)} pairs -> condition '{CONDITION}'")

    reals = [_load_real(t.explanations_dir, p["drug"], p["mutation"]) for p in pairs]
    corrupted = stratified_derange(reals, [p["mutation"] for p in pairs], seed=args.seed)

    # sanity: contact flag must be untouched everywhere
    bad = [i for i, ((a, _), (b, _)) in enumerate(zip(reals, corrupted))
           if a.get(CONTACT_FIELD) != b.get(CONTACT_FIELD)]
    assert not bad, f"contact flag changed for {len(bad)} pairs"
    changed = sum(1 for (a, _), (b, _) in zip(reals, corrupted)
                  if a.get("min_distance_to_ligand_angstrom")
                  != b.get("min_distance_to_ligand_angstrom"))
    print(f"contact flag preserved for all {len(pairs)} pairs; "
          f"distance changed for {changed}")

    base = config.DATA_DIR / ("ablation" if t.name == "HIV1_PR" else f"ablation_{t.name}")
    cond_dir = base / CONDITION
    n_new = 0
    for i, p in enumerate(pairs):
        ctx, ddg = corrupted[i]
        existed = (cond_dir / f"{p['drug']}_{p['mutation']}.json").exists()
        try:
            generate_for_condition(p["drug"], p["mutation"], ddg, ctx,
                                   "corrupted", cond_dir, target=t)
        except Exception as exc:  # noqa: BLE001
            print(f"  FAILED {p['drug']} {p['mutation']}: {exc}")
            continue
        if not existed:
            n_new += 1
            time.sleep(args.sleep)
    print(f"generated {n_new} new -> {cond_dir}")

    res = evaluate_faithfulness(explanations_dir=cond_dir,
                                ground_truth_path=t.ground_truth_path, target=t)
    res["condition"] = CONDITION
    out = config.VALIDATION_DIR / "faithfulness_ablation_contactpreserved.parquet"
    res.to_parquet(out, index=False)

    # --- compare against the existing full / corrupted arms -------------------
    prev = pd.read_parquet(config.VALIDATION_DIR / "faithfulness_ablation.parquet")
    allres = pd.concat([prev, res], ignore_index=True)
    key = ["drug", "mutation"]
    piv = allres.pivot_table(index=key, columns="condition", values="score",
                             aggfunc="first").dropna(
                                 subset=["full", "corrupted", CONDITION])

    print("\n" + "=" * 68)
    print(f"n paired = {len(piv)}")
    for c in ("full", "minimal", "corrupted", CONDITION):
        if c in piv:
            print(f"  {c:16s} mean={piv[c].mean():.3f}  "
                  f"%correct={100*(piv[c]==2).mean():.0f}%")
    for a, b in (("full", CONDITION), ("full", "corrupted"),
                 (CONDITION, "corrupted")):
        d = piv[a] - piv[b]
        if d.abs().sum() == 0:
            print(f"\n{a} - {b}: identical scores")
            continue
        w, p = stats.wilcoxon(piv[a], piv[b], alternative="greater")
        print(f"\n{a} - {b}: Δ={d.mean():+.3f}  Wilcoxon one-sided p={p:.4g} "
              f"({int((d>0).sum())} up / {int((d<0).sum())} down / {int((d==0).sum())} tied)")

    piv.to_csv(config.VALIDATION_DIR / "contactpreserved_paired.csv")
    print(f"\nwrote {out.name} and contactpreserved_paired.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
