"""04 (GPU): Batched Uni-Dock benchmark — all drugs vs each receptor per call.

Uni-Dock's cost is dominated by per-receptor grid construction + CUDA init
(~3.4 s), while each additional ligand in a ``--gpu_batch`` call adds only
~0.15 s. So instead of one call per (drug, mutant) pair, this docks ALL benchmark
drugs against each receptor in a SINGLE call and iterates over receptors. That
turns ~1,600 calls into ~286, cutting the full sweep from hours to ~20 min.

Also: always uses ``--search_mode`` (never raw ``--exhaustiveness``, which on
Uni-Dock triggers a pathological ~90 s/dock step count). ``detail`` recovers
affinities closest to CPU Vina (WT ~ -9.1 vs -9.35) at the same speed as balance.

Usage (on the GPU box, env activated)::

    python scripts/04_gpu_batch.py --subset all --search-mode detail --replicates 3
    python scripts/04_gpu_batch.py --subset primary          # quick
    python scripts/04_gpu_batch.py --drugs DRV NFV --subset prevalence --min-isolates 20

Output matches scripts/04: data/docking_results/benchmark_docking.parquet with
columns drug, mutation, delta_g, delta_delta_g, n_poses, n_ok.
"""

import argparse
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

import config  # noqa: E402
from services.docking import get_benchmark_smiles, smiles_to_pdbqt  # noqa: E402
from services.mutation_panel import load_panel  # noqa: E402


def _select(panel: pd.DataFrame, subset: str, min_isolates: int) -> pd.DataFrame:
    if subset == "all":
        return panel
    if subset == "primary":
        return panel[panel["is_primary"]].reset_index(drop=True)
    if subset == "prevalence":
        return panel[panel["n_isolates"] >= min_isolates].reset_index(drop=True)
    raise ValueError(subset)


def _parse_best_affinity(pdbqt_path: Path):
    """Best (most negative) REMARK VINA RESULT affinity from a Uni-Dock output."""
    if not pdbqt_path.exists():
        return None, 0
    energies = []
    for line in pdbqt_path.read_text().splitlines():
        if "VINA RESULT" in line:
            try:
                energies.append(float(line.split(":", 1)[1].split()[0]))
            except (ValueError, IndexError):
                pass
    if not energies:
        return None, 0
    return min(energies), len(energies)


def dock_receptor_batch(receptor_pdbqt: Path, ligand_files: dict, search_mode: str,
                        n_poses: int, unidock_bin: str) -> dict:
    """Dock all ligands against one receptor in a single Uni-Dock call.

    ``ligand_files`` maps drug -> ligand pdbqt Path. Returns drug -> (delta_g,
    n_poses). Missing receptor or failed call yields all-None.
    """
    if not Path(receptor_pdbqt).exists():
        return {d: (None, 0) for d in ligand_files}

    with tempfile.TemporaryDirectory() as td:
        outdir = Path(td)
        cmd = [
            unidock_bin,
            "--receptor", str(receptor_pdbqt),
            "--gpu_batch", *[str(p) for p in ligand_files.values()],
            "--center_x", str(config.DOCKING_CENTER[0]),
            "--center_y", str(config.DOCKING_CENTER[1]),
            "--center_z", str(config.DOCKING_CENTER[2]),
            "--size_x", str(config.DOCKING_BOX_SIZE[0]),
            "--size_y", str(config.DOCKING_BOX_SIZE[1]),
            "--size_z", str(config.DOCKING_BOX_SIZE[2]),
            "--scoring", "vina",
            "--num_modes", str(n_poses),
            "--search_mode", search_mode,
            "--dir", str(outdir),
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            print(f"    unidock failed on {Path(receptor_pdbqt).stem}: "
                  f"{exc.stderr.strip().splitlines()[-1] if exc.stderr else exc}")
            return {d: (None, 0) for d in ligand_files}

        out = {}
        for drug, lig in ligand_files.items():
            # Uni-Dock writes "<ligand_stem>_out.pdbqt" into --dir.
            cand = outdir / f"{Path(lig).stem}_out.pdbqt"
            if not cand.exists():
                matches = sorted(outdir.glob(f"{Path(lig).stem}*.pdbqt"))
                cand = matches[0] if matches else cand
            out[drug] = _parse_best_affinity(cand)
        return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", default="HIV1_PR",
                        help="Docking target: HIV1_PR / pr (default) or HIV1_RT / rt.")
    parser.add_argument("--subset", choices=["primary", "prevalence", "all"],
                        default="all")
    parser.add_argument("--min-isolates", type=int, default=20)
    parser.add_argument("--drugs", nargs="+", default=None)
    parser.add_argument("--search-mode", choices=["fast", "balance", "detail"],
                        default="detail")
    parser.add_argument("--replicates", type=int, default=1)
    parser.add_argument("--n-poses", type=int, default=config.VINA_NUM_POSES)
    parser.add_argument("--unidock-bin", default="unidock")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    t = config.set_active_target(args.target)
    if args.out is None:
        args.out = config.DOCKING_DIR / "benchmark_docking.parquet"
    print(f"Target: {t.name} ({t.label})")

    # Target-scoped benchmark-SMILES cache so RT never reuses PI's cached SMILES.
    ligands_cache = (config.DATA_DIR / t.subdir if t.subdir else config.DATA_DIR) \
        / "ligands" / "benchmark_drugs.json"
    smiles = get_benchmark_smiles(cache_path=ligands_cache)
    drugs = args.drugs or [
        d for d in config.PI_DRUGS
        if (config.PANELS_DIR / f"{d}.parquet").exists() and d in smiles
    ]

    # Per-drug target mutations (subset filter), and the union of receptors to dock.
    targets = {}          # drug -> set of mutation names
    union = set()
    for drug in drugs:
        sub = _select(load_panel(drug, config.PANELS_DIR), args.subset, args.min_isolates)
        targets[drug] = set(sub["mutation"])
        union |= targets[drug]
    # Only dock mutants whose receptor PDBQT actually exists.
    receptors = {"WT": config.STRUCTURES_DIR / "wildtype.pdbqt"}
    for mut in sorted(union):
        p = config.MUTANTS_DIR / f"{mut}.pdbqt"
        if p.exists():
            receptors[mut] = p
    missing = sorted(m for m in union if m not in receptors)

    print(f"Drugs: {drugs}")
    print(f"Subset: {args.subset} | search_mode: {args.search_mode} | "
          f"replicates: {args.replicates}")
    print(f"Receptors to dock: {len(receptors)} (WT + {len(receptors)-1} mutants)")
    if missing:
        print(f"Skipping {len(missing)} mutant(s) with no receptor structure: {missing}")

    # Convert each drug SMILES to a ligand PDBQT file once (named by drug so the
    # Uni-Dock output "<drug>_out.pdbqt" is identifiable).
    lig_dir = Path(tempfile.mkdtemp(prefix="resistscope_ligs_"))
    ligand_files = {}
    for drug in drugs:
        lp = lig_dir / f"{drug}.pdbqt"
        lp.write_text(smiles_to_pdbqt(smiles[drug]))
        ligand_files[drug] = lp

    # Dock: iterate receptors, batch all drugs, average replicates.
    t0 = time.time()
    # accum[(drug, mut)] = [list of delta_g across replicates], and n_poses
    accum = {}
    n_poses_seen = {}
    rec_items = list(receptors.items())
    for i, (mut, recpath) in enumerate(rec_items, 1):
        reps = []
        for _ in range(args.replicates):
            reps.append(dock_receptor_batch(recpath, ligand_files, args.search_mode,
                                            args.n_poses, args.unidock_bin))
        for drug in drugs:
            vals = [r[drug][0] for r in reps if r[drug][0] is not None]
            npose = max((r[drug][1] for r in reps), default=0)
            accum[(drug, mut)] = (sum(vals) / len(vals)) if vals else None
            n_poses_seen[(drug, mut)] = npose
        if i % 10 == 0 or i == len(rec_items):
            print(f"  {i}/{len(rec_items)} receptors  "
                  f"({(time.time()-t0)/60:.1f} min elapsed)")

    # Assemble rows: WT per drug + each drug's target mutations.
    rows = []
    for drug in drugs:
        wt_dg = accum.get((drug, "WT"))
        rows.append({"drug": drug, "mutation": "WT", "delta_g": wt_dg,
                     "delta_delta_g": 0.0 if wt_dg is not None else None,
                     "n_poses": n_poses_seen.get((drug, "WT"), 0), "n_ok": 1})
        for mut in sorted(targets[drug]):
            dg = accum.get((drug, mut))
            ddg = (dg - wt_dg) if (dg is not None and wt_dg is not None) else None
            rows.append({"drug": drug, "mutation": mut, "delta_g": dg,
                         "delta_delta_g": ddg,
                         "n_poses": n_poses_seen.get((drug, mut), 0),
                         "n_ok": 1 if dg is not None else 0})

    results = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    results.to_parquet(args.out, index=False)

    print("\n" + "=" * 64)
    print(f"{'drug':<6}{'docked':>10}{'failed':>10}{'mean_ddG':>12}{'max_ddG':>12}")
    print("-" * 64)
    for drug in drugs:
        r = results[(results.drug == drug) & (results.mutation != "WT")]
        docked = int(r.delta_g.notna().sum())
        failed = int(r.delta_g.isna().sum())
        print(f"{drug:<6}{docked:>10}{failed:>10}"
              f"{r.delta_delta_g.mean():>12.3f}{r.delta_delta_g.max():>12.3f}")
    print("=" * 64)
    print(f"Total wall time: {(time.time()-t0)/60:.1f} min")
    print(f"Saved: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
