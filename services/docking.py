"""SMILES -> 3D -> PDBQT -> AutoDock Vina docking pipeline.

Converts candidate compounds (SMILES) to Vina-ready PDBQT with RDKit + meeko,
then docks them against prepared HIV-1 protease receptors (wildtype + mutants).

Key facts enforced here:
- RDKit does NOT emit PDBQT; meeko does. ``AddHs`` must precede ``EmbedMolecule``.
- meeko 0.7's ``PDBQTWriterLegacy.write_string`` returns a 3-tuple
  ``(pdbqt_string, success, error_msg)`` and must be unpacked.
- Vina affinity is kcal/mol; MORE NEGATIVE = stronger binding.
"""

import json
import os
import shutil
import subprocess
import tempfile
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

import config

# Deterministic conformer generation for reproducible docking inputs.
_EMBED_SEED = 42

# Default Uni-Dock (GPU) executable name; override via dock_against_panel.
_UNIDOCK_BIN = "unidock"


def smiles_to_pdbqt(smiles: str) -> str:
    """Convert a SMILES string to a Vina-ready PDBQT string, in memory.

    Adds explicit hydrogens, embeds a 3D conformer (ETKDGv3), MMFF-optimizes the
    geometry, then writes PDBQT via meeko. Raises ``ValueError`` on an invalid
    SMILES, failed embedding, or a meeko write failure.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

    mol = Chem.AddHs(mol)  # explicit H required by meeko and for docking

    params = AllChem.ETKDGv3()
    params.randomSeed = _EMBED_SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        # Retry with random coordinates for awkward geometries.
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError(f"RDKit could not embed a 3D conformer for {smiles!r}")

    # Geometry cleanup (best-effort; some ligands lack MMFF parameters).
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:  # noqa: BLE001
        pass

    mk_prep = MoleculePreparation()
    molsetup_list = mk_prep(mol)  # returns a LIST
    pdbqt_string, success, error_msg = PDBQTWriterLegacy.write_string(molsetup_list[0])
    if not success:
        raise ValueError(f"meeko failed to write PDBQT for {smiles!r}: {error_msg}")
    return pdbqt_string


def dock_single(
    ligand_pdbqt: str,
    receptor_pdbqt_path: Path,
    center: tuple = config.DOCKING_CENTER,
    box_size: tuple = config.DOCKING_BOX_SIZE,
    exhaustiveness: int = config.VINA_EXHAUSTIVENESS,
    n_poses: int = config.VINA_NUM_POSES,
    cpu: int = config.VINA_NUM_CPUS,
    seed: int = 0,
) -> dict:
    """Dock one ligand (PDBQT string) against one receptor (PDBQT file).

    Returns ``{"delta_g": float, "n_poses": int, "all_energies": list[float]}``
    where ``delta_g`` is the best (most negative) Vina affinity in kcal/mol.
    ``seed=0`` lets Vina pick a random seed; pass a nonzero value for a
    reproducible run. On ANY failure returns
    ``{"delta_g": None, "n_poses": 0, "error": str}`` so a single bad dock never
    crashes a batch run.
    """
    try:
        from vina import Vina

        v = Vina(sf_name="vina", cpu=cpu, seed=seed, verbosity=0)
        v.set_receptor(str(receptor_pdbqt_path))
        v.set_ligand_from_string(ligand_pdbqt)
        v.compute_vina_maps(center=list(center), box_size=list(box_size))
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

        energies = v.energies()  # numpy array, shape (n_poses, >=1)
        all_energies = [float(row[0]) for row in energies]
        return {
            "delta_g": float(energies[0][0]),
            "n_poses": len(energies),
            "all_energies": all_energies,
        }
    except Exception as exc:  # noqa: BLE001 - keep batches alive
        return {"delta_g": None, "n_poses": 0, "error": str(exc)}


def dock_single_unidock(
    ligand_pdbqt_path: Path,
    receptor_pdbqt_path: Path,
    center: tuple = config.DOCKING_CENTER,
    box_size: tuple = config.DOCKING_BOX_SIZE,
    exhaustiveness: int = config.VINA_EXHAUSTIVENESS,
    n_poses: int = config.VINA_NUM_POSES,
    unidock_bin: str = _UNIDOCK_BIN,
    search_mode: str = None,
    seed: int = 0,
) -> dict:
    """GPU dock of one ligand against one receptor via Uni-Dock.

    Uni-Dock implements the SAME Vina scoring function on GPU, so affinities are
    directly comparable to :func:`dock_single`. Both inputs must be PDBQT files
    (Uni-Dock reads ligands from disk via ``--gpu_batch``, not from a string).

    ``search_mode`` (``fast``/``balance``/``detail``) overrides
    ``exhaustiveness`` when set — that is Uni-Dock's preferred effort knob.

    Returns the same dict shape as :func:`dock_single`. Requires ``unidock`` on
    PATH and a CUDA GPU; on any failure returns ``delta_g=None``.

    NOTE: exercised only where a GPU + Uni-Dock are installed. Run the validation
    in docs/UNIDOCK_SETUP.md before trusting a full batch.
    """
    ligand_pdbqt_path = Path(ligand_pdbqt_path)
    receptor_pdbqt_path = Path(receptor_pdbqt_path)
    if not receptor_pdbqt_path.exists():
        return {"delta_g": None, "n_poses": 0, "error": "missing receptor structure"}

    try:
        with tempfile.TemporaryDirectory() as td:
            outdir = Path(td) / "out"
            outdir.mkdir(parents=True, exist_ok=True)
            cmd = [
                unidock_bin,
                "--receptor", str(receptor_pdbqt_path),
                "--gpu_batch", str(ligand_pdbqt_path),
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(box_size[0]),
                "--size_y", str(box_size[1]),
                "--size_z", str(box_size[2]),
                "--num_modes", str(n_poses),
                "--scoring", "vina",
                "--dir", str(outdir),
            ]
            if search_mode:
                cmd += ["--search_mode", search_mode]
            else:
                cmd += ["--exhaustiveness", str(exhaustiveness)]
            if seed:
                cmd += ["--seed", str(seed)]

            subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Uni-Dock writes "<ligand_stem>_out.pdbqt" into --dir.
            stem = ligand_pdbqt_path.stem
            out = outdir / f"{stem}_out.pdbqt"
            if not out.exists():
                cands = sorted(outdir.glob("*_out.pdbqt")) or sorted(outdir.glob("*.pdbqt"))
                out = cands[0] if cands else None
            if out is None:
                return {"delta_g": None, "n_poses": 0, "error": "no Uni-Dock output"}

            energies = []
            for line in out.read_text().splitlines():
                if "VINA RESULT" in line:
                    energies.append(float(line.split(":", 1)[1].split()[0]))
            if not energies:
                return {"delta_g": None, "n_poses": 0, "error": "no VINA RESULT lines"}
            return {
                "delta_g": min(energies),  # most negative = best
                "n_poses": len(energies),
                "all_energies": energies,
            }
    except subprocess.CalledProcessError as exc:
        return {"delta_g": None, "n_poses": 0, "error": f"unidock failed: {exc.stderr}"}
    except Exception as exc:  # noqa: BLE001 - keep batches alive
        return {"delta_g": None, "n_poses": 0, "error": str(exc)}


def _dock_job(job: dict) -> dict:
    """Worker: dock one (mutation, receptor) pair, averaging over replicates.

    ``job`` carries the ligand PDBQT, receptor path, mutation label, docking box,
    and replicate count. Runs ``dock_single`` ``replicates`` times with distinct
    seeds and averages the best-pose affinity (more stable than a single
    stochastic run). Returns a per-mutation result row. Missing receptor files
    yield ``delta_g=None`` rather than raising.
    """
    receptor = Path(job["receptor_path"])
    label = job["mutation"]
    if not receptor.exists():
        return {"mutation": label, "delta_g": None, "n_poses": 0,
                "n_ok": 0, "error": "missing receptor structure"}

    replicates = job.get("replicates", 1)
    base_seed = job.get("seed", _EMBED_SEED)
    backend = job.get("backend", "vina")
    energies = []
    n_poses = 0
    last_err = None
    for r in range(replicates):
        if backend == "unidock":
            res = dock_single_unidock(
                job["ligand_pdbqt_path"], receptor,
                center=job["center"], box_size=job["box_size"],
                exhaustiveness=job["exhaustiveness"], n_poses=job["n_poses"],
                unidock_bin=job.get("unidock_bin", _UNIDOCK_BIN),
                search_mode=job.get("search_mode"), seed=base_seed + r,
            )
        else:
            res = dock_single(
                job["ligand_pdbqt"], receptor,
                center=job["center"], box_size=job["box_size"],
                exhaustiveness=job["exhaustiveness"], n_poses=job["n_poses"],
                cpu=job["cpu"], seed=base_seed + r,
            )
        if res["delta_g"] is not None:
            energies.append(res["delta_g"])
            n_poses = max(n_poses, res["n_poses"])
        else:
            last_err = res.get("error")

    if not energies:
        return {"mutation": label, "delta_g": None, "n_poses": 0,
                "n_ok": 0, "error": last_err}
    return {
        "mutation": label,
        "delta_g": sum(energies) / len(energies),
        "n_poses": n_poses,
        "n_ok": len(energies),
        "error": None,
    }


def dock_against_panel(
    smiles: str,
    drug_name: str,
    panel_df: pd.DataFrame,
    structures_dir: Path = config.STRUCTURES_DIR,
    mutants_dir: Path = config.MUTANTS_DIR,
    exhaustiveness: int = config.VINA_EXHAUSTIVENESS,
    n_poses: int = config.VINA_NUM_POSES,
    replicates: int = 1,
    n_workers: int = 1,
    cpu_per_worker: int = None,
    seed: int = _EMBED_SEED,
    backend: str = "vina",
    unidock_bin: str = _UNIDOCK_BIN,
    search_mode: str = None,
    progress: bool = True,
) -> pd.DataFrame:
    """Dock a compound against wildtype + every mutation in ``panel_df``.

    The SMILES is converted to PDBQT once and reused for all receptors.

    ``backend="vina"`` (default) uses CPU AutoDock Vina; ``backend="unidock"``
    uses the GPU Uni-Dock engine (same Vina scoring function, ~10-50x faster) and
    forces sequential execution since a GPU is not shared across processes.

    Threading (Vina only): by default (``n_workers=1``) docks run sequentially
    with each Vina using ALL cores — robust for long unattended runs.
    ``n_workers>1`` runs that many concurrent docks, each Vina limited to
    ``cpu_per_worker`` threads (default ``cores // n_workers``); parallel can
    thrash on memory-bandwidth-bound machines, so benchmark first.
    ``replicates`` averages multiple seeded runs per pair.

    Returns a DataFrame with columns ``drug, mutation, delta_g, delta_delta_g,
    n_poses, n_ok`` where ``mutation == "WT"`` is the wildtype reference and
    ``delta_delta_g = delta_g_mutant - delta_g_wildtype`` (positive = weaker
    binding = resistance).
    """
    ligand_pdbqt = smiles_to_pdbqt(smiles)

    # Resolve per-dock thread count: sequential uses all cores (cpu=0), parallel
    # splits the machine so workers do not oversubscribe.
    if n_workers <= 1:
        n_workers = 1
        cpu = 0 if cpu_per_worker is None else cpu_per_worker
    else:
        cores = os.cpu_count() or 2
        cpu = cpu_per_worker if cpu_per_worker else max(1, cores // n_workers)

    # Uni-Dock reads ligands from disk; write the ligand PDBQT once for the whole
    # panel and force sequential execution (single GPU).
    ligand_tmp_dir = None
    ligand_pdbqt_path = None
    if backend == "unidock":
        n_workers = 1
        ligand_tmp_dir = tempfile.mkdtemp(prefix=f"resistscope_{drug_name}_")
        ligand_pdbqt_path = str(Path(ligand_tmp_dir) / f"{drug_name}.pdbqt")
        Path(ligand_pdbqt_path).write_text(ligand_pdbqt)

    wildtype_pdbqt = Path(structures_dir) / "wildtype.pdbqt"

    def _make_job(mutation: str, receptor_path: str) -> dict:
        return {
            "mutation": mutation, "receptor_path": receptor_path,
            "ligand_pdbqt": ligand_pdbqt, "ligand_pdbqt_path": ligand_pdbqt_path,
            "center": config.DOCKING_CENTER, "box_size": config.DOCKING_BOX_SIZE,
            "exhaustiveness": exhaustiveness, "n_poses": n_poses, "cpu": cpu,
            "replicates": replicates, "seed": seed, "backend": backend,
            "unidock_bin": unidock_bin, "search_mode": search_mode,
        }

    jobs = [_make_job("WT", str(wildtype_pdbqt))]
    for mut in panel_df["mutation"]:
        jobs.append(_make_job(mut, str(Path(mutants_dir) / f"{mut}.pdbqt")))

    results: list[dict] = []
    try:
        if n_workers == 1:
            for i, job in enumerate(jobs, 1):
                results.append(_dock_job(job))
                if progress:
                    print(f"  [{drug_name}] {i}/{len(jobs)} {job['mutation']}")
        else:
            with Pool(processes=n_workers) as pool:
                for i, res in enumerate(pool.imap_unordered(_dock_job, jobs), 1):
                    results.append(res)
                    if progress and (i % 10 == 0 or i == len(jobs)):
                        print(f"  [{drug_name}] {i}/{len(jobs)} docked")
    finally:
        if ligand_tmp_dir:
            shutil.rmtree(ligand_tmp_dir, ignore_errors=True)

    by_mut = {r["mutation"]: r for r in results}
    wt = by_mut.get("WT", {})
    wt_dg = wt.get("delta_g")

    rows = []
    for r in results:
        dg = r["delta_g"]
        ddg = (dg - wt_dg) if (dg is not None and wt_dg is not None) else None
        rows.append({
            "drug": drug_name,
            "mutation": r["mutation"],
            "delta_g": dg,
            "delta_delta_g": ddg,
            "n_poses": r["n_poses"],
            "n_ok": r.get("n_ok", 0),
        })
    out = pd.DataFrame(rows)
    # Sort with WT first, then by mutation.
    out["_wt"] = (out["mutation"] != "WT").astype(int)
    out = out.sort_values(["_wt", "mutation"]).drop(columns="_wt").reset_index(drop=True)
    return out


def get_benchmark_smiles(
    cache_path: Path = None,
    refresh: bool = False,
) -> dict[str, str]:
    """Fetch SMILES for all benchmark PIs from PubChem, caching to JSON.

    Uses ``config.PUBCHEM_CIDS``. Returns ``{drug: smiles}``. Results are cached
    to ``data/ligands/benchmark_drugs.json`` and reused unless ``refresh``.
    """
    if cache_path is None:
        cache_path = config.DATA_DIR / "ligands" / "benchmark_drugs.json"
    cache_path = Path(cache_path)

    if cache_path.exists() and not refresh:
        return json.loads(cache_path.read_text())

    smiles: dict[str, str] = {}
    for drug, cid in config.PUBCHEM_CIDS.items():
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
            f"/property/SMILES,IsomericSMILES,CanonicalSMILES/JSON"
        )
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        props = resp.json()["PropertyTable"]["Properties"][0]
        # Prefer isomeric (stereochemistry matters for these chiral drugs).
        smiles[drug] = (
            props.get("IsomericSMILES")
            or props.get("SMILES")
            or props.get("CanonicalSMILES")
        )

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(smiles, indent=2))
    return smiles
