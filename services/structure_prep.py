"""PDB download, cleaning, point mutagenesis, and PDBQT conversion.

Prepares a receptor (the active :class:`~targets.Target`'s wildtype PDB) and its
resistance mutants for AutoDock Vina docking. The scientific validity of every
downstream score depends on a handful of non-obvious, *target-specific* details,
all sourced from the Target so protease and RT are handled by one code path:

- The co-crystallized ligand's occupancy-weighted centroid (``target.
  ligand_hetcodes``) defines the docking box center (e.g. ``ROC`` = saquinavir
  in protease's 3OXC).
- ``target.chains`` are kept; a point mutation is applied to ``target.
  mutate_chains`` — both monomers for the protease homodimer (A+B), only p66
  (A) for the RT heterodimer.
- ``target.protonation`` (protease only): the catalytic Asp25 dyad is
  asymmetrically protonated — chain A neutral (``ASH``), chain B charged
  (``ASP``) — fixed after adding hydrogens for wildtype and every mutant. RT
  sets ``protonation=None`` and uses standard pH-7 hydrogens.

Dependencies: pdbfixer (from OpenMM), biopython, meeko.
"""

import shutil
import subprocess
import sys
from pathlib import Path

from Bio.PDB import PDBIO, PDBParser, Select
from openmm.app import PDBFile
from pdbfixer import PDBFixer

import config
from targets import Protonation, Target

# One-letter -> three-letter amino acid code, for PDBFixer mutation strings.
AA_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}

# Cached OpenMM force field for the clash-relief minimization of mutant side
# chains (created lazily on first use).
_FORCEFIELD = None


def _resolve(target: Target | None) -> Target:
    """Default an optional target to the active one."""
    return target if target is not None else config.ACTIVE_TARGET


def _get_forcefield():
    """Return a cached ``amber14-all`` force field for mutant minimization."""
    global _FORCEFIELD
    if _FORCEFIELD is None:
        from openmm.app import ForceField
        _FORCEFIELD = ForceField("amber14-all.xml")
    return _FORCEFIELD


# --- Structure cleaning ------------------------------------------------------

def _ligand_centroid(structure, het_codes: list[str]) -> tuple[float, float, float]:
    """Occupancy-weighted centroid of all atoms belonging to ``het_codes``.

    Saquinavir (ROC) is modeled in two altLoc conformations; BioPython exposes
    these as disordered atoms. Every altLoc is included, weighted by occupancy,
    so the result reproduces the deposited-frame centroid regardless of which
    conformer is "selected".
    """
    codes = {c.strip().upper() for c in het_codes}
    sx = sy = sz = sw = 0.0
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip().upper() not in codes:
                    continue
                for atom in residue:
                    # Expand disordered atoms into their individual altLocs.
                    variants = (
                        atom.disordered_get_list()
                        if atom.is_disordered()
                        else [atom]
                    )
                    for var in variants:
                        occ = var.get_occupancy()
                        occ = 1.0 if occ is None else float(occ)
                        x, y, z = var.get_coord()
                        sx += x * occ
                        sy += y * occ
                        sz += z * occ
                        sw += occ
        break  # first model only
    if sw == 0.0:
        raise ValueError(
            f"No atoms found for het codes {sorted(codes)}; cannot compute "
            f"docking box center."
        )
    return (sx / sw, sy / sw, sz / sw)


class _ReceptorSelect(Select):
    """Keep only standard protein residues in the requested chains.

    Rejects the ligand, sulfate/formate ions, waters, and every other
    heteroatom (their residue id hetero-flag is non-blank), and drops any
    chain not in ``chains_to_keep``.
    """

    def __init__(self, chains_to_keep, strip_waters=True):
        self.chains = {c.upper() for c in chains_to_keep}
        self.strip_waters = strip_waters

    def accept_chain(self, chain):
        return chain.id.upper() in self.chains

    def accept_residue(self, residue):
        hetflag, _, _ = residue.get_id()
        if hetflag != " ":
            # 'W' = water, 'H_XXX' = heteroatom (ligand, SO4, FMT, ...).
            return False
        return True

    def accept_atom(self, atom):
        # For disordered protein side chains, keep a single altLoc so the
        # receptor has one position per atom.
        if atom.is_disordered() and atom.get_altloc() not in (" ", "", "A"):
            return False
        return True


def clean_structure(
    pdb_path: Path,
    chains_to_keep: list[str] = None,
    het_codes_to_strip: list[str] = None,
    strip_waters: bool = True,
    output_path: Path = None,
    target: Target | None = None,
) -> Path:
    """Strip ligand/ions/waters, keep the target's chains, save ``wildtype.pdb``.

    BEFORE stripping, the occupancy-weighted centroid of the ligand
    (``het_codes_to_strip``, default = the target's ``ligand_hetcodes``) is
    computed and PRINTED — this is the value to set as the target's
    ``docking_center`` in targets.py. The printout also reports the offset from
    the currently configured center (for RT this reveals the placeholder box).
    """
    t = _resolve(target)
    if chains_to_keep is None:
        chains_to_keep = list(t.chains)
    if het_codes_to_strip is None:
        het_codes_to_strip = list(t.ligand_hetcodes)
    pdb_path = Path(pdb_path)
    if output_path is None:
        output_path = t.structures_dir / "wildtype.pdb"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("wt", str(pdb_path))

    # --- Docking box center: ligand centroid BEFORE stripping ---
    center = _ligand_centroid(structure, het_codes_to_strip)
    cx, cy, cz = center
    print(
        f"[{t.name}] Ligand ({','.join(het_codes_to_strip)}) occupancy-weighted "
        f"centroid: ({cx:.3f}, {cy:.3f}, {cz:.3f})"
    )
    dcx, dcy, dcz = t.docking_center
    print(
        f"  {t.name}.docking_center = ({dcx:.3f}, {dcy:.3f}, {dcz:.3f}); "
        f"offset = ({cx - dcx:+.3f}, {cy - dcy:+.3f}, {cz - dcz:+.3f}) A"
    )
    print("  -> Set this as docking_center for this target in targets.py if it differs.")

    # --- Strip and write the receptor ---
    io = PDBIO()
    io.set_structure(structure)
    io.save(
        str(output_path),
        select=_ReceptorSelect(chains_to_keep, strip_waters=strip_waters),
    )
    print(f"Cleaned receptor written to {output_path}")
    return output_path


# --- Hydrogens & catalytic-dyad protonation ----------------------------------

def _dyad_variants(topology, protonation: Protonation) -> list:
    """Build a per-residue variants list for asymmetric catalytic-dyad protonation.

    Returns a list aligned to ``topology.residues()`` (the shape
    ``Modeller.addHydrogens`` expects): the catalytic residue (the
    ``protonation.resnum``-th residue) of ``protonation.protonated_chain`` is set
    to ``"ASH"`` (protonated, neutral) and of ``protonation.deprotonated_chain``
    to ``"ASP"`` (deprotonated, charged); every other residue is ``None`` (pH
    default).

    Locating the catalytic residue by ordinal position is robust to PDBFixer's
    chain-B renumbering (chain B -> residues 101-199).
    """
    prot_chain = protonation.protonated_chain.upper()
    deprot_chain = protonation.deprotonated_chain.upper()
    ordinal = protonation.resnum
    variants: list = []
    for chain in topology.chains():
        cid = chain.id.upper()
        for i, res in enumerate(chain.residues()):
            is_catalytic = (i == ordinal - 1) and res.name in ("ASP", "ASH")
            if is_catalytic and cid == prot_chain:
                variants.append(protonation.resname_protonated)
            elif is_catalytic and cid == deprot_chain:
                variants.append(protonation.resname_deprotonated)
            else:
                variants.append(None)
    return variants


def _add_hydrogens(fixer: PDBFixer, protonation: Protonation | None,
                   ph: float = 7.0) -> None:
    """Add hydrogens to a PDBFixer structure.

    When ``protonation`` is given (HIV protease), uses ``Modeller.addHydrogens``
    with a per-residue ``variants`` list so the catalytic dyad is built
    asymmetrically (chain A ASH / chain B ASP). When ``None`` (e.g. RT), adds
    standard pH-``ph`` hydrogens. Updates ``fixer.topology``/``positions`` in place.
    """
    from openmm.app import Modeller

    modeller = Modeller(fixer.topology, fixer.positions)
    variants = _dyad_variants(modeller.topology, protonation) if protonation else None
    modeller.addHydrogens(pH=ph, variants=variants)
    fixer.topology = modeller.topology
    fixer.positions = modeller.positions


def add_hydrogens_and_fix_protonation(
    pdb_path: Path, ph: float = 7.0, target: Target | None = None,
) -> Path:
    """Add missing atoms/hydrogens, enforcing the target's dyad protonation.

    For a target with ``protonation`` (protease) this builds the asymmetric
    Asp25 dyad; for ``protonation=None`` (RT) it adds standard hydrogens.
    Overwrites ``pdb_path`` in place and returns it.
    """
    t = _resolve(target)
    pdb_path = Path(pdb_path)
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    _add_hydrogens(fixer, t.protonation, ph)

    with open(pdb_path, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    return pdb_path


# --- Point mutagenesis -------------------------------------------------------

def _minimize_mutated_residues(fixer: PDBFixer, mutated_atom_indices: set) -> bool:
    """Relieve steric clashes introduced by a mutation via local minimization.

    Builds an OpenMM system, freezes every atom except those in
    ``mutated_atom_indices`` (mass 0), and runs a short local minimization so
    only the newly built side chain(s) relax. Freezing the rest keeps the pocket
    in the crystal frame, so ``config.DOCKING_CENTER`` stays valid.

    Requires the topology to have a complete standard-protonation hydrogen set
    (so the force field can build a system). Returns True on success; on any
    force-field/minimization error it leaves positions unchanged and returns
    False (the caller falls back to the un-minimized structure).
    """
    try:
        from openmm import Context, LangevinIntegrator, LocalEnergyMinimizer, unit
        from openmm.app import NoCutoff

        forcefield = _get_forcefield()
        system = forcefield.createSystem(
            fixer.topology, nonbondedMethod=NoCutoff, constraints=None
        )
        for i in range(system.getNumParticles()):
            if i not in mutated_atom_indices:
                system.setParticleMass(i, 0.0)  # freeze
        integrator = LangevinIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
        )
        context = Context(system, integrator)
        context.setPositions(fixer.positions)
        LocalEnergyMinimizer.minimize(context, maxIterations=200)
        fixer.positions = context.getState(getPositions=True).getPositions()
        return True
    except Exception as exc:  # noqa: BLE001 - minimization is best-effort
        print(f"    (minimization skipped: {type(exc).__name__}: {exc})")
        return False


def generate_mutant(
    wildtype_pdb: Path,
    chain_id: str,          # kept for signature compatibility; target drives chains
    position: int,
    wildtype_aa: str,
    mutant_aa: str,
    output_dir: Path = None,
    target: Target | None = None,
) -> Path:
    """Generate a point-mutant receptor PDB using PDBFixer.

    The mutation is applied to every chain in ``target.mutate_chains`` — both
    monomers for the protease homodimer (A+B), only p66 (A) for the RT
    heterodimer. Missing atoms/hydrogens are rebuilt and the target's dyad
    protonation (protease only) is re-applied. Saves to
    ``{output_dir}/{wt}{position}{mut}.pdb`` and returns that path.
    """
    t = _resolve(target)
    wildtype_pdb = Path(wildtype_pdb)
    output_dir = Path(output_dir) if output_dir is not None else t.mutants_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    from openmm.app import Modeller
    from openmm.app.element import hydrogen

    mut3 = AA_ONE_TO_THREE[mutant_aa.upper()]
    mutation_name = f"{wildtype_aa.upper()}{position}{mutant_aa.upper()}"

    fixer = PDBFixer(filename=str(wildtype_pdb))

    # Start from a clean heavy-atom topology: strip any pre-existing hydrogens so
    # the mutated residue and (protease) the catalytic dyad are rehydrogenated
    # consistently.
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.delete([a for a in modeller.topology.atoms() if a.element == hydrogen])
    fixer.topology, fixer.positions = modeller.topology, modeller.positions

    # Map the position to each mutated chain's actual residue by ordinal position.
    # This is robust to PDBFixer's chain renumbering (protease chain B -> 101-199,
    # where residue 82 is stored as 182). We mutate FROM the residue actually
    # present in the structure, not the consensus wildtype: e.g. 3OXC carries a
    # few natural/engineered variants (I33, A67, A95), so a consensus-named
    # mutation like "L33F" must be applied as ILE-33-PHE on this structure.
    # NOTE: ordinal mapping assumes no unresolved residues before `position` in
    # the mutated chain; verify for RT crystal structures with N-terminal gaps.
    mutate = {c.upper() for c in t.mutate_chains}
    sites: list[tuple[str, int, str]] = []  # (chain, resseq, current_resname)
    for chain in fixer.topology.chains():
        cid = chain.id.upper()
        if cid not in mutate:
            continue
        residues = list(chain.residues())
        if position - 1 >= len(residues):
            continue
        res = residues[position - 1]
        sites.append((cid, int(res.id), res.name))

    if not sites:
        raise ValueError(f"{mutation_name}: no target chains found for mutation.")

    # Apply to each mutated chain, skipping chains that already carry the target
    # residue (a structural polymorphism matching the mutation target).
    for cid, resseq, current in sites:
        if current == mut3:
            continue
        fixer.applyMutations([f"{current}-{resseq}-{mut3}"], cid)

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Relieve clashes from the rigidly-placed side chain: add a standard-
    # protonation H set, minimize only the mutated residues (rest frozen), then
    # re-hydrogenate with the target's dyad protonation (protease) or standard
    # hydrogens (RT).
    fixer.addMissingHydrogens(7.0)
    mutated_atom_indices: set = set()
    for chain in fixer.topology.chains():
        if chain.id.upper() in mutate:
            res = list(chain.residues())[position - 1]
            mutated_atom_indices.update(a.index for a in res.atoms())
    _minimize_mutated_residues(fixer, mutated_atom_indices)

    strip = Modeller(fixer.topology, fixer.positions)
    strip.delete([a for a in strip.topology.atoms() if a.element == hydrogen])
    fixer.topology, fixer.positions = strip.topology, strip.positions
    _add_hydrogens(fixer, t.protonation, 7.0)

    out_path = output_dir / f"{mutation_name}.pdb"
    with open(out_path, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    return out_path


# --- Receptor PDBQT preparation ----------------------------------------------

def prepare_receptor_pdbqt(pdb_path: Path) -> Path:
    """Convert a receptor PDB to PDBQT for Vina using meeko.

    Uses ``--read_pdb`` (meeko's gemmi-based reader) rather than ``-i``, which
    in meeko 0.7.x requires the optional ProDy dependency. Tries the
    ``mk_prepare_receptor.py`` CLI first, then ``python -m
    meeko.cli.mk_prepare_receptor``. Returns the ``.pdbqt`` path.
    """
    pdb_path = Path(pdb_path)
    out_stem = pdb_path.with_suffix("")           # meeko appends .pdbqt
    pdbqt_path = pdb_path.with_suffix(".pdbqt")

    # --read_pdb uses gemmi (installed); -o sets the basename; -p writes PDBQT.
    base_args = ["--read_pdb", str(pdb_path), "-o", str(out_stem), "-p"]

    # Candidate invocations, in order of preference.
    invocations = []
    cli = shutil.which("mk_prepare_receptor.py") or shutil.which("mk_prepare_receptor")
    if cli:
        invocations.append([cli, *base_args])
    invocations.append([sys.executable, "-m", "meeko.cli.mk_prepare_receptor", *base_args])

    last_err = None
    for cmd in invocations:
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except FileNotFoundError as exc:
            last_err = exc
            continue
        except subprocess.CalledProcessError as exc:
            last_err = exc
            # A real failure (bad input) — surface stderr but keep trying
            # alternative entry points in case this one is just missing.
            continue
        if pdbqt_path.exists():
            return pdbqt_path

    # If meeko wrote to a slightly different name, try to locate it.
    if pdbqt_path.exists():
        return pdbqt_path
    msg = f"Failed to prepare receptor PDBQT for {pdb_path}."
    if isinstance(last_err, subprocess.CalledProcessError) and last_err.stderr:
        msg += f"\nLast error:\n{last_err.stderr}"
    raise RuntimeError(msg)


# --- Wildtype preparation ----------------------------------------------------

def prepare_wildtype(
    raw_pdb: Path = None,
    force: bool = False,
    target: Target | None = None,
) -> dict[str, Path]:
    """Run the full wildtype receptor pipeline: clean -> H+protonation -> PDBQT.

    Returns ``{"pdb": wildtype.pdb, "pdbqt": wildtype.pdbqt}`` under the target's
    structures dir. If both outputs already exist and ``force`` is False, the
    (cheap) clean + hydrogen steps are still re-run only when the PDB is missing;
    a present PDBQT is left as-is.
    """
    t = _resolve(target)
    if raw_pdb is None:
        raw_pdb = config.RAW_DIR / f"{t.pdb_id}.pdb"
    raw_pdb = Path(raw_pdb)

    wt_pdb = t.structures_dir / "wildtype.pdb"
    wt_pdbqt = t.structures_dir / "wildtype.pdbqt"

    if wt_pdbqt.exists() and wt_pdb.exists() and not force:
        print(f"[{t.name}] Wildtype already prepared: {wt_pdb.name}, {wt_pdbqt.name}")
        return {"pdb": wt_pdb, "pdbqt": wt_pdbqt}

    if not raw_pdb.exists():
        raise FileNotFoundError(
            f"Raw PDB not found at {raw_pdb}. Run scripts/01_download_data.py."
        )

    clean_structure(raw_pdb, target=t)
    add_hydrogens_and_fix_protonation(wt_pdb, target=t)
    pdbqt = prepare_receptor_pdbqt(wt_pdb)
    print(f"[{t.name}] Wildtype receptor ready: {wt_pdb.name}, {pdbqt.name}")
    return {"pdb": wt_pdb, "pdbqt": pdbqt}


# --- Mutant cache builder ----------------------------------------------------

def collect_unique_mutations(
    panels_dir: Path = None,
    target: Target | None = None,
) -> list[tuple[str, int, str, str]]:
    """Return the union of unique mutations across all of the target's panels.

    Each item is ``(mutation_name, position, wildtype_aa, mutant_aa)``, sorted
    by position then mutant residue. ``panels_dir`` defaults to the target's.
    """
    from services.mutation_panel import load_panel

    t = _resolve(target)
    panels_dir = Path(panels_dir) if panels_dir is not None else t.panels_dir

    mutations: dict[str, tuple[int, str, str]] = {}
    for parquet in sorted(Path(panels_dir).glob("*.parquet")):
        panel = load_panel(parquet.stem, panels_dir)
        for _, row in panel.iterrows():
            mutations[row["mutation"]] = (
                int(row["position"]),
                str(row["wildtype_aa"]),
                str(row["mutant_aa"]),
            )
    return [
        (name, pos, wt, mut)
        for name, (pos, wt, mut) in sorted(
            mutations.items(), key=lambda kv: (kv[1][0], kv[1][2])
        )
    ]


def build_mutant_cache(
    panels_dir: Path = None,
    limit: int = None,
    target: Target | None = None,
) -> list[Path]:
    """Generate mutant PDB + PDBQT for every unique mutation across the panels.

    Skips mutations whose PDBQT already exists. ``limit`` caps the number of
    mutations processed (useful for a quick smoke test). Returns the list of
    PDBQT paths that exist after the run. All paths default to the target's dirs.
    """
    t = _resolve(target)
    panels_dir = Path(panels_dir) if panels_dir is not None else t.panels_dir
    wildtype_pdb = t.structures_dir / "wildtype.pdb"
    if not wildtype_pdb.exists():
        raise FileNotFoundError(
            f"Wildtype receptor not found at {wildtype_pdb}. "
            f"Run structure_prep.prepare_wildtype() first."
        )

    ordered = collect_unique_mutations(panels_dir, target=t)
    if limit is not None:
        ordered = ordered[:limit]
    total = len(ordered)

    pdbqt_paths: list[Path] = []
    for i, (name, pos, wt, mut) in enumerate(ordered, start=1):
        pdbqt_path = t.mutants_dir / f"{name}.pdbqt"
        if pdbqt_path.exists():
            print(f"Skipping {name} (cached) ({i}/{total})")
            pdbqt_paths.append(pdbqt_path)
            continue
        try:
            mutant_pdb = generate_mutant(wildtype_pdb, "A", pos, wt, mut, target=t)
            pdbqt_paths.append(prepare_receptor_pdbqt(mutant_pdb))
            print(f"Generated {name} ({i}/{total})")
        except Exception as exc:  # noqa: BLE001 - keep the batch alive
            print(f"FAILED {name} ({i}/{total}): {type(exc).__name__}: {exc}")
    return pdbqt_paths
