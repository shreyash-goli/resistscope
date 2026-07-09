"""PDB download, cleaning, point mutagenesis, and PDBQT conversion.

Prepares the HIV-1 protease receptor (PDB 3OXC) and its resistance mutants for
AutoDock Vina docking. The scientific validity of every downstream score
depends on a handful of non-obvious details enforced here:

- The co-crystallized ligand in 3OXC is **saquinavir with het code ``ROC``**
  (not ``SQV`` / ``938``). Its occupancy-weighted centroid defines the docking
  box center.
- HIV-1 protease is a **C2-symmetric homodimer**; both chains A and B are kept
  and every point mutation is applied to both chains.
- The catalytic **Asp25 dyad is asymmetrically protonated** — chain A neutral
  (protonated), chain B charged (deprotonated). This is fixed after adding
  hydrogens, for the wildtype and every mutant.

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

# One-letter -> three-letter amino acid code, for PDBFixer mutation strings.
AA_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}

# The catalytic aspartate residue number in each protease monomer.
ASP25_RESNUM = 25

# Cached OpenMM force field for the clash-relief minimization of mutant side
# chains (created lazily on first use).
_FORCEFIELD = None


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
    chains_to_keep: list[str] = config.CHAINS_TO_KEEP,
    het_codes_to_strip: list[str] = config.LIGAND_HETCODES_TO_STRIP,
    strip_waters: bool = True,
    output_path: Path = None,
) -> Path:
    """Strip ligand/ions/waters, keep chains A+B, save as ``wildtype.pdb``.

    BEFORE stripping, the occupancy-weighted centroid of the ligand
    (``het_codes_to_strip``, default ``["ROC"]`` = saquinavir) is computed and
    PRINTED — this is the value to use for ``config.DOCKING_CENTER``. The
    printout also reports the offset from the current configured center.

    Note: the co-crystal ligand het code in 3OXC is ``ROC``, not ``938``.
    """
    pdb_path = Path(pdb_path)
    if output_path is None:
        output_path = config.STRUCTURES_DIR / "wildtype.pdb"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("wt", str(pdb_path))

    # --- Docking box center: ligand centroid BEFORE stripping ---
    center = _ligand_centroid(structure, het_codes_to_strip)
    cx, cy, cz = center
    print(
        f"Ligand ({','.join(het_codes_to_strip)}) occupancy-weighted centroid: "
        f"({cx:.3f}, {cy:.3f}, {cz:.3f})"
    )
    dcx, dcy, dcz = config.DOCKING_CENTER
    print(
        f"  config.DOCKING_CENTER = ({dcx:.3f}, {dcy:.3f}, {dcz:.3f}); "
        f"offset = ({cx - dcx:+.3f}, {cy - dcy:+.3f}, {cz - dcz:+.3f}) A"
    )
    print("  -> Set this as DOCKING_CENTER in config.py if it differs.")

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

def _asp25_variants(topology) -> list:
    """Build a per-residue variants list for asymmetric Asp25 protonation.

    Returns a list aligned to ``topology.residues()`` (the shape
    ``Modeller.addHydrogens`` expects): the catalytic Asp25 (the 25th residue)
    of ``config.ASP25_PROTONATED_CHAIN`` is set to ``"ASH"`` (protonated,
    neutral) and of ``config.ASP25_DEPROTONATED_CHAIN`` to ``"ASP"``
    (deprotonated, charged); every other residue is ``None`` (pH default).

    Locating Asp25 by ordinal position is robust to PDBFixer's chain-B
    renumbering (chain B -> residues 101-199).
    """
    prot_chain = config.ASP25_PROTONATED_CHAIN.upper()
    deprot_chain = config.ASP25_DEPROTONATED_CHAIN.upper()
    variants: list = []
    for chain in topology.chains():
        cid = chain.id.upper()
        for i, res in enumerate(chain.residues()):
            is_catalytic = (i == ASP25_RESNUM - 1) and res.name in ("ASP", "ASH")
            if is_catalytic and cid == prot_chain:
                variants.append("ASH")
            elif is_catalytic and cid == deprot_chain:
                variants.append("ASP")
            else:
                variants.append(None)
    return variants


def _add_hydrogens_with_asp25_fix(fixer: PDBFixer, ph: float = 7.0) -> None:
    """Add hydrogens to a PDBFixer structure with correct Asp25 protonation.

    Uses ``openmm.app.Modeller.addHydrogens`` with a per-residue ``variants``
    list so the catalytic dyad is built asymmetrically (chain A ASH / chain B
    ASP) with proper hydrogen geometry, instead of hand-editing atoms
    afterwards. Updates ``fixer.topology`` and ``fixer.positions`` in place.
    """
    from openmm.app import Modeller

    modeller = Modeller(fixer.topology, fixer.positions)
    variants = _asp25_variants(modeller.topology)
    modeller.addHydrogens(pH=ph, variants=variants)
    fixer.topology = modeller.topology
    fixer.positions = modeller.positions


def add_hydrogens_and_fix_protonation(pdb_path: Path, ph: float = 7.0) -> Path:
    """Add missing atoms/hydrogens and enforce the Asp25 protonation asymmetry.

    Overwrites ``pdb_path`` in place with the hydrogenated, protonation-fixed
    structure and returns it.
    """
    pdb_path = Path(pdb_path)
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    _add_hydrogens_with_asp25_fix(fixer, ph)

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
    chain_id: str,          # kept for signature compatibility; both chains used
    position: int,
    wildtype_aa: str,
    mutant_aa: str,
    output_dir: Path = config.MUTANTS_DIR,
) -> Path:
    """Generate a point-mutant receptor PDB using PDBFixer.

    The mutation is applied to BOTH chains A and B (HIV protease is a homodimer,
    and the major DRMs occur in both monomers). Missing atoms/hydrogens are
    rebuilt and the Asp25 protonation asymmetry is re-applied. Saves to
    ``{output_dir}/{wt}{position}{mut}.pdb`` and returns that path.
    """
    wildtype_pdb = Path(wildtype_pdb)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    from openmm.app import Modeller
    from openmm.app.element import hydrogen

    mut3 = AA_ONE_TO_THREE[mutant_aa.upper()]
    mutation_name = f"{wildtype_aa.upper()}{position}{mutant_aa.upper()}"

    fixer = PDBFixer(filename=str(wildtype_pdb))

    # Start from a clean heavy-atom topology: strip any pre-existing hydrogens so
    # the mutated residue and the Asp25 dyad are rehydrogenated consistently.
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.delete([a for a in modeller.topology.atoms() if a.element == hydrogen])
    fixer.topology, fixer.positions = modeller.topology, modeller.positions

    # Map the protease position to each chain's actual residue by ordinal
    # position. This is robust to PDBFixer's chain-B renumbering (101-199), where
    # protease residue 82 is stored as residue 182. We mutate FROM the residue
    # actually present in the structure, not the consensus wildtype: 3OXC carries
    # a few natural/engineered variants (e.g. I33, A67, A95), so a consensus-named
    # mutation like "L33F" must be applied as ILE-33-PHE on this structure.
    keep = {c.upper() for c in config.CHAINS_TO_KEEP}
    targets: list[tuple[str, int, str]] = []  # (chain, resseq, current_resname)
    for chain in fixer.topology.chains():
        cid = chain.id.upper()
        if cid not in keep:
            continue
        residues = list(chain.residues())
        if position - 1 >= len(residues):
            continue
        target = residues[position - 1]
        targets.append((cid, int(target.id), target.name))

    if not targets:
        raise ValueError(f"{mutation_name}: no target chains found for mutation.")

    # Apply to BOTH chains of the homodimer, skipping chains that already carry
    # the target residue (a 3OXC polymorphism matching the mutation target).
    for cid, resseq, current in targets:
        if current == mut3:
            continue
        fixer.applyMutations([f"{current}-{resseq}-{mut3}"], cid)

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Relieve clashes from the rigidly-placed side chain: add a standard-
    # protonation H set, minimize only the mutated residues (rest frozen), then
    # re-hydrogenate with the correct asymmetric Asp25 protonation.
    fixer.addMissingHydrogens(7.0)
    mutated_atom_indices: set = set()
    for chain in fixer.topology.chains():
        if chain.id.upper() in keep:
            res = list(chain.residues())[position - 1]
            mutated_atom_indices.update(a.index for a in res.atoms())
    _minimize_mutated_residues(fixer, mutated_atom_indices)

    strip = Modeller(fixer.topology, fixer.positions)
    strip.delete([a for a in strip.topology.atoms() if a.element == hydrogen])
    fixer.topology, fixer.positions = strip.topology, strip.positions
    _add_hydrogens_with_asp25_fix(fixer, 7.0)

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
) -> dict[str, Path]:
    """Run the full wildtype receptor pipeline: clean -> H+protonation -> PDBQT.

    Returns ``{"pdb": wildtype.pdb, "pdbqt": wildtype.pdbqt}``. If both outputs
    already exist and ``force`` is False, the (cheap) clean + hydrogen steps are
    still re-run only when the PDB is missing; a present PDBQT is left as-is.
    """
    if raw_pdb is None:
        raw_pdb = config.RAW_DIR / f"{config.WILDTYPE_PDB_ID}.pdb"
    raw_pdb = Path(raw_pdb)

    wt_pdb = config.STRUCTURES_DIR / "wildtype.pdb"
    wt_pdbqt = config.STRUCTURES_DIR / "wildtype.pdbqt"

    if wt_pdbqt.exists() and wt_pdb.exists() and not force:
        print(f"Wildtype already prepared: {wt_pdb.name}, {wt_pdbqt.name}")
        return {"pdb": wt_pdb, "pdbqt": wt_pdbqt}

    if not raw_pdb.exists():
        raise FileNotFoundError(
            f"Raw PDB not found at {raw_pdb}. Run scripts/01_download_data.py."
        )

    clean_structure(raw_pdb)
    add_hydrogens_and_fix_protonation(wt_pdb)
    pdbqt = prepare_receptor_pdbqt(wt_pdb)
    print(f"Wildtype receptor ready: {wt_pdb.name}, {pdbqt.name}")
    return {"pdb": wt_pdb, "pdbqt": pdbqt}


# --- Mutant cache builder ----------------------------------------------------

def collect_unique_mutations(
    panels_dir: Path = config.PANELS_DIR,
) -> list[tuple[str, int, str, str]]:
    """Return the union of unique mutations across all panels.

    Each item is ``(mutation_name, position, wildtype_aa, mutant_aa)``, sorted
    by position then mutant residue.
    """
    from services.mutation_panel import load_panel

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
    panels_dir: Path = config.PANELS_DIR,
    limit: int = None,
) -> list[Path]:
    """Generate mutant PDB + PDBQT for every unique mutation across all panels.

    Skips mutations whose PDBQT already exists. ``limit`` caps the number of
    mutations processed (useful for a quick smoke test). Returns the list of
    PDBQT paths that exist after the run.
    """
    wildtype_pdb = config.STRUCTURES_DIR / "wildtype.pdb"
    if not wildtype_pdb.exists():
        raise FileNotFoundError(
            f"Wildtype receptor not found at {wildtype_pdb}. "
            f"Run structure_prep.prepare_wildtype() first."
        )

    ordered = collect_unique_mutations(panels_dir)
    if limit is not None:
        ordered = ordered[:limit]
    total = len(ordered)

    pdbqt_paths: list[Path] = []
    for i, (name, pos, wt, mut) in enumerate(ordered, start=1):
        pdbqt_path = config.MUTANTS_DIR / f"{name}.pdbqt"
        if pdbqt_path.exists():
            print(f"Skipping {name} (cached) ({i}/{total})")
            pdbqt_paths.append(pdbqt_path)
            continue
        try:
            mutant_pdb = generate_mutant(wildtype_pdb, "A", pos, wt, mut)
            pdbqt_paths.append(prepare_receptor_pdbqt(mutant_pdb))
            print(f"Generated {name} ({i}/{total})")
        except Exception as exc:  # noqa: BLE001 - keep the batch alive
            print(f"FAILED {name} ({i}/{total}): {type(exc).__name__}: {exc}")
    return pdbqt_paths
