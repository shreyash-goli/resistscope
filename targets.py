"""Docking targets for ResistScope.

A :class:`Target` bundles everything specific to one drug-resistance system —
receptor structure, docking box, resistance dataset, drug panel, and reference
sequence — so the pipeline can be pointed at HIV-1 protease *or* reverse
transcriptase without touching the docking / scoring / explanation code.

``config.py`` re-exports one active target's values under the legacy
module-level names (``HIV1_PR`` by default), so all existing code and committed
data keep working unchanged. RT is purely additive: its artifacts live under
``data/rt/`` (``subdir="rt"``), so building it never touches the PI data.

The split is: **target-specific** data lives here (structure, box, dataset,
drugs, mutations); **method** constants (Vina effort, ΔΔG thresholds, Claude
model) stay in ``config.py`` because they are shared across targets.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

PROJECT_ROOT = Path(__file__).parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"  # shared across targets (raw downloads)


@dataclass(frozen=True)
class Protonation:
    """Asymmetric catalytic-residue protonation (the HIV-1 protease Asp25 dyad).

    HIV-1 protease requires one Asp25 protonated (neutral, ``ASH``) and one
    deprotonated (charged, ``ASP``). Targets without such a requirement (e.g.
    reverse transcriptase) set ``protonation=None`` and skip the fix-up.
    """

    resnum: int
    protonated_chain: str
    deprotonated_chain: str
    resname_protonated: str = "ASH"
    resname_deprotonated: str = "ASP"


@dataclass(frozen=True)
class Target:
    """Everything the pipeline needs to know about one resistance system."""

    name: str          # stable id, e.g. "HIV1_PR"
    label: str         # human label, e.g. "HIV-1 protease"
    subdir: str        # "" -> legacy flat data dirs; else data/<subdir>/...

    # --- receptor structure ---
    pdb_id: str
    pdb_url: str
    ligand_hetcodes: tuple      # het codes stripped before docking (box center)
    chains: tuple               # chains kept in the receptor
    mutate_chains: tuple        # chains a point mutation is applied to
    reference_seq: str          # consensus-B AA sequence over the genotyped region
    n_positions: int            # genotype columns P1..P{n_positions}
    protonation: Optional[Protonation]

    # --- docking box (Angstroms) ---
    docking_center: tuple
    docking_box_size: tuple

    # --- resistance dataset (Stanford HIVdb genotype-phenotype) ---
    dataset_filename: str
    dataset_urls: tuple         # tried in order
    drug_columns: tuple         # fold-resistance columns present in the TSV

    # --- drug panel ---
    drugs: dict                 # abbrev -> full name (dropdown / demo)
    pubchem_cids: dict          # abbrev -> PubChem CID (SMILES lookup)
    primary_mutations: frozenset

    # --- per-target data directories (PR keeps the legacy flat layout) ---
    def _dir(self, name: str) -> Path:
        base = DATA_DIR / self.subdir if self.subdir else DATA_DIR
        return base / name

    @property
    def structures_dir(self) -> Path:
        return self._dir("structures")

    @property
    def mutants_dir(self) -> Path:
        return self.structures_dir / "mutants"

    @property
    def panels_dir(self) -> Path:
        return self._dir("panels")

    @property
    def docking_dir(self) -> Path:
        return self._dir("docking_results")

    @property
    def explanations_dir(self) -> Path:
        return self._dir("explanations")

    @property
    def validation_dir(self) -> Path:
        return self._dir("validation")


# =============================================================================
# HIV-1 protease (the original, validated target)
# =============================================================================

# HIV-1 protease consensus-B sequence, residues 1-99 (1-indexed).
_PR_REFERENCE_SEQ = (
    "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQ"
    "ILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
)

HIV1_PR = Target(
    name="HIV1_PR",
    label="HIV-1 protease",
    subdir="",  # legacy flat dirs: data/structures, data/panels, ... (committed)
    pdb_id="3OXC",
    pdb_url="https://files.rcsb.org/download/3OXC.pdb",
    ligand_hetcodes=("ROC",),      # saquinavir in 3OXC (not SQV / 938)
    chains=("A", "B"),             # C2-symmetric homodimer
    mutate_chains=("A", "B"),      # DRMs applied to both monomers
    reference_seq=_PR_REFERENCE_SEQ,
    n_positions=99,
    protonation=Protonation(resnum=25, protonated_chain="A", deprotonated_chain="B"),
    docking_center=(5.341, -1.893, 14.179),   # occupancy-weighted ROC centroid
    docking_box_size=(22, 22, 22),
    dataset_filename="PI_DataSet.txt",
    dataset_urls=(
        "https://hivdb.stanford.edu/download/GenoPhenoDatasets/PI_DataSet.txt",
        "https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/PI_DataSet.txt",
        "https://hivdb.stanford.edu/download/GenoPhenoDatasets/PI_DataSet.Full.txt",
    ),
    drug_columns=("FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"),
    drugs={
        "ATV": "atazanavir",
        "DRV": "darunavir",
        "LPV": "lopinavir",
        "SQV": "saquinavir",
        "IDV": "indinavir",
        "NFV": "nelfinavir",
        "RTV": "ritonavir",
    },
    pubchem_cids={
        "ATV": 148192,
        "DRV": 213039,
        "LPV": 92727,
        "SQV": 441243,
        "IDV": 5362440,
        "NFV": 64143,
        "RTV": 392622,
    },
    # IAS-USA 2019/2022 major PI mutations + Stanford HIVdb major PI list.
    primary_mutations=frozenset({
        "L23I", "L24I", "D30N", "V32I", "L33F",
        "M46I", "M46L", "I47V", "I47A", "G48V", "G48M",
        "I50L", "I50V", "I54V", "I54L", "I54M", "I54A", "I54T", "I54S",
        "L76V",
        "V82A", "V82T", "V82F", "V82S", "V82L", "V82M", "V82C",
        "I84V", "I84A", "I84C", "N88S", "N88D", "L90M",
    }),
)


# =============================================================================
# HIV-1 reverse transcriptase (NNRTI-focused extension)
# =============================================================================
#
# SCIENTIFIC SCOPE — read before trusting RT scores:
#   The docking-ΔΔG method models loss of *binding affinity*. That is a sound
#   proxy for NNRTIs, which bind an allosteric hydrophobic pocket (~10 A from
#   the polymerase active site) competitively, exactly like PIs. It is NOT a
#   sound proxy for NRTIs: NRTI resistance is about nucleotide *incorporation*
#   and ATP-mediated *excision* (e.g. TAMs), not pocket affinity — rigid
#   docking of a nucleoside analog will not capture it. So this target is
#   scoped to the NNRTI pocket and the NNRTI drugs only.
#
# TWO VALUES TO CONFIRM BEFORE A PRODUCTION RUN (both self-verify, see below):
#   1. docking_center — placeholder here. structure_prep.clean_structure() prints
#      the occupancy-weighted centroid of `ligand_hetcodes` for THIS PDB; set the
#      printed value here (same workflow the PR box used). Until then, RT docking
#      is not physically meaningful.
#   2. reference_seq — consensus-B RT over positions 1-240. Anchored below with
#      assertions on canonical DRM residues so a transcription error fails loudly
#      at import; still cross-check against Stanford HIVdb consensus B before a
#      real run.

# HIV-1 RT consensus-B sequence, residues 1-240 (1-indexed). The NNRTI-pocket
# and polymerase DRMs (100,101,103,106,108,138,179,181,188,190,221,227,230,...)
# all fall inside this window.
_RT_REFERENCE_SEQ = (
    "PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPV"  # 1-60
    "FAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPL"  # 61-120
    "DEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVI"  # 121-180
    "YQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHPDKWT"  # 181-240
)

# 3V81: HIV-1 RT in complex with rilpivirine (an NNRTI) in the allosteric
# pocket. p66 (chain A) carries the pocket; p51 (chain B) is structural. We keep
# both chains but only mutate p66 (chain A), where the RT DRMs act.
# NOTE: het code + docking box must be confirmed after downloading this PDB.
HIV1_RT = Target(
    name="HIV1_RT",
    label="HIV-1 reverse transcriptase (NNRTI pocket)",
    subdir="rt",  # data/rt/structures, data/rt/panels, ... (never touches PI data)
    pdb_id="3V81",
    pdb_url="https://files.rcsb.org/download/3V81.pdb",
    ligand_hetcodes=("RIL",),      # rilpivirine (TMC278) in 3V81 — CONFIRM het code
    chains=("A", "B"),             # p66 (A) + p51 (B) heterodimer
    mutate_chains=("A",),          # RT DRMs act on the p66 subunit only
    reference_seq=_RT_REFERENCE_SEQ,
    n_positions=240,
    protonation=None,              # no catalytic-dyad protonation fix-up for RT
    docking_center=(0.0, 0.0, 0.0),           # PLACEHOLDER — set from clean_structure() print
    docking_box_size=(24, 24, 24),            # NNRTI pocket is roomy; refine after box check
    dataset_filename="NNRTI_DataSet.txt",
    dataset_urls=(
        "https://hivdb.stanford.edu/download/GenoPhenoDatasets/NNRTI_DataSet.txt",
        "https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/NNRTI_DataSet.txt",
    ),
    # NNRTI fold-resistance columns in the Stanford NNRTI dataset.
    drug_columns=("NVP", "EFV", "ETR", "RPV"),
    drugs={
        "NVP": "nevirapine",
        "EFV": "efavirenz",
        "ETR": "etravirine",
        "RPV": "rilpivirine",
        "DOR": "doravirine",
    },
    pubchem_cids={
        "NVP": 4463,
        "EFV": 64139,
        "ETR": 193962,
        "RPV": 6451164,
        "DOR": 58460047,
    },
    # Major NNRTI resistance mutations (IAS-USA 2022 / Stanford HIVdb NNRTI).
    primary_mutations=frozenset({
        "A98G", "L100I", "K101E", "K101P", "K103N", "K103S",
        "V106A", "V106M", "V108I", "E138K", "E138A", "E138G", "E138Q",
        "V179D", "V179F", "V179L", "Y181C", "Y181I", "Y181V",
        "Y188L", "Y188C", "Y188H", "G190A", "G190S", "G190E",
        "H221Y", "P225H", "F227C", "F227L", "M230L", "M230I",
    }),
)


# =============================================================================
# Registry
# =============================================================================

TARGETS: dict[str, Target] = {t.name: t for t in (HIV1_PR, HIV1_RT)}

# Aliases for CLI convenience (--target rt / pr / protease / ...).
_ALIASES = {
    "pr": "HIV1_PR", "protease": "HIV1_PR", "hiv1_pr": "HIV1_PR", "pi": "HIV1_PR",
    "rt": "HIV1_RT", "reverse_transcriptase": "HIV1_RT", "hiv1_rt": "HIV1_RT",
    "nnrti": "HIV1_RT",
}


def get_target(name: str) -> Target:
    """Resolve a target by canonical name or alias (case-insensitive)."""
    if name in TARGETS:
        return TARGETS[name]
    canonical = _ALIASES.get(name.strip().lower())
    if canonical:
        return TARGETS[canonical]
    raise KeyError(
        f"Unknown target {name!r}. Known: {sorted(TARGETS)} "
        f"(aliases: {sorted(_ALIASES)})."
    )


# --- Integrity checks: fail loudly at import if a reference sequence drifted ---
assert len(HIV1_PR.reference_seq) == HIV1_PR.n_positions, "PR reference must be 99 aa"
assert len(HIV1_RT.reference_seq) == HIV1_RT.n_positions, "RT reference must be 240 aa"

# Anchor a handful of canonical DRM wildtype residues so a transcription slip in
# the reference sequence (off-by-one, wrong letter) is caught immediately rather
# than silently mis-naming mutations downstream.
_PR_ANCHORS = {30: "D", 32: "V", 46: "M", 50: "I", 54: "I", 82: "V", 84: "I", 90: "L"}
for _pos, _aa in _PR_ANCHORS.items():
    assert HIV1_PR.reference_seq[_pos - 1] == _aa, (
        f"PR reference position {_pos} should be {_aa}, got "
        f"{HIV1_PR.reference_seq[_pos - 1]}"
    )

_RT_ANCHORS = {
    41: "M", 65: "K", 67: "D", 70: "K", 100: "L", 103: "K", 106: "V",
    108: "V", 138: "E", 179: "V", 181: "Y", 184: "M", 188: "Y", 190: "G",
    215: "T", 219: "K", 225: "P", 230: "M",
}
for _pos, _aa in _RT_ANCHORS.items():
    assert HIV1_RT.reference_seq[_pos - 1] == _aa, (
        f"RT reference position {_pos} should be {_aa}, got "
        f"{HIV1_RT.reference_seq[_pos - 1]} — verify against Stanford consensus B"
    )
