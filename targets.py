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

import json
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

    @property
    def ground_truth_path(self) -> Path:
        """Curated/literature-derived resistance-mechanism reference for this target.

        Target-scoped so RT mechanisms never collide with the committed PI file:
        PI -> ``data/mechanism_ground_truth.json``; RT -> ``data/rt/...``.
        """
        base = DATA_DIR / self.subdir if self.subdir else DATA_DIR
        return base / "mechanism_ground_truth.json"

    @property
    def is_user(self) -> bool:
        """True for a bring-your-own target loaded from ``data/user_targets/``."""
        return self.subdir.startswith("user/")


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
# STATUS OF THE TWO VALUES THAT NEEDED CONFIRMATION (both now resolved):
#   1. docking_center — CONFIRMED. Set to the occupancy-weighted centroid of the
#      chain-A NVP (nevirapine) atoms in 3V81 = (41.105, 52.332, 49.098), the
#      NNRTI pocket of the A/B heterodimer. structure_prep.clean_structure() will
#      re-print this centroid during step 03 as a cross-check.
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

# 3V81: HIV-1 RT with an NNRTI (nevirapine, het NVP) bound in the allosteric
# pocket. CONFIRMED from the deposited coordinates: the ligand present in the
# NNRTI pocket is NEVIRAPINE (het NVP), not rilpivirine — the "TMC278/RIL" in the
# paper title refers to a related structure in the same study, not to 3V81's
# coordinates. The asymmetric unit holds two heterodimers (A/B and C/D) plus DNA
# (T/P, E/F); we keep the A/B copy: chain A = p66 (555 res, carries the pocket),
# chain B = p51 (412 res, structural). Only p66 (chain A) is mutated.
HIV1_RT = Target(
    name="HIV1_RT",
    label="HIV-1 reverse transcriptase (NNRTI pocket)",
    subdir="rt",  # data/rt/structures, data/rt/panels, ... (never touches PI data)
    pdb_id="3V81",
    pdb_url="https://files.rcsb.org/download/3V81.pdb",
    ligand_hetcodes=("NVP",),      # nevirapine in the NNRTI pocket of 3V81 (CONFIRMED)
    chains=("A", "B"),             # p66 (A) + p51 (B) heterodimer
    mutate_chains=("A",),          # RT DRMs act on the p66 subunit only
    reference_seq=_RT_REFERENCE_SEQ,
    n_positions=240,
    protonation=None,              # no catalytic-dyad protonation fix-up for RT
    docking_center=(41.105, 52.332, 49.098),  # occ-weighted centroid of chain-A NVP
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


# =============================================================================
# Bring-your-own (user) targets — triage-only, loaded from data/user_targets/
# =============================================================================
# A user target is assembled in the app (upload a PDB → confirm the pocket →
# a Claude agent proposes the drugs + resistance mutations) and persisted as a
# draft JSON. It becomes a first-class Target here (so the selector, structure
# viewer, and triage path all work) but with no genotype dataset — hence no
# reference sequence / panels / validation. Its "panel" is the agent's mutation
# list (see the API's resistance-panel fallback).

USER_TARGETS_DIR = DATA_DIR / "user_targets"


def _resolve_draft(d: dict) -> dict:
    """Normalise a draft (scripts/10 nested form OR a flat confirmed form)."""
    struct = d.get("structure", {})
    sugg = struct.get("suggested", {})
    bio = d.get("biology", {})
    drugs = d.get("drugs") or bio.get("drugs", [])
    muts = d.get("mutations") or [m["mutation"] for m in bio.get("resistance_mutations", [])]
    pdb_id = (d.get("pdb_id") or Path(d.get("source_pdb", d["target_id"])).stem).upper()
    center = d.get("docking_center") or sugg.get("docking_center") or (0.0, 0.0, 0.0)
    return {
        "target_id": d["target_id"].upper(),
        "label": d.get("label", d["target_id"]),
        "pdb_id": pdb_id,
        "ligand_hetcode": d.get("ligand_hetcode") or sugg.get("ligand_hetcode"),
        "docking_center": tuple(float(x) for x in center),
        "chains": tuple(d.get("chains") or sugg.get("chains") or ("A",)),
        "mutate_chains": tuple(d.get("mutate_chains") or sugg.get("mutate_chains") or ("A",)),
        "drugs": {x["abbrev"]: x["name"] for x in drugs},
        "mutations": [m for m in muts if m],
    }


def target_from_draft(d: dict) -> Target:
    """Build a triage-only :class:`Target` from a user-target draft dict."""
    r = _resolve_draft(d)
    het = r["ligand_hetcode"]
    return Target(
        name=r["target_id"], label=r["label"], subdir=f"user/{r['target_id']}",
        pdb_id=r["pdb_id"], pdb_url="",
        ligand_hetcodes=((het,) if het else ()),
        chains=r["chains"], mutate_chains=r["mutate_chains"],
        reference_seq="", n_positions=0, protonation=None,
        docking_center=r["docking_center"], docking_box_size=(24, 24, 24),
        dataset_filename="", dataset_urls=(), drug_columns=tuple(r["drugs"].keys()),
        drugs=r["drugs"], pubchem_cids={}, primary_mutations=frozenset(r["mutations"]),
    )


def load_user_targets() -> list[str]:
    """Load every draft in ``data/user_targets/`` into the registry; return names."""
    added = []
    if not USER_TARGETS_DIR.exists():
        return added
    for f in sorted(USER_TARGETS_DIR.glob("*.json")):
        try:
            t = target_from_draft(json.loads(f.read_text()))
        except Exception:  # noqa: BLE001 - a malformed draft must not break startup
            continue
        TARGETS[t.name] = t
        added.append(t.name)
    return added


load_user_targets()


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
