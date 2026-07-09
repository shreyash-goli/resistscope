from pathlib import Path

# === Paths ===
PROJECT_ROOT = Path(__file__).parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
STRUCTURES_DIR = DATA_DIR / "structures"
MUTANTS_DIR = STRUCTURES_DIR / "mutants"
PANELS_DIR = DATA_DIR / "panels"
DOCKING_DIR = DATA_DIR / "docking_results"
EXPLANATIONS_DIR = DATA_DIR / "explanations"
VALIDATION_DIR = DATA_DIR / "validation"

# === PDB ===
# 3OXC: wildtype HIV-1 protease complexed with saquinavir, 1.16 Angstrom
# We strip the ligand. Using a ligand-bound structure (then stripping) gives
# a more realistic closed-flap pocket conformation than apo structures like
# 2PC0 where flaps are open and not representative of drug binding.
WILDTYPE_PDB_ID = "3OXC"
WILDTYPE_PDB_URL = "https://files.rcsb.org/download/3OXC.pdb"
LIGAND_HETCODES_TO_STRIP = ["ROC"]   # Saquinavir in 3OXC (not SQV or 938)
CHAINS_TO_KEEP = ["A", "B"]  # HIV-1 protease is a homodimer

# === Docking box ===
# Centered on the catalytic Asp25/Asp25' dyad region.
# These coordinates are for cleaned 3OXC after ligand stripping.
# IMPORTANT: Verify these after structure prep by checking that the box
# encompasses the active site cavity. Use the co-crystallized saquinavir
# position as a reference before stripping it.
DOCKING_CENTER = (5.341, -1.893, 14.179)  # Occupancy-weighted centroid of ROC
                                           # altLoc A: (5.100, -3.435, 15.307)
                                           # altLoc B: (5.592, -0.287, 13.006)
DOCKING_BOX_SIZE = (22, 22, 22)      # Angstroms, generous for PI binding
VINA_EXHAUSTIVENESS = 16             # Balance speed vs accuracy; 32 for final
VINA_NUM_POSES = 5
VINA_NUM_CPUS = 0                    # 0 = use all available

#Also worth noting for the clean_structure() function: because saquinavir has two altLoc conformations (A and B), the PDB parser may return duplicate atom records. When you're computing the ligand centroid programmatically (before stripping) to confirm these coordinates match, use BioPython's is_disordered() check and take only altLoc A atoms, or average across both weighted by occupancy. The occupancy-weighted number you just computed manually (5.341, -1.893, 14.179) is the ground truth — the script just needs to reproduce it consistently.


# === Protonation ===
# HIV protease catalytic mechanism requires one Asp25 protonated, one not.
# Chain A Asp25: protonated (neutral, ASH)
# Chain B Asp25: deprotonated (charged, ASP)
# PDBFixer addMissingHydrogens at pH 7.0 will deprotonate both by default.
# We must manually fix this after hydrogen addition.
ASP25_PROTONATED_CHAIN = "A"
ASP25_DEPROTONATED_CHAIN = "B"

# === Dataset ===
RHEE_DATASET_URL = (
    "https://hivdb.stanford.edu/pages/published_analysis/"
    "genophenoPNAS2006/"
)
# The PI dataset is a TSV. Format:
# Row = one isolate. Columns = drug fold-resistance values + binary mutation
# indicators. Response variable is log fold-resistance (continuous).
# Mutation columns are binary (0/1) indicating presence of amino acid change.
PI_DRUGS = {
    "ATV": "atazanavir",
    "DRV": "darunavir",
    "LPV": "lopinavir",
    "SQV": "saquinavir",
    "IDV": "indinavir",
    "NFV": "nelfinavir",
    "RTV": "ritonavir",
}

# PubChem CIDs for SMILES lookup
PUBCHEM_CIDS = {
    "ATV": 148192,
    "DRV": 213039,
    "LPV": 92727,
    "SQV": 441243,
    "IDV": 5362440,
    "NFV": 64143,
    "RTV": 392622,
}

# === Claude API ===
# Opus 4.8 is Anthropic's most capable current model — the right default for
# structural-mechanism reasoning. Swap to "claude-sonnet-5" to cut cost/latency
# if generating many explanations. (The original "claude-sonnet-4-6" is also a
# valid current model but a tier below Opus for scientific reasoning.)
CLAUDE_MODEL = "claude-opus-4-8"
CLAUDE_MAX_TOKENS = 1024

# === Scoring ===
# Delta-delta-G thresholds (kcal/mol) for flagging concerning mutations
DDG_WARNING_THRESHOLD = 0.5    # Mild concern
DDG_DANGER_THRESHOLD = 1.5     # Likely resistance

# === Primary (major) PI resistance mutations ===
# Curated set of MAJOR protease-inhibitor resistance mutations, used to flag
# panel rows with is_primary=True and to drive the "count known DRMs" baseline.
# Source: IAS-USA 2019/2022 "Update of the Drug Resistance Mutations in HIV-1"
# (Wensing et al.) major PI mutations + Stanford HIVdb major PI list. Each entry
# is a specific wildtype-position-mutant substitution (accessory/minor mutations
# and polymorphisms are intentionally excluded).
PRIMARY_PI_MUTATIONS = frozenset({
    "L23I",
    "L24I",
    "D30N",
    "V32I",
    "L33F",
    "M46I", "M46L",
    "I47V", "I47A",
    "G48V", "G48M",
    "I50L", "I50V",
    "I54V", "I54L", "I54M", "I54A", "I54T", "I54S",
    "L76V",
    "V82A", "V82T", "V82F", "V82S", "V82L", "V82M", "V82C",
    "I84V", "I84A", "I84C",
    "N88S", "N88D",
    "L90M",
})
