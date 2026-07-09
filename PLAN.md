# ResistScope — Comprehensive Build Document
### Built with Claude: Life Sciences Hackathon, Build Track

---

## 0. What This Is

A resistance-aware antiviral compound triage tool. A drug discovery scientist
(the named user) pastes a SMILES string for a candidate compound targeting HIV-1
protease, and the tool:

1. Docks it against wildtype and a panel of ~25 known clinical resistance mutants
2. Computes a robustness score (how much binding degrades across the panel)
3. Generates per-mutation mechanistic explanations via Claude
4. Validates the whole pipeline against real measured fold-resistance data

**MVP scope:** HIV-1 protease inhibitors only. RT is a post-hackathon extension.

**Why Gladstone cares:** Gladstone's Infectious Disease Institute, HOPE
Collaboratory, and Center for HIV Cure Research are all actively working on HIV.
Drug resistance is a central problem in their research. This tool addresses it
directly.

---

## 1. Repository Structure

```
resistscope/
├── config.py                     # All constants, paths, parameters
├── requirements.txt
├── environment.yml               # Conda env (preferred for rdkit/vina/openmm)
├── README.md
│
├── data/
│   ├── raw/                      # Downloaded datasets (gitignored)
│   │   └── PI_DataSet.txt        # Rhee 2006 PI genotype-phenotype data
│   ├── structures/
│   │   ├── wildtype.pdb          # Cleaned 3OXC chain A+B, ligand stripped
│   │   ├── wildtype.pdbqt        # Receptor prepared for Vina
│   │   └── mutants/              # One PDB + PDBQT per mutation
│   │       ├── V82A.pdb
│   │       ├── V82A.pdbqt
│   │       └── ...
│   ├── panels/                   # Precomputed mutation panels per drug
│   │   ├── ATV.parquet
│   │   ├── DRV.parquet
│   │   └── ...
│   ├── ligands/                  # SMILES + 3D conformers for benchmark drugs
│   │   └── benchmark_drugs.json
│   ├── docking_results/
│   │   └── benchmark_docking.parquet
│   ├── explanations/             # Cached Claude explanations
│   │   └── {drug}_{mutation}.json
│   ├── mechanism_ground_truth.json
│   └── validation/
│       ├── scores_vs_fold_resistance.parquet
│       └── validation_plot.png
│
├── services/
│   ├── __init__.py
│   ├── mutation_panel.py         # Parse Rhee dataset, build per-drug panels
│   ├── structure_prep.py         # PDB download, clean, mutate, PDBQT convert
│   ├── docking.py                # SMILES→3D→PDBQT→Vina→score
│   ├── scoring.py                # Delta-delta-G, aggregation, baselines
│   ├── explanation.py            # Claude API for mechanistic rationales
│   └── validation.py             # Predicted vs real fold-resistance
│
├── scripts/
│   ├── 01_download_data.py       # Fetch Rhee dataset + PDB structure
│   ├── 02_build_panels.py        # Parse dataset into per-drug panels
│   ├── 03_build_mutant_cache.py  # Generate all mutant structures
│   ├── 04_dock_benchmark.py      # Dock 7 approved PIs against full panel
│   ├── 05_validate.py            # Correlation analysis + plots
│   ├── 06_generate_explanations.py
│   └── 07_faithfulness_eval.py
│
├── api/
│   ├── __init__.py
│   └── main.py                   # FastAPI app
│
├── frontend/                     # React + Vite (stretch goal)
│   └── ...
│
├── notebooks/
│   └── demo.ipynb                # Fallback demo if frontend is cut
│
├── Dockerfile
└── demo.py                       # Scripted CLI demo for video
```

---

## 2. Configuration (config.py)

Every magic number, path, and parameter lives here. Nothing hardcoded in
service modules.

```python
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
CLAUDE_MODEL = "claude-sonnet-4-6"
CLAUDE_MAX_TOKENS = 1024

# === Scoring ===
# Delta-delta-G thresholds (kcal/mol) for flagging concerning mutations
DDG_WARNING_THRESHOLD = 0.5    # Mild concern
DDG_DANGER_THRESHOLD = 1.5     # Likely resistance
```

---

## 3. Data Schemas

### 3a. Mutation Panel (per drug, parquet)

| Column | Type | Description |
|---|---|---|
| mutation | str | e.g. "V82A" |
| position | int | residue number, e.g. 82 |
| wildtype_aa | str | one-letter, e.g. "V" |
| mutant_aa | str | one-letter, e.g. "A" |
| mean_log_fold_resistance | float | mean log(fold-resistance) from Rhee data |
| n_isolates | int | number of isolates with this mutation |
| is_primary | bool | whether this is a primary/major DRM |

### 3b. Docking Results (parquet)

| Column | Type | Description |
|---|---|---|
| drug | str | drug abbreviation, e.g. "DRV" |
| mutation | str | "WT" for wildtype, else e.g. "V82A" |
| delta_g | float | best Vina binding affinity (kcal/mol, negative = better) |
| delta_delta_g | float | delta_g_mutant - delta_g_wildtype (positive = worse binding) |
| n_poses | int | number of poses returned |

### 3c. Explanation Cache (JSON per drug-mutation pair)

```json
{
  "drug": "DRV",
  "mutation": "V82A",
  "delta_delta_g": 1.2,
  "structural_context": {
    "position": 82,
    "wt_aa": "V",
    "mut_aa": "A",
    "volume_change": "decrease",
    "charge_change": "none",
    "distance_from_ligand_centroid_angstrom": 4.3,
    "contacts_ligand_directly": true,
    "region": "active_site_S1_subpocket"
  },
  "explanation": "V82A removes a methyl group from ...",
  "model": "claude-sonnet-4-6",
  "timestamp": "2026-07-11T14:30:00Z"
}
```

### 3d. Mechanism Ground Truth (manually curated JSON)

```json
{
  "V82A": {
    "mechanism": "Reduces van der Waals contact with P1/P1' groups of PIs by removing hydrophobic bulk in the S1 subpocket",
    "source": "IAS-USA 2024 drug resistance mutations figure",
    "affects_drugs": ["ATV", "LPV", "IDV", "NFV"]
  },
  "D30N": {
    "mechanism": "Disrupts hydrogen bond between Asp30 side chain and the P2 aniline NH of nelfinavir specifically",
    "source": "IAS-USA 2024; Wlodawer & Vondrasek 1998",
    "affects_drugs": ["NFV"]
  }
}
```

### 3e. Validation Output (parquet)

| Column | Type | Description |
|---|---|---|
| drug | str | drug abbreviation |
| scoring_method | str | "weighted_mean", "worst_case", "simple_mean", "mutation_count", "wt_only" |
| spearman_rho | float | Spearman correlation with real fold-resistance |
| spearman_pvalue | float | |
| pearson_r | float | Pearson correlation |
| n_mutations | int | number of mutations in comparison |

---

## 4. Function Signatures (every public function in services/)

### services/mutation_panel.py

```python
def download_rhee_dataset(output_dir: Path = RAW_DIR) -> Path:
    """Download the PI genotype-phenotype TSV from Stanford HIVdb.
    Returns path to the downloaded file.
    Raises ConnectionError if download fails."""

def parse_pi_dataset(filepath: Path) -> pd.DataFrame:
    """Parse the raw TSV into a clean DataFrame.
    Returns DataFrame with columns: isolate_id, drug, log_fold_resistance,
    and one binary column per mutation position.
    Handles: missing values (drop rows with NaN resistance for a given drug),
    duplicate columns (merge), zero-variance columns (drop)."""

def build_mutation_panel(
    df: pd.DataFrame,
    drug: str,
    min_isolates: int = 3
) -> pd.DataFrame:
    """For a given drug, compute per-mutation summary statistics.
    Only includes mutations present in >= min_isolates isolates.
    Returns DataFrame matching the Mutation Panel schema above."""

def build_all_panels(output_dir: Path = PANELS_DIR) -> dict[str, Path]:
    """Build and save panels for all 7 PI drugs. Returns {drug: path}."""

def load_panel(drug: str, panels_dir: Path = PANELS_DIR) -> pd.DataFrame:
    """Load a precomputed panel from parquet."""
```

### services/structure_prep.py

```python
def download_pdb(
    pdb_id: str = WILDTYPE_PDB_ID,
    output_dir: Path = STRUCTURES_DIR
) -> Path:
    """Download PDB file from RCSB. Returns path to downloaded file."""

def clean_structure(
    pdb_path: Path,
    chains_to_keep: list[str] = CHAINS_TO_KEEP,
    het_codes_to_strip: list[str] = LIGAND_HETCODES_TO_STRIP,
    strip_waters: bool = True
) -> Path:
    """Strip ligands, waters, non-protein atoms. Keep only specified chains.
    Save as wildtype.pdb. Returns path.

    IMPORTANT: Before stripping the ligand, record its centroid coordinates
    and save to config or a sidecar JSON. This centroid defines the docking
    box center. Print it to stdout so the developer can update DOCKING_CENTER
    in config.py."""

def add_hydrogens_and_fix_protonation(
    pdb_path: Path,
    ph: float = 7.0
) -> Path:
    """Use PDBFixer to:
    1. Add missing heavy atoms
    2. Add hydrogens at specified pH
    3. FIX Asp25 protonation: after PDBFixer adds H, manually ensure
       chain A Asp25 has the OD2 hydrogen (protonated/neutral) and
       chain B Asp25 does not (charged). This is critical for correct
       docking scores. See config.ASP25_PROTONATED_CHAIN.
    Save result. Returns path.

    Implementation note on protonation fix:
    PDBFixer will protonate both Asp25 residues at pH 7.0 (their pKa in
    the free enzyme is ~3.5, but in the dimer active site the effective
    pKa shifts to ~5-6, so one should be protonated). After
    addMissingHydrogens(), iterate over atoms in chain B residue 25,
    find the HD2 atom on the OD2 oxygen, and delete it. This gives the
    asymmetric protonation state the enzyme actually uses."""

def generate_mutant(
    wildtype_pdb: Path,
    chain_id: str,
    position: int,
    wildtype_aa: str,
    mutant_aa: str,
    output_dir: Path = MUTANTS_DIR
) -> Path:
    """Generate a point mutant using PDBFixer.

    Steps:
    1. Load wildtype PDB into PDBFixer
    2. Call fixer.applyMutations([f"{wildtype_aa}-{position}-{mutant_aa}"],
                                 chain_id)
    3. fixer.findMissingResidues(), fixer.findMissingAtoms()
    4. fixer.addMissingAtoms()
    5. Save to output_dir/{mutation_name}.pdb

    HIV protease is a homodimer. Many resistance mutations affect both
    chains. Apply the mutation to BOTH chains A and B unless the mutation
    is chain-specific (none of the major DRMs are).

    Returns path to mutant PDB."""

def prepare_receptor_pdbqt(pdb_path: Path) -> Path:
    """Convert a PDB to PDBQT for Vina using meeko's mk_prepare_receptor.

    from meeko import Polymer
    with open(pdb_path) as f:
        pdb_string = f.read()
    polymer = Polymer.from_pdb_string(pdb_string)
    # Write PDBQT
    pdbqt_path = pdb_path.with_suffix('.pdbqt')
    polymer.write_pdbqt_file(str(pdbqt_path))

    If meeko's receptor prep is problematic, fall back to:
    subprocess.run(['mk_prepare_receptor.py', '-i', str(pdb_path),
                    '-o', str(pdbqt_path.stem), '-p'])

    Returns path to PDBQT file."""

def build_mutant_cache(panels_dir: Path = PANELS_DIR) -> list[Path]:
    """Load all panels, collect unique mutations across all drugs,
    generate mutant PDB + PDBQT for each. Returns list of paths.
    Skip mutations that already exist in the cache."""
```

### services/docking.py

```python
def smiles_to_pdbqt(smiles: str) -> str:
    """Convert SMILES to PDBQT string for Vina, in memory (no files).

    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    mk_prep = MoleculePreparation()
    molsetup_list = mk_prep(mol)
    pdbqt_string = PDBQTWriterLegacy.write_string(molsetup_list[0])
    return pdbqt_string

    IMPORTANT: meeko does NOT assign protonation states. The input mol
    must have explicit hydrogens. RDKit's AddHs() handles this but does
    not consider pH. For most PI-like molecules at physiological pH this
    is fine. If a specific protonation state matters, use scrub.py from
    the molscrub package before meeko."""

def dock_single(
    ligand_pdbqt: str,
    receptor_pdbqt_path: Path,
    center: tuple[float, float, float] = DOCKING_CENTER,
    box_size: tuple[int, int, int] = DOCKING_BOX_SIZE,
    exhaustiveness: int = VINA_EXHAUSTIVENESS,
    n_poses: int = VINA_NUM_POSES,
    cpu: int = VINA_NUM_CPUS
) -> dict:
    """Dock one ligand against one receptor. Returns:
    {"delta_g": float, "n_poses": int, "all_energies": list[float]}

    from vina import Vina
    v = Vina(sf_name='vina')
    v.set_receptor(str(receptor_pdbqt_path))
    v.set_ligand_from_string(ligand_pdbqt)
    v.compute_vina_maps(center=list(center), box_size=list(box_size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

    # Parse energies from output
    # v.energies() returns numpy array, first column is total energy
    energies = v.energies()
    best_energy = float(energies[0][0])  # kcal/mol

    Return best_energy as delta_g.

    Wrap in try/except. On failure (timeout, bad geometry, no poses found),
    return {"delta_g": None, "n_poses": 0, "error": str(e)}.
    Do NOT let one failed dock crash the whole pipeline."""

def dock_against_panel(
    smiles: str,
    drug_name: str,
    panel: pd.DataFrame,
    structures_dir: Path = STRUCTURES_DIR,
    mutants_dir: Path = MUTANTS_DIR
) -> pd.DataFrame:
    """Dock a compound against wildtype + all mutations in a panel.
    Returns DataFrame matching the Docking Results schema.

    Steps:
    1. Convert SMILES to PDBQT string once (reuse for all receptors)
    2. Dock against wildtype.pdbqt → get WT delta_g
    3. For each mutation in panel:
       a. Load mutants/{mutation}.pdbqt
       b. Dock → get mutant delta_g
       c. Compute delta_delta_g = mutant_delta_g - wt_delta_g
    4. Return results DataFrame

    Parallelize with multiprocessing.Pool if >10 mutations. Each dock
    is independent. Use Pool(cpu_count() - 1) to leave one core free."""

def get_benchmark_smiles() -> dict[str, str]:
    """Fetch SMILES for all 7 benchmark PIs from PubChem PUG REST API.

    import requests
    smiles = {}
    for drug, cid in PUBCHEM_CIDS.items():
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
        resp = requests.get(url)
        data = resp.json()
        smiles[drug] = data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    return smiles

    Cache result to data/ligands/benchmark_drugs.json so we don't
    re-fetch every time."""
```

### services/scoring.py

```python
def compute_robustness_scores(
    docking_results: pd.DataFrame,
    panel: pd.DataFrame
) -> dict[str, float]:
    """Compute three robustness scores from docking results for one drug.

    Join docking_results with panel on mutation to get prevalence weights
    and ground-truth fold-resistance.

    Returns {
        "simple_mean_ddg": float,        # mean of all delta_delta_g values
        "prevalence_weighted_ddg": float, # weighted by mutation prevalence
        "worst_case_ddg": float,          # max delta_delta_g (worst mutation)
        "robustness_0_100": float,        # normalized 0-100, 100=perfectly robust
        "n_mutations_scored": int,
        "n_mutations_failed": int         # docking failures
    }

    Normalization for robustness_0_100:
    score = max(0, 100 - (prevalence_weighted_ddg * 33.3))
    This maps 0 kcal/mol → 100, 3 kcal/mol → 0. Adjust scale if needed."""

def compute_baseline_mutation_count(panel: pd.DataFrame) -> pd.DataFrame:
    """Baseline 1: simple count of known resistance mutations.
    For each drug, the 'predicted resistance' is just the number of
    primary DRMs in the panel. No docking involved.
    Returns DataFrame with drug, mutation, baseline_score columns."""

def compute_baseline_wt_only(
    docking_results: pd.DataFrame
) -> pd.DataFrame:
    """Baseline 2: wildtype docking score only (mutation-blind).
    Uses only the WT delta_g as the prediction. This tests whether
    just knowing how well something binds WT is sufficient.
    Returns DataFrame with drug, baseline_score columns."""

def compare_all_methods(
    docking_results: pd.DataFrame,
    panel: pd.DataFrame
) -> pd.DataFrame:
    """Compute Spearman and Pearson correlations between each scoring
    method and real mean_log_fold_resistance from the Rhee dataset.
    Returns DataFrame matching Validation Output schema."""
```

### services/explanation.py

```python
def build_structural_context(
    mutation: str,
    docking_result: dict,
    wildtype_pdb: Path = STRUCTURES_DIR / "wildtype.pdb"
) -> dict:
    """Extract structural context for a mutation to feed to Claude.

    Using BioPython:
    1. Parse wildtype PDB
    2. Find the mutated residue
    3. Compute distance from residue CA to the docking box center
       (proxy for distance to ligand)
    4. Look up amino acid properties (volume, charge, hydrophobicity)
       for both WT and mutant residues
    5. Determine which subpocket/region the residue is in based on
       position number (S1, S1', S2, S2', flap, etc.)

    Returns dict matching the structural_context schema above."""

def generate_explanation(
    drug: str,
    mutation: str,
    delta_delta_g: float,
    structural_context: dict,
    cache_dir: Path = EXPLANATIONS_DIR
) -> str:
    """Call Claude API to generate a mechanistic explanation.

    Check cache first: if {drug}_{mutation}.json exists, return cached.

    System prompt:
    'You are a structural biologist specializing in HIV-1 protease
    drug resistance. Given structural context about a mutation and its
    effect on drug binding, provide a 2-3 sentence mechanistic
    hypothesis for WHY this mutation causes resistance to this specific
    drug. Ground your explanation in the structural facts provided.
    Do not speculate beyond the given data. Be specific about
    molecular interactions (hydrogen bonds, van der Waals contacts,
    steric clashes, electrostatic changes).'

    User prompt:
    f'Drug: {drug_full_name}
    Mutation: {mutation}
    Delta-delta-G: {delta_delta_g:.2f} kcal/mol (positive = worse binding)
    Structural context: {json.dumps(structural_context, indent=2)}

    Explain the structural mechanism by which this mutation likely
    causes resistance to this drug.'

    Cache the result. Return the explanation text."""

def evaluate_faithfulness(
    explanations_dir: Path = EXPLANATIONS_DIR,
    ground_truth_path: Path = DATA_DIR / "mechanism_ground_truth.json"
) -> pd.DataFrame:
    """Use Claude as a judge to score explanation faithfulness.

    For each mutation in ground_truth:
    1. Load the generated explanation
    2. Load the ground truth mechanism
    3. Ask Claude to score on a 0-2 scale:
       0 = contradicts known mechanism
       1 = consistent but vague or incomplete
       2 = correctly identifies the key structural mechanism

    System prompt for judge:
    'You are evaluating whether an AI-generated explanation of HIV
    protease drug resistance correctly identifies the known structural
    mechanism. Score 0 if the explanation contradicts the known
    mechanism, 1 if it is consistent but vague or misses the key
    interaction, 2 if it correctly identifies the primary structural
    basis. Respond with only the score and a one-sentence justification.'

    Returns DataFrame with mutation, drug, score, justification columns."""
```

### services/validation.py

```python
def run_full_validation(
    docking_dir: Path = DOCKING_DIR,
    panels_dir: Path = PANELS_DIR,
    output_dir: Path = VALIDATION_DIR
) -> dict:
    """Run the complete validation pipeline.

    1. Load all docking results and panels
    2. For each drug, compute all scoring methods + baselines
    3. Compute correlations with real fold-resistance
    4. Generate scatter plot (predicted vs real)
    5. Save everything

    Returns summary dict with headline numbers."""

def plot_validation(
    results: pd.DataFrame,
    output_path: Path = VALIDATION_DIR / "validation_plot.png"
) -> Path:
    """Generate a multi-panel figure:
    - One scatter plot per scoring method
    - X axis: real mean log fold-resistance (from Rhee)
    - Y axis: predicted delta-delta-G
    - Annotate with Spearman rho and p-value
    - Color points by drug

    Use matplotlib, not plotly. Keep it publication-quality but simple.
    Returns path to saved figure."""
```

---

## 5. Key Technical Gotchas

### 5a. PDBQT Conversion

RDKit does NOT output PDBQT. You need **meeko** (pip install meeko) for both
ligand and receptor preparation. This is the single most common point of
confusion in Vina tutorials.

Ligand: `MoleculePreparation` + `PDBQTWriterLegacy` (in-memory, no temp files)
Receptor: `mk_prepare_receptor.py` CLI or `Polymer.from_pdb_string()` Python API

### 5b. Asp25/Asp25' Protonation

This is the #1 source of bad docking scores for HIV protease. The catalytic
dyad has an asymmetric protonation state in the active enzyme. PDBFixer at
pH 7.0 will deprotonate both (both charged), which is wrong. After adding
hydrogens, you must manually delete the HD2 atom from one of the two Asp25
residues. See the function signature for add_hydrogens_and_fix_protonation.

### 5c. Crystal Waters

Conserved "flap water" mediates binding in HIV protease. For the hackathon,
strip all waters and accept the noise. Document this as a known limitation.
Adding explicit waters is a post-hackathon improvement.

### 5d. Docking Box Placement

The docking box center MUST be determined empirically from the structure, not
guessed. The clean_structure function should print the ligand centroid
coordinates BEFORE stripping the ligand. Use those coordinates as
DOCKING_CENTER in config.py. Verify by checking that the box fully encloses
the active site cavity (Asp25 dyad + S1/S1'/S2/S2' subpockets).

### 5e. PDB 3OXC Is NOT Apo

3OXC contains saquinavir (het code 938). This is actually desirable: the
protein is in a drug-bound (closed-flap) conformation, which is more
appropriate for docking PIs than an apo open-flap structure like 2PC0.
Strip the ligand but use the protein coordinates as-is.

### 5f. Homodimer Mutations

HIV protease is a C2-symmetric homodimer. Most resistance mutations (V82A,
I84V, etc.) occur in both chains. When generating mutants, apply the mutation
to BOTH chain A and chain B.

### 5g. PDBFixer vs PyMOL for Mutagenesis

Use PDBFixer (from OpenMM) as the primary tool. It handles point mutations
programmatically with `applyMutations()`, adds missing atoms, and runs
energy minimization to resolve clashes. It requires conda install
(pdbfixer + openmm). PyMOL's mutagenesis wizard is finicky when scripted
headlessly. Only fall back to PyMOL if PDBFixer doesn't handle a specific
mutation well.

PDBFixer mutation API:
```python
from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename=str(pdb_path))
fixer.applyMutations(["VAL-82-ALA"], "A")  # Chain A
fixer.applyMutations(["VAL-82-ALA"], "B")  # Chain B
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)
# Fix Asp25 protonation here (see 5b)
with open(output_path, 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)
```

### 5h. Rhee Dataset Format Quirks

The PI dataset from Stanford HIVdb is tab-separated. Column headers may
include special characters. The resistance values are log10(fold-resistance).
Some isolates have missing values for certain drugs (encoded as NA or empty).
Drop these per-drug, not globally. Mutation columns are binary (0/1) with
names like "10F" meaning position 10 mutated to Phe.

### 5i. Vina Python API Pattern

```python
from vina import Vina

v = Vina(sf_name='vina')
v.set_receptor(str(receptor_pdbqt_path))
v.set_ligand_from_string(ligand_pdbqt_string)  # NOT from_file
v.compute_vina_maps(
    center=[15.0, 0.5, -2.0],
    box_size=[22, 22, 22]
)
v.dock(exhaustiveness=16, n_poses=5)
energies = v.energies()  # numpy array, shape (n_poses, 4+)
best_affinity = float(energies[0][0])  # kcal/mol, negative = good
```

---

## 6. Environment Setup

```yaml
# environment.yml
name: resistscope
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - rdkit
  - openmm
  - pdbfixer
  - vina
  - numpy
  - scipy
  - pandas
  - matplotlib
  - biopython
  - pip
  - pip:
    - meeko
    - fastapi
    - uvicorn
    - requests
    - anthropic
    - pyarrow
```

Install: `conda env create -f environment.yml && conda activate resistscope`

Test that everything works:
```bash
python -c "from rdkit import Chem; print('RDKit OK')"
python -c "from vina import Vina; print('Vina OK')"
python -c "from pdbfixer import PDBFixer; print('PDBFixer OK')"
python -c "from meeko import MoleculePreparation; print('Meeko OK')"
```

---

## 7. Timeline (Jul 7-12)

### Jul 7 (today, remaining hours)
- [ ] Create conda env, verify all imports
- [ ] Create repo with full directory structure
- [ ] Run `scripts/01_download_data.py`: download Rhee PI dataset + 3OXC PDB
- [ ] Eyeball the Rhee TSV: check headers, delimiters, NA encoding
- [ ] Eyeball 3OXC in any PDB viewer: identify saquinavir, note its centroid

### Jul 8 — Data + Structures
- [ ] `scripts/02_build_panels.py`: parse Rhee into per-drug mutation panels
- [ ] Sanity check: how many mutations per drug? Are counts reasonable?
- [ ] `services/structure_prep.py`: clean 3OXC, strip saquinavir, record
      centroid, update DOCKING_CENTER in config.py
- [ ] Add hydrogens with Asp25 protonation fix
- [ ] Generate wildtype PDBQT
- [ ] `scripts/03_build_mutant_cache.py`: generate all mutant PDB + PDBQT files
- [ ] Spot-check 2-3 mutants by superimposing on wildtype

### Jul 9 — Docking Pipeline (highest risk day)
- [ ] `services/docking.py`: get SMILES→PDBQT working for one drug (darunavir)
- [ ] Dock darunavir against wildtype. Sanity check: known binding affinity
      for DRV is roughly -12 to -14 kcal/mol. If you get -5 or -20, something
      is wrong (box placement, protonation, or PDBQT conversion).
- [ ] Dock darunavir against 3 mutants. Check that delta_delta_g has the
      right sign (positive = worse binding for resistant mutants)
- [ ] If sanity checks pass: scale to full panel for darunavir
- [ ] If sanity checks fail: debug box placement, reprotonation, receptor prep

### Jul 10 — Scale + Validate (load-bearing day)
- [ ] `scripts/04_dock_benchmark.py`: dock all 7 drugs × full panel
      (expect ~210 runs, maybe 1-3 hours depending on CPU and exhaustiveness)
- [ ] `services/scoring.py`: compute all scoring methods + baselines
- [ ] `scripts/05_validate.py`: correlation analysis
- [ ] **CHECKPOINT: look at the Spearman rho.** If ≥0.3, you have a result.
      If <0.2, check for systematic errors before proceeding. If negative,
      something is fundamentally wrong.

### Jul 11 — Explanations + Interface
- [ ] Build `data/mechanism_ground_truth.json` by hand (~15 mutations,
      ~2 hours of literature work with IAS-USA figures)
- [ ] `services/explanation.py`: generate explanations for all drug-mutation
      pairs with |delta_delta_g| > 0.5
- [ ] `scripts/07_faithfulness_eval.py`: run the judge
- [ ] Build the interface: either React frontend OR Jupyter notebook demo
      (decide based on how much time is left)
- [ ] Wire up `api/main.py` if doing React

### Jul 12 — Polish + Submit
- [ ] `demo.py` or `demo.ipynb`: scripted end-to-end walkthrough on darunavir
- [ ] README with setup, architecture diagram, headline validation numbers
- [ ] Record demo video if required
- [ ] Submit

**If behind schedule, cut in this order:**
1. Frontend → deliver as CLI + notebook
2. Faithfulness eval → paper extension only
3. Claude explanations → paper extension only
4. NEVER cut: docking pipeline + validation (this IS the project)

---

## 8. Claude Code Prompts

Paste these sequentially. Each is self-contained and testable before moving on.
They reference this document's schemas and patterns.

---

### Prompt 1 — Project Scaffold

```
Create a Python project called "resistscope" with the following exact structure:

resistscope/
├── config.py
├── requirements.txt
├── environment.yml
├── README.md
├── data/raw/ data/structures/ data/structures/mutants/ data/panels/
│   data/docking_results/ data/explanations/ data/validation/
├── services/__init__.py services/mutation_panel.py services/structure_prep.py
│   services/docking.py services/scoring.py services/explanation.py
│   services/validation.py
├── scripts/01_download_data.py through 07_faithfulness_eval.py (empty stubs)
├── api/__init__.py api/main.py
├── notebooks/demo.ipynb (empty)
└── demo.py (empty)

For config.py, use this exact content:
[paste the config.py section from section 2 of this document]

For environment.yml, use this exact content:
[paste the environment.yml from section 6]

Write a README.md stub with: project name "ResistScope", one-paragraph
description ("Resistance-aware triage tool for antiviral compounds targeting
HIV-1 protease"), setup instructions pointing to environment.yml, and a
"Quick Start" section that says "Coming soon."

Do not install anything yet, just create the files.
```

---

### Prompt 2 — Download Script

```
Write scripts/01_download_data.py for the resistscope project.

This script does two things:

1. Download the HIV-1 protease genotype-phenotype dataset from Stanford HIVdb.
   URL: https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/
   This page hosts several TSV files. We need the PI (protease inhibitor)
   dataset. Download it to data/raw/PI_DataSet.txt.
   
   The page may serve the files directly or via links. Try fetching the page
   first, parse for .txt or .tsv links containing "PI", and download. If that
   fails, try common filenames like PI_DataSet.txt or PI.txt.

2. Download PDB structure 3OXC from RCSB.
   URL: https://files.rcsb.org/download/3OXC.pdb
   Save to data/raw/3OXC.pdb

Use requests library. Print progress. Verify files are non-empty after download.
Import paths from config.py.
```

---

### Prompt 3 — Mutation Panel Builder

```
Write services/mutation_panel.py and scripts/02_build_panels.py for resistscope.

The Rhee PI dataset (data/raw/PI_DataSet.txt) is a tab-separated file where:
- Each row is one viral isolate
- There are columns for drug resistance (log fold-resistance) for each of the
  7 protease inhibitors: ATV, DRV, LPV, SQV, IDV, NFV, RTV (though column
  names may use full drug names or other abbreviations)
- There are binary (0/1) columns for each mutation, named like "10F" meaning
  position 10, mutated amino acid F (phenylalanine)
- Some resistance values may be missing (NA or empty)

Implement these functions in services/mutation_panel.py with these exact
signatures and behaviors:

def parse_pi_dataset(filepath: Path) -> pd.DataFrame:
    """Parse raw TSV. Handle NA values. Return clean DataFrame."""

def build_mutation_panel(df, drug, min_isolates=3) -> pd.DataFrame:
    """For one drug, compute per-mutation stats. Return DataFrame with columns:
    mutation (str like "V82A"), position (int), wildtype_aa (str),
    mutant_aa (str), mean_log_fold_resistance (float), n_isolates (int)"""

def build_all_panels(output_dir) -> dict:
    """Build panels for all 7 drugs, save as parquet."""

def load_panel(drug, panels_dir) -> pd.DataFrame:
    """Load a panel from parquet."""

The mutation column names in the dataset use a compact format. Position 82
mutated to Alanine would be column "82A" or similar. The wildtype amino acid
at each position can be inferred from the HIV-1 protease reference sequence:
PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIG
TVLVGPTPVNIIGRNLLTQIGCTLNF

Use this reference sequence to determine what the wildtype AA is at each
position. Then a column named "82A" means wildtype V (position 82 in the
reference is V) mutated to A, so mutation = "V82A".

scripts/02_build_panels.py should call build_all_panels() and print a summary:
drug name, number of mutations, number of isolates.

Use config.PI_DRUGS for drug name mappings.
```

---

### Prompt 4 — Structure Preparation

```
Write services/structure_prep.py and scripts/03_build_mutant_cache.py.

This module prepares the HIV-1 protease receptor for docking.

Dependencies: pdbfixer (from openmm), biopython, meeko

Key functions:

def clean_structure(pdb_path, chains_to_keep=["A","B"],
                    het_codes_to_strip=["938"], strip_waters=True) -> Path:
    BEFORE stripping the ligand (het code 938 = saquinavir), compute and
    PRINT the centroid (mean x,y,z) of the ligand atoms. This will be used
    as DOCKING_CENTER. Then strip the ligand, any other heterogens, and
    waters. Keep only chains A and B. Save as data/structures/wildtype.pdb.

def add_hydrogens_and_fix_protonation(pdb_path, ph=7.0) -> Path:
    Use PDBFixer:
    1. fixer = PDBFixer(filename=str(pdb_path))
    2. fixer.findMissingResidues()
    3. fixer.findMissingAtoms()
    4. fixer.addMissingAtoms()
    5. fixer.addMissingHydrogens(ph)
    6. CRITICAL: Fix Asp25 protonation. After step 5, iterate through
       chain B residue 25 atoms. Find the hydrogen on OD2 (atom name HD2).
       Delete it to leave chain B Asp25 deprotonated (charged) while chain A
       Asp25 remains protonated (neutral). This asymmetric protonation is
       required for correct docking scores.
    7. Save with PDBFile.writeFile()

def generate_mutant(wildtype_pdb, chain_id, position, wildtype_aa,
                    mutant_aa, output_dir) -> Path:
    Use PDBFixer.applyMutations(). Apply the mutation to BOTH chains A and B
    (HIV protease is a homodimer). Then findMissingAtoms, addMissingAtoms,
    addMissingHydrogens, fix Asp25 protonation, save.

def prepare_receptor_pdbqt(pdb_path) -> Path:
    Use meeko's mk_prepare_receptor.py via subprocess:
    subprocess.run(['mk_prepare_receptor.py', '-i', str(pdb_path),
                    '-o', str(pdb_path.stem), '-p'], check=True)
    The -p flag writes a .pdbqt file. Return the PDBQT path.
    
    If mk_prepare_receptor.py is not on PATH, try:
    python -m meeko.scripts.mk_prepare_receptor

scripts/03_build_mutant_cache.py:
1. Load all mutation panels from data/panels/
2. Collect the union of unique mutations across all drugs
3. For each unique mutation, call generate_mutant() then prepare_receptor_pdbqt()
4. Skip if the PDBQT already exists in data/structures/mutants/
5. Print progress: "Generated V82A (3/27)"
```

---

### Prompt 5 — Docking Pipeline

```
Write services/docking.py and scripts/04_dock_benchmark.py.

This is the core computational module. It docks small molecules against
HIV-1 protease structures using AutoDock Vina.

Key functions:

def smiles_to_pdbqt(smiles: str) -> str:
    Convert SMILES to PDBQT string in memory using RDKit + meeko.
    Steps:
    1. mol = Chem.MolFromSmiles(smiles)
    2. mol = Chem.AddHs(mol)
    3. AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    4. AllChem.MMFFOptimizeMolecule(mol)
    5. mk_prep = MoleculePreparation()
    6. molsetup_list = mk_prep(mol)
    7. pdbqt_string = PDBQTWriterLegacy.write_string(molsetup_list[0])
    Return pdbqt_string.

def dock_single(ligand_pdbqt_string, receptor_pdbqt_path,
                center, box_size, exhaustiveness=16, n_poses=5) -> dict:
    Use Vina Python API:
    v = Vina(sf_name='vina')
    v.set_receptor(str(receptor_pdbqt_path))
    v.set_ligand_from_string(ligand_pdbqt_string)
    v.compute_vina_maps(center=list(center), box_size=list(box_size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    energies = v.energies()
    Return {"delta_g": float(energies[0][0]), "n_poses": len(energies)}
    
    Wrap in try/except. On ANY exception, return
    {"delta_g": None, "n_poses": 0, "error": str(e)}

def dock_against_panel(smiles, drug_name, panel_df,
                       structures_dir, mutants_dir) -> pd.DataFrame:
    1. Convert SMILES to PDBQT once
    2. Dock against wildtype → wt_delta_g
    3. For each mutation in panel_df:
       - dock against mutants/{mutation}.pdbqt
       - compute delta_delta_g = mutant_delta_g - wt_delta_g
    4. Return DataFrame: drug, mutation, delta_g, delta_delta_g, n_poses

def get_benchmark_smiles() -> dict:
    Fetch SMILES from PubChem for all 7 PIs using config.PUBCHEM_CIDS.
    Cache to data/ligands/benchmark_drugs.json.

scripts/04_dock_benchmark.py:
1. Load or fetch benchmark SMILES
2. For each drug in config.PI_DRUGS:
   a. Load its mutation panel
   b. Call dock_against_panel()
   c. Print progress per drug
3. Concatenate all results, save to data/docking_results/benchmark_docking.parquet
4. Print summary: drug, n_mutations_docked, n_failures, mean_delta_delta_g

Use config values for DOCKING_CENTER, DOCKING_BOX_SIZE, VINA_EXHAUSTIVENESS.
Use multiprocessing.Pool for parallelism within each drug's panel.
Expected runtime: 1-3 hours for all 7 drugs on a modern laptop.
```

---

### Prompt 6 — Scoring + Validation

```
Write services/scoring.py, services/validation.py, and scripts/05_validate.py.

scoring.py implements:

def compute_robustness_scores(docking_results_df, panel_df) -> dict:
    Join on mutation column. Compute:
    - simple_mean_ddg: mean of delta_delta_g (excluding None/NaN)
    - prevalence_weighted_ddg: weighted mean using n_isolates as weights
    - worst_case_ddg: max delta_delta_g
    - robustness_0_100: max(0, 100 - prevalence_weighted_ddg * 33.3)
    Return as dict.

def compute_baseline_mutation_count(panel_df) -> float:
    Return count of mutations with mean_log_fold_resistance > 1.0
    (i.e., >10x fold resistance). This is the naive "count known DRMs" baseline.

def compute_baseline_wt_only(docking_results_df) -> float:
    Return the wildtype delta_g value only. No mutation information used.

validation.py implements:

def run_full_validation(docking_parquet, panels_dir, output_dir) -> pd.DataFrame:
    For each drug:
    1. Load docking results and panel
    2. For each mutation, pair: predicted delta_delta_g vs real
       mean_log_fold_resistance from the panel
    3. Compute Spearman and Pearson correlations for each scoring method
    4. Also compute correlations for the two baselines
    Return DataFrame with columns: drug, scoring_method, spearman_rho,
    spearman_pvalue, pearson_r, n_mutations

def plot_validation(results_df, docking_parquet, panels_dir, output_path):
    Create a figure with subplots:
    - Main plot: scatter of delta_delta_g (y) vs mean_log_fold_resistance (x)
      across all drugs, colored by drug
    - Annotate with overall Spearman rho
    - Secondary panel: bar chart comparing Spearman rho across scoring methods
      and baselines
    Use matplotlib. Save as PNG at 300 DPI.

scripts/05_validate.py:
1. Run full validation
2. Print a formatted table of results
3. Generate and save the plot
4. Print the headline number: "Overall Spearman rho = X.XX (p = X.XX)"
```

---

### Prompt 7 — Explanation Service

```
Write services/explanation.py and scripts/06_generate_explanations.py.

explanation.py uses the Anthropic Python SDK to generate per-mutation
mechanistic explanations.

def build_structural_context(mutation, docking_result_dict,
                             wildtype_pdb_path) -> dict:
    Using BioPython's PDB parser:
    1. Parse the wildtype PDB
    2. Extract the residue at the mutation position
    3. Compute distance from residue CA to DOCKING_CENTER (proxy for
       distance to active site center)
    4. Use a hardcoded dict of amino acid properties:
       AA_PROPERTIES = {
           'A': {'volume': 88.6, 'charge': 'neutral', 'hydrophobicity': 1.8},
           'V': {'volume': 140.0, 'charge': 'neutral', 'hydrophobicity': 4.2},
           ... (all 20 standard AAs)
       }
    5. Determine the subpocket region based on position:
       SUBPOCKET_MAP = {
          range(23,33): "catalytic_triad",
          range(43,58): "flap_region",
          range(80,85): "S1_subpocket",
          range(25,30): "active_site_floor",
          ... (fill in from HIV protease structural knowledge)
       }
    Return dict matching the structural_context schema from the build doc.

def generate_explanation(drug, mutation, delta_delta_g,
                         structural_context, cache_dir) -> str:
    Check cache: if data/explanations/{drug}_{mutation}.json exists, load and
    return the explanation field.

    Otherwise, call the Anthropic API:
    import anthropic
    client = anthropic.Anthropic()  # uses ANTHROPIC_API_KEY env var
    
    response = client.messages.create(
        model="claude-sonnet-4-6",
        max_tokens=1024,
        system="You are a structural biologist specializing in HIV-1 protease
                drug resistance. Given structural context about a mutation and
                its effect on drug binding, provide a 2-3 sentence mechanistic
                hypothesis for WHY this mutation causes resistance to this
                specific drug. Ground your explanation in the structural facts
                provided. Be specific about molecular interactions.",
        messages=[{
            "role": "user",
            "content": f"Drug: {drug}\nMutation: {mutation}\n..."
        }]
    )
    explanation = response.content[0].text
    
    Save full result (including context, timestamp, model) to cache as JSON.
    Return explanation.

scripts/06_generate_explanations.py:
1. Load benchmark docking results
2. For each (drug, mutation) pair where |delta_delta_g| > 0.5 kcal/mol:
   a. Build structural context
   b. Generate explanation (cache-aware)
3. Print: "Generated X explanations, Y from cache"
4. Add a sleep(0.5) between API calls to stay under rate limits
```

---

### Prompt 8 — API + Frontend

```
Write api/main.py (FastAPI) and a React frontend in frontend/.

api/main.py:

POST /triage
  Request body: {"smiles": str, "target": "HIV1_PR"}
  Response: {
    "robustness_score": float,
    "wildtype_binding": float,
    "mutations": [
      {"mutation": "V82A", "delta_g": float, "delta_delta_g": float,
       "severity": "low"|"medium"|"high", "explanation": str},
      ...
    ],
    "scoring_methods": {"simple_mean": float, "weighted": float, "worst_case": float}
  }
  
  Implementation: call dock_against_panel for the input SMILES, compute scores,
  generate explanations for top 5 worst mutations only (to limit API calls and
  latency). Use threading to run explanation generation in parallel.

GET /benchmark
  Returns the precomputed validation data: scatter plot data points + Spearman rho.

GET /health
  Returns {"status": "ok", "n_mutants_cached": int}

Frontend (React + Vite + Tailwind):
- Input: SMILES text field with a "Triage" button, plus a dropdown to select
  a benchmark drug instead of typing SMILES
- Loading state while docking runs (warn user it takes 2-5 minutes)
- Results: large robustness score (0-100) with color coding (green/yellow/red),
  sortable mutation table with delta_delta_g bars and explanation expandable rows,
  and a "Validation" tab showing the benchmark scatter plot

Keep the design clean and clinical, not flashy. Use a monospace font for
SMILES display. White background, minimal color.

PLEASE use the frontend-design skill of Claude to design this

If time is short, skip the frontend entirely. Write notebooks/demo.ipynb
instead: a Jupyter notebook that walks through one drug (darunavir) end to end,
showing the docking, scoring, explanation, and validation steps with inline
plots. This is a perfectly acceptable demo format.
```

---

### Prompt 9 — Demo + Polish

```
Write demo.py as a command-line demo script and finalize README.md.

demo.py:
1. Print a banner: "ResistScope — HIV-1 Protease Resistance Triage"
2. Take an optional --smiles argument, default to darunavir SMILES
3. Run the full triage pipeline:
   a. Convert SMILES, dock against wildtype, print binding affinity
   b. Dock against top 10 most clinically relevant mutations
   c. Print a formatted table: mutation, delta_g, delta_delta_g, severity
   d. Print robustness score
   e. Generate and print explanation for the worst mutation
   f. Print validation summary (if benchmark data exists)
4. Total runtime should be ~5-10 minutes for 10 mutations

README.md should include:
- Project title and one-paragraph description
- Architecture diagram (ASCII art matching the one in section 1 of the build doc)
- Setup instructions (conda env create, download scripts)
- Quick start (python demo.py)
- Validation results (fill in after running: "ResistScope achieves Spearman
  rho = X.XX against real clinical fold-resistance data")
- Known limitations (crystal waters stripped, single protonation state for
  ligands, no flexible receptor docking, HIV-1 protease only)
- Future work (RT target, flexible docking, ensemble methods, wet-lab validation)
- Citation: Rhee et al. 2006 PNAS, Stanford HIVdb
```

---

## 9. Paper Path (Post-Hackathon)

**Target venue:** NeurIPS 2026 AI4D3 workshop (Drug Discovery and Development)
or MLSB workshop. Both in Sydney, Dec 6-12. Submission deadlines likely August
2026 based on prior years.

**Paper title idea:** "ResistScope: Structure-Based Resistance Prediction for
Antiviral Compounds with Interpretable Mechanistic Explanations"

**What to add beyond the hackathon:**
1. Extend to RT inhibitors (6 NRTIs + 3 NNRTIs in the Rhee dataset)
2. Ablate scoring methods properly (mean vs weighted vs worst-case)
3. Formal faithfulness evaluation of explanations against IAS-USA ground truth
4. Compare against Stanford HIVdb's own penalty scoring system as a third baseline
5. Cross-validation: leave-one-drug-out to test generalization

**What makes it publishable:**
- Empirical claim: "automated structure-based rescoring predicts real clinical
  fold-resistance with Spearman rho = X" — testable, reproducible
- Interpretability claim: "LLM-generated mechanistic explanations agree with
  expert-annotated mechanisms Y% of the time" — novel evaluation angle
- Practical claim: "this reduces days of manual triage to minutes" — tool paper

**Non-archival workshops welcome concurrent submissions and work in progress.**