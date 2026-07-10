import os
from pathlib import Path

import targets
from targets import Target, get_target, TARGETS  # noqa: F401 (re-exported)

# === Paths ===
PROJECT_ROOT = targets.PROJECT_ROOT


# === .env loading (zero-dependency) ===
def _load_dotenv(path: Path) -> None:
    """Populate os.environ from a .env file, without overwriting existing vars.

    A tiny, dependency-free parser (no python-dotenv): ``KEY=VALUE`` per line,
    ``#`` comments and blank lines ignored, optional surrounding quotes and a
    leading ``export`` stripped. Existing environment variables win, so an
    explicitly exported key still overrides the file. Runs once on import, so
    every entrypoint (CLI, scripts, API) sees ANTHROPIC_API_KEY from .env.
    """
    if not path.exists():
        return
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[len("export "):].lstrip()
        key, sep, value = line.partition("=")
        if not sep:
            continue
        key = key.strip()
        value = value.strip().strip('"').strip("'")
        if key and key not in os.environ:
            os.environ[key] = value


_load_dotenv(PROJECT_ROOT / ".env")

# Shared across all targets.
DATA_DIR = targets.DATA_DIR
RAW_DIR = targets.RAW_DIR


# =============================================================================
# Active target
# =============================================================================
# All *target-specific* constants below (paths, PDB, docking box, drugs,
# mutations) are re-exported from ACTIVE_TARGET under their legacy names so the
# existing pipeline and committed data keep working unchanged. HIV-1 protease is
# the default; the data-generation scripts switch it per run with
# ``config.set_active_target("rt")`` (see targets.py).
#
# Because these are module-level names read at *call time*, a switch propagates
# to any code that does ``config.X``. Functions that captured a value as a
# default argument (e.g. ``dock_single(center=config.DOCKING_CENTER)``) keep the
# import-time target — the target-aware entrypoints accept an explicit target or
# directory instead of relying on the mutated global.

ACTIVE_TARGET: Target = targets.HIV1_PR


def _bind(t: Target) -> None:
    """Re-export ``t``'s fields under the legacy module-level names."""
    g = globals()
    g["ACTIVE_TARGET"] = t

    # --- per-target data directories ---
    g["STRUCTURES_DIR"] = t.structures_dir
    g["MUTANTS_DIR"] = t.mutants_dir
    g["PANELS_DIR"] = t.panels_dir
    g["DOCKING_DIR"] = t.docking_dir
    g["EXPLANATIONS_DIR"] = t.explanations_dir
    g["VALIDATION_DIR"] = t.validation_dir

    # --- receptor structure ---
    g["WILDTYPE_PDB_ID"] = t.pdb_id
    g["WILDTYPE_PDB_URL"] = t.pdb_url
    g["LIGAND_HETCODES_TO_STRIP"] = list(t.ligand_hetcodes)
    g["CHAINS_TO_KEEP"] = list(t.chains)
    g["MUTATE_CHAINS"] = list(t.mutate_chains)
    g["REFERENCE_SEQUENCE"] = t.reference_seq
    g["N_POSITIONS"] = t.n_positions

    # --- docking box ---
    g["DOCKING_CENTER"] = t.docking_center
    g["DOCKING_BOX_SIZE"] = t.docking_box_size

    # --- catalytic-dyad protonation (HIV protease only; None for RT) ---
    g["PROTONATION"] = t.protonation
    g["ASP25_PROTONATED_CHAIN"] = t.protonation.protonated_chain if t.protonation else None
    g["ASP25_DEPROTONATED_CHAIN"] = t.protonation.deprotonated_chain if t.protonation else None

    # --- resistance dataset + drug panel ---
    g["DATASET_FILENAME"] = t.dataset_filename
    g["DATASET_URLS"] = list(t.dataset_urls)
    g["DRUG_COLUMNS"] = list(t.drug_columns)
    g["DRUGS"] = t.drugs
    g["PUBCHEM_CIDS"] = t.pubchem_cids
    g["PRIMARY_MUTATIONS"] = t.primary_mutations
    # Legacy aliases (protease-era names), kept so existing code/imports work.
    g["PI_DRUGS"] = t.drugs
    g["PRIMARY_PI_MUTATIONS"] = t.primary_mutations


def set_active_target(name) -> Target:
    """Switch the active target (by name/alias or Target) and re-bind globals.

    Call this once, early, in a data-generation script (before other modules
    read ``config.X``). Returns the resolved Target.
    """
    t = name if isinstance(name, Target) else get_target(name)
    _bind(t)
    return t


_bind(ACTIVE_TARGET)

# Landing page for the Rhee 2006 PNAS analysis (best-effort scrape in 01).
RHEE_DATASET_URL = (
    "https://hivdb.stanford.edu/pages/published_analysis/"
    "genophenoPNAS2006/"
)


# =============================================================================
# Method constants (shared across all targets)
# =============================================================================

# === Docking (AutoDock Vina / Uni-Dock) ===
VINA_EXHAUSTIVENESS = 16             # Balance speed vs accuracy; 32 for final
VINA_NUM_POSES = 5
VINA_NUM_CPUS = 0                    # 0 = use all available

# === Claude API ===
# IMPORTANT: the larger models (Opus 4.8/4.6, Sonnet 4.6) run a bio-safety
# classifier that FALSE-POSITIVE REFUSES HIV drug-resistance content
# (stop_reason="refusal", empty output) even though it is mainstream, published
# science. Haiku 4.5 does not refuse and produces accurate, specific structural
# explanations — and it is well-suited to this short, heavily-grounded task.
# Used for both explanation generation and the faithfulness judge.
CLAUDE_MODEL = "claude-haiku-4-5"
CLAUDE_MAX_TOKENS = 1024

# === Scoring ===
# Delta-delta-G thresholds (kcal/mol) for flagging concerning mutations.
DDG_WARNING_THRESHOLD = 0.5    # Mild concern
DDG_DANGER_THRESHOLD = 1.5     # Likely resistance
