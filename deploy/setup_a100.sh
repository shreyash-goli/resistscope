#!/usr/bin/env bash
# Set up a fresh NVIDIA GPU instance (e.g. Prime Intellect A100) to run the
# ResistScope Uni-Dock backend. Run this ON the instance, from the project root
# (the directory containing config.py), AFTER extracting the deploy tarball.
#
#   bash deploy/setup_a100.sh
#
# It installs Miniconda (if missing), creates a lean docking env, installs
# Uni-Dock, and verifies the GPU. Structures are already bundled, so no
# openmm/pdbfixer/vina/structure-prep is needed here.
set -euo pipefail

echo "==================================================================="
echo " ResistScope GPU setup"
echo "==================================================================="

# --- 0. GPU visible? ------------------------------------------------------
echo "== nvidia-smi =="
if ! nvidia-smi; then
  echo "ERROR: no NVIDIA GPU visible. Uni-Dock requires CUDA." >&2
  exit 1
fi

# --- 1. conda -------------------------------------------------------------
if ! command -v conda >/dev/null 2>&1; then
  echo "== Installing Miniconda =="
  curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
  bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
fi
# shellcheck disable=SC1091
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"

# --- 2. lean docking env --------------------------------------------------
# Only what the Uni-Dock path needs: RDKit + meeko (SMILES->PDBQT ligand prep),
# pandas/pyarrow (panels + results), requests (PubChem, optional since cached).
ENV=resistscope_gpu
if ! conda env list | grep -q "^${ENV}\b"; then
  echo "== Creating conda env '${ENV}' =="
  conda create -y -n "$ENV" -c conda-forge \
      python=3.11 rdkit pandas numpy scipy pyarrow requests pip
fi
conda activate "$ENV"
python -m pip install --quiet meeko gemmi

# --- 3. Uni-Dock ----------------------------------------------------------
if ! command -v unidock >/dev/null 2>&1; then
  echo "== Installing Uni-Dock =="
  if ! conda install -y -n "$ENV" -c conda-forge unidock; then
    echo "conda-forge unidock failed; trying pip unidock_tools ..."
    python -m pip install unidock_tools || {
      echo "ERROR: could not install Uni-Dock. Build from source:" >&2
      echo "  https://github.com/dptech-corp/Uni-Dock" >&2
      exit 1
    }
  fi
fi

echo "== Uni-Dock check =="
unidock --version || echo "WARN: 'unidock --version' failed; check the install."

# --- 4. quick import check ------------------------------------------------
python - <<'PY'
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy
import pandas, gemmi
print("Python deps OK (rdkit, meeko, pandas, gemmi)")
PY

echo
echo "==================================================================="
echo " Setup complete. Next:"
echo "   conda activate ${ENV}"
echo "   # 1) VALIDATE (required first run):"
echo "   python scripts/04_dock_benchmark.py --backend unidock --drugs DRV \\"
echo "       --subset primary --out data/docking_results/_validation_drv.parquet"
echo "   #    -> confirm DRV vs WT delta_g is about -9 to -11 kcal/mol"
echo "   # 2) FULL sweep:"
echo "   python scripts/04_dock_benchmark.py --backend unidock --subset all \\"
echo "       --search-mode balance --replicates 3"
echo "==================================================================="
