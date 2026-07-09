# GPU docking with Uni-Dock

ResistScope's default docking backend is CPU AutoDock Vina (~20–30 s/dock,
so ~9–11 h for the full ~1,600-dock benchmark). The **Uni-Dock** backend runs
the *same Vina scoring function* on an NVIDIA GPU, cutting the full sweep to
~20–30 min. This guide sets it up on a GPU box and validates it before you
trust a full run.

> **Scores are directly comparable.** Uni-Dock implements Vina's scoring
> function, so affinities should match the CPU numbers (darunavir vs wildtype
> ≈ −9 to −10 kcal/mol) up to search stochasticity. If your validation run lands
> far outside that, stop and debug before scaling — same rule as the CPU sanity
> dock.

> **Status:** the `unidock` backend in `services/docking.py` is written against
> Uni-Dock's documented CLI but has NOT been run on a GPU from this machine.
> The validation step below is mandatory the first time.

---

## 1. Requirements

- NVIDIA GPU (Uni-Dock needs CUDA; it does **not** run on Apple/AMD GPUs).
- CUDA toolkit ≥ 11.
- Verify the GPU is visible:
  ```bash
  nvidia-smi
  ```

## 2. Install Uni-Dock

Pick one:

**A. conda-forge (simplest, if available for your CUDA):**
```bash
conda activate resistscope
conda install -c conda-forge unidock
unidock --version        # confirm it runs and sees the GPU
```

**B. Uni-Dock Tools (pip wrapper + CLI):**
```bash
pip install unidock_tools
# provides the `unidock` engine; see https://github.com/dptech-corp/Uni-Dock
```

**C. Build from source** (if no prebuilt binary matches your CUDA):
follow https://github.com/dptech-corp/Uni-Dock — CMake + CUDA build producing a
`unidock` executable, then put it on `PATH` or pass `--unidock-bin /path/to/unidock`.

## 3. Get the data onto the GPU box

The receptor cache (`data/structures/`) and panels (`data/panels/`) are
gitignored, so either **regenerate** them (self-contained) or **copy** them.

**Regenerate (recommended — reproducible):**
```bash
git clone <repo> && cd gladstone-hackathon
conda env create -f environment.yml && conda activate resistscope
python -m pip install gemmi                     # meeko needs it (see README)
python scripts/01_download_data.py              # Rhee dataset + 3OXC
python scripts/02_build_panels.py               # per-drug panels
python scripts/03_build_mutant_cache.py         # wildtype + 285 mutant PDBQTs (~15 min)
```

**Or copy** the prebuilt artifacts from your laptop:
```bash
rsync -av data/structures/ gpubox:.../data/structures/
rsync -av data/panels/ gpubox:.../data/panels/
rsync -av data/ligands/ gpubox:.../data/ligands/   # cached benchmark SMILES
```

## 4. Validate before the full run  ⚠️ required

Dock darunavir against wildtype + its primary DRMs on the GPU and check the
numbers match the CPU reference:

```bash
python scripts/04_dock_benchmark.py \
    --backend unidock --unidock-bin unidock \
    --drugs DRV --subset primary \
    --out data/docking_results/_validation_drv.parquet
```

Then confirm in Python:
```python
import pandas as pd
r = pd.read_parquet("data/docking_results/_validation_drv.parquet")
print(r[["mutation","delta_g","delta_delta_g"]])
```

Pass criteria (same as the CPU sanity dock):
- **DRV vs WT `delta_g` ≈ −9 to −11 kcal/mol** (matches CPU ≈ −9.3). Far outside
  ⇒ box/receptor/scoring problem — stop and debug.
- Known DRV major mutations (I50V, I47V, V32I, I54M, I84V) should give
  **`delta_delta_g` ≥ 0** (weaker binding). Darunavir's ddGs are legitimately
  small — it is the highest-barrier PI.

Optionally compare a handful of pairs head-to-head against the CPU parquet; the
per-pair `delta_g` should agree within ~0.5–1 kcal/mol (search noise).

## 5. Run the full benchmark

Once validation passes:
```bash
python scripts/04_dock_benchmark.py \
    --backend unidock --subset all \
    --search-mode balance \
    --out data/docking_results/benchmark_docking.parquet
```
- `--search-mode {fast,balance,detail}` is Uni-Dock's effort knob (overrides
  `--exhaustiveness`); `balance` is a good default, `detail` for the final run.
- `--replicates 3` averages seeded runs for noise-robust ddG (cheap on GPU;
  recommended given how small darunavir's ddGs are).
- Output parquet feeds directly into `scripts/05_validate.py` (Step 5).

## 6. How the backend maps onto Uni-Dock

`services/docking.py :: dock_single_unidock` invokes:
```
unidock --receptor <mutant>.pdbqt --gpu_batch <drug>.pdbqt \
        --center_x/y/z <DOCKING_CENTER> --size_x/y/z 22 22 22 \
        --scoring vina --num_modes 5 --dir <tmp> [--search_mode M | --exhaustiveness E]
```
and parses the best `REMARK VINA RESULT` affinity from `<tmp>/<drug>_out.pdbqt`.
Box center/size come from `config.DOCKING_CENTER` / `config.DOCKING_BOX_SIZE`.

**Throughput note.** ResistScope docks *one drug against many receptors*, whereas
Uni-Dock's `--gpu_batch` batches *many ligands against one receptor*. The current
backend calls Uni-Dock once per (drug, receptor) — already GPU-fast, but a future
optimization is to invert the loop (iterate receptors, batch all 6 drugs per
receptor via `--gpu_batch`) to exploit Uni-Dock's ligand batching. Left as an
enhancement; not required for the hackathon.
