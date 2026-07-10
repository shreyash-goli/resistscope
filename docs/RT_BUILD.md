# Reverse Transcriptase (RT) extension — build runbook

ResistScope was built for HIV-1 **protease**. This doc is the runbook for the
**reverse transcriptase (NNRTI)** extension, built on the `Target` abstraction
in [`targets.py`](../targets.py).

## Scientific scope (read first)

The docking-ΔΔG method models **loss of binding affinity**. That is a sound
proxy for **NNRTIs** — they bind an allosteric hydrophobic pocket ~10 Å from the
polymerase active site, competitively, exactly like PIs bind the protease active
site. It is **not** a sound proxy for **NRTIs**: NRTI resistance is about
nucleotide *incorporation* and ATP-mediated *excision* (TAMs), not pocket
affinity, so rigid docking of a nucleoside analog will not capture it.

**→ This target is scoped to the NNRTI pocket and NNRTI drugs only**
(nevirapine, efavirenz, etravirine, rilpivirine; doravirine has no column in the
Rhee dataset and is skipped). Do not add NRTIs to this target expecting docking
to explain their resistance.

## What's done and verified

| Layer | File | Status |
|---|---|---|
| Target abstraction | [`targets.py`](../targets.py) | ✅ `HIV1_PR` + `HIV1_RT`, reference sequences anchor-checked |
| Config re-export + `set_active_target()` | [`config.py`](../config.py) | ✅ PI values bit-identical; RT switch verified |
| Panel parsing (ref seq, drug cols, DRMs, 240 positions) | [`services/mutation_panel.py`](../services/mutation_panel.py) | ✅ built on **real** NNRTI data: NVP 453 / EFV 455 / ETR 323 / RPV 135 muts |
| Structure prep (heterodimer, **no** Asp25 dyad, mutate p66 only) | [`services/structure_prep.py`](../services/structure_prep.py) | ✅ compiles; target-driven |
| Download (RT dataset + 3V81 PDB) | [`scripts/01_download_data.py`](../scripts/01_download_data.py) | ✅ downloaded + committed under `data/raw/` |
| Build panels | [`scripts/02_build_panels.py`](../scripts/02_build_panels.py) | ✅ run; `data/rt/panels/*.parquet` committed |
| Docking box + ligand het code | [`targets.py`](../targets.py) | ✅ **confirmed**: center `(41.105,52.332,49.098)`, het `NVP` (was wrongly `RIL`) |
| Build mutant cache | [`scripts/03_build_mutant_cache.py`](../scripts/03_build_mutant_cache.py) | ⏳ needs GPU box (`--target rt`) |
| Dock (CPU + GPU) | [`scripts/04_dock_benchmark.py`](../scripts/04_dock_benchmark.py) / [`04_gpu_batch.py`](../scripts/04_gpu_batch.py) | ⏳ needs GPU box (`--target rt`) |
| Validation + rigorous benchmark | [`scripts/05_validate.py`](../scripts/05_validate.py) + [`scripts/08_benchmark.py`](../scripts/08_benchmark.py) | ✅ `--target rt`-ready (runs after docking) |
| Ground truth (mechanisms) | [`scripts/09_build_ground_truth.py`](../scripts/09_build_ground_truth.py) | ✅ agent-built + verified: K103N / Y181C / Y188L (extend `--primary`) |
| Explanations (NNRTI pocket + prompts) | [`services/explanation.py`](../services/explanation.py) + [`scripts/06`](../scripts/06_generate_explanations.py) | ⚠️ wired; RT region map needs expert review |
| Faithfulness eval | [`scripts/07_faithfulness_eval.py`](../scripts/07_faithfulness_eval.py) | ✅ `--target rt` (SDK bug fixed; scores vs agent ground truth) |
| API + UI (target selector) | [`api/main.py`](../api/main.py) + [`frontend/`](../frontend) | ✅ per-request `?target=`; header target switcher |

RT artifacts live under `data/rt/` (`subdir="rt"`), so **building RT never
touches the committed PI data**. `data/rt/structures/` is gitignored (large,
regenerable), matching the PI layout.

## Two data values — now CONFIRMED

Both were previously placeholders; both are now resolved (2026-07, from the
deposited 3V81 coordinates — see [`targets.py`](../targets.py)):

1. **Docking box** — ✅ `HIV1_RT.docking_center = (41.105, 52.332, 49.098)`, the
   occupancy-weighted centroid of the **chain-A NVP** atoms (the NNRTI pocket of
   the A/B heterodimer). `structure_prep.clean_structure()` re-prints this during
   step 03 as a cross-check.
2. **Ligand het code** — ✅ **CORRECTED**: the earlier `ligand_hetcodes=("RIL",)`
   was **wrong**. 3V81's deposited coordinates contain **nevirapine (het `NVP`)**
   in the NNRTI pocket, not rilpivirine — the "TMC278/RIL" in the paper title is a
   related structure in the same study. It is now `ligand_hetcodes=("NVP",)`.
   (Also present but *not* used for the box: `MRG`/`ATM` nucleotides at the
   polymerase active site ~10 Å away — do not center on those.) The reference
   sequence is anchored at 18 canonical DRM positions (import fails loudly on a
   transcription slip) and the built panels put the wildtype residue at Y181/Y188/
   G190/V179 exactly as expected.

## Remaining work

The pipeline is `--target`-aware, the two data blockers are resolved, panels are
built from real data, the API is target-aware, and the UI has a target selector.
What's left is the **GPU docking run** (below) plus one scientific review item:

- **Expert review** — [`services/explanation.py`](../services/explanation.py)
  `_RT_SUBPOCKET_REGIONS` (rim 100–108, core 179–190 incl. Y181/Y188/G190,
  primer-grip wall 227–236) is a **coarse contiguous-range approximation** of a
  discontiguous pocket; have a structural expert review the labels before
  trusting RT explanations.
- **Ground truth — now built agentically.** `data/rt/mechanism_ground_truth.json`
  is seeded by the literature agent (`scripts/09`, Claude + PubMed, self-verified;
  K103N / Y181C / Y188L to start). Extend with `--primary` for the full NNRTI DRM
  set. Faithfulness (07) can score RT the moment RT explanations exist.

## Runbook — GPU docking hand-off

Data prep (01, 02, 09) is **already done locally and committed**. The remaining
steps need the docking stack; run 03–07 on the GPU box (see
[`deploy/setup_a100.sh`](../deploy/setup_a100.sh)).

```bash
# --- already done locally (no GPU) — committed under data/rt/ ---
# python scripts/01_download_data.py --target rt      # NNRTI dataset + 3V81 PDB  ✅
# python scripts/02_build_panels.py  --target rt      # NVP/EFV/ETR/RPV panels    ✅
# python scripts/09_build_ground_truth.py --target rt --primary   # agent ground truth ✅ (seeded)

# --- on the GPU box (conda env from deploy/setup_a100.sh) ---
python scripts/03_build_mutant_cache.py --target rt  # p66 mutant receptors -> data/rt/structures
#   ^ box is already set; step 03 re-prints the NVP centroid — confirm it matches (41.105, 52.332, 49.098)
python scripts/04_gpu_batch.py --target rt --subset all --search-mode detail   # (or 04_dock_benchmark.py on CPU)
python scripts/05_validate.py  --target rt           # correlations + figure -> data/rt/validation
python scripts/08_benchmark.py --target rt           # enrichment CIs / ROC-AUC / de-confounding
python scripts/06_generate_explanations.py --target rt --cite   # NNRTI-pocket explanations + PubMed
python scripts/07_faithfulness_eval.py --target rt   # scores RT explanations vs the agent ground truth
```

After the run, `data/rt/{docking_results,validation,explanations}` are populated
and the app's RT target lights up automatically — no code changes needed.

## API / frontend — DONE (target-aware)

`api/main.py` now resolves `config.get_target(...)` **per request**: `/targets`,
`/health`, `/drugs`, `/drug/{abbrev}`, `/benchmark`, `/validation/plot`,
`/structure/*`, and `/triage` all take `?target=` (default `HIV1_PR`) and read
that target's dirs — no global mutation. The React app has a **target selector**
in the header; selecting RT loads its NNRTI drug list + agent-built mechanisms
and shows a "docking pending" banner until `data/rt/docking_results/` exists,
after which precomputed ΔΔG + the Validation tab populate automatically.
