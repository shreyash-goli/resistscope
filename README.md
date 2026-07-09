# ResistScope

**Resistance-aware triage for antiviral compounds targeting HIV-1 protease.**
Paste a SMILES string for a candidate protease inhibitor and ResistScope docks it
against wildtype HIV-1 protease and a panel of clinical resistance mutants,
computes a robustness score (how much binding degrades across the panel),
generates a per-mutation mechanistic explanation with Claude (grounded in
PubMed literature), and validates the whole pipeline against real measured
fold-resistance data (Rhee et al. 2006, Stanford HIVdb). Built for the
*Built with Claude: Life Sciences* hackathon (Gladstone Institutes).

MVP scope: HIV-1 protease inhibitors. Reverse transcriptase is a post-hackathon
extension.

---

## Architecture

```
                 SMILES (candidate PI)
                        │
       ┌────────────────▼───────────────────┐
       │  services/docking.py               │   RDKit + meeko -> PDBQT
       │  SMILES -> 3D -> AutoDock Vina      │   (GPU: Uni-Dock backend)
       │  vs wildtype + 285 mutant receptors │
       └────────────────┬───────────────────┘
                        │ ΔΔG per mutation
       ┌────────────────▼───────────────────┐
       │  services/scoring.py               │   robustness score (0-100)
       │  services/validation.py            │   vs Rhee fold-resistance
       └────────────────┬───────────────────┘
                        │
       ┌────────────────▼───────────────────┐
       │  services/explanation.py           │   Claude (Haiku 4.5) +
       │  structural context -> mechanism   │   PubMed citations +
       │  + faithfulness eval (Claude judge) │   0-2 faithfulness score
       └────────────────┬───────────────────┘
                        │
          api/main.py (FastAPI)  ──►  frontend/ (React + Vite + Tailwind)
          demo.py (CLI)          ──►  notebooks/demo.ipynb

  data pipeline: scripts/01 download → 02 panels → 03 mutant cache →
                 04 dock → 05 validate → 06 explain → 07 faithfulness
```

Key data:
- **Receptor:** PDB `3OXC` (HIV-1 protease + saquinavir, 1.16 Å), ligand stripped,
  Asp25 dyad asymmetrically protonated, homodimer mutations on both chains.
- **Mutation panels:** parsed from the Stanford HIVdb PI genotype-phenotype
  dataset; 285 unique mutant receptors built with PDBFixer + meeko.

---

## Setup

The scientific stack (RDKit, OpenMM, PDBFixer, AutoDock Vina, meeko) is installed
via conda:

```bash
conda env create -f environment.yml
conda activate resistscope
python -m pip install gemmi          # required by meeko 0.7.x (not auto-pulled)
```

Then reproduce the data + docking (run once):

```bash
python scripts/01_download_data.py        # Rhee dataset + 3OXC PDB
python scripts/02_build_panels.py         # per-drug mutation panels
python scripts/03_build_mutant_cache.py   # wildtype + 285 mutant receptors (~15 min)
python scripts/04_gpu_batch.py --subset all --search-mode detail   # docking (GPU: docs/UNIDOCK_SETUP.md)
python scripts/05_validate.py             # correlation analysis + figure
```

Explanations + faithfulness need an Anthropic key (uses **Claude Haiku 4.5** —
the larger models refuse HIV drug-resistance content via a bio-safety classifier):

```bash
export ANTHROPIC_API_KEY=sk-ant-...
python scripts/06_generate_explanations.py --ground-truth --cite   # + PubMed citations
python scripts/07_faithfulness_eval.py
```

---

## Quick Start

**CLI demo** (docks darunavir against the top clinical mutations end-to-end):

```bash
python demo.py                 # darunavir, top 10 mutations (~5-10 min)
python demo.py --mutations 5   # faster
python demo.py --smiles "CC(C)..."   # a custom compound
```

**Notebook demo** (self-contained, runs off cached results in <1 min):

```bash
jupyter notebook notebooks/demo.ipynb
```

**Interactive web app** (FastAPI + React):

```bash
# Terminal 1 — backend. Use `python -m uvicorn`, NOT bare `uvicorn`
# (bare uvicorn resolves to system Python without the deps). Set ANTHROPIC_API_KEY
# for live-triage explanations.
python -m uvicorn api.main:app --port 8000

# Terminal 2 — frontend (proxies /api -> :8000)
cd frontend && npm install && npm run dev     # http://localhost:5173
```

Benchmark drugs load instantly (precomputed); a custom SMILES runs live docking
(~2-5 min on CPU).

---

## Validation results

Validated against measured Rhee fold-resistance across 6 PIs (1,591 drug-mutation
pairs). Honest summary:

| metric | result |
|---|---|
| Per-mutation Spearman ρ (pooled) | **≈ 0.00** — weak; the measured target is confounded by co-occurring mutations |
| Darunavir major-DRM Spearman ρ | **≈ 0.40** (p = 0.03) — strong per-drug signal |
| Top-40 ΔΔG **DRM enrichment** | **≈ 2.9×** (35% vs 12% base rate) — top predictions recover known DRMs |
| Explanation **faithfulness** | **72%** correctly identify the expert mechanism (mean 1.61 / 2, n = 46) |

The headline is not a single correlation: rigid single-mutation docking does *not*
quantitatively predict pooled clinical fold-resistance, **but** its top-ranked
predictions are meaningfully enriched for real resistance mutations, darunavir
validates well, and the LLM explanations agree with expert-annotated mechanisms
72% of the time.

---

## Known limitations

- **Crystal waters stripped**, including the conserved flap water that bridges
  ligand to flaps — systematically weakens absolute affinities (documented; ΔΔG
  differences partially cancel it).
- **Single ligand protonation state** (RDKit `AddHs`, no pH enumeration).
- **Rigid-receptor docking** — only the mutated side chain is relaxed; no flexible
  side chains or backbone.
- **Confounded validation target** — per-mutation fold-resistance is averaged over
  isolates carrying *other* mutations, capping any single-mutation method.
- **HIV-1 protease only.**
- **Model refusals** — Opus/Sonnet decline HIV-resistance prompts; explanations
  run on Haiku 4.5.

## Future work

- Reverse transcriptase target (NRTIs/NNRTIs are in the same Rhee dataset).
- Flexible-receptor / ensemble docking; explicit flap water.
- Deconfounded target (e.g. Stanford HIVdb penalty scores) and leave-one-drug-out
  cross-validation.
- Native *Claude for Life Sciences* connectors (PubMed/ChEMBL MCP) in an
  enterprise context; wet-lab validation.

## Citations

- Rhee S-Y, et al. *HIV-1 protease and reverse-transcriptase mutations for drug
  resistance surveillance.* PNAS 2006. Stanford HIV Drug Resistance Database
  (https://hivdb.stanford.edu).
- PDB `3OXC`; IAS-USA drug resistance mutations figures.
- Docking: AutoDock Vina (Trott & Olson 2010; Eberhardt et al. 2021),
  Uni-Dock (Yu et al. 2023), meeko.
