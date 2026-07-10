# ResistScope

**Resistance-aware triage for antiviral compounds targeting HIV-1 protease.**
Paste a SMILES string for a candidate protease inhibitor and ResistScope docks it
against wildtype HIV-1 protease and a panel of clinical resistance mutants,
computes a robustness score (how much binding degrades across the panel),
generates a per-mutation mechanistic explanation with Claude (grounded in
PubMed literature), and validates the whole pipeline against real measured
fold-resistance data (Rhee et al. 2006, Stanford HIVdb). Built for the
*Built with Claude: Life Sciences* hackathon (Gladstone Institutes).

**Highlights**
- **Interactive 3D structure viewer** (Mol\*/NGL) — every mutation renders in the
  receptor with the inhibitor in the pocket, the mutated residue highlighted and
  colored by ΔΔG, right next to Claude's explanation.
- **Agentic literature grounding** — a Claude tool-use agent (`scripts/09`)
  researches each mutation's mechanism over live PubMed, cites the papers it
  actually read, and self-verifies grounding before an entry enters the
  benchmark. This is how the resistance-mechanism ground truth is built.
- **Rigorous, honest validation** (`scripts/08`) — every claim carries a
  significance test or CI: permutation p-values + bootstrap CIs on DRM
  enrichment, ROC/PR-AUC for DRM recovery, and an explicit de-confounding
  analysis that reports what the method *cannot* do.
- **Two targets** — HIV-1 protease (validated) and HIV-1 reverse transcriptase
  (NNRTI pocket); the whole pipeline + API + UI are target-aware.
- **Bring your own target** — in the app, upload a receptor PDB and name a
  protein; the structure is parsed for the pocket and a Claude agent researches
  its inhibitors + resistance mutations from PubMed, producing a new triage
  target in seconds (demonstrated on influenza neuraminidase / Tamiflu — the
  agent recovers H275Y, E119V, R292K, …). Triage generalizes to any target;
  validation attaches wherever a resistance dataset exists.
- **Live docking, no precompute** — a custom SMILES docks live through a
  pluggable backend (local CPU, or a remote **GPU worker** you connect via
  `RESISTSCOPE_DOCKING_URL`); benchmark drugs stay a precomputed cache.

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
       │  services/literature_agent.py      │   agentic PubMed research
       │  structural context -> mechanism   │   + self-verified ground truth
       └────────────────┬───────────────────┘
                        │
          api/main.py (FastAPI, target-aware)  ──►  frontend/ (React + Vite + Tailwind
          demo.py (CLI)                        ──►    + NGL 3D viewer, target selector)
                                               ──►  notebooks/demo.ipynb

  data pipeline: scripts/01 download → 02 panels → 03 mutant cache → 04 dock →
                 05 validate → 06 explain → 07 faithfulness →
                 08 rigorous benchmark → 09 agentic ground truth
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
export ANTHROPIC_API_KEY=sk-ant-...   # or put it in .env (auto-loaded)
python scripts/09_build_ground_truth.py --primary --cite   # agentic PubMed ground truth
python scripts/06_generate_explanations.py --ground-truth --cite   # + PubMed citations
python scripts/07_faithfulness_eval.py    # Claude-judge faithfulness (curated + agent split)
python scripts/08_benchmark.py            # permutation p / bootstrap CI / ROC-AUC / de-confounding
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

## Live docking (custom SMILES, no precompute)

Benchmark drugs are a precomputed cache, but any custom SMILES docks **live** —
`/triage` runs it against the resistance panel through a pluggable backend
([services/docking_backend.py](services/docking_backend.py)), resolved once at
startup:

| Backend | When | How |
|---|---|---|
| **Remote GPU** | `RESISTSCOPE_DOCKING_URL` set | API forwards to a Uni-Dock worker you run — **this is how you "connect a GPU"** |
| **Local** | docking stack installed on the API host | in-process CPU Vina (or GPU Uni-Dock) |
| **None** | neither | clean **503** with how-to-enable guidance (never an opaque 500) |

**Connect a GPU** (fast live docking, seconds instead of minutes):

```bash
# on a GPU box (e.g. A100) — one-time
bash deploy/setup_a100.sh && conda activate resistscope_gpu
python scripts/03_build_mutant_cache.py --target pr    # build receptors on this box
python docking_worker.py                                # serves POST /dock on :9000

# point the API at it (any host)
RESISTSCOPE_DOCKING_URL=http://<gpu-host>:9000 python -m uvicorn api.main:app --port 8000
```

The web app shows a live-docking status dot (green = a backend is ready) so users
know a custom SMILES will actually run before they submit it.

## Validation results

Validated against measured Rhee fold-resistance across 6 PIs (1,591 drug-mutation
pairs). Honest summary:

Every number below carries a significance test or confidence interval
(`scripts/08_benchmark.py`: permutation p-values, bootstrap CIs, ROC/PR-AUC):

| metric | result |
|---|---|
| Top-40 ΔΔG **DRM enrichment** | **2.86×** (95% CI 1.63–4.08, permutation p < 0.001) — top predictions recover known DRMs |
| Pooled **DRM-recovery ROC-AUC** | **0.51** (95% CI 0.46–0.56) — ≈ chance as a *global* ranker |
| Darunavir major-DRM Spearman ρ | **≈ 0.40** (p = 0.03) — strong per-drug signal |
| Per-mutation Spearman ρ (pooled) | **≈ 0.00** — the measured target is confounded by co-occurring mutations |
| **De-confounding** (single-/≤2-mutation isolates) | does **not** rescue the magnitude correlation (ρ ≤ 0) — an honest bound |
| Explanation **faithfulness** | **72%** correctly identify the expert mechanism (mean 1.65 / 2, n = 46) |

The honest headline: rigid single-mutation docking ΔΔG is a **coarse DRM-triage
flag, not a quantitative resistance predictor** — its extreme predictions are
significantly enriched for real resistance mutations (≈3×, p < 0.001) even though
it ranks at chance overall, and de-confounding the target does not rescue the
magnitude correlation. Darunavir validates well per-drug, and the LLM
explanations agree with expert-annotated mechanisms 72% of the time.

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

- **Reverse transcriptase (NNRTI) target** — now wired end-to-end and its two
  data blockers are resolved (docking box + het code confirmed from 3V81; note the
  co-crystal ligand is nevirapine/`NVP`, *not* rilpivirine as first assumed). Real
  NNRTI panels are built and agent ground truth is seeded; the remaining step is
  the GPU docking run (`docs/RT_BUILD.md`).
- **Flexible-receptor / ensemble docking; explicit flap water** — the most likely
  path to a docking signal that beats the DRM-enrichment flag.
- **Stanford HIVdb penalty-score target** — single-mutation-isolate
  de-confounding was tested and does *not* rescue the correlation
  (`scripts/08`); the expert-derived penalty scores are the next de-confounded
  target to try, with leave-one-drug-out CV.
- **Wet-lab validation** and native *Claude for Life Sciences* MCP connectors
  (the literature agent already grounds in live PubMed via tool use).

## Citations

- Rhee S-Y, et al. *HIV-1 protease and reverse-transcriptase mutations for drug
  resistance surveillance.* PNAS 2006. Stanford HIV Drug Resistance Database
  (https://hivdb.stanford.edu).
- PDB `3OXC`; IAS-USA drug resistance mutations figures.
- Docking: AutoDock Vina (Trott & Olson 2010; Eberhardt et al. 2021),
  Uni-Dock (Yu et al. 2023), meeko.
