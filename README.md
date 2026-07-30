# ResistScope

**Resistance-aware triage for antiviral compounds — for any target.**
Paste a candidate inhibitor's SMILES and ResistScope docks it against a wildtype
receptor plus a panel of clinical resistance mutants, scores how well its binding
*survives* the panel (a 0–100 robustness score), and — grounded in PubMed —
explains, per mutation, *why* each one threatens the compound. It's built to
answer two questions that are usually asked separately:

1. **Does the prediction hold up?** A benchmark against measured
   fold-resistance ([scripts/08](scripts/08_benchmark.py)): top-ΔΔG predictions
   are **~3× enriched** for real resistance mutations (p < 0.001), reported with
   every confidence interval and caveat — including where the method is only
   chance-level.
2. **Does the *explanation* hold up?** A faithfulness check
   ([scripts/07](scripts/07_faithfulness_eval.py)) that scores each LLM-written
   mechanism against expert-annotated biology (IAS-USA mutation lists), Claude as
   judge — **72%** recover the correct primary mechanism. Most tools in this space
   score; almost none ask whether the model's *reasoning* matches known biology.

The method is **target-agnostic**: it ships validated on **HIV-1 protease**
(Rhee et al. 2006, Stanford HIVdb), is wired end-to-end for **HIV-1 reverse
transcriptase**, and a Claude agent stands up a brand-new target — influenza
neuraminidase, SARS-CoV-2 Mpro, whatever PDB/mmCIF you upload — in seconds. Built
for the *Built with Claude: Life Sciences* hackathon (Gladstone Institutes).

**Who it's for.** An antiviral medicinal chemist has 20 candidate inhibitors
(for HIV protease, SARS-CoV-2 Mpro, influenza NA — pick a target) and needs to
decide which ones to push into (slow, expensive) phenotypic resistance assays.
They can't run all 20. ResistScope ranks the candidates by how well their
predicted binding holds up across the panel of known clinical resistance
mutations, flags *which* mutations break each compound, and — grounded in
PubMed — explains *why*, so the chemist prioritizes the compounds least
vulnerable to escape before committing wet-lab time.

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
- **Bring your own target** — in the app, upload a receptor PDB/mmCIF and name a
  protein; the structure is parsed for the pocket and a Claude agent researches
  its inhibitors + resistance mutations from PubMed, producing a new triage
  target in seconds. Two are committed as worked examples: **influenza
  neuraminidase / Tamiflu** (the agent recovers H275Y, E119V, R292K, …) and
  **SARS-CoV-2 main protease (Mpro) / Paxlovid** (recovers the nirmatrelvir
  escape set E166V, E166A, S144, …). Triage generalizes to any target;
  validation attaches wherever a resistance dataset exists. Building a BYO
  target's mutant receptors for live docking is GPU-gated — see
  [docs/BYO_TARGET_GPU.md](docs/BYO_TARGET_GPU.md).
- **Live docking, no precompute** — a custom SMILES docks live through a
  pluggable backend (local CPU, or a remote **GPU worker** you connect via
  `RESISTSCOPE_DOCKING_URL`); benchmark drugs stay a precomputed cache.

---

## The stack — what's in the box

ResistScope is a full pipeline, HTTP API, and web app — not a notebook. Every
layer is target-aware.

**Structural / docking engine** (`services/`, Python 3.11 + conda)
- [structure_prep.py](services/structure_prep.py) — PDB/mmCIF download, cleaning
  (PDBFixer), point mutagenesis, and PDBQT conversion (meeko).
- [docking.py](services/docking.py) — SMILES → 3D (RDKit) → PDBQT → **AutoDock
  Vina** (CPU) or **Uni-Dock** (GPU) against wildtype + every mutant receptor.
- [docking_backend.py](services/docking_backend.py) — pluggable backend resolved
  once at startup: **local** in-process, **remote GPU** worker, or a clean
  **503** if neither is configured. [docking_worker.py](docking_worker.py) is the
  standalone GPU service (`POST /dock`).
- [scoring.py](services/scoring.py) — per-mutation ΔΔG → 0–100 robustness score.

**Claude intelligence layer**
- [explanation.py](services/explanation.py) — Claude (**Haiku 4.5**; larger models
  refuse HIV-resistance content) writes a structural mechanism per mutation, plus
  a **Claude-as-judge** faithfulness scorer (0–2 vs. expert ground truth).
- [literature_agent.py](services/literature_agent.py) — a **Claude tool-use
  agent** that researches mechanisms over **live PubMed**, cites the papers it
  read, and self-verifies grounding.
- [target_builder.py](services/target_builder.py) / [pdb_intake.py](services/pdb_intake.py)
  — the **bring-your-own-target** agent: from a name + structure, parse the
  pocket and agent-propose the inhibitor + resistance-mutation set.

**Validation / benchmark** (`services/`, scipy + scikit-learn)
- [validation.py](services/validation.py) — predicted vs. measured fold-resistance,
  correlation + figures.
- [benchmark.py](services/benchmark.py) — permutation p-values, bootstrap CIs,
  ROC/PR-AUC, and de-confounding.

**API** — [api/main.py](api/main.py), **FastAPI**, target-aware. Endpoints:
`/targets`, `/targets/intake|assemble|save` (BYO flow), `/benchmark`, `/drug/{abbrev}`,
`/triage` (live docking), `/structure/receptor|ligand` (3D), `/validation/plot`, `/health`.

**Frontend** — [frontend/](frontend/), **React 18 + Vite + Tailwind**, with an
**NGL** WebGL 3D viewer. Components: `StructureViewer`, `MutationTable`, `ScoreCard`,
`ValidationTab`, `AddTargetModal` (upload → pocket → agent → new target).

**Reproducible pipeline** — `scripts/01`→`10`: download → panels → mutant cache →
dock → validate → explain → faithfulness → rigorous benchmark → agentic ground
truth → add BYO target. Plus [demo.py](demo.py) (CLI), `notebooks/demo.ipynb`, and
[scripts/smoke_test.py](scripts/smoke_test.py).

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
                 08 rigorous benchmark → 09 agentic ground truth →
                 10 add BYO target (agent-researched receptor + mutations)
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

### Demo path (no docking stack required)

The reliable walkthrough — everything below is served from the committed cache,
so it works with **just the API + frontend running**, no conda docking env, no
GPU, no Anthropic key:

1. **Pick a benchmark drug** — start with **darunavir (DRV)** or **saquinavir
   (SQV)**. It loads instantly and returns a full scored panel.
2. **Read the robustness score** and open the worst mutations (highest ΔΔG).
3. **Click a mutation** — the 3D viewer highlights the mutated residue in the
   pocket and shows Claude's PubMed-cited mechanistic explanation.
4. **Open the Validation tab** — the honest headline: top ΔΔG predictions are
   ~3× enriched for real resistance mutations (p < 0.001), but this is a *coarse
   triage flag, not a quantitative fold-resistance predictor*.

A **custom SMILES** is the optional live path — it only runs when a docking
backend is configured (local stack or a `RESISTSCOPE_DOCKING_URL` GPU worker);
otherwise it returns a clean 503 with how-to-enable guidance, never an error.

**Smoke test** (confirms the app boots and serves cached results, in seconds):

```bash
python -m uvicorn api.main:app --port 8000    # in one terminal
python scripts/smoke_test.py                  # checks /health, /targets, /benchmark, /drug/DRV
```

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

## Results — two findings

ResistScope reports two independent validations: does the *prediction* hold up,
and does the *explanation* hold up? They are separate questions with separate
ground truth, and we report both — including where each falls short.

### Finding 1 — Predictive benchmark (does the score hold up?)

Validated against measured Rhee fold-resistance across 6 PIs (1,591 drug-mutation
pairs). Every number carries a significance test or confidence interval
([scripts/08_benchmark.py](scripts/08_benchmark.py): permutation p-values,
bootstrap CIs, ROC/PR-AUC):

| metric | result |
|---|---|
| Top-40 ΔΔG **DRM enrichment** | **2.86×** (95% CI 1.63–4.08, permutation p < 0.001) — top predictions recover known DRMs |
| Pooled **DRM-recovery ROC-AUC** | **0.51** (95% CI 0.46–0.56) — ≈ chance as a *global* ranker |
| Darunavir major-DRM Spearman ρ | **≈ 0.40** (p = 0.03) — strong per-drug signal |
| Per-mutation Spearman ρ (pooled) | **≈ 0.00** — the measured target is confounded by co-occurring mutations |
| **De-confounding** (single-/≤2-mutation isolates) | does **not** rescue the magnitude correlation (ρ ≤ 0) — an honest bound |

The honest headline: rigid single-mutation docking ΔΔG is a **coarse DRM-triage
flag, not a quantitative resistance predictor** — its extreme predictions are
significantly enriched for real resistance mutations (≈3×, p < 0.001) even though
it ranks at chance overall, and de-confounding the target does not rescue the
magnitude correlation. Darunavir validates well per-drug.

### Finding 2 — Explanation faithfulness (does the *reasoning* hold up?)

The tool doesn't just score — it explains, per mutation, *why* a residue change
threatens the compound. That explanation is itself something you can validate. We
check whether Claude's mechanistic rationale matches the biology experts already
documented, rather than just sounding plausible. This is the finding that's
distinctly ours; most predictive-docking work never asks it.

- **Expert ground truth.** [data/mechanism_ground_truth.json](data/mechanism_ground_truth.json)
  — curated resistance mechanisms sourced from **IAS-USA mutation lists** and the
  primary structural literature (e.g. V82A → *"IAS-USA 2022; King et al. 2004"*).
  It is extended agentically: [scripts/09](scripts/09_build_ground_truth.py) +
  [literature_agent.py](services/literature_agent.py) research each mechanism over
  **live PubMed**, cite the papers, and self-verify grounding before an entry is
  admitted.
- **The check.** [scripts/07](scripts/07_faithfulness_eval.py) scores each cached
  explanation against that ground truth with **Claude as judge** on a 0–2 rubric:
  **0** = contradicts the known mechanism, **1** = consistent but vague, **2** =
  identifies the correct primary mechanism. Results split by provenance
  (hand-curated vs. agent-built) and land in
  [data/validation/faithfulness_scores.parquet](data/validation/faithfulness_scores.parquet).

| judge score | count (n = 46) | share |
|---|---|---|
| **2** — correct primary mechanism | 33 | **72%** |
| **1** — consistent but vague | 10 | 22% |
| **0** — contradicts known biology | 3 | 7% |
| Mean faithfulness | | **1.65 / 2** |

The takeaway: the LLM explanations recover the expert-annotated primary mechanism
**72%** of the time, and most of the remainder (10 of 13 misses) under-specify
rather than contradict — but **3 explanations do get the mechanism wrong**, which
is exactly why the faithfulness gate exists and why the score is reported openly
rather than assumed.

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
- **Quantitative validation is HIV-1 protease only.** Triage runs on any target
  (RT, BYO influenza NA / SARS-CoV-2 Mpro), but rigorous fold-resistance
  validation only exists where a measured resistance dataset does — so far, the
  Rhee/Stanford PI data.
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
