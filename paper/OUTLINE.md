# Paper outline + abstract (working draft)

**Venue:** NeurIPS 2026 AI4D3 workshop (non-archival, WIP-friendly). **Deadline: Aug 29, 2026 (AoE).**
**Format:** ~4 pages + refs. Concurrent submission OK.

---

## Title (options)

1. **When Structure Helps: Target-Dependent Docking and Prior-Dependent Explanation for HIV Drug-Resistance Triage** *(lead)*
2. Knowing When to Trust Structure: Two Regimes for Structure-Based Antiviral Resistance Prediction and Interpretation
3. ResistScope: Structure-Based HIV Resistance Triage — and an Honest Account of When It Works

---

## Abstract (~230 words)

Antiviral drug resistance forces medicinal chemists to triage candidate inhibitors by how well their
binding survives known resistance mutations — a slow, manual judgment. We present **ResistScope**, a
target-agnostic pipeline that docks a candidate against a wildtype receptor and a panel of clinical
resistance mutants, scores binding robustness (ΔΔG), and generates a per-mutation mechanistic
explanation with an LLM. Rather than report a single headline number, we ask two questions that have
real ground truth. **First, does the docking score predict measured fold-resistance?** Across HIV-1
**protease** (6 PIs) and **reverse transcriptase** (4 NNRTIs) using Stanford HIVdb data, the answer is
**target-dependent**: rigid single-structure ΔΔG is a *chance-level* ranker for the flexible,
network-driven protease active site (ROC-AUC 0.51 — below a zero-cost mutation-prevalence baseline of
0.72), but a *strong* one for the rigid, contact-driven NNRTI pocket (DRM enrichment 8.35×, ROC-AUC
0.70; de-confounded single-mutation Spearman ρ≈0.4). **Second, do the explanations reflect real
biology?** Judged against IAS-USA and agentically-curated PubMed mechanisms, 71–72% recover the
correct primary mechanism. An **attribution ablation** (real vs absent vs corrupted structural context)
shows the explanations genuinely condition on the pipeline; across both targets, the structural
context's benefit **scales inversely with the model's prior knowledge** — large for obscure, distal
mutations, negligible for textbook ones. The unifying finding: structure-based methods, predictive and
interpretive, deliver value precisely where non-structural shortcuts fail.

---

## Contributions

1. **An honest, two-target benchmark** of single-structure docking-ΔΔG for HIV resistance, with
   significance-tested baselines — showing the method's value is **target-dependent** (works for the
   NNRTI pocket, fails for protease) rather than uniform.
2. **A faithfulness evaluation of the LLM's mechanistic explanations** against expert + agentic ground
   truth — most docking work never asks whether the *reasoning* matches known biology.
3. **An attribution ablation** isolating whether faithfulness comes from the structural pipeline or the
   model's prior, revealing a **prior-strength dose-response**: structural grounding helps exactly where
   the model's memorized knowledge is weak.
4. A reproducible, **target-agnostic** implementation (protease, RT, bring-your-own-target).

---

## 1. Introduction
- Resistance is the central failure mode for antivirals; triaging N candidates into slow phenotypic
  assays is a real bottleneck (Gladstone framing).
- Two questions usually asked separately, both with ground truth: *does the prediction hold?* (measured
  fold-resistance), *does the explanation hold?* (documented mechanisms).
- Thesis: **structure helps where non-structural shortcuts fail** — for prediction, the shortcut is
  prevalence/expert lists; for explanation, the shortcut is the LLM's memorized biology.
- Contributions (above).

## 2. Related work
- Structure-based resistance / docking rescoring; substrate-envelope hypothesis (protease).
- Rule-based expert systems: Stanford HIVdb penalty scores; IAS-USA mutation lists. *(baseline)*
- LLMs in molecular biology; **faithfulness / explanation evaluation** (the gap we fill).
- Rhee et al. 2006 genotype–phenotype datasets.

## 3. Methods
- **Pipeline:** SMILES → 3D (RDKit/meeko) → AutoDock Vina / Uni-Dock vs wildtype + mutant receptors;
  per-mutation ΔΔG; 0–100 robustness score. Receptor prep (PDBFixer point mutagenesis + clash-relief
  minimization); protease homodimer + Asp25 dyad; RT p66 heterodimer, NNRTI pocket box.
- **Targets & data:** PR (3OXC, 6 PIs, 285 mutants); RT (3V81, 4 NNRTIs, 424 mutants). Measured
  fold-resistance = Rhee/Stanford HIVdb. DRM labels = IAS-USA.
- **Explanations:** structural context (Δvolume/charge, ligand distance, subpocket) + ΔΔG → LLM
  (Haiku 4.5; larger models refuse HIV-resistance prompts) → 2–3-sentence mechanism.
- **Ground truth for mechanisms:** curated (IAS-USA + primary lit) + an **agentic PubMed builder** that
  researches, cites, and self-verifies each mechanism.
- **Faithfulness eval:** LLM-as-judge, 0/1/2 rubric vs ground-truth mechanism.
- **Attribution ablation:** regenerate explanations under `full` (real context+ΔΔG) / `minimal`
  (drug+mutation only) / `corrupted` (deranged pipeline fields, intrinsic chemistry kept); same judge.
- **Baselines (predictive):** mutation prevalence, physicochemical (|Δvolume|), **Stanford HIVdb
  penalty score** *[TODO: fetch + integrate]*.

## 4. Results

### 4.1 Predictive benchmark is target-dependent  (Fig 2, Table 1)
- Every claim carries permutation p / bootstrap CI.
- **Protease:** top-40 DRM enrichment 2.86× (p<0.001) but pooled ROC-AUC **0.510** (CI 0.465–0.557,
  chance); per-mutation Spearman ≈0; de-confounding does **not** rescue (ρ≤0). **Loses to baselines:**
  mutation-prevalence DRM-recovery AUC **0.719**; |Δvolume| magnitude ρ 0.14 vs docking 0.00.
- **RT:** top-40 enrichment **8.35×** (CI 6.7–9.7); pooled ROC-AUC **0.695** (CI 0.62–0.76, above
  chance); **de-confounding rescues** — confounded ρ −0.11 → single-/≤2-mut ρ +0.36/+0.41 (p<0.001).
  WT ΔG −9 to −12 kcal/mol (box validated).
- **Read:** rigid ΔΔG tracks resistance for the compact, contact-driven NNRTI pocket; fails for the
  flexible, co-evolving, envelope-constrained protease active site. Pocket biophysics predicts
  trustworthiness.
- **Table 1 — baseline head-to-head** (DRM-recovery ROC-AUC / magnitude ρ vs measured fold-resistance):

  | predictor | PR ROC / ρ | RT ROC / ρ |
  |---|---|---|
  | docking ΔΔG | 0.51 / 0.00 | **0.70** / −0.11 (de-conf +0.4) |
  | prevalence | 0.72 / −0.08 | 0.65 / 0.13 |
  | \|Δvolume\| | 0.49 / 0.14 | 0.65 / −0.01 |
  | **HIVdb penalty** | **0.86 / 0.36** | **0.93 / 0.27** |

  The expert system (HIVdb) is strongest on both — **docking never beats it retrospectively**; on PR
  docking loses even to zero-cost prevalence, on RT it beats the docking-free baselines. **Docking's
  value is prospective:** HIVdb/prevalence need clinical/expert data a novel compound or target lacks;
  docking is the only method applicable with zero clinical history. *(HIVdb-recovery ROC is partly
  circular vs the IAS-USA DRM labels; the magnitude-ρ column is the fair comparison.)*
- **Bootstrap 95% CIs** (2000 resamples): PR docking ROC 0.51 **[0.46,0.56]** (includes chance); RT
  docking 0.70 **[0.62,0.77]** (excludes chance); HIVdb tight/high on both (PR 0.86 [0.83,0.90], RT 0.93
  [0.90,0.96]); docking magnitude ρ PR [−0.05,+0.05], RT [−0.16,−0.05].
- **Leave-one-drug-out jackknife:** dropping any single drug leaves PR at chance (ROC [0.50,0.53]) and
  RT above chance (ROC [0.67,0.74]) — the target-dependent finding is not driven by any one drug.

### 4.2 Explanation faithfulness  (Fig 3a)
- Recovers the correct primary mechanism **72%** (protease, n=46) / **71%** (RT, n=55) *under the
  Claude-Haiku judge*; most misses under-specify rather than contradict.
- RT misses cluster on **distal/allosteric** DRMs (E138–K101 inter-subunit salt bridge, A98G) where
  single-structure docking is genuinely blind — the model faithfully reports "no direct effect."
- **Judge robustness (cross-model), two checks.** *(a) Absolute level, hard sample:* GPT-5.6 vs
  Claude-Haiku on a 21-item stress sample (over-weighted with 0/1 cases) — 95% within-1 but 48% exact
  (κ 0.23/0.38 unwtd/lin-wtd), Claude systematically stricter (all disagreements GPT-higher). *(b) The
  contrasts, full ablation set:* re-grading all **138 protease ablation items** with a second model
  (Sonnet 5) **reproduces the §4.3 contrasts almost exactly** — full−minimal Δ+0.35 (p=0.001),
  full−corrupted Δ+0.65 (p<1e-4), vs Haiku's +0.43/+0.63 — with strong Sonnet–Haiku agreement
  (exact 78%, within-1 100%, κ 0.64/0.72). **⇒ the paired ablation contrasts are judge-independent;
  only the *absolute* faithfulness % is judge-dependent** (strict Haiku ~72% ↔ lenient GPT). Cross-vendor
  disagreement concentrates on subtle accessory/allosteric mutations. *(Human-expert anchor: future work
  — 11-case adjudication sheet prepared.)*

### 4.3 Attribution ablation & the prior-strength dose-response  (Fig 3b — headline)
- **Corruption degrades faithfulness on both targets** (PR full-vs-corrupted Δ+0.63 p<0.001; RT Δ+0.31
  p=0.006) → explanations genuinely condition on the structural input, not just the mutation name.
- **Dose-response (both targets, split by ligand distance = proxy for how textbook the DRM is):**
  the structural context's benefit scales **inversely with the model's prior** (minimal-condition score):

  | target · subset | minimal | full | full − minimal |
  |---|---|---|---|
  | PR · distal (n=19) | 0.74 | 1.47 | **+0.74, p=0.004** |
  | RT · distal (n=30) | 1.13 | 1.50 | **+0.37, p=0.036** |
  | PR · contact (n=27) | 1.78 | 2.00 | +0.22, p=0.014 |
  | RT · contact (n=25) | 1.92 | 1.76 | −0.16, n.s. |

- **Read:** structural grounding helps LLM explanations *exactly where the model's memorized knowledge
  is weakest* — the prospective setting (novel/obscure mutations & targets).
- **Formalized, ceiling-controlled (pooled n=101):** full>minimal Δ+0.27 (p=7e-4), full>corrupted
  Δ+0.46 (p<1e-4). Gain by prior: minimal=0 → **+1.23 (85% improved)**, =1 → +0.62 (66%), =2 → −0.16
  (ceiling). Among pairs with room to improve (minimal<2, n=45): Δ+0.80 (p<1e-4), gradient persists
  (Spearman ρ=−0.40, p=0.006) — not merely headroom. **Headline stat: when the model's prior alone
  gets the mechanism wrong, structural context corrects it 85% of the time.**
- **Which field carries it? (field-level ablation, protease n=46, drop one context group):** the lift is
  **geometric, not energetic or chemical.** Dropping **ligand distance** costs the most (Δ+0.22,
  p=0.051), **subpocket** next (Δ+0.15, p=0.09); dropping the **docking ΔΔG** (Δ−0.02, p=0.78) or the
  physicochemistry (Δ−0.02, p=0.76) does nothing. The pipeline's explanatory value is its *geometric
  localization* of the mutation vs the ligand — consistent with ΔΔG being a chance-level predictor (its
  **energy is noise, its geometry is the signal**). Unifies both findings: across prediction and
  explanation, structural *geometry* carries value, not the docking *energy*. *(Marginal p at n=46;
  RT / low-prior subset would confirm.)*

## 5. Discussion
- **Unifying thesis:** structure-based methods pay off where the non-structural shortcut fails —
  network/allosteric resistance for prediction, unfamiliar mutations for explanation.
- **When to trust structure-based resistance triage:** rigid contact-driven pockets, single-mutation
  regimes; distrust for flexible/allosteric/co-evolving sites.
- **Limitations:** rigid single-structure docking (no ensemble/flap-water/flexible receptor);
  crystal waters stripped; RT subpocket map is a coarse approximation; RT mechanism ground truth is
  agent-built (needs human spot-check); **LLM-as-judge is only fairly reproducible across models**
  (GPT-5.6 vs Claude κ≈0.23–0.38, Claude stricter) so the absolute faithfulness % is judge-dependent —
  the paper leans on paired within-judge contrasts, which are offset-invariant, and still needs a
  human-expert anchor; modest per-bin n; compound-level scoring aggregation
  (mean/weighted/worst robustness score) is underpowered (6 PI + 4 NNRTI compounds) — reported
  descriptively only, so all quantitative validation is at the per-mutation level.

## 6. Conclusion
- Two honest findings + a unifying account of when structure helps; a reproducible target-agnostic tool.

---

## Figure / table plan
- **F1** Pipeline schematic (SMILES → dock vs panel → ΔΔG + robustness + explanation).
- **F2** Predictive benchmark, PR vs RT: enrichment bars + ROC + de-confounding, with baselines.
- **Table 1** Baseline head-to-head (docking / prevalence / physicochem / HIVdb), both tasks, both targets.
- **F3a** Faithfulness distribution (0/1/2), PR & RT.
- **F3b** Prior-strength dose-response (the headline interpretability figure).

## Open items (map to rigor todos)
- [x] HIVdb penalty-score baseline → Table 1 (done — PR 0.86/0.36, RT 0.93/0.27; via Sierra GraphQL, cached)
- [x] bootstrap CIs on Table 1 + leave-one-drug-out jackknife (done)
- [ ] Dose-response regression (ceiling-aware) → F3b
- [ ] Field-level ablation (drop ΔΔG-only / distance-only / subpocket-only)
- [ ] Human spot-check of the judge + 2nd judge / IRR (esp. RT agent-built GT)
- [ ] Leave-one-drug-out CV; scoring ablation (mean/weighted/worst)
- [ ] Fix README DRV metric inconsistency
