# When Structure Helps: Target-Dependent Docking and Prior-Dependent Explanation for HIV Drug-Resistance Triage

*Working draft for NeurIPS 2026 AI4D3 (non-archival). ~4 pages. All numbers are current as of the analysis in `paper/OUTLINE.md`.*

## Abstract

Antiviral drug resistance forces medicinal chemists to triage candidate inhibitors by how well their binding survives known resistance mutations — a slow, manual judgment. We present **ResistScope**, a target-agnostic pipeline that docks a candidate against a wildtype receptor plus a panel of clinical resistance mutants, scores binding robustness (ΔΔG), and generates a per-mutation mechanistic explanation with an LLM. Rather than report a single headline number, we ask two questions that have real ground truth. **First, does the docking score predict measured fold-resistance?** Across HIV-1 protease (6 protease inhibitors) and reverse transcriptase (4 NNRTIs) using Stanford HIVdb data, the answer is **target-dependent**: rigid single-structure ΔΔG is a chance-level ranker for the flexible, network-driven protease active site (ROC-AUC 0.51 [0.46, 0.56], below a zero-cost mutation-prevalence baseline of 0.72), but a strong one for the rigid, contact-driven NNRTI pocket (DRM enrichment 8.35×, ROC-AUC 0.70 [0.62, 0.77]; de-confounded single-mutation Spearman ρ ≈ 0.4). **Second, do the explanations reflect real biology?** An LLM judge, scoring against IAS-USA and agentically-curated PubMed mechanisms, marks 71–72% as recovering the correct primary mechanism (a judge-dependent absolute; the ablation contrasts below are judge-robust across three models). An **attribution ablation** (real vs. absent vs. corrupted structural context) shows the explanations genuinely condition on the pipeline; across both targets the structural context's benefit **scales inversely with the model's prior knowledge** — large for obscure, distal mutations, negligible for textbook ones — and a field-level ablation localizes the benefit to the mutation's *geometry* (ligand distance, subpocket), not the docking energy. The unifying finding: structure-based methods, predictive and interpretive, deliver value precisely where non-structural shortcuts fail.

## 1. Introduction

Resistance is the central failure mode of antiviral therapy. A medicinal chemist with twenty candidate inhibitors must decide which few to advance into slow, expensive phenotypic resistance assays — a decision currently made by hand, mutation by mutation. A natural computational aid is to dock each candidate against a panel of known clinical resistance mutants and rank candidates by how well predicted binding *survives* the panel.

Such a tool raises two questions that are usually asked separately, and that both have ground truth:

1. **Does the prediction hold up?** Stanford HIVdb records *measured* fold-resistance — how many-fold less potent a drug becomes against a given mutant — for essentially every approved HIV drug. This turns "is the score reasonable?" into "does the score predict a measured number, and does it beat the trivial baselines a chemist could use instead?"
2. **Does the *explanation* hold up?** If the tool also explains *why* a mutation threatens a compound, that explanation can be checked against mechanisms experts have already documented. Most structure-based resistance work never asks whether the model's *reasoning* — not just its score — matches known biology.

Our central finding, across both questions, is that **structure-based methods pay off exactly where the non-structural shortcut fails.** For prediction, the shortcut is mutation prevalence or an expert rule-base; structure adds value only for targets whose resistance is structurally local. For explanation, the shortcut is the LLM's memorized biology; structure adds value only for mutations the model hasn't memorized.

**Contributions.** (i) An honest, two-target benchmark of single-structure docking-ΔΔG for HIV resistance with significance-tested baselines, showing the method's value is *target-dependent*. (ii) A faithfulness evaluation of the LLM's mechanistic explanations against expert and agentic ground truth. (iii) An attribution ablation isolating whether faithfulness comes from the structural pipeline or the model's prior, revealing a *prior-strength dose-response* and localizing the effect to geometry. (iv) A reproducible, target-agnostic implementation (protease, RT, bring-your-own-target).

## 2. Related work

**Structure-based resistance prediction and docking rescoring** have a long history; the substrate-envelope hypothesis in particular explains why protease resistance resists simple pocket-affinity models. **Rule-based expert systems** — the Stanford HIVdb penalty algorithm and IAS-USA mutation lists — are the clinical standard and our strongest baseline. Separately, **LLMs are increasingly used for molecular reasoning**, but their *faithfulness* — whether generated mechanisms match established biology — is rarely evaluated in this domain; that gap is our second contribution. Genotype–phenotype resistance data come from Rhee et al. (2006) and the Stanford HIVdb.

## 3. Methods

**Pipeline.** A candidate SMILES is embedded (RDKit/meeko) and docked (AutoDock Vina on CPU, or Uni-Dock on GPU) against a cleaned wildtype receptor and a per-drug panel of point-mutant receptors; the per-mutation ΔΔG (mutant minus wildtype binding) feeds a 0–100 robustness score. Mutant receptors are built by PDBFixer point mutagenesis with a frozen-environment clash-relief minimization of the mutated side chain.

**Targets and data.** *Protease* (PR): PDB 3OXC, homodimer, catalytic Asp25 dyad asymmetrically protonated; 6 PIs (ATV, DRV, IDV, LPV, NFV, SQV); 285 mutant receptors. *Reverse transcriptase* (RT): PDB 3V81, p66/p51 heterodimer, NNRTI-pocket docking box centered on the chain-A nevirapine ligand (41.1, 52.3, 49.1); 4 NNRTIs (NVP, EFV, ETR, RPV); 424 mutant receptors. We scope RT to NNRTIs (allosteric pocket binders) — rigid docking of the pocket is a sound proxy for NNRTI affinity but not for NRTI mechanisms. Measured fold-resistance and DRM labels come from Stanford HIVdb / IAS-USA.

**Explanations.** For each mutation we assemble a structural-context record (residue identity, Δvolume/charge/hydrophobicity, distance to the ligand, direct-contact flag, subpocket) plus the docking ΔΔG, and prompt an LLM (Claude Haiku 4.5; larger models decline HIV-resistance prompts via a safety classifier) for a 2–3 sentence mechanistic hypothesis.

**Ground truth for mechanisms.** A curated set (IAS-USA + primary literature) plus an *agentic* builder that researches each mechanism over live PubMed, cites the papers, and self-verifies grounding before admitting an entry (RT ground truth: 22 mechanisms).

**Faithfulness evaluation.** An LLM-as-judge scores each explanation against the ground-truth mechanism on a 0/1/2 rubric (0 = contradicts, 1 = consistent but vague, 2 = correct primary mechanism). We assess judge robustness with a second, independent model (§4.2).

**Attribution ablation.** We regenerate the ground-truth explanations under three context conditions — *full* (real context + ΔΔG), *minimal* (drug + mutation identity only), *corrupted* (production-shaped context but with the pipeline-derived fields deranged across pairs, intrinsic chemistry kept) — with the same model and judge. A *field-level* variant drops one context group at a time.

**Predictive benchmark.** Every claim carries a permutation p-value or bootstrap 95% CI: top-N DRM enrichment, threshold-free ROC/PR-AUC for DRM recovery, a de-confounding analysis (rebuild the fold-resistance target from progressively less co-mutated isolate subsets), a leave-one-drug-out jackknife, and baselines (mutation prevalence, physicochemical |Δvolume|, and the Stanford HIVdb penalty score fetched via the Sierra GraphQL API).

## 4. Results

### 4.1 The predictive benchmark is target-dependent

Rigid single-structure docking-ΔΔG behaves oppositely on the two enzymes.

**Protease — chance-level, loses to trivial baselines.** Top-40 DRM enrichment is 2.86× (95% CI 1.63–4.08, permutation p < 0.001), but the pooled DRM-recovery ROC-AUC is only 0.51 [0.46, 0.56] and per-mutation Spearman correlation with measured fold-resistance is ≈ 0. De-confounding does *not* rescue it (ρ ≤ 0 on single-/few-mutation isolate subsets).

**RT — a genuine ranker, and de-confounding rescues.** Top-40 enrichment is 8.35× (95% CI 6.7–9.7, p < 0.001); pooled ROC-AUC is 0.70 [0.62, 0.77], comfortably above chance; and the pooled magnitude correlation goes from −0.11 (confounded) to +0.36…+0.41 (p < 0.001) on single-/≤2-mutation isolates — the co-mutation confound was masking a real single-mutation signal. Wildtype absolute affinities (−9 to −12 kcal/mol) confirm the docking box is on the pocket.

**Table 1 (baseline head-to-head; DRM-recovery ROC-AUC / magnitude ρ vs. measured fold-resistance).**

| predictor | protease | RT |
|---|---|---|
| docking ΔΔG | 0.51 / 0.00 | **0.70** / −0.11 (de-conf +0.4) |
| mutation prevalence | 0.72 / −0.08 | 0.65 / 0.13 |
| \|Δ side-chain volume\| | 0.49 / 0.14 | 0.65 / −0.01 |
| **Stanford HIVdb penalty** | **0.86 / 0.36** | **0.93 / 0.27** |

The expert system (HIVdb) is strongest on both, and **docking never beats it retrospectively**: on protease docking loses even to zero-cost prevalence; on RT it beats the docking-free baselines but not HIVdb. A leave-one-drug-out jackknife confirms the pattern is not driven by any single drug (protease ROC range [0.50, 0.53]; RT [0.67, 0.74]). *(HIVdb's DRM-recovery ROC is partly circular — its penalties and the IAS-USA labels are both expert-derived — so the magnitude-ρ column is the fair comparison.)*

**Why the difference, and why docking is still worth running.** NNRTIs bind a compact, rigid, hydrophobic pocket where resistance is direct steric/contact disruption — exactly what rigid ΔΔG captures. Protease resistance is flexible-active-site, co-evolving, and substrate-envelope-constrained — what rigid docking misses. So the pocket's biophysics predicts where structure-based prediction is trustworthy. Crucially, the baselines that beat docking (prevalence, HIVdb) require *clinical/expert data that a novel compound or novel target does not have*; docking is the only method applicable to a candidate with zero clinical history. Its value is prospective, not retrospective.

### 4.2 Explanation faithfulness

Under the Claude-Haiku judge, explanations recover the correct primary mechanism 72% of the time on protease (n = 46) and 71% on RT (n = 55); most misses under-specify rather than contradict. On RT the outright misses cluster on distal/allosteric DRMs — e.g. the E138–K101 inter-subunit salt bridge, A98G — where single-structure docking is genuinely blind, and the model faithfully reports "no direct effect."

**Judge robustness (cross-model).** On a 21-item sample stratified to over-represent the judge's hard (0/1) cases, an independent judge (GPT-5.6) agrees with Claude-Haiku 95% within one point but only 48% exactly (Cohen's κ 0.23 unweighted, 0.38 linear-weighted); Claude is systematically *stricter* (every disagreement is the other model scoring higher). We therefore treat the *absolute* faithfulness percentage as **judge-dependent** and rest our headline claims on the *paired, within-judge* contrasts of §4.3, where a uniform strictness offset cancels. Disagreements concentrate on subtle accessory/allosteric mutations; both judges agree on canonical DRMs. Critically, re-grading all 138 protease ablation items with a *second* model (Sonnet 5) reproduces the §4.3 contrasts almost exactly (full−minimal Δ = +0.35, p = 0.001; full−corrupted Δ = +0.65, p < 1e-4, vs. Haiku's +0.43/+0.63; Sonnet–Haiku agreement κ = 0.64/0.72, within-1 100%) — the paired contrasts are judge-independent even though the absolute level is not. A human-expert anchor is future work.

### 4.3 Attribution ablation: the pipeline drives faithfulness, geometrically, where the prior is weak

**Context is genuinely used.** Pooled across both targets (n = 101 pairs), full context beats minimal (drug + mutation only) by Δ = +0.27 (p = 7e-4) and beats corrupted context by Δ = +0.46 (p < 1e-4). Feeding *wrong* structural facts is worse than feeding none — the explanations condition on the pipeline, they do not merely pattern-match the mutation name.

**Prior-strength dose-response.** The benefit of correct context scales inversely with the model's prior knowledge (proxied by the minimal-condition score). Where the prior alone is *wrong* (minimal = 0), context corrects the explanation **85%** of the time (gain +1.23); where the prior is vague (minimal = 1), gain +0.62; where the prior is already correct (minimal = 2), context is at ceiling. Restricting to pairs with room to improve (minimal < 2, n = 45), context still helps Δ = +0.80 (p < 1e-4) and the gradient persists (Spearman ρ = −0.40, p = 0.006) — it is a behavioral effect, not the ceiling. **Structural grounding helps exactly where the model's memorized knowledge is weakest** — the prospective setting of novel or obscure mutations.

**It's the geometry, not the energy.** A field-level ablation (protease, n = 46) drops one context group at a time. Removing the *ligand distance* costs the most (Δ = +0.22, p = 0.051) and the *subpocket* next (Δ = +0.15, p = 0.09); removing the *docking ΔΔG* (Δ = −0.02, p = 0.78) or the physicochemistry (Δ = −0.02, p = 0.76) does nothing. The pipeline's explanatory value is its geometric localization of the mutation relative to the ligand — consistent with the docking *energy* being a chance-level predictor (its energy is noise; its geometry is the signal). The chemistry is redundant with the mutation string itself.

## 5. Discussion

The two findings share one account: **structure-based methods deliver value where the non-structural shortcut fails.** For prediction, the shortcut is prevalence/expert rules, and structure helps only when resistance is structurally local (the NNRTI pocket, not the protease active site). For explanation, the shortcut is the LLM's parametric biology, and structure helps only for mutations the model has not memorized — and specifically through the mutation's *geometry*. This yields a practical rule for practitioners: trust structure-based resistance triage for rigid, contact-driven pockets and single-mutation regimes; distrust it for flexible, allosteric, co-evolving sites — and use it prospectively, where clinical and expert baselines do not exist.

**Limitations.** Rigid single-structure docking (no ensemble, flexible receptor, or explicit flap water); crystal waters stripped; the RT subpocket map is a coarse approximation of a discontiguous pocket; RT mechanism ground truth is agent-built and needs a human-expert anchor; the LLM judge's *absolute* scores are judge-dependent (cross-vendor κ ≈ 0.23–0.38 on hard cases, cross-Claude-model κ ≈ 0.64–0.72) though the paired contrasts reproduce across judges — a human-expert anchor is still needed; per-bin n is modest; compound-level scoring aggregation is underpowered (6 + 4 compounds), so all quantitative validation is per-mutation.

## 6. Conclusion

An honest two-target study of structure-based HIV resistance triage yields two findings and one account of when structure helps: predictive value is target-dependent (rigid NNRTI pocket yes, flexible protease no), and interpretive value is prior-dependent and geometric (structure grounds explanations exactly where the model's memory is weak). Both say the same thing — structure earns its keep where non-structural shortcuts run out — which is the prospective setting the tool is built for.

## References

*(to fill: Rhee et al. 2006 PNAS; Stanford HIVdb / Sierra; IAS-USA mutation lists; AutoDock Vina — Trott & Olson 2010, Eberhardt et al. 2021; Uni-Dock — Yu et al. 2023; meeko; PDBFixer/OpenMM; RDKit; PDB 3OXC, 3V81.)*
