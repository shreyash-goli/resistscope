# Red-team pass — ResistScope paper (self-adversarial review)

Reviewing as a skeptical AI4D3/MLSB reviewer. Ordered by severity. Each item:
the attack → the fix. "Quick" = textual/honesty edit; "Work" = new experiment/analysis.

## CRITICAL — a reviewer could reject or heavily criticize

**C1. "Docking energy is noise" is confounded with *crude* docking.** [Work/Quick]
Both the protease chance-level ΔΔG and the field-level "drop_ΔΔG does nothing" use rigid,
single-pose, waters-stripped Vina scoring. A reviewer will say: your ΔΔG is noise *because the
protocol is crude*, not because binding energy is uninformative in principle. The claim "its energy
is noise, its geometry is the signal" over-generalizes from one weak protocol.
→ Fix (quick): scope every such claim to "rigid single-structure Vina ΔΔG (this protocol)"; state
explicitly the null may be protocol-limited; the ensemble/flexible-docking future-work becomes the
direct test of this. (Work: the predictive-lift experiment would resolve it — deferred by choice.)

**C2. LLM-judge circularity — the #1 vulnerability, worst for RT (triple-Claude).** [Work]
RT faithfulness = Claude explanations, judged by Claude, against *Claude-agent-built* ground truth.
Cross-model checks (GPT, Sonnet) establish *consistency*, not *correctness* — all are LLMs that may
share biases, and RT's ground truth is itself LLM-generated. Without any non-LLM signal, a reviewer
may treat §4.2 absolute numbers as uninterpretable.
→ Fix: (a) lead the faithfulness headline with the *curated* protease ground truth; flag RT GT as
agent-built and de-emphasize its absolute %. (b) **Reconsider the skipped human anchor** — even a
~15-item expert spot-check is the one thing that breaks the loop; it directly answers the biggest
criticism. Given this is the top weakness, the cost/benefit of the human anchor has shifted.

**C3. The central "prospective value" claim is asserted, never tested.** [Work/Quick]
"Docking is the only method applicable to a novel compound with no clinical data" is the paper's
answer to "why use docking when HIVdb wins?" — but no experiment shows prospective generalization to a
novel scaffold. It's the load-bearing justification and it's an argument, not a result.
→ Fix (quick): explicitly label it motivation/hypothesis, not a finding. (Better, Work): a
leave-one-compound-out / scaffold holdout where docking recovers a held-out drug's DRMs using *zero*
clinical data for that drug — the per-drug ROC (RT 0.63–0.77) is close; reframe it as the
prospective-per-drug evidence and say so.

**C4. Ignored the supervised-ML-on-resistance literature.** [Quick]
A substantial line trains ML on the Stanford/Rhee genotype→phenotype data and achieves high accuracy
(e.g., work from the Shafer/Rhee group and others). A reviewer will ask why these strong, obvious
baselines aren't cited or compared.
→ Fix: cite them in Related Work and position docking as *training-free / prospective* — it needs no
resistance dataset, which supervised ML fundamentally does. This turns the omission into the paper's
differentiator and reinforces C3. (Add specific citations — do not fabricate.)

## SHOULD-FIX — reviewers will ding these

**S1. Overclaim "strong" for RT ROC-AUC 0.70.** [Quick] 0.70 is moderate. → "genuinely
above-chance / moderate," in abstract + §4.1.

**S2. De-confounding ρ has no CI; the CI is wide.** [Quick] Measured: single-mut ρ=0.36,
**95% CI [0.12, 0.56]**, p=6e-4, leave-one-out [0.34, 0.40] over 30 mutations. → Report it as
"ρ=0.36 [0.12–0.56], robust to leave-one-out" — discloses the width *and* the robustness. Turns an
attackable claim into a rigor point.

**S3. Dose-response magnitude is headroom-confounded.** [Quick] Mean gain +1.23 (minimal=0) vs
+0.62 (minimal=1) is partly mechanical — minimal=0 has 2× the headroom. → Lead with the
ceiling-robust statistic, P(context improves) = **85% vs 66%**, and cite the minimal<2 Spearman
(−0.40) as support; de-emphasize the magnitude numbers.

**S4. "85% rescue" headline rests on n=13.** [Quick] minimal=0 bin is n=13 (11/13). → State the n
explicitly; don't let the headline float without it.

**S5. Field-level ablation is underpowered.** [Quick] distance p=0.051, subpocket p=0.09 (not
<0.05); n=46, protease only; "drop_ΔΔG does nothing" is absence-of-evidence. → Soften "the lift is
geometric" to "suggests the lift is geometric (marginal, n=46, protease only)"; flag replication.

**S6. Universal thesis from 2 HIV targets.** [Quick] "Structure helps where shortcuts fail" is
generalized from n=2, both HIV. → Frame as a hypothesis supported by two contrasting targets, not a
law; note more enzymes/targets needed.

**S7. Thesis softness — structure loses to HIVdb everywhere on known data.** [Quick] The value is
*applicability* (prospective, training-free), not superior accuracy. → Align abstract/title/discussion
to that honest claim (ties to C3/C4).

## MINOR / polish

- **M1.** Protease 2.86× enrichment (p<0.001) vs ROC 0.51 looks contradictory → one sentence:
  extreme-tail enrichment coexists with chance-level *global* ranking.
- **M2.** Related work thin (4 refs) → add substrate-envelope (Schiffer lab), LLM-faithfulness eval,
  ML-docking positioning (DiffDock/Boltz-style), and the C4 supervised-ML refs.
- **M3.** Vina score ≠ binding free energy → one caveat sentence.
- **M4.** Corrupted-context derangement can coincide on low-cardinality fields (region); directed
  (antipodal) corruption is a stronger variant — noted TODO.
- **M5.** No multiple-comparison correction mentioned across many p-values → one acknowledging line.
- **M6.** The tool's actual output (0–100 robustness score) is never validated (underpowered, 6+4
  compounds) — already in limitations; a reviewer may still note the headline product is unvalidated.

## Net assessment
Findings are real and honestly framed; the paper's risk is **under-hedging** (S1–S5 are quick fixes)
and **two structural gaps**: the LLM-judge loop (C2) and the untested-but-load-bearing prospective
claim (C3/C4). None require the deferred predictive-lift. The highest-leverage single action is a
small human-expert faithfulness spot-check (C2).
