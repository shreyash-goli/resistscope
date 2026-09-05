# Domain primer for the faithfulness spot-check

Everything you (or an LLM) need to know to grade whether an explanation correctly
describes *why a mutation causes drug resistance*. Prepend this to the grading sheet.

---

## 1. The pathway — HIV replication, and where our two enzymes act

HIV is a **retrovirus**. To reproduce inside a human cell it runs a lifecycle, and
two steps are enzymes we can drug:

1. **Reverse transcription** — HIV's genome is RNA. Early in infection, the enzyme
   **reverse transcriptase (RT)** copies that RNA into DNA so it can hijack the cell.
   *(Target #1 = RT.)*
2. …integration, and the cell mass-produces new viral proteins as one long
   "polyprotein" chain…
3. **Maturation** — right as new virus particles bud off, the enzyme **protease (PR)**
   cuts that long polyprotein into the individual working proteins. Without this cut,
   the new virus is a dud (non-infectious). *(Target #2 = PR.)*

**Drugs** block one of these enzymes → HIV can't replicate. **Resistance** = the virus
mutates the enzyme so the *drug* no longer binds well, but the *enzyme still works*.
That's the whole game: a resistance mutation weakens drug binding while keeping the
enzyme functional.

---

## 2. The two targets

| code | enzyme | what it does | drug class | our structure |
|---|---|---|---|---|
| **PR** | HIV-1 **protease** | cuts the polyprotein during **maturation** | **PIs** (protease inhibitors) | 3OXC |
| **RT** | HIV-1 **reverse transcriptase** | copies viral **RNA → DNA** | **NNRTIs** (here) | 3V81 |

- **PR** is a small **homodimer** (two identical 99-amino-acid halves). PIs sit in its
  **active site** and physically block the groove where the polyprotein would be cut.
  Key spots you'll see named: the **active site**, the **S1/S1′ subpockets** (where the
  drug's bulky groups sit), the **flaps** (two flexible loops that close over the drug),
  and the catalytic **Asp25** residues.

- **RT** is a **heterodimer** (a big subunit p66 + a smaller p51). Two different drug
  classes hit RT, but **we only work with NNRTIs**:
  - **NNRTIs** bind an **allosteric pocket** (the "NNRTI-binding pocket", NNIBP) about
    **10 Å away** from where the copying actually happens. They don't block the active
    site directly — they wedge into this pocket and **lock the enzyme in a shape that
    can't move/catalyze**. Think of jamming a stick in the hinge rather than plugging the
    mouth. Key residues lining this pocket: **L100, K101, K103, V106, V108, V179, Y181,
    Y188, G190, F227, W229, M230** (on p66) and **E138** (contributed by p51).
  - *(The other RT class, NRTIs, works by a totally different chemistry at the active
    site — we deliberately excluded them, because our "does the drug fit the pocket"
    method only makes sense for pocket-binders like NNRTIs and PIs.)*

---

## 3. The drug codes

**PIs (protease inhibitors — used with PR):**

| code | drug |
|---|---|
| ATV | atazanavir |
| DRV | darunavir |
| IDV | indinavir |
| LPV | lopinavir |
| NFV | nelfinavir |
| SQV | saquinavir |

All six bind the **protease active site**. (In the clinic most are "boosted" with a
low dose of ritonavir, written e.g. "ATV/r" — that's a pharmacokinetic booster, not
relevant to the binding mechanism.)

**NNRTIs (non-nucleoside RT inhibitors — used with RT):**

| code | drug | note |
|---|---|---|
| NVP | nevirapine | 1st-generation |
| EFV | efavirenz | 1st-generation |
| ETR | etravirine | 2nd-generation (harder to resist) |
| RPV | rilpivirine | 2nd-generation (harder to resist) |

1st-gen (NVP, EFV) are rigid and easily defeated by a *single* mutation (e.g. K103N,
Y181C). 2nd-gen (ETR, RPV) are flexible and "wiggle" to stay bound, so they resist
resistance better — it takes more mutations to break them.

---

## 4. How to read a mutation

Format: **[normal amino acid] [position] [new amino acid]**.

- **V82A** = at position **82**, the normal residue **V**aline was replaced by **A**lanine.
- **K103N** = position **103**, **K** (lysine) → **N** (asparagine).
- **Y181C** = position **181**, **Y** (tyrosine) → **C** (cysteine).

The **number** is just which residue in the enzyme's chain. Amino acids are the 20
building blocks of proteins; each has a different **size, charge, and chemistry**, so
swapping one for another changes the local shape/charge of the pocket — and *that* is
what can weaken drug binding.

**One-letter amino-acid cheat sheet** (with the traits that matter for binding):

| letter | amino acid | rough character |
|---|---|---|
| A | alanine | tiny, greasy |
| V, L, I | valine / leucine / isoleucine | medium, greasy (hydrophobic) |
| F, Y, W | phenylalanine / tyrosine / tryptophan | big, flat **aromatic rings** (good at stacking against drugs) |
| M | methionine | medium, greasy, flexible |
| G | glycine | smallest, very flexible |
| P | proline | rigid, kinks the backbone |
| S, T | serine / threonine | small, polar (can H-bond) |
| C | cysteine | small, can form special bonds |
| N, Q | asparagine / glutamine | polar (H-bond) |
| D, E | aspartate / glutamate | **negatively charged** |
| K, R, H | lysine / arginine / histidine | **positively charged** |

---

## 5. The mechanism vocabulary (why a mutation → resistance)

Almost every resistance explanation is one of these physical stories. When you grade,
you're checking the explanation names the *right one*:

- **Lost bulk / cavity** — a big residue → small one (e.g. **V82A**, Val→Ala) removes
  material the drug was leaning on, so the drug loses **van der Waals contact** and
  binds more weakly. Also **Y181C / Y188L** (big aromatic ring → small residue) — the
  drug loses the **aromatic stacking** it relied on.
- **Steric clash** — a small residue → bigger one (e.g. **G190A**, Gly→Ala adds a methyl)
  bumps into the drug, pushing it out.
- **Lost hydrogen bond / charge change** — swapping a charged/polar residue (e.g.
  **K103N**, or **E138K** flipping negative→positive) removes an electrostatic
  interaction or a salt bridge that was holding things in place.
- **Allosteric / indirect** — the residue doesn't touch the drug directly (it's far
  away, e.g. **E138** on the *other* subunit, ~40 Å from an NNRTI) but its change
  reshapes the pocket or the enzyme's motion, weakening binding indirectly. These are
  the subtle ones.

For PIs, the drug's parts that sit in the pockets are called **P1/P1′, P2/P2′** groups —
you'll see explanations mention the drug losing contact with a subpocket.

---

## 6. What you're actually grading

Each item gives you:
- **Known mechanism** (ground truth): the established, literature-backed reason this
  mutation causes resistance.
- **AI explanation**: what our model wrote.

Score how well the **AI explanation matches the known primary mechanism**:
- **2** = names the correct primary structural mechanism (same physical story as above).
- **1** = consistent/plausible but vague, or misses the key interaction.
- **0** = contradicts the known mechanism (e.g. calls a real resistance mutation
  "negligible", or invents a different mechanism).

Judge the **mechanism**, not the writing quality or the exact numbers. If the known
mechanism is "loses aromatic stacking" and the AI says "loses van der Waals contact from
the aromatic ring," that's a **2** (same story). If the AI says "no effect, too far
away" but the known mechanism is a real salt-bridge disruption, that's a **0**.
