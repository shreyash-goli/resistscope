# Adjudicate the 11 judge-disagreement cases

The two LLM judges (Claude-Haiku vs GPT-5.6) disagreed on these. Read the **known mechanism** vs the **AI explanation** and give YOUR score (0/1/2) — this anchors which judge is closer to truth.

Rubric: **0**=contradicts known mechanism, **1**=consistent but vague, **2**=correct primary mechanism.

---

## #1  RT · ETR · L100I
*(Claude judge said 1, GPT-5.6 said 2)*

**Known mechanism:** L100I destabilizes NNRTI binding primarily through increased drug dissociation (koff increased ~12-fold) without changing the association rate, indicating loss of stabilizing interactions in the NNRTI pocket rim. The isoleucine substitution introduces steric and/or thermodynamic barriers that weaken van der Waals contacts between the residue and the bound inhibitor, reducing overall nevirapine affinity by ~17-fold.

**AI explanation:** The L100I substitution increases local hydrophobicity at the NNRTI pocket rim while maintaining essentially identical van der Waals volume, which subtly alters the electrostatic and hydrophobic balance of the binding pocket microenvironment despite being 7.6 Å from etravirine. This conformational perturbation at the pocket periphery likely propagates inward to destabilize the precise geometry required for optimal etravirine positioning, reducing favorable aromatic and hydrophobic contacts between the inhibitor and the binding site. The modest energetic penalty (+0.28 kcal/mol) reflects a cumulative effect of this indirect allosteric modulation rather than direct steric clash, consistent with L100I's location on the pocket rim where it can influence drug accommodation without catastrophic structural disruption.

**YOUR score (0 / 1 / 2):** ____

---

## #2  RT · EFV · L100I
*(Claude judge said 0, GPT-5.6 said 2)*

**Known mechanism:** L100I destabilizes NNRTI binding primarily through increased drug dissociation (koff increased ~12-fold) without changing the association rate, indicating loss of stabilizing interactions in the NNRTI pocket rim. The isoleucine substitution introduces steric and/or thermodynamic barriers that weaken van der Waals contacts between the residue and the bound inhibitor, reducing overall nevirapine affinity by ~17-fold.

**AI explanation:** The L100I substitution reduces efavirenz binding through a subtle conformational mechanism rather than direct contact: leucine's longer, more flexible side chain permits accommodation of the NNRTI pocket geometry that positions efavirenz optimally, while isoleucine's branched-at-γ structure creates geometric constraints that shift the local backbone or adjacent side-chain orientations at the pocket rim. This reorientation propagates through the binding pocket architecture, causing efavirenz to adopt a suboptimal binding pose approximately 7.6 Å away, reducing favorable van der Waals interactions across the extended contact surface without introducing steric clashes. The 0.7 kcal/mol hydrophobicity increase and minimal volume change indicate the resistance mechanism is driven by altered spatial fit rather than energetic destabilization, consistent with the modest energetic penalty observed.

**YOUR score (0 / 1 / 2):** ____

---

## #3  PR · ATV · I50L
*(Claude judge said 1, GPT-5.6 said 2)*

**Known mechanism:** Ile50 at the flap tip contacts the inhibitor. I50L reshapes flap-tip packing in a way that selectively reduces atazanavir binding (while often increasing susceptibility to other PIs); an atazanavir signature mutation.

**AI explanation:** The I50L substitution reduces atazanavir binding through a subtle loss of favorable van der Waals contacts in the flap region. Although the volume remains essentially unchanged, leucine's slightly shorter branched side chain (reduced hydrophobicity) creates a small cavity at the 2.2 Å interface with the ligand, weakening the close-range hydrophobic packing interactions that normally stabilize atazanavir's positioning. This modest energetic penalty accumulated at a direct contact point 5.7 Å from the ligand centroid destabilizes the inhibitor binding without imposing steric clashes, allowing the protease to maintain catalytic competence while evading drug inhibition.

**YOUR score (0 / 1 / 2):** ____

---

## #5  RT · RPV · E138A
*(Claude judge said 1, GPT-5.6 said 2)*

**Known mechanism:** E138A disrupts the E138-K101 salt bridge between the p51 and p66 subunits of HIV-1 RT, reducing NNRTI-induced opening of the fingers and thumb subdomains and maintaining a catalytically competent RT-template/primer-dNTP complex despite NNRTI binding. The loss of negative charge at position 138 (glutamate to alanine) eliminates a key electrostatic interaction required for stabilizing the NNRTI-susceptible conformational state.

**AI explanation:** The E138A mutation reduces rilpivirine binding through an indirect allosteric mechanism rather than direct contact. The negatively charged glutamate at position 138 is located ~43 Å from the NNRTI pocket, and its replacement with small, neutral alanine eliminates a critical electrostatic anchor point that stabilizes the overall three-dimensional architecture of the palm domain; loss of this negative charge removes an organizing electrostatic influence that maintains the precise geometric presentation of rilpivirine-binding residues. The 50 Ų volume reduction coupled with the loss of negative charge likely permits localized conformational relaxation in the NNRTI-binding pocket, destabilizing favorable van der Waals and electrostatic complementarity with the inhibitor despite the direct contact points remaining geometrically intact.

**YOUR score (0 / 1 / 2):** ____

---

## #6  RT · EFV · F227L
*(Claude judge said 1, GPT-5.6 said 2)*

**Known mechanism:** F227L is a secondary NNRTI-resistance mutation that confers minimal resistance as a single mutation but enhances resistance when combined with primary mutations such as K103N or Y188L. The phenylalanine-to-leucine substitution reduces the aromatic bulk at position 227 in the primer-grip NNRTI-binding wall, allowing conformational flexibility in the NNRTI pocket that partially accommodates inhibitor binding despite primary resistance mutations present.

**AI explanation:** The F227L mutation reduces efavirenz binding through a loss of aromatic stabilization rather than direct contact—the bulky phenylalanine's π-electron system participates in favorable van der Waals interactions with efavirenz's aromatic rings from the primer grip wall, and its larger hydrophobic surface contributes to the optimal shape complementarity of the NNRTI-binding pocket at this distance. Substitution with leucine's smaller, saturated side chain eliminates these aromatic interactions and creates a subtle conformational void that destabilizes ligand positioning, weakening overall binding affinity despite remaining on the pocket perimeter. This represents a cumulative loss of hydrophobic stabilization in a region critical for maintaining the pocket geometry that accommodates efavirenz's tricyclic scaffold.

**YOUR score (0 / 1 / 2):** ____

---

## #7  RT · EFV · Y188L
*(Claude judge said 1, GPT-5.6 said 2)*

**Known mechanism:** Y188L replaces a tyrosine with leucine in the NNRTI binding pocket core, eliminating the aromatic hydroxyl group's capacity for hydrogen bonding and introducing steric constraints. This mutation increases both the dissociation rate (koff) and reduces the binding rate (kon) for NNRTIs like nevirapine and efavirenz, thereby conferring 17–90-fold resistance through loss of stabilizing interactions and increased thermodynamic/steric barriers to inhibitor binding.

**AI explanation:** The Y188L mutation replaces a bulky aromatic tyrosine with a smaller leucine within the NNRTI pocket core, creating a cavity that reduces favorable van der Waals contacts with efavirenz despite the residue not directly contacting the ligand. This loss of complementary packing at 6.7 Å minimum distance disrupts the precise shape complementarity of the binding pocket, destabilizing the inhibitor's orientation and positioning. Additionally, the substitution of tyrosine's hydroxyl group with leucine's nonpolar methyl side chain diminishes local hydrophobic interactions that normally contribute to binding affinity in this deep pocket region.

**YOUR score (0 / 1 / 2):** ____

---

## #9  PR · ATV · L90M
*(Claude judge said 0, GPT-5.6 said 1)*

**Known mechanism:** Leu90 does not contact the inhibitor directly. L90M repositions the side chain against the catalytic Asp25 loop (residues 23-25), perturbing active-site geometry and the catalytic aspartate environment, which indirectly weakens inhibitor binding.

**AI explanation:** The L90M substitution reduces atazanavir binding through an indirect mechanism involving disruption of the protease dimer interface stability rather than direct ligand contact. Leucine's larger, branched hydrophobic side chain provides optimal van der Waals packing within the C-terminal dimer region at this distance from the active site, and its replacement with methionine's smaller, linear side chain creates a subtle structural perturbation that destabilizes the dimeric architecture necessary for proper inhibitor accommodation. This allosteric effect compromises the precise conformational state required for optimal atazanavir positioning and binding geometry, despite the minimal absolute volume change.

**YOUR score (0 / 1 / 2):** ____

---

## #10  PR · ATV · L33F
*(Claude judge said 1, GPT-5.6 said 2)*

**Known mechanism:** Leu33 lies near the flap/substrate cleft. L33F introduces an aromatic side chain that alters hydrophobic packing against the inhibitor; an accessory mutation in the darunavir/fosamprenavir/atazanavir pathways.

**AI explanation:** The L33F substitution introduces a bulky aromatic ring at position 33 in the flap hinge region, which despite being 9.1 Å from the ligand, likely creates unfavorable steric or conformational effects that destabilize atazanavir binding. The 23.2 Å³ volume increase combined with the loss of hydrophobicity upon adding the polar aromatic ring may disrupt the precise positioning of the flap dynamics or indirectly compress the binding pocket geometry, reducing optimal ligand contacts. This indirect allosteric effect at the flap hinge—a region critical for protease conformational changes during inhibitor binding—weakens the overall stabilization of the atazanavir-bound state despite the mutation not directly contacting the drug.

**YOUR score (0 / 1 / 2):** ____

---

## #12  PR · SQV · L90M
*(Claude judge said 0, GPT-5.6 said 1)*

**Known mechanism:** Leu90 does not contact the inhibitor directly. L90M repositions the side chain against the catalytic Asp25 loop (residues 23-25), perturbing active-site geometry and the catalytic aspartate environment, which indirectly weakens inhibitor binding.

**AI explanation:** The L90M substitution reduces saquinavir binding primarily through loss of optimal hydrophobic packing in the C-terminal dimer region. Although the volume decrease is modest (~3.8 Å³), the reduction in hydrophobicity (-1.9 units) weakens van der Waals contacts that stabilize the inhibitor-binding pocket architecture, since leucine's branched aliphatic chain provides superior hydrophobic complementarity compared to methionine's linear thioether moiety. At 8.2 Å minimum distance from the ligand, position 90 indirectly maintains binding pocket geometry through cooperative hydrophobic interactions; the L→M change subtly destabilizes this network, destabilizing saquinavir's binding pose without directly clashing, explaining the modest but meaningful resistance (~0.38 kcal/mol).

**YOUR score (0 / 1 / 2):** ____

---

## #14  PR · IDV · L90M
*(Claude judge said 0, GPT-5.6 said 1)*

**Known mechanism:** Leu90 does not contact the inhibitor directly. L90M repositions the side chain against the catalytic Asp25 loop (residues 23-25), perturbing active-site geometry and the catalytic aspartate environment, which indirectly weakens inhibitor binding.

**AI explanation:** The L90M substitution reduces indinavir binding through a subtle loss of hydrophobic complementarity in the protease binding pocket. Although leucine and methionine are both hydrophobic residues with minimal volume change, methionine's lower hydrophobicity and slightly smaller aliphatic footprint create a small cavity or reduce optimal van der Waals packing contacts with indinavir's nonpolar surfaces at ~8 Å distance. This modest packing deficiency is amplified by methionine's longer, more flexible side chain, which may adopt suboptimal conformations that fail to maintain tight hydrophobic interactions that stabilize the inhibitor in the active site, cumulatively destabilizing the complex by approximately 0.34 kcal/mol.

**YOUR score (0 / 1 / 2):** ____

---

## #20  RT · NVP · E138G
*(Claude judge said 0, GPT-5.6 said 1)*

**Known mechanism:** E138G removes the negatively charged glutamic acid residue at position 138 on the p51 subunit, eliminating its role in the hydrogen-bonding network with NNRTIs and destabilizing critical water-mediated interactions that stabilize inhibitor binding in the NNRTI pocket. The glycine substitution, lacking the charge and bulky side chain, disrupts the charged interaction network at the p66/p51 subunit interface that normally positions Glu138 to maintain indirect contacts with NNRTIs through bridging water molecules, thereby reducing nevirapine binding affinity.

**AI explanation:** The E138G mutation causes minimal direct resistance (ΔΔG ≈ -0.01 kcal/mol) because position 138 is structurally remote from the NNRTI-binding pocket, located >43 Å away in the palm domain. The loss of the bulky, negatively charged glutamate side chain reduces steric constraints in the catalytic core, likely permitting subtle conformational relaxation that propagates to distant allosteric effects on the NNRTI-binding site rather than through direct contact. This long-range destabilization of the inhibitor-bound state probably arises from altered domain dynamics or electrostatic topology, where the charge neutralization and volume loss disrupt the precise structural scaffold that optimally positions the binding pocket for nevirapine accommodation.

**YOUR score (0 / 1 / 2):** ____

---
