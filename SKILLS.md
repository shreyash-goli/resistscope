---
name: hiv-protease-docking
description: >
  Domain knowledge for building the ResistScope HIV-1 protease resistance
  triage pipeline. Use this skill whenever writing, editing, or debugging code
  that touches HIV-1 protease structure preparation, point mutagenesis,
  SMILES-to-PDBQT conversion, AutoDock Vina docking, or delta-delta-G scoring
  in this project. This encodes environment-specific constraints (correct
  library APIs, protonation requirements, the specific PDB structure and its
  quirks, docking box coordinates) that are NOT obvious from general knowledge
  and that, if gotten wrong, produce silently incorrect docking scores.
  Read this before writing any structure_prep, docking, or scoring code.
---

# HIV-1 Protease Docking — Project Skill

This project docks candidate compounds against a panel of HIV-1 protease
resistance mutants and predicts drug resistance. The scientific validity of
every result depends on getting a handful of non-obvious details right. Getting
any of them wrong produces docking scores that *look* plausible but are
scientifically meaningless. This skill enforces those details.

## The structure: PDB 3OXC

- **Use PDB 3OXC** as the wildtype receptor. It is HIV-1 protease at 1.16 Å
  resolution, co-crystallized with saquinavir.
- **Saquinavir's het code in 3OXC is `ROC`**, NOT `SQV` and NOT `938`. Code
  that strips ligands by het code MUST use `ROC` or it will silently strip
  nothing and leave the ligand in the pocket, breaking docking.
- **3OXC is ligand-bound, not apo.** This is intentional and correct: the
  protein is in a closed-flap, drug-binding-competent conformation. Do NOT
  substitute an apo structure like 2PC0 (open flaps, wrong conformation for
  docking inhibitors). Strip the ligand but keep the protein coordinates.
- **HIV-1 protease is a C2-symmetric homodimer.** Keep BOTH chains A and B.
  The active site sits at the dimer interface. A single chain is not a valid
  receptor.
- The structure also contains **sulfate ions** (star-shaped in viewers). These
  are crystallization artifacts, not part of the ligand. Strip them along with
  waters and other heterogens.

## The docking box

- **Center: `(5.341, -1.893, 14.179)`** in the deposited 3OXC Cartesian frame.
  This is the occupancy-weighted centroid of the saquinavir (ROC) ligand,
  which occupies the active site. It was computed from the two modeled altLoc
  conformations (altLoc A occ 0.51 at `(5.100, -3.435, 15.307)`; altLoc B
  occ 0.49 at `(5.592, -0.287, 13.006)`).
- **Box size: `(22, 22, 22)` Å.** Large enough to cover the full active site
  cavity including the S1/S1'/S2/S2' subpockets.
- When verifying the centroid programmatically, saquinavir's two altLoc
  conformations produce duplicate atom records. Use BioPython's
  `atom.is_disordered()` and either take a single altLoc or occupancy-weight
  across both. The occupancy-weighted target is `(5.341, -1.893, 14.179)`.

## Asp25/Asp25' protonation — the #1 silent failure

HIV-1 protease has a catalytic aspartate dyad: Asp25 in chain A and Asp25 in
chain B. **In the active enzyme these are asymmetrically protonated** — one is
protonated (neutral, ASH) and one is deprotonated (charged, ASP). This
asymmetry is essential to the electrostatics of the binding pocket.

PDBFixer's `addMissingHydrogens(7.0)` will protonate BOTH Asp25 residues (or
deprotonate both, depending on version) symmetrically. This is wrong and will
skew every docking score in the project.

**Required fix, applied after `addMissingHydrogens()` and to BOTH wildtype and
every mutant:**
- Keep chain A Asp25 protonated (neutral).
- Deprotonate chain B Asp25 by finding and deleting the hydrogen atom on its
  OD2 oxygen (atom name `HD2`).

Every structure that goes into docking — wildtype AND all mutants — must pass
through this fix. Do not skip it for mutants.

## Point mutagenesis

- **Use PDBFixer, not PyMOL.** PDBFixer's `applyMutations()` works cleanly in
  headless/scripted mode; PyMOL's mutagenesis wizard is finicky when scripted.
- **Apply each mutation to BOTH chains A and B.** Because the enzyme is a
  homodimer, resistance mutations like V82A, I84V, L90M occur in both
  monomers. Mutating only one chain models a heterodimer that does not reflect
  the biology.
- Mutation string format for PDBFixer is `RESNAME-POSITION-NEWRESNAME`, e.g.
  `VAL-82-ALA`, applied per chain:
  ```python
  fixer.applyMutations(["VAL-82-ALA"], "A")
  fixer.applyMutations(["VAL-82-ALA"], "B")
  ```
- After mutating: `findMissingResidues()`, `findMissingAtoms()`,
  `addMissingAtoms()`, `addMissingHydrogens(7.0)`, then the Asp25 protonation
  fix, then save.

## SMILES → PDBQT (ligand prep)

- **RDKit does NOT output PDBQT.** You need `meeko` for the conversion.
- Correct pattern (in memory, no temp files):
  ```python
  from rdkit import Chem
  from rdkit.Chem import AllChem
  from meeko import MoleculePreparation, PDBQTWriterLegacy

  mol = Chem.MolFromSmiles(smiles)
  mol = Chem.AddHs(mol)                          # explicit H required by meeko
  AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())  # 3D coords
  AllChem.MMFFOptimizeMolecule(mol)              # geometry cleanup

  mk_prep = MoleculePreparation()
  molsetup_list = mk_prep(mol)                   # returns a LIST
  pdbqt_string = PDBQTWriterLegacy.write_string(molsetup_list[0])
  ```
- `mk_prep(mol)` returns a list; take element `[0]` unless doing reactive
  docking.
- Meeko does not assign protonation states or generate 3D coordinates — RDKit
  must do both first. `AddHs()` before `EmbedMolecule()`, in that order.

## Receptor prep (PDB → PDBQT)

- Use meeko's `mk_prepare_receptor.py` (CLI via subprocess) or the `Polymer`
  Python API. The receptor PDBQT must be regenerated for every mutant, not
  just wildtype.

## AutoDock Vina API

- Correct pattern:
  ```python
  from vina import Vina
  v = Vina(sf_name='vina')
  v.set_receptor(str(receptor_pdbqt_path))
  v.set_ligand_from_string(ligand_pdbqt_string)   # from_string, not from_file
  v.compute_vina_maps(center=[5.341, -1.893, 14.179], box_size=[22, 22, 22])
  v.dock(exhaustiveness=16, n_poses=5)
  energies = v.energies()          # numpy array, shape (n_poses, 4+)
  best_affinity = float(energies[0][0])   # kcal/mol, MORE NEGATIVE = better
  ```
- Binding affinity is in kcal/mol and **more negative means stronger binding.**
- **delta_delta_g = delta_g_mutant − delta_g_wildtype.** A POSITIVE
  delta_delta_g means the mutant binds WORSE than wildtype, i.e. resistance.
  If resistant mutants are coming out with negative delta_delta_g, something
  is wrong.

## Sanity checks — run these before trusting any batch of scores

- **Darunavir against wildtype should score roughly −12 to −14 kcal/mol.**
  If you get −5 (way too weak) or −20 (unphysically strong), the box
  placement, protonation, or PDBQT conversion is broken. Stop and debug before
  scaling up. Do NOT run the full 7-drug batch on top of a broken single-dock.
- Known resistance mutations (V82A for several PIs, I84V, L90M) should produce
  positive delta_delta_g for the drugs they're known to affect. If they don't,
  investigate before proceeding.

## Failure handling

- Wrap every individual dock in try/except. A single ligand-receptor pair can
  fail (bad geometry, no poses, timeout) without meaning the batch is broken.
  On failure, record `delta_g = None` and continue; never let one failed dock
  crash a multi-hour batch run.

## Known limitations to document, not fix (for the hackathon)

- Crystal waters (including the conserved flap water) are stripped. This adds
  noise but avoids the complexity of placing waters in mutant structures.
- Ligands use a single RDKit-assigned protonation state, not pH-enumerated
  protomers.
- Rigid receptor docking (no flexible side chains). These are all reasonable
  post-hackathon improvements and honest limitations for a one-week build.