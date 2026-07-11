# Building a bring-your-own target's receptors on the GPU worker

Takes a user-added target (e.g. `FLU_NA`, influenza neuraminidase) from
**"docking pending"** to **live ΔΔG triage**: build the wildtype + mutant
receptors, then serve them through the docking worker. Everything up to this
point (add the target in the app, confirm the pocket, agent-propose the
mutations) needs no GPU; **this** is the GPU-gated step.

The pipeline is already wired for it:
- `scripts/03` derives the build list from the target's mutation set when there
  are no genotype panels (`collect_unique_mutations` fallback).
- `generate_mutant` maps mutations by **author residue number (resSeq)**, so a
  structure numbered 83..468 (like NA) mutates residue 275 correctly — not the
  275th residue in the chain.

---

## 0. One-time setup on the GPU box

The **lean** `deploy/setup_a100.sh` env (RDKit + meeko + Uni-Dock) can *dock* but
not *build* receptors — building needs OpenMM + PDBFixer. Create the full env so
one box does both:

```bash
git clone <this repo> && cd resistscope
conda env create -f environment.yml        # RDKit, OpenMM, PDBFixer, Vina, meeko
conda activate resistscope
python -m pip install gemmi                 # meeko 0.7.x needs it
# GPU docking engine (skip to use CPU Vina):
bash deploy/setup_a100.sh                   # or follow docs/UNIDOCK_SETUP.md for Uni-Dock
export ANTHROPIC_API_KEY=sk-ant-...         # optional, for live explanations (or use .env)
```

> Alternative split: build receptors on any machine with the full env, then
> `rsync data/user/<ID>/structures/` to the GPU box and run only the worker there.

---

## 1. Get the target + its structure onto the box

The target draft (`data/user_targets/<ID>.json`) is committed, so it's already in
the repo. Its **raw structure is gitignored** (`data/raw/`), so fetch/copy it:

```bash
# added by PDB id (e.g. 2HU4) — just download it:
curl -sL https://files.rcsb.org/download/2HU4.pdb -o data/raw/2HU4.pdb
# added by upload — copy the PDB the app saved (data/raw/<ID>.pdb, mmCIF already
# converted to PDB on save) from the API host to the same path here.
```

Confirm the target resolves, with the pocket/box the app saved:

```bash
python -c "import config; t=config.get_target('FLU_NA'); \
  print(t.label, '| het', t.ligand_hetcodes, '| box', t.docking_center, \
        '| mutate', t.mutate_chains, '| muts', sorted(t.primary_mutations))"
```

---

## 2. Build the wildtype + mutant receptors

```bash
# smoke-test two mutants first (~1 min):
python scripts/03_build_mutant_cache.py --target FLU_NA --limit 2
# then the full set:
python scripts/03_build_mutant_cache.py --target FLU_NA
```

This writes `data/user/FLU_NA/structures/wildtype.{pdb,pdbqt}` and
`.../mutants/<MUT>.pdbqt` for every mutation. The wildtype step prints the
occupancy-weighted ligand centroid — **confirm it matches** `t.docking_center`
(the pocket the app picked).

**Check the residue numbering matched.** The agent may propose a mutation in more
than one convention (NA has N1 *and* N2 numbering, e.g. `H274Y` vs `H275Y`) —
keep the one whose residue is actually present in your structure, and verify a
built mutant changed the intended residue:

```bash
grep -E "^ATOM.{9}CA .(TYR) A 275" data/user/FLU_NA/structures/mutants/H275Y.pdb
# (a hit means residue 275 is now TYR — the mutation landed on the right residue)
```

If a mutant failed to build, it's logged (`FAILED …`) and simply omitted from the
panel — the run is resumable (existing PDBQTs are skipped).

---

## 3. Serve the worker + point the API at it

```bash
# on the GPU box — serves POST /dock on :9000 using the receptors just built
python docking_worker.py
#   GET /health shows targets_ready: {"FLU_NA": <n mutants>, ...}

# on the API host (any machine)
RESISTSCOPE_DOCKING_URL=http://<gpu-host>:9000 python -m uvicorn api.main:app --port 8000
```

In the app: pick **FLU_NA**, paste a SMILES (e.g. oseltamivir
`CCC(CC)O[C@@H]1C=C(C[C@H]([C@H]1NC(C)=O)N)C(=O)O`), **Triage** → live ΔΔG against
H275Y / E119V / R292K / N294S, with Claude explanations. The "docking pending"
banner disappears automatically once `data/user/FLU_NA/structures/mutants/` has
receptors — no code change, no redeploy.

---

## Notes & limits

- **Triage-only.** A BYO target has no genotype–phenotype dataset, so there's no
  Validation tab for it (the method is validated on protease/RT). If the agent
  flags that a public dataset exists, wiring it (fetch → panels → `scripts/05/08`)
  is the path to validating the new target — future work.
- **Numbering.** Mutations map by author resSeq, falling back to ordinal for
  chains PDBFixer renumbers. If your structure uses a numbering that doesn't match
  the DRM convention, rename the mutations to the structure's numbering before
  step 2.
- **mmCIF.** Uploaded `.cif` is converted to PDB on save; a very large assembly
  converts lossily (PDB caps at 99,999 atoms / single-char chain ids) — fine for
  a single receptor + ligand.
- **Pocket confirmation.** Auto-detection picks the largest non-solvent het group;
  always confirm in the Add-target dialog that it's the inhibitor site, not a
  glycan/ion (e.g. NA structures carry NAG glycans).
```
