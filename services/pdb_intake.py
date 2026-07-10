"""Parse an uploaded/downloaded PDB into the structural half of a new target.

The "bring your own target" flow has two halves: the **biology** (drugs +
resistance mutations), which the Claude agent assembles
(``services/target_builder.py``), and the **structure** (which chains form the
receptor, which ligand marks the pocket, where the docking box goes), which is
read straight from the deposited coordinates here — exactly the workflow used to
confirm the HIV-RT box (occupancy-weighted centroid of the co-crystal ligand).

No docking stack needed (BioPython only), so this runs on the thin API host.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

# Het codes that are solvent / ions / cryo-additives, never the pocket ligand.
_NON_LIGAND = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "ACT", "FMT", "CL",
    "NA", "K", "MG", "CA", "ZN", "MN", "CD", "NAG", "MAN", "BMA", "FUC", "DMS",
    "IOD", "BR", "NO3", "TRS", "MES", "EPE", "CO3", "PG4", "1PE",
}


def _centroid(atoms) -> tuple:
    coords = np.array([a.coord for a in atoms], dtype=float)
    occ = np.array([a.get_occupancy() or 1.0 for a in atoms], dtype=float)
    c = (coords * occ[:, None]).sum(0) / occ.sum()
    return tuple(round(float(x), 3) for x in c)


def parse_pdb(source: str | Path) -> dict:
    """Parse a PDB (path or raw text) into a structural draft for a target.

    Returns a dict with:
      ``chains``          – [{id, n_residues, kind}]   (kind: protein/nucleic/other)
      ``ligands``         – [{hetcode, chain, n_atoms, centroid}] candidate pocket ligands
      ``suggested``       – {chains, mutate_chains, ligand_hetcode, docking_center}
      ``title``           – PDB header title (feeds the agent as a hint)
    Best-effort; a caller should let the user confirm the suggestion.
    """
    # Distinguish a filesystem path from raw PDB text: PDB content has newlines
    # (and is long), a path does not — so never feed a 2 MB string to Path().
    if isinstance(source, Path):
        text = source.read_text()
    elif isinstance(source, str) and "\n" not in source and len(source) < 4096 and Path(source).exists():
        text = Path(source).read_text()
    else:
        text = str(source)
    # Title from the header (helps the agent identify the protein).
    title = " ".join(l[10:80].strip() for l in text.splitlines()
                     if l.startswith("TITLE")).strip()

    import io
    struct = PDBParser(QUIET=True).get_structure("up", io.StringIO(text))
    model = next(iter(struct))

    chains, ligands = [], {}
    aa = set("ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
    nuc = set("DA DT DG DC DU A U G C".split())
    for ch in model:
        prot = sum(1 for r in ch if r.get_resname().strip() in aa)
        nucl = sum(1 for r in ch if r.get_resname().strip() in nuc)
        kind = "protein" if prot >= max(nucl, 1) else ("nucleic" if nucl else "other")
        chains.append({"id": ch.id, "n_residues": prot or nucl or len(list(ch)), "kind": kind})
        for r in ch:
            code = r.get_resname().strip().upper()
            if r.id[0] == " " or code in _NON_LIGAND:
                continue
            atoms = [a for a in r]
            key = (code, ch.id)
            ligands.setdefault(key, []).extend(atoms)

    lig_list = [{"hetcode": code, "chain": cid, "n_atoms": len(atoms),
                 "centroid": _centroid(atoms)}
                for (code, cid), atoms in ligands.items()]
    lig_list.sort(key=lambda x: x["n_atoms"], reverse=True)  # biggest = likely inhibitor

    protein_chains = [c["id"] for c in chains if c["kind"] == "protein"]
    best = lig_list[0] if lig_list else None
    suggested = {
        "chains": protein_chains[:2] or [c["id"] for c in chains][:1],
        "mutate_chains": ([best["chain"]] if best and best["chain"] in protein_chains
                          else protein_chains[:1]),
        "ligand_hetcode": best["hetcode"] if best else None,
        "docking_center": best["centroid"] if best else None,
    }
    return {"title": title, "chains": chains, "ligands": lig_list, "suggested": suggested}
