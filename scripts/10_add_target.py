"""10: Bring your own target — assemble a new resistance target from a PDB + name.

Ties the two halves together:
  • structure  → services/pdb_intake.py parses the PDB (chains, ligands, pocket)
  • biology    → services/target_builder.py agentically assembles the inhibitor
                 panel + resistance mutations + dataset status from PubMed

and writes a reviewable draft to ``data/user_targets/<id>.json``. This is the
"upload a PDB, we find the drugs/mutations, then build mutant receptors" flow —
everything up to the (GPU) receptor build, which runs on the docking worker.

Usage::

    python scripts/10_add_target.py --pdb 2HU4 --protein "influenza A neuraminidase"
    python scripts/10_add_target.py --pdb ./my_receptor.pdb --protein "SARS-CoV-2 Mpro" \\
        --hint "nirmatrelvir resistance"
"""

import argparse
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402
from services.pdb_intake import parse_pdb  # noqa: E402
from services.target_builder import assemble_target  # noqa: E402


def _load_pdb(pdb: str) -> tuple[str, str]:
    """Return (pdb_text, source_label). ``pdb`` is a local path or a 4-char id."""
    p = Path(pdb)
    if p.exists():
        return p.read_text(), p.name
    pid = pdb.upper()
    url = f"https://files.rcsb.org/download/{pid}.pdb"
    print(f"Fetching {url} …")
    with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310
        return r.read().decode(), pid


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True, help="PDB id (e.g. 2HU4) or path to a .pdb file.")
    ap.add_argument("--protein", required=True, help="Protein/target name for the agent.")
    ap.add_argument("--hint", default="", help="Optional context (drug class, resistance focus).")
    ap.add_argument("--no-agent", action="store_true", help="Structural intake only (skip Claude).")
    args = ap.parse_args()

    text, source = _load_pdb(args.pdb)
    # Persist the raw PDB so the API's structure viewer + ligand extraction work.
    raw_id = Path(source).stem.upper()
    config.RAW_DIR.mkdir(parents=True, exist_ok=True)
    (config.RAW_DIR / f"{raw_id}.pdb").write_text(text)
    intake = parse_pdb(text)
    print(f"\nStructure ({source}): {intake['title'][:70]}")
    print(f"  protein chains: {[c['id'] for c in intake['chains'] if c['kind']=='protein']}")
    print(f"  suggested pocket ligand: {intake['suggested']['ligand_hetcode']} "
          f"@ {intake['suggested']['docking_center']}")
    print(f"  → CONFIRM this is the inhibitor pocket (not a glycan/ion) before docking.")

    spec = {}
    citations = []
    if not args.no_agent:
        print(f"\nAgent researching '{args.protein}' …")
        res = assemble_target(args.protein, hint=(args.hint or intake["title"]))
        if res.error:
            print(f"  agent error: {res.error}")
        else:
            spec, citations = res.spec, res.citations
            print(f"  drugs: {[d['abbrev'] for d in spec.get('drugs', [])]}")
            print(f"  resistance mutations: {[m['mutation'] for m in spec.get('resistance_mutations', [])]}")
            print(f"  dataset: {spec.get('dataset', {}).get('exists')} "
                  f"({spec.get('dataset', {}).get('name', '—')})  "
                  f"→ {'validatable' if spec.get('dataset',{}).get('exists') else 'triage-only'}")

    target_id = (spec.get("target_id") or Path(source).stem).upper()
    draft = {
        "target_id": target_id,
        "label": spec.get("label", args.protein),
        "source_pdb": source,
        "structure": intake,
        "biology": spec,
        "citations": citations,
        "provenance": "byo-agent",
        "status": "draft — confirm pocket + mutations, then build receptors on the GPU worker",
    }
    out_dir = config.DATA_DIR / "user_targets"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{target_id}.json"
    out.write_text(json.dumps(draft, indent=2))

    print(f"\nSaved target draft → {out}")
    print("Next: confirm the pocket + mutation list, then on the GPU worker run "
          "structure prep + scripts/03 to build the mutant receptors, and triage "
          "any SMILES against it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
