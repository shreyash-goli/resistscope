"""Claude-generated mechanistic explanations for resistance mutations.

For a (drug, mutation) pair, we build a structural-context dict from the
wildtype receptor (residue identity, size/charge change, distance to the active
site, subpocket) and ask Claude for a short, grounded mechanistic hypothesis of
*why* the mutation reduces binding of that specific drug. Explanations are
cached to ``data/explanations/{drug}_{mutation}.json``.

Uses the Anthropic SDK (``anthropic``); requires ``ANTHROPIC_API_KEY`` (or an
``ant auth login`` profile). Model/params come from ``config``.
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

import config

# --- Amino-acid physicochemical properties ---
# volume (A^3, Zamyatnin 1972), charge at pH 7, Kyte-Doolittle hydrophobicity.
AA_PROPERTIES = {
    "A": {"name": "Ala", "volume": 88.6,  "charge": "neutral",  "hydrophobicity": 1.8},
    "R": {"name": "Arg", "volume": 173.4, "charge": "positive", "hydrophobicity": -4.5},
    "N": {"name": "Asn", "volume": 114.1, "charge": "neutral",  "hydrophobicity": -3.5},
    "D": {"name": "Asp", "volume": 111.1, "charge": "negative", "hydrophobicity": -3.5},
    "C": {"name": "Cys", "volume": 108.5, "charge": "neutral",  "hydrophobicity": 2.5},
    "E": {"name": "Glu", "volume": 138.4, "charge": "negative", "hydrophobicity": -3.5},
    "Q": {"name": "Gln", "volume": 143.8, "charge": "neutral",  "hydrophobicity": -3.5},
    "G": {"name": "Gly", "volume": 60.1,  "charge": "neutral",  "hydrophobicity": -0.4},
    "H": {"name": "His", "volume": 153.2, "charge": "positive", "hydrophobicity": -3.2},
    "I": {"name": "Ile", "volume": 166.7, "charge": "neutral",  "hydrophobicity": 4.5},
    "L": {"name": "Leu", "volume": 166.7, "charge": "neutral",  "hydrophobicity": 3.8},
    "K": {"name": "Lys", "volume": 168.6, "charge": "positive", "hydrophobicity": -3.9},
    "M": {"name": "Met", "volume": 162.9, "charge": "neutral",  "hydrophobicity": 1.9},
    "F": {"name": "Phe", "volume": 189.9, "charge": "neutral",  "hydrophobicity": 2.8},
    "P": {"name": "Pro", "volume": 112.7, "charge": "neutral",  "hydrophobicity": -1.6},
    "S": {"name": "Ser", "volume": 89.0,  "charge": "neutral",  "hydrophobicity": -0.8},
    "T": {"name": "Thr", "volume": 116.1, "charge": "neutral",  "hydrophobicity": -0.7},
    "W": {"name": "Trp", "volume": 227.8, "charge": "neutral",  "hydrophobicity": -0.9},
    "Y": {"name": "Tyr", "volume": 193.6, "charge": "neutral",  "hydrophobicity": -1.3},
    "V": {"name": "Val", "volume": 140.0, "charge": "neutral",  "hydrophobicity": 4.2},
}

# --- HIV-1 protease subpocket / structural region by residue position ---
# Coarse, single-region-per-position assignment for interpretability. Ranges are
# grounded in HIV-1 protease structural biology (catalytic dyad ~25, flaps
# 43-58, S1 subpocket lining ~78-85 incl. V82/I84). ``(start, end)`` is inclusive
# on start, exclusive on end.
SUBPOCKET_REGIONS = [
    (range(1, 10), "N_terminal_region"),
    (range(10, 23), "beta_sheet_scaffold"),
    (range(23, 28), "catalytic_dyad_region"),        # 23,24,25(Asp),26,27
    (range(28, 33), "S2_S2prime_subpocket"),          # D29,D30,V32
    (range(33, 43), "flap_hinge_region"),
    (range(43, 59), "flap_region"),                   # flaps incl. M46,I47,G48,I50,I54
    (range(59, 78), "core_scaffold"),                 # incl. L76
    (range(78, 86), "S1_S1prime_subpocket"),          # direct P1/P1' contact: V82,I84
    (range(86, 100), "C_terminal_dimer_region"),      # incl. N88,L90
]

# A residue atom within this distance of any ligand (ROC) atom is a direct
# van der Waals contact.
LIGAND_CONTACT_THRESHOLD_A = 4.5

# Lazily-loaded coordinates of the co-crystal ligand (saquinavir / ROC) atoms
# from the raw 3OXC structure — used for accurate residue-ligand contact
# distances (the cleaned wildtype has the ligand stripped).
_LIGAND_COORDS = None


def _ligand_atom_coords():
    """Return an (N, 3) array of ROC ligand atom coords from raw 3OXC, or None."""
    global _LIGAND_COORDS
    if _LIGAND_COORDS is None:
        raw = config.RAW_DIR / f"{config.WILDTYPE_PDB_ID}.pdb"
        if not raw.exists():
            _LIGAND_COORDS = False  # sentinel: unavailable
        else:
            codes = {c.strip().upper() for c in config.LIGAND_HETCODES_TO_STRIP}
            struct = PDBParser(QUIET=True).get_structure("raw", str(raw))
            coords = [
                atom.coord
                for atom in struct[0].get_atoms()
                if atom.get_parent().get_resname().strip().upper() in codes
            ]
            _LIGAND_COORDS = np.array(coords, dtype=float) if coords else False
    return None if _LIGAND_COORDS is False else _LIGAND_COORDS


def _region_for(position: int) -> str:
    for rng, label in SUBPOCKET_REGIONS:
        if position in rng:
            return label
    return "other_region"


def _charge_change(wt: str, mut: str) -> str:
    wc = AA_PROPERTIES[wt]["charge"]
    mc = AA_PROPERTIES[mut]["charge"]
    if wc == mc:
        return "none"
    return f"{wc}_to_{mc}"


def _volume_change(wt: str, mut: str) -> tuple[str, float]:
    delta = AA_PROPERTIES[mut]["volume"] - AA_PROPERTIES[wt]["volume"]
    if delta > 10:
        label = "increase"
    elif delta < -10:
        label = "decrease"
    else:
        label = "minimal"
    return label, round(delta, 1)


def _parse_mutation(mutation: str) -> tuple[str, int, str]:
    """Split e.g. 'V82A' -> ('V', 82, 'A')."""
    wt = mutation[0].upper()
    mut = mutation[-1].upper()
    position = int(mutation[1:-1])
    return wt, position, mut


def build_structural_context(
    mutation: str,
    docking_result: dict = None,
    wildtype_pdb: Path = None,
) -> dict:
    """Assemble the structural context Claude will reason over.

    ``docking_result`` may carry ``delta_g``/``delta_delta_g`` (optional; passed
    through for convenience). Distances are measured from the mutated residue's
    CA (chain A) to ``config.DOCKING_CENTER`` (the ligand centroid).

    Returns a dict matching the explanation-cache ``structural_context`` schema.
    """
    if wildtype_pdb is None:
        wildtype_pdb = config.STRUCTURES_DIR / "wildtype.pdb"
    wt_aa, position, mut_aa = _parse_mutation(mutation)

    if wt_aa not in AA_PROPERTIES or mut_aa not in AA_PROPERTIES:
        raise ValueError(f"Non-standard residue in mutation {mutation!r}")

    # Distance from the residue CA to the docking-box center (ligand centroid),
    # plus the accurate minimum distance from any residue atom to any ligand atom.
    structure = PDBParser(QUIET=True).get_structure("wt", str(wildtype_pdb))
    center = np.array(config.DOCKING_CENTER, dtype=float)
    ligand = _ligand_atom_coords()
    distance = None          # CA -> ligand centroid
    min_lig_dist = None      # nearest residue atom -> nearest ligand atom
    chain_a = next((c for c in structure[0] if c.id == "A"), None)
    if chain_a is not None:
        residues = list(chain_a)
        if 1 <= position <= len(residues):
            res = residues[position - 1]  # chain A is protease residues 1-99
            ca = res["CA"] if "CA" in res else next(iter(res), None)
            if ca is not None:
                distance = float(np.linalg.norm(np.asarray(ca.coord) - center))
            if ligand is not None:
                res_coords = np.array([a.coord for a in res], dtype=float)
                if len(res_coords):
                    # pairwise min distance between residue atoms and ligand atoms
                    d = np.linalg.norm(
                        res_coords[:, None, :] - ligand[None, :, :], axis=-1
                    )
                    min_lig_dist = float(d.min())

    vol_label, vol_delta = _volume_change(wt_aa, mut_aa)
    hydro_delta = round(
        AA_PROPERTIES[mut_aa]["hydrophobicity"] - AA_PROPERTIES[wt_aa]["hydrophobicity"], 1
    )

    context = {
        "position": position,
        "wt_aa": wt_aa,
        "mut_aa": mut_aa,
        "wt_residue": AA_PROPERTIES[wt_aa]["name"],
        "mut_residue": AA_PROPERTIES[mut_aa]["name"],
        "volume_change": vol_label,
        "volume_change_A3": vol_delta,
        "charge_change": _charge_change(wt_aa, mut_aa),
        "hydrophobicity_change": hydro_delta,
        "distance_from_ligand_centroid_angstrom": (
            round(distance, 1) if distance is not None else None
        ),
        "min_distance_to_ligand_angstrom": (
            round(min_lig_dist, 1) if min_lig_dist is not None else None
        ),
        "contacts_ligand_directly": (
            bool(min_lig_dist is not None and min_lig_dist < LIGAND_CONTACT_THRESHOLD_A)
        ),
        "region": _region_for(position),
    }
    if docking_result:
        if docking_result.get("delta_delta_g") is not None:
            context["delta_delta_g"] = float(docking_result["delta_delta_g"])
        if docking_result.get("delta_g") is not None:
            context["delta_g"] = float(docking_result["delta_g"])
    return context


SYSTEM_PROMPT = (
    "You are a structural biologist specializing in HIV-1 protease drug "
    "resistance. Given structural context about a point mutation and its effect "
    "on inhibitor binding, provide a 2-3 sentence mechanistic hypothesis for WHY "
    "this mutation causes resistance to this specific drug. Ground your "
    "explanation in the structural facts provided (residue size/charge change, "
    "distance to the active site, subpocket). Be specific about molecular "
    "interactions (van der Waals contacts, hydrogen bonds, steric clashes, "
    "electrostatic changes). Do not speculate beyond the given data, and do not "
    "restate the numbers back — interpret them. No preamble."
)


def generate_explanation(
    drug: str,
    mutation: str,
    delta_delta_g: float,
    structural_context: dict,
    cache_dir: Path = None,
    model: str = None,
) -> str:
    """Return a cached or freshly-generated mechanistic explanation.

    Cache key is ``{drug}_{mutation}.json`` in ``cache_dir``. On a cache hit the
    stored explanation is returned without an API call.
    """
    if cache_dir is None:
        cache_dir = config.EXPLANATIONS_DIR
    cache_dir = Path(cache_dir)
    cache_path = cache_dir / f"{drug}_{mutation}.json"

    if cache_path.exists():
        return json.loads(cache_path.read_text())["explanation"]

    import anthropic

    model = model or config.CLAUDE_MODEL
    drug_full = config.PI_DRUGS.get(drug, drug)
    ddg_str = (
        f"{delta_delta_g:+.2f} kcal/mol (positive = weaker binding = resistance)"
        if delta_delta_g is not None
        else "not available"
    )
    user_prompt = (
        f"Drug: {drug_full} ({drug}), an HIV-1 protease inhibitor\n"
        f"Mutation: {mutation}\n"
        f"Predicted delta-delta-G: {ddg_str}\n"
        f"Structural context:\n{json.dumps(structural_context, indent=2)}\n\n"
        f"Explain the structural mechanism by which {mutation} likely reduces "
        f"{drug_full} binding."
    )

    client = anthropic.Anthropic()  # ANTHROPIC_API_KEY or ant-auth profile
    response = client.messages.create(
        model=model,
        max_tokens=config.CLAUDE_MAX_TOKENS,
        system=SYSTEM_PROMPT,
        messages=[{"role": "user", "content": user_prompt}],
    )
    explanation = next(
        (b.text for b in response.content if b.type == "text"), ""
    ).strip()

    record = {
        "drug": drug,
        "mutation": mutation,
        "delta_delta_g": (float(delta_delta_g) if delta_delta_g is not None else None),
        "structural_context": structural_context,
        "explanation": explanation,
        "model": model,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(record, indent=2))
    return explanation


# --- Faithfulness evaluation (Claude-as-judge) -------------------------------

JUDGE_SYSTEM_PROMPT = (
    "You are evaluating whether an AI-generated explanation of HIV-1 protease "
    "drug resistance correctly identifies the known structural mechanism. Score "
    "0 if the explanation contradicts the known mechanism, 1 if it is consistent "
    "but vague or misses the key interaction, 2 if it correctly identifies the "
    "primary structural basis. Judge only structural/mechanistic agreement, not "
    "writing style, and give a one-sentence justification."
)

# Structured-output schema so the 0/1/2 score is always machine-readable.
FAITHFULNESS_SCHEMA = {
    "type": "object",
    "properties": {
        "score": {"type": "integer", "enum": [0, 1, 2]},
        "justification": {"type": "string"},
    },
    "required": ["score", "justification"],
    "additionalProperties": False,
}


def _faithfulness_pairs(explanations_dir: Path, ground_truth: dict) -> list[dict]:
    """Enumerate cached explanations that have a ground-truth mechanism.

    For each mutation in ``ground_truth``, finds cached ``{drug}_{mutation}.json``
    explanations for the drugs it affects. Returns dicts with keys
    ``mutation, drug, explanation, mechanism`` — no API calls.
    """
    explanations_dir = Path(explanations_dir)
    pairs = []
    for mutation, gt in ground_truth.items():
        drugs = gt.get("affects_drugs") or list(config.PI_DRUGS)
        for drug in drugs:
            cache = explanations_dir / f"{drug}_{mutation}.json"
            if not cache.exists():
                continue
            expl = json.loads(cache.read_text()).get("explanation", "")
            if expl:
                pairs.append({
                    "mutation": mutation, "drug": drug,
                    "explanation": expl, "mechanism": gt["mechanism"],
                })
    return pairs


def evaluate_faithfulness(
    explanations_dir: Path = None,
    ground_truth_path: Path = None,
    model: str = None,
) -> pd.DataFrame:
    """Score cached explanations against curated ground-truth mechanisms.

    Uses Claude as a judge (0=contradicts, 1=consistent-but-vague, 2=correct)
    with a structured-output schema so scores parse reliably. Returns a
    DataFrame with columns ``mutation, drug, score, justification`` (empty if no
    cached explanation overlaps the ground truth).
    """
    if explanations_dir is None:
        explanations_dir = config.EXPLANATIONS_DIR
    if ground_truth_path is None:
        ground_truth_path = config.DATA_DIR / "mechanism_ground_truth.json"

    ground_truth = json.loads(Path(ground_truth_path).read_text())
    pairs = _faithfulness_pairs(explanations_dir, ground_truth)
    if not pairs:
        return pd.DataFrame(
            columns=["mutation", "drug", "score", "justification"]
        )

    import anthropic

    model = model or config.CLAUDE_MODEL
    client = anthropic.Anthropic()

    rows = []
    for p in pairs:
        user_prompt = (
            f"Mutation: {p['mutation']}\n"
            f"Drug: {config.PI_DRUGS.get(p['drug'], p['drug'])} ({p['drug']})\n\n"
            f"Known structural mechanism (ground truth):\n{p['mechanism']}\n\n"
            f"AI-generated explanation to evaluate:\n{p['explanation']}\n\n"
            f"Score how faithfully the explanation captures the known mechanism."
        )
        response = client.messages.create(
            model=model,
            max_tokens=512,
            system=JUDGE_SYSTEM_PROMPT,
            messages=[{"role": "user", "content": user_prompt}],
            output_config={"format": {"type": "json_schema", "schema": FAITHFULNESS_SCHEMA}},
        )
        text = next((b.text for b in response.content if b.type == "text"), "{}")
        verdict = json.loads(text)
        rows.append({
            "mutation": p["mutation"],
            "drug": p["drug"],
            "score": int(verdict["score"]),
            "justification": verdict.get("justification", ""),
        })

    return pd.DataFrame(rows)
