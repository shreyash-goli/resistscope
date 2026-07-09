"""FastAPI application exposing the ResistScope triage pipeline.

Endpoints:
  POST /triage     dock a SMILES against the resistance panel, score it, and
                   explain the worst mutations.
  GET  /benchmark  precomputed validation data (scatter points + correlations).
  GET  /health     liveness + cached-mutant count.

Run::

    uvicorn api.main:app --reload --port 8000

NOTE: /triage runs live AutoDock Vina docking against the resistance panel —
that is minutes on CPU. For an interactive UI, wire it to the GPU Uni-Dock
backend (see docs/UNIDOCK_SETUP.md) or dock a smaller panel.
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402
from fastapi import FastAPI, HTTPException  # noqa: E402
from fastapi.middleware.cors import CORSMiddleware  # noqa: E402
from pydantic import BaseModel  # noqa: E402

import config  # noqa: E402
from services.docking import dock_against_panel  # noqa: E402
from services.explanation import build_structural_context, generate_explanation  # noqa: E402
from services.mutation_panel import load_panel  # noqa: E402
from services.scoring import compute_robustness_scores  # noqa: E402

app = FastAPI(title="ResistScope", version="0.1.0",
              description="Resistance-aware triage for HIV-1 protease inhibitors")
app.add_middleware(
    CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"],
)

# Benchmark drug -> canonical SMILES (for the dropdown / demo).
BENCHMARK_SMILES_PATH = config.DATA_DIR / "ligands" / "benchmark_drugs.json"


class TriageRequest(BaseModel):
    smiles: str
    target: str = "HIV1_PR"
    subset: str = "primary"   # "primary" | "prevalence" | "all"


def _resistance_panel(subset: str = "primary", min_isolates: int = 20) -> pd.DataFrame:
    """Union resistance panel across all drugs (one row per unique mutation).

    Aggregates the per-drug panels into a single deduplicated panel used to
    triage a novel compound: max prevalence and any is_primary across drugs,
    mean measured fold-resistance. Filtered to the requested subset.
    """
    frames = []
    for parquet in sorted(config.PANELS_DIR.glob("*.parquet")):
        frames.append(load_panel(parquet.stem, config.PANELS_DIR))
    allp = pd.concat(frames, ignore_index=True)
    agg = allp.groupby("mutation").agg(
        position=("position", "first"),
        wildtype_aa=("wildtype_aa", "first"),
        mutant_aa=("mutant_aa", "first"),
        mean_log_fold_resistance=("mean_log_fold_resistance", "mean"),
        n_isolates=("n_isolates", "max"),
        is_primary=("is_primary", "any"),
    ).reset_index()
    if subset == "primary":
        agg = agg[agg["is_primary"]]
    elif subset == "prevalence":
        agg = agg[agg["n_isolates"] >= min_isolates]
    return agg.reset_index(drop=True)


def _severity(ddg: float) -> str:
    if ddg is None:
        return "unknown"
    if ddg >= config.DDG_DANGER_THRESHOLD:
        return "high"
    if ddg >= config.DDG_WARNING_THRESHOLD:
        return "medium"
    return "low"


@app.get("/health")
def health() -> dict:
    n_mutants = len(list(config.MUTANTS_DIR.glob("*.pdbqt")))
    wt_ready = (config.STRUCTURES_DIR / "wildtype.pdbqt").exists()
    return {"status": "ok", "n_mutants_cached": n_mutants, "wildtype_ready": wt_ready}


@app.get("/drugs")
def drugs() -> dict:
    """Benchmark PIs available for instant (precomputed) triage."""
    import json
    smiles = {}
    if BENCHMARK_SMILES_PATH.exists():
        smiles = json.loads(BENCHMARK_SMILES_PATH.read_text())
    out = []
    for d, name in config.PI_DRUGS.items():
        if (config.PANELS_DIR / f"{d}.parquet").exists():
            out.append({"abbrev": d, "name": name, "smiles": smiles.get(d)})
    return {"drugs": out}


@app.get("/drug/{abbrev}")
def drug_precomputed(abbrev: str) -> dict:
    """Instant triage result for a benchmark PI from precomputed docking."""
    abbrev = abbrev.upper()
    dock_path = config.DOCKING_DIR / "benchmark_docking.parquet"
    if abbrev not in config.PI_DRUGS or not dock_path.exists():
        raise HTTPException(404, f"No precomputed results for {abbrev}.")
    dock = pd.read_parquet(dock_path)
    drug_dock = dock[dock["drug"] == abbrev]
    if drug_dock.empty:
        raise HTTPException(404, f"No docking rows for {abbrev}.")
    panel = load_panel(abbrev, config.PANELS_DIR)
    return _assemble_result(abbrev, drug_dock, panel, precomputed=True)


@app.get("/benchmark")
def benchmark() -> dict:
    """Precomputed validation: per-mutation scatter points + overall correlation."""
    merged_p = config.VALIDATION_DIR / "scores_vs_fold_resistance.parquet"
    corr_p = config.VALIDATION_DIR / "validation_correlations.parquet"
    if not merged_p.exists() or not corr_p.exists():
        raise HTTPException(404, "Validation data not found. Run scripts/05_validate.py.")
    merged = pd.read_parquet(merged_p)
    corr = pd.read_parquet(corr_p)
    overall = corr[(corr.drug == "OVERALL") & (corr.scoring_method == "docking_ddg")]
    points = merged[["drug", "mutation", "delta_delta_g",
                     "mean_log_fold_resistance", "is_primary"]].to_dict("records")
    faithfulness = None
    faith_p = config.VALIDATION_DIR / "faithfulness_scores.parquet"
    if faith_p.exists():
        f = pd.read_parquet(faith_p)
        faithfulness = {
            "n": int(len(f)),
            "mean": float(f["score"].mean()),
            "pct_correct": float(100 * (f["score"] == 2).mean()),
            "pct_ok": float(100 * (f["score"] >= 1).mean()),
            "dist": {str(k): int(v) for k, v in
                     f["score"].value_counts().reindex([0, 1, 2], fill_value=0).items()},
        }
    return {
        "points": points,
        "overall_spearman_rho": float(overall["spearman_rho"].iloc[0]) if len(overall) else None,
        "overall_spearman_pvalue": float(overall["spearman_pvalue"].iloc[0]) if len(overall) else None,
        "n_mutations": int(len(merged)),
        "faithfulness": faithfulness,
    }


@app.get("/validation/plot")
def validation_plot():
    """Serve the precomputed validation figure (PNG)."""
    from fastapi.responses import FileResponse
    p = config.VALIDATION_DIR / "validation_plot.png"
    if not p.exists():
        raise HTTPException(404, "validation_plot.png not found. Run scripts/05_validate.py.")
    return FileResponse(str(p), media_type="image/png")


def _assemble_result(drug_name, results, panel, precomputed: bool) -> dict:
    """Build the triage response from a docking-results DataFrame.

    ``precomputed=True`` loads cached explanations ({drug}_{mut}.json); otherwise
    it generates explanations for the top-5 worst mutations in parallel (live).
    """
    import json

    wt_row = results[results["mutation"] == "WT"]
    wt_dg = (float(wt_row["delta_g"].iloc[0])
             if len(wt_row) and wt_row["delta_g"].notna().any() else None)
    scores = compute_robustness_scores(results, panel)

    muts = results[results["mutation"] != "WT"].dropna(subset=["delta_delta_g"]).copy()
    muts = muts.sort_values("delta_delta_g", ascending=False)
    primary = set(panel[panel["is_primary"]]["mutation"]) if "is_primary" in panel else set()

    explanations, citations = {}, {}
    if precomputed:
        for mut in muts["mutation"]:
            p = config.EXPLANATIONS_DIR / f"{drug_name}_{mut}.json"
            if p.exists():
                rec = json.loads(p.read_text())
                explanations[mut] = rec.get("explanation")
                citations[mut] = rec.get("citations") or []
    else:
        def _explain(row) -> tuple:
            ctx = build_structural_context(
                row["mutation"], {"delta_delta_g": float(row["delta_delta_g"]),
                                  "delta_g": float(row["delta_g"])})
            try:
                expl = generate_explanation(drug_name, row["mutation"],
                                            float(row["delta_delta_g"]), ctx)
            except Exception as exc:  # noqa: BLE001
                expl = f"(explanation unavailable: {exc})"
            return row["mutation"], expl
        top5 = muts.head(5)
        if len(top5):
            with ThreadPoolExecutor(max_workers=5) as pool:
                for mut, expl in pool.map(_explain, [r for _, r in top5.iterrows()]):
                    explanations[mut] = expl

    mutation_rows = [{
        "mutation": row["mutation"],
        "delta_g": float(row["delta_g"]) if pd.notna(row["delta_g"]) else None,
        "delta_delta_g": float(row["delta_delta_g"]),
        "severity": _severity(float(row["delta_delta_g"])),
        "is_primary": row["mutation"] in primary,
        "explanation": explanations.get(row["mutation"]),
        "citations": citations.get(row["mutation"], []),
    } for _, row in muts.iterrows()]

    return {
        "drug": drug_name,
        "robustness_score": scores.get("robustness_0_100"),
        "wildtype_binding": wt_dg,
        "mutations": mutation_rows,
        "scoring_methods": {
            "simple_mean": scores.get("simple_mean_ddg"),
            "weighted": scores.get("prevalence_weighted_ddg"),
            "worst_case": scores.get("worst_case_ddg"),
        },
        "n_mutations_scored": scores.get("n_mutations_scored"),
        "n_mutations_failed": scores.get("n_mutations_failed"),
        "precomputed": precomputed,
    }


@app.post("/triage")
def triage(req: TriageRequest) -> dict:
    """Dock a compound against the resistance panel, score it, explain the worst."""
    panel = _resistance_panel(req.subset)
    try:
        results = dock_against_panel(req.smiles, "QUERY", panel)
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(400, f"Docking failed: {exc}")
    return _assemble_result("QUERY", results, panel, precomputed=False)
