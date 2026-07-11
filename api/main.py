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
from services.explanation import build_structural_context, generate_explanation  # noqa: E402
from services.mutation_panel import load_panel  # noqa: E402
from services.scoring import compute_robustness_scores  # noqa: E402

# NOTE: services.docking imports the heavy scientific stack (meeko / RDKit /
# AutoDock Vina) at module load. Only live /triage needs it, so it is imported
# lazily inside that handler — the read-only endpoints (/benchmark, /drug,
# precomputed results) then boot without the docking env installed.

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


class IntakeRequest(BaseModel):
    pdb_id: str | None = None      # fetch from RCSB
    pdb_text: str | None = None    # or an uploaded PDB


class AssembleRequest(BaseModel):
    protein: str
    hint: str = ""


class SaveTargetRequest(BaseModel):
    target_id: str
    label: str
    pdb_id: str | None = None
    pdb_text: str | None = None
    ligand_hetcode: str | None = None
    docking_center: list[float]
    chains: list[str]
    mutate_chains: list[str]
    drugs: list[dict]              # [{abbrev, name}]
    mutations: list[str]
    binding_site: str = ""
    validatable: bool = False
    citations: list = []


def _target(name: str):
    """Resolve a Target by name/alias per request, or 400 on an unknown name."""
    try:
        return config.get_target(name)
    except KeyError as exc:
        raise HTTPException(400, str(exc))


def _panel_from_mutations(mutations) -> pd.DataFrame:
    """Build a minimal triage panel from a bare mutation list (BYO targets).

    User targets have no genotype-phenotype dataset, so their "panel" is the
    agent-proposed resistance-mutation set. No measured fold-resistance is
    available, so ``mean_log_fold_resistance`` is NaN and every row is primary.
    """
    rows = []
    for m in sorted(mutations):
        try:
            wt, pos, mut = m[0], int(m[1:-1]), m[-1]
        except (ValueError, IndexError):
            continue
        rows.append({"mutation": m, "position": pos, "wildtype_aa": wt, "mutant_aa": mut,
                     "mean_log_fold_resistance": float("nan"), "n_isolates": 0,
                     "is_primary": True})
    return pd.DataFrame(rows)


def _resistance_panel(subset: str = "primary", min_isolates: int = 20,
                      panels_dir=None, target=None) -> pd.DataFrame:
    """Union resistance panel across all drugs (one row per unique mutation).

    Aggregates the per-drug panels into a single deduplicated panel used to
    triage a novel compound: max prevalence and any is_primary across drugs,
    mean measured fold-resistance. Filtered to the requested subset. For a
    user/BYO target with no panels, falls back to the target's mutation list.
    """
    if panels_dir is None:
        panels_dir = target.panels_dir if target is not None else config.PANELS_DIR
    parquets = sorted(panels_dir.glob("*.parquet")) if panels_dir.exists() else []
    if not parquets:
        if target is not None and target.primary_mutations:
            return _panel_from_mutations(target.primary_mutations)
        raise HTTPException(404, "No resistance panel for this target.")
    frames = []
    for parquet in parquets:
        frames.append(load_panel(parquet.stem, panels_dir))
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


@app.get("/targets")
def targets_available() -> dict:
    """List configured targets: curated (PI/RT) and user-added (BYO) alike."""
    out = []
    for name, t in config.TARGETS.items():
        mutants = t.mutants_dir.exists() and any(t.mutants_dir.glob("*.pdbqt"))
        dock_ready = (t.docking_dir / "benchmark_docking.parquet").exists() or mutants
        out.append({
            "name": name, "label": t.label,
            "n_panels": len(list(t.panels_dir.glob("*.parquet"))) if t.panels_dir.exists() else 0,
            "docking_ready": dock_ready,
            "user": t.is_user,
            "n_mutations": len(t.primary_mutations),
            "validatable": bool(t.dataset_urls),  # curated dataset wired
        })
    return {"targets": out, "default": config.ACTIVE_TARGET.name}


def _fetch_pdb(pdb_id: str) -> str:
    """Download a PDB from RCSB (raises HTTPException on failure)."""
    import urllib.request
    try:
        with urllib.request.urlopen(  # noqa: S310
                f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", timeout=60) as r:
            return r.read().decode()
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(400, f"Could not fetch PDB {pdb_id}: {exc}")


@app.post("/targets/intake")
def targets_intake(req: IntakeRequest) -> dict:
    """Step 1 of add-target: parse a PDB into a structural draft (chains, pocket)."""
    from services.pdb_intake import parse_pdb
    text = req.pdb_text or (_fetch_pdb(req.pdb_id) if req.pdb_id else None)
    if not text:
        raise HTTPException(400, "Provide pdb_id or pdb_text.")
    return {"intake": parse_pdb(text)}


@app.post("/targets/assemble")
def targets_assemble(req: AssembleRequest) -> dict:
    """Step 2 of add-target: Claude agent proposes drugs + resistance mutations."""
    from services.target_builder import assemble_target
    res = assemble_target(req.protein, hint=req.hint)
    if res.error:
        raise HTTPException(502, f"Agent failed: {res.error}")
    return {"spec": res.spec, "citations": res.citations, "n_searches": res.n_searches}


@app.post("/targets/save")
def targets_save(req: SaveTargetRequest) -> dict:
    """Step 3 of add-target: persist the confirmed draft + hot-add to the registry."""
    import json
    import targets as T
    tid = req.target_id.upper().replace(" ", "_").replace("-", "_")
    raw_id = (req.pdb_id or tid).upper()
    raw_path = config.RAW_DIR / f"{raw_id}.pdb"
    if req.pdb_text:
        # Uploaded content may be mmCIF; store as PDB so the viewer + ligand
        # extraction + docking pipeline stay PDB-only.
        from services.pdb_intake import to_pdb_text
        raw_path.write_text(to_pdb_text(req.pdb_text))
    elif req.pdb_id and not raw_path.exists():
        raw_path.write_text(_fetch_pdb(req.pdb_id))

    draft = {
        "target_id": tid, "label": req.label, "pdb_id": raw_id, "source_pdb": raw_id,
        "ligand_hetcode": req.ligand_hetcode, "docking_center": req.docking_center,
        "chains": req.chains, "mutate_chains": req.mutate_chains,
        "drugs": req.drugs, "mutations": req.mutations,
        "biology": {"binding_site": req.binding_site,
                    "dataset": {"exists": req.validatable}},
        "citations": req.citations, "provenance": "byo-app",
        "status": "draft — build receptors on the GPU worker to enable docking",
    }
    T.USER_TARGETS_DIR.mkdir(parents=True, exist_ok=True)
    (T.USER_TARGETS_DIR / f"{tid}.json").write_text(json.dumps(draft, indent=2))
    t = T.target_from_draft(draft)
    T.TARGETS[t.name] = t  # hot-add so it appears immediately in /targets
    return {"saved": t.name, "label": t.label, "n_mutations": len(t.primary_mutations)}


@app.get("/health")
def health(target: str = "HIV1_PR") -> dict:
    from services.docking_backend import get_backend
    t = _target(target)
    n_mutants = len(list(t.mutants_dir.glob("*.pdbqt"))) if t.mutants_dir.exists() else 0
    wt_ready = (t.structures_dir / "wildtype.pdbqt").exists()
    backend = get_backend()
    return {"status": "ok", "target": t.name,
            "n_mutants_cached": n_mutants, "wildtype_ready": wt_ready,
            "docking_backend": backend.label, "live_docking": backend.available}


@app.get("/drugs")
def drugs(target: str = "HIV1_PR") -> dict:
    """Benchmark drugs available for instant (precomputed) triage, for a target."""
    import json
    t = _target(target)
    smiles = {}
    if BENCHMARK_SMILES_PATH.exists():
        smiles = json.loads(BENCHMARK_SMILES_PATH.read_text())
    out = []
    for d, name in t.drugs.items():
        # Curated targets list drugs that have a built panel; user targets have
        # no genotype panels, so list all inhibitors the agent proposed.
        if t.is_user or (t.panels_dir / f"{d}.parquet").exists():
            out.append({"abbrev": d, "name": name, "smiles": smiles.get(d)})
    return {"target": t.name, "drugs": out}


@app.get("/drug/{abbrev}")
def drug_precomputed(abbrev: str, target: str = "HIV1_PR") -> dict:
    """Instant triage result for a benchmark drug from precomputed docking."""
    t = _target(target)
    abbrev = abbrev.upper()
    dock_path = t.docking_dir / "benchmark_docking.parquet"
    if abbrev not in t.drugs or not dock_path.exists():
        raise HTTPException(404, f"No precomputed results for {abbrev} ({t.name}).")
    dock = pd.read_parquet(dock_path)
    drug_dock = dock[dock["drug"] == abbrev]
    if drug_dock.empty:
        raise HTTPException(404, f"No docking rows for {abbrev}.")
    panel = load_panel(abbrev, t.panels_dir)
    return _assemble_result(abbrev, drug_dock, panel, precomputed=True, target=t)


@app.get("/benchmark")
def benchmark(target: str = "HIV1_PR") -> dict:
    """Precomputed validation: per-mutation scatter points + overall correlation."""
    t = _target(target)
    merged_p = t.validation_dir / "scores_vs_fold_resistance.parquet"
    corr_p = t.validation_dir / "validation_correlations.parquet"
    if not merged_p.exists() or not corr_p.exists():
        raise HTTPException(404, f"No validation data for {t.name}. Run scripts/05_validate.py --target {target}.")
    merged = pd.read_parquet(merged_p)
    corr = pd.read_parquet(corr_p)
    overall = corr[(corr.drug == "OVERALL") & (corr.scoring_method == "docking_ddg")]
    points = merged[["drug", "mutation", "delta_delta_g",
                     "mean_log_fold_resistance", "is_primary"]].to_dict("records")
    faithfulness = None
    faith_p = t.validation_dir / "faithfulness_scores.parquet"
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
    # Rigorous benchmark (permutation p-values, bootstrap CIs, ROC/PR-AUC,
    # de-confounding analysis) from scripts/08_benchmark.py, if present.
    import json
    rigorous = None
    rig_p = t.validation_dir / "benchmark_metrics.json"
    if rig_p.exists():
        rigorous = json.loads(rig_p.read_text())

    return {
        "target": t.name,
        "points": points,
        "overall_spearman_rho": float(overall["spearman_rho"].iloc[0]) if len(overall) else None,
        "overall_spearman_pvalue": float(overall["spearman_pvalue"].iloc[0]) if len(overall) else None,
        "n_mutations": int(len(merged)),
        "faithfulness": faithfulness,
        "rigorous": rigorous,
    }


@app.get("/validation/plot")
def validation_plot(target: str = "HIV1_PR"):
    """Serve the precomputed validation figure (PNG) for a target."""
    from fastapi.responses import FileResponse
    t = _target(target)
    p = t.validation_dir / "validation_plot.png"
    if not p.exists():
        raise HTTPException(404, f"validation_plot.png not found for {t.name}. "
                                 f"Run scripts/05_validate.py --target {target}.")
    return FileResponse(str(p), media_type="image/png")


def _assemble_result(drug_name, results, panel, precomputed: bool, target=None) -> dict:
    """Build the triage response from a docking-results DataFrame.

    ``precomputed=True`` loads cached explanations ({drug}_{mut}.json); otherwise
    it generates explanations for the top-5 worst mutations in parallel (live).
    Reads/writes explanations from ``target``'s dir (defaults to the active one).
    """
    import json
    t = target if target is not None else config.ACTIVE_TARGET

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
            p = t.explanations_dir / f"{drug_name}_{mut}.json"
            if p.exists():
                rec = json.loads(p.read_text())
                explanations[mut] = rec.get("explanation")
                citations[mut] = rec.get("citations") or []
    else:
        def _explain(row) -> tuple:
            ctx = build_structural_context(
                row["mutation"], {"delta_delta_g": float(row["delta_delta_g"]),
                                  "delta_g": float(row["delta_g"])}, target=t)
            try:
                expl = generate_explanation(drug_name, row["mutation"],
                                            float(row["delta_delta_g"]), ctx, target=t)
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


@app.get("/structure/receptor")
def structure_receptor(target: str = "HIV1_PR"):
    """Serve the cleaned wildtype receptor PDB (ligand stripped) for the viewer."""
    from fastapi.responses import PlainTextResponse
    t = config.get_target(target)
    pdb = t.structures_dir / "wildtype.pdb"
    if pdb.exists():
        return PlainTextResponse(pdb.read_text(), media_type="chemical/x-pdb")
    # Fallback (e.g. a user target before its receptors are built): serve the
    # raw deposited PDB so the viewer still shows the structure + pocket.
    raw = config.RAW_DIR / f"{t.pdb_id}.pdb"
    if raw.exists():
        return PlainTextResponse(raw.read_text(), media_type="chemical/x-pdb")
    raise HTTPException(404, f"No structure for {t.name}. Run scripts/03 or upload a PDB.")


@app.get("/structure/ligand")
def structure_ligand(target: str = "HIV1_PR"):
    """Serve the co-crystal ligand as a minimal PDB (one altloc), from the raw PDB.

    The cleaned receptor has the ligand stripped, so the pocket-filling inhibitor
    is extracted here from the deposited structure (protease saquinavir/ROC,
    RT nevirapine/NVP) for overlay in the 3D viewer.
    """
    from fastapi.responses import PlainTextResponse
    t = config.get_target(target)
    raw = config.RAW_DIR / f"{t.pdb_id}.pdb"
    if not raw.exists():
        raise HTTPException(404, f"Raw {t.pdb_id}.pdb not found. Run scripts/01.")
    codes = {c.strip().upper() for c in t.ligand_hetcodes}
    keep = []
    for line in raw.read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper()
        altloc = line[16]
        if resname in codes and altloc in (" ", "A"):
            keep.append(line)
    if not keep:
        raise HTTPException(404, f"No ligand {sorted(codes)} atoms in {t.pdb_id}.pdb.")
    body = "\n".join(keep) + "\nEND\n"
    return PlainTextResponse(body, media_type="chemical/x-pdb")


@app.post("/triage")
def triage(req: TriageRequest) -> dict:
    """Dock a compound against the resistance panel, score it, explain the worst.

    Docking runs through the configured backend (local CPU, remote GPU worker, or
    none). If no backend is available the response is a clean 503 with how-to-
    enable guidance rather than an opaque 500.
    """
    from services.docking_backend import get_backend, DockingUnavailable
    t = _target(req.target)
    if not (t.mutants_dir.exists() and any(t.mutants_dir.glob("*.pdbqt"))) \
            and get_backend().label.startswith("local"):
        raise HTTPException(503, f"{t.label} has no built mutant receptors yet "
                                 f"(run scripts/03 --target {req.target}).")
    panel = _resistance_panel(req.subset, panels_dir=t.panels_dir, target=t)
    try:
        results = get_backend().dock(req.smiles, panel, t)
    except DockingUnavailable as exc:
        raise HTTPException(503, str(exc))
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(502, f"Docking failed: {exc}")
    return _assemble_result("QUERY", results, panel, precomputed=False, target=t)
