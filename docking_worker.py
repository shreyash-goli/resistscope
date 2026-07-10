"""ResistScope GPU docking worker — run this ON a machine with a GPU.

The main API (`api/main.py`) stays thin and forwards live-docking requests here.
This worker holds the heavy docking stack (RDKit / meeko / Uni-Dock) and the
built mutant receptors, and does the actual GPU docking.

Setup on the GPU box (e.g. an A100)::

    bash deploy/setup_a100.sh          # builds the Uni-Dock conda env
    conda activate resistscope_gpu
    # make sure the target's receptors exist (scripts/03 --target ...)
    python docking_worker.py           # serves on :9000

Then point the API at it::

    RESISTSCOPE_DOCKING_URL=http://<gpu-host>:9000 python -m uvicorn api.main:app

Contract (must match services/docking_backend.RemoteDockingBackend):
  POST /dock   {smiles, target, panel:[{mutation,...}], engine?} -> {results:[...]}
  GET  /health -> {status, engine, gpu, targets_ready}
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd  # noqa: E402
from fastapi import FastAPI, HTTPException  # noqa: E402
from pydantic import BaseModel  # noqa: E402

import config  # noqa: E402

app = FastAPI(title="ResistScope docking worker", version="0.1.0")


class DockRequest(BaseModel):
    smiles: str
    target: str = "HIV1_PR"
    panel: list[dict]
    engine: str = "unidock"   # GPU by default; "vina" for CPU fallback


def _gpu_visible() -> bool:
    """Best-effort check that an NVIDIA GPU is present (for /health only)."""
    import shutil
    import subprocess
    if not shutil.which("nvidia-smi"):
        return False
    try:
        return subprocess.run(["nvidia-smi"], capture_output=True, timeout=10).returncode == 0
    except Exception:  # noqa: BLE001
        return False


@app.get("/health")
def health() -> dict:
    ready = {}
    for name, t in config.TARGETS.items():
        n = len(list(t.mutants_dir.glob("*.pdbqt"))) if t.mutants_dir.exists() else 0
        ready[name] = n
    return {"status": "ok", "engine": "unidock", "gpu": _gpu_visible(),
            "targets_ready": ready}


@app.post("/dock")
def dock(req: DockRequest) -> dict:
    """Dock a SMILES against a target's resistance panel on the GPU."""
    from services.docking import dock_against_panel  # heavy stack lives here

    try:
        t = config.get_target(req.target)
    except KeyError as exc:
        raise HTTPException(400, str(exc))
    if not (t.mutants_dir.exists() and any(t.mutants_dir.glob("*.pdbqt"))):
        raise HTTPException(409, f"No built receptors for {t.name} on this worker "
                                 f"(run scripts/03_build_mutant_cache.py --target {req.target}).")

    panel = pd.DataFrame(req.panel)
    if "mutation" not in panel.columns or panel.empty:
        raise HTTPException(400, "panel must be a non-empty list of {mutation, ...} rows.")

    try:
        results = dock_against_panel(
            req.smiles, "QUERY", panel,
            structures_dir=t.structures_dir, mutants_dir=t.mutants_dir,
            backend=req.engine,
        )
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(500, f"docking failed: {exc}")
    return {"target": t.name, "n_receptors": int(len(results)),
            "results": results.to_dict("records")}


if __name__ == "__main__":
    import uvicorn
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--host", default="0.0.0.0")
    ap.add_argument("--port", type=int, default=9000)
    args = ap.parse_args()
    uvicorn.run(app, host=args.host, port=args.port)
