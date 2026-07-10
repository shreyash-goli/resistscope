"""Pluggable docking backend so the API can dock locally *or* on a remote GPU.

The API process is intentionally thin — it should not require the heavy docking
stack (RDKit / meeko / AutoDock Vina / Uni-Dock). Live ``/triage`` therefore
dispatches the actual docking through one of three backends, resolved once at
startup (env wins):

  1. ``RESISTSCOPE_DOCKING_URL`` set  → :class:`RemoteDockingBackend`
     (POST the SMILES + panel to a GPU worker the user runs; see
     ``docking_worker.py`` + ``deploy/setup_a100.sh``). **This is how a user
     "connects a GPU".**
  2. else, if the docking stack imports here → :class:`LocalDockingBackend`
     (in-process CPU Vina, or GPU Uni-Dock if ``RESISTSCOPE_DOCKING_ENGINE=unidock``).
  3. else → :class:`NullDockingBackend`, whose ``dock`` raises
     :class:`DockingUnavailable` with actionable guidance — so a misconfigured
     server returns a clean 503, not an opaque 500.

Every backend takes a SMILES + a resistance-panel DataFrame + a Target and
returns the standard docking-results DataFrame (``drug, mutation, delta_g,
delta_delta_g, n_poses, n_ok``), so the rest of the pipeline is backend-agnostic.
"""

from __future__ import annotations

import os

import pandas as pd


class DockingUnavailable(RuntimeError):
    """No usable docking backend is configured (raised by the Null backend)."""


_GUIDANCE = (
    "Live docking needs a docking backend, and none is configured on this server. "
    "Either (1) run the GPU worker on a machine with a GPU "
    "(bash deploy/setup_a100.sh, then `python docking_worker.py`) and set "
    "RESISTSCOPE_DOCKING_URL=http://<that-host>:9000 for the API, or (2) install "
    "the local docking env (conda env create -f environment.yml) and restart the "
    "API inside it. Benchmark drugs still work instantly — they are precomputed."
)


class DockingBackend:
    """Interface: dock a SMILES against a resistance panel for a target."""

    label = "abstract"

    def dock(self, smiles: str, panel_df: pd.DataFrame, target) -> pd.DataFrame:
        raise NotImplementedError

    @property
    def available(self) -> bool:
        return True


class LocalDockingBackend(DockingBackend):
    """Dock in-process using the local docking stack (needs meeko/RDKit/Vina)."""

    def __init__(self, engine: str = "vina"):
        self.engine = engine  # "vina" (CPU) or "unidock" (GPU)
        self.label = f"local:{engine}"

    def dock(self, smiles: str, panel_df: pd.DataFrame, target) -> pd.DataFrame:
        from services.docking import dock_against_panel
        return dock_against_panel(
            smiles, "QUERY", panel_df,
            structures_dir=target.structures_dir,
            mutants_dir=target.mutants_dir,
            backend=self.engine,
        )


class RemoteDockingBackend(DockingBackend):
    """Dispatch docking to a remote GPU worker over HTTP (the "connect a GPU" path)."""

    def __init__(self, url: str, timeout: int = 1800):
        self.url = url.rstrip("/")
        self.timeout = timeout
        self.label = f"remote:{self.url}"

    def dock(self, smiles: str, panel_df: pd.DataFrame, target) -> pd.DataFrame:
        import requests
        payload = {
            "smiles": smiles,
            "target": target.name,
            "panel": panel_df.to_dict("records"),
        }
        try:
            r = requests.post(f"{self.url}/dock", json=payload, timeout=self.timeout)
        except requests.RequestException as exc:
            raise DockingUnavailable(
                f"GPU docking worker at {self.url} is unreachable ({exc}). "
                f"Is it running and RESISTSCOPE_DOCKING_URL correct?"
            )
        if r.status_code >= 400:
            raise RuntimeError(f"docking worker error {r.status_code}: {r.text[:300]}")
        return pd.DataFrame(r.json()["results"])

    def health(self) -> dict:
        import requests
        r = requests.get(f"{self.url}/health", timeout=10)
        r.raise_for_status()
        return r.json()


class NullDockingBackend(DockingBackend):
    """No backend available — dock() fails loudly and helpfully."""

    label = "none"

    @property
    def available(self) -> bool:
        return False

    def dock(self, smiles: str, panel_df: pd.DataFrame, target) -> pd.DataFrame:
        raise DockingUnavailable(_GUIDANCE)


_BACKEND: DockingBackend | None = None


def get_backend(refresh: bool = False) -> DockingBackend:
    """Resolve the docking backend once (cached). See module docstring for order."""
    global _BACKEND
    if _BACKEND is not None and not refresh:
        return _BACKEND

    url = os.environ.get("RESISTSCOPE_DOCKING_URL", "").strip()
    if url:
        _BACKEND = RemoteDockingBackend(url)
        return _BACKEND

    try:  # local stack present?
        import meeko  # noqa: F401
        import rdkit  # noqa: F401
        engine = os.environ.get("RESISTSCOPE_DOCKING_ENGINE", "vina").strip() or "vina"
        _BACKEND = LocalDockingBackend(engine=engine)
    except Exception:  # noqa: BLE001 - any import failure means "no local stack"
        _BACKEND = NullDockingBackend()
    return _BACKEND
