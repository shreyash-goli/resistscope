#!/usr/bin/env python3
"""Smoke test for the ResistScope API — no docking stack, no GPU, no API key.

Confirms the read-only, precomputed surface a reviewer actually touches boots
and serves cached results:

    GET /health          backend is up (reports the docking backend it resolved)
    GET /targets         registry lists HIV-1 protease (+ any BYO targets)
    GET /benchmark       rigorous validation metrics are present (n_pairs, enrichment)
    GET /drug/DRV        one precomputed drug returns a scored resistance panel

It deliberately does NOT hit /triage (that runs live docking, minutes on CPU) or
anything that needs meeko/RDKit/Vina or an Anthropic key. Start the server, then
run this:

    python -m uvicorn api.main:app --port 8000      # terminal 1
    python scripts/smoke_test.py                     # terminal 2
    python scripts/smoke_test.py --url http://host:8000

Exit code 0 = all checks passed, 1 = something failed (CI-friendly).
"""

import argparse
import json
import sys
import urllib.error
import urllib.request

DEFAULT_URL = "http://localhost:8000"
TIMEOUT = 30

# Green-check / red-x without pulling in a color lib.
OK, BAD = "\033[32m✓\033[0m", "\033[31m✗\033[0m"


class CheckFailed(Exception):
    """A single check's assertion did not hold."""


def _get(base: str, path: str) -> dict:
    """GET base+path and return parsed JSON, or raise CheckFailed with context."""
    url = base.rstrip("/") + path
    try:
        with urllib.request.urlopen(url, timeout=TIMEOUT) as resp:  # noqa: S310
            body = resp.read().decode()
    except urllib.error.HTTPError as exc:
        raise CheckFailed(f"{path} -> HTTP {exc.code} {exc.reason}")
    except urllib.error.URLError as exc:
        raise CheckFailed(
            f"cannot reach {url} ({exc.reason}). "
            f"Is the server running? `python -m uvicorn api.main:app --port 8000`"
        )
    try:
        return json.loads(body)
    except json.JSONDecodeError:
        raise CheckFailed(f"{path} -> response was not JSON: {body[:120]!r}")


def check_health(base: str) -> str:
    d = _get(base, "/health")
    if d.get("status") != "ok":
        raise CheckFailed(f"status was {d.get('status')!r}, expected 'ok'")
    return f"backend={d.get('docking_backend')} live_docking={d.get('live_docking')}"


def check_targets(base: str) -> str:
    d = _get(base, "/targets")
    names = {t["name"] for t in d.get("targets", [])}
    if "HIV1_PR" not in names:
        raise CheckFailed(f"HIV1_PR missing from targets: {sorted(names)}")
    return f"{len(names)} targets ({', '.join(sorted(names))}), default={d.get('default')}"


def check_benchmark(base: str) -> str:
    d = _get(base, "/benchmark")
    rig = d.get("rigorous")
    if not rig or "n_pairs" not in rig:
        raise CheckFailed("no rigorous benchmark metrics (run scripts/08_benchmark.py)")
    enr = {e["top_n"]: e["enrichment"] for e in rig.get("enrichment", [])}
    if not enr:
        raise CheckFailed("benchmark metrics present but enrichment table is empty")
    top40 = enr.get(40)
    detail = f"n_pairs={rig['n_pairs']}"
    if top40 is not None:
        detail += f", top-40 DRM enrichment={top40:.2f}x"
    return detail


def check_drug(base: str, abbrev: str) -> str:
    d = _get(base, f"/drug/{abbrev}")
    muts = d.get("mutations")
    if not muts:
        raise CheckFailed(f"{abbrev} returned no mutations")
    if d.get("robustness_score") is None:
        raise CheckFailed(f"{abbrev} has no robustness_score")
    return f"{abbrev} robustness={d['robustness_score']:.0f}/100 over {len(muts)} mutations"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--url", default=DEFAULT_URL, help=f"API base URL (default {DEFAULT_URL})")
    ap.add_argument("--drug", default="DRV", help="precomputed benchmark drug to check (default DRV)")
    args = ap.parse_args()

    print(f"ResistScope smoke test -> {args.url}\n")
    checks = [
        ("GET /health", lambda: check_health(args.url)),
        ("GET /targets", lambda: check_targets(args.url)),
        ("GET /benchmark", lambda: check_benchmark(args.url)),
        (f"GET /drug/{args.drug}", lambda: check_drug(args.url, args.drug)),
    ]

    failed = 0
    for label, fn in checks:
        try:
            detail = fn()
            print(f"  {OK} {label:22s} {detail}")
        except CheckFailed as exc:
            failed += 1
            print(f"  {BAD} {label:22s} {exc}")

    print()
    if failed:
        print(f"{BAD} {failed}/{len(checks)} checks failed")
        return 1
    print(f"{OK} all {len(checks)} checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
