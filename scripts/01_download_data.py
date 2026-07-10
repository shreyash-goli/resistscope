"""01: Download the target's genotype-phenotype dataset and its PDB structure.

Downloads two inputs into ``data/raw/`` for the chosen target:

1. The Stanford HIVdb genotype-phenotype dataset (protease PI, or RT NNRTI).
2. The wildtype receptor PDB (protease 3OXC, or RT 3V81) from RCSB.

The Stanford landing page is a JavaScript app that does not embed the file
links in its HTML, so we cannot reliably scrape it. We still *try* to parse it
for a matching .txt/.tsv link, then fall back to the target's canonical URLs.

Run from anywhere::

    python scripts/01_download_data.py              # protease (default)
    python scripts/01_download_data.py --target rt  # reverse transcriptase (NNRTIs)
"""

import argparse
import re
import sys
from pathlib import Path

import requests

# Make the project root importable so ``import config`` works regardless of the
# current working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config  # noqa: E402

# --- Sources -----------------------------------------------------------------

# Landing page for the Rhee 2006 PNAS analysis (mostly JS-rendered).
RHEE_LANDING_PAGE = config.RHEE_DATASET_URL

REQUEST_TIMEOUT = 60  # seconds
HEADERS = {"User-Agent": "ResistScope/0.1 (HIV drug-resistance triage)"}


# --- Helpers -----------------------------------------------------------------

def _human(num_bytes: int) -> str:
    """Format a byte count as a short human-readable string."""
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} GB"


def download_file(url: str, dest: Path, min_bytes: int = 1) -> bool:
    """Stream ``url`` to ``dest``, printing progress.

    Returns True on success (HTTP 200 and >= ``min_bytes`` written), False
    otherwise. Never raises for an ordinary HTTP/connection failure — it prints
    the problem and returns False so callers can try the next candidate.
    """
    print(f"  GET {url}")
    try:
        with requests.get(
            url, stream=True, timeout=REQUEST_TIMEOUT, headers=HEADERS
        ) as resp:
            if resp.status_code != 200:
                print(f"    -> HTTP {resp.status_code}, skipping")
                return False

            total = int(resp.headers.get("Content-Length", 0))
            dest.parent.mkdir(parents=True, exist_ok=True)
            written = 0
            with open(dest, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=64 * 1024):
                    if not chunk:
                        continue
                    fh.write(chunk)
                    written += len(chunk)
                    if total:
                        pct = 100 * written / total
                        print(
                            f"\r    {_human(written)} / {_human(total)} "
                            f"({pct:5.1f}%)",
                            end="",
                            flush=True,
                        )
                    else:
                        print(
                            f"\r    {_human(written)} downloaded",
                            end="",
                            flush=True,
                        )
            print()  # newline after the progress line
    except requests.RequestException as exc:
        print(f"    -> request failed: {exc}")
        return False

    if written < min_bytes:
        print(f"    -> only {written} bytes written (expected >= {min_bytes})")
        return False

    print(f"    -> saved {_human(written)} to {dest}")
    return True


def find_dataset_links_on_page(url: str, match_token: str) -> list[str]:
    """Fetch the landing page and return absolute URLs to matching .txt/.tsv files.

    Best-effort: returns an empty list if the page can't be fetched or has no
    links whose href contains ``match_token`` (the page is JS-rendered, so this
    usually finds nothing and we fall back to the target's canonical URLs).
    """
    try:
        resp = requests.get(url, timeout=REQUEST_TIMEOUT, headers=HEADERS)
        resp.raise_for_status()
    except requests.RequestException as exc:
        print(f"  Could not fetch landing page ({exc}); using fallback URLs.")
        return []

    token = match_token.lower()
    hrefs = re.findall(r'href=["\']([^"\']+)["\']', resp.text, flags=re.I)
    links: list[str] = []
    for href in hrefs:
        low = href.lower()
        if (token in low) and (low.endswith(".txt") or low.endswith(".tsv")):
            links.append(requests.compat.urljoin(url, href))
    # De-duplicate while preserving order.
    seen: set[str] = set()
    unique = [x for x in links if not (x in seen or seen.add(x))]
    if unique:
        print(f"  Found {len(unique)} candidate link(s) on the page.")
    return unique


# --- Downloads ---------------------------------------------------------------

def download_dataset(raw_dir: Path, target) -> Path:
    """Download the target's genotype-phenotype dataset to ``raw_dir``."""
    print(f"[1/2] {target.label} genotype-phenotype dataset")
    dest = raw_dir / target.dataset_filename

    # These datasets are hundreds of KB; guard against a truncated/HTML error page.
    min_bytes = 50_000

    # Prefer any link discovered on the page, then the target's canonical URLs.
    token = target.dataset_filename.split("_")[0]  # "PI" / "NNRTI"
    candidates = find_dataset_links_on_page(RHEE_LANDING_PAGE, token) + list(target.dataset_urls)
    for url in candidates:
        if download_file(url, dest, min_bytes=min_bytes):
            return dest

    raise ConnectionError(
        f"Failed to download the {target.name} dataset from all known sources. "
        f"Check network access or update {target.name}.dataset_urls in targets.py."
    )


def download_pdb_structure(raw_dir: Path, target) -> Path:
    """Download the target's wildtype PDB structure to ``raw_dir``."""
    print(f"[2/2] PDB structure {target.pdb_id}")
    dest = raw_dir / f"{target.pdb_id}.pdb"
    # A real PDB file is tens/hundreds of KB; a 404 stub would be far smaller.
    if download_file(target.pdb_url, dest, min_bytes=10_000):
        return dest
    raise ConnectionError(
        f"Failed to download PDB {target.pdb_id} from {target.pdb_url}"
    )


def verify_nonempty(path: Path) -> None:
    """Raise if ``path`` is missing or empty."""
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"Downloaded file is missing or empty: {path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target", default="HIV1_PR",
        help="Docking target: HIV1_PR / pr (default) or HIV1_RT / rt.",
    )
    args = parser.parse_args()
    t = config.set_active_target(args.target)

    raw_dir = config.RAW_DIR
    raw_dir.mkdir(parents=True, exist_ok=True)
    print(f"Target: {t.name} ({t.label})")
    print(f"Download dir: {raw_dir}\n")

    data_path = download_dataset(raw_dir, t)
    print()
    pdb_path = download_pdb_structure(raw_dir, t)
    print()

    print("Verifying downloads...")
    for path in (data_path, pdb_path):
        verify_nonempty(path)
        print(f"  OK  {path}  ({_human(path.stat().st_size)})")

    print("\nDone. Both files downloaded and non-empty.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
