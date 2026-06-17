from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.accession_linkage import load_accession_linkage_table


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download GenBank records for linked phage accessions."
    )
    parser.add_argument(
        "--linkage",
        type=Path,
        default=Path("data/curation/accession_linkage.tsv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/curation/downloads/phages_genbank"),
    )
    parser.add_argument("--sleep", type=float, default=0.35)
    args = parser.parse_args()

    linkage = load_accession_linkage_table(args.linkage)
    phages = linkage[
        (linkage["entity_type"] == "phage")
        & (linkage["accession"].astype(str).str.strip() != "")
    ].copy()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = 0
    existing = 0
    for accession in sorted(set(phages["accession"].astype(str))):
        path = args.output_dir / f"{_safe_filename(accession)}.gb"
        if path.exists() and path.stat().st_size > 0:
            existing += 1
            continue
        text = _download_text(_efetch_url(accession))
        if "LOCUS" not in text[:200]:
            raise ValueError(f"Downloaded record for {accession} does not look like GenBank")
        path.write_text(text, encoding="utf-8")
        downloaded += 1
        time.sleep(args.sleep)
    print(f"downloaded\t{downloaded}")
    print(f"existing\t{existing}")


def _efetch_url(accession: str) -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gb",
        "retmode": "text",
    }
    return f"{EUTILS_BASE}?{urllib.parse.urlencode(params)}"


def _download_text(url: str) -> str:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "SABR-curation/0.1 (phage-genbank-download)"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read().decode("utf-8")


def _safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")


if __name__ == "__main__":
    main()
