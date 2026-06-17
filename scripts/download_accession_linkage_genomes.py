from __future__ import annotations

import argparse
from datetime import date
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

from crispr_phage_predictor.accession_linkage import (
    ACCESSION_LINKAGE_COLUMNS,
    load_accession_linkage_table,
    validate_accession_linkage_table,
)


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DOWNLOADED_RECORD_COLUMNS = [
    "record_type",
    "name",
    "accession",
    "local_path",
    "record_count",
    "total_bp",
    "source_url",
    "download_date",
    "notes",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download FASTA records referenced by accession_linkage.tsv."
    )
    parser.add_argument(
        "--linkage",
        type=Path,
        default=Path("data/curation/accession_linkage.tsv"),
    )
    parser.add_argument(
        "--downloaded-records",
        type=Path,
        default=Path("data/curation/downloaded_records.tsv"),
    )
    parser.add_argument(
        "--downloads-root",
        type=Path,
        default=Path("data/curation/downloads"),
    )
    parser.add_argument("--sleep", type=float, default=0.35)
    args = parser.parse_args()

    linkage = load_accession_linkage_table(args.linkage)
    downloaded = _load_downloaded_records(args.downloaded_records)
    existing_accessions = set(downloaded["accession"].astype(str))
    new_records = []

    for index, row in linkage.iterrows():
        accession = str(row.get("accession", "")).strip()
        entity_type = str(row.get("entity_type", "")).strip()
        if not accession or entity_type == "cocktail":
            continue
        if str(row.get("local_path", "")).strip():
            continue
        target_dir = args.downloads_root / ("phages" if entity_type == "phage" else "bacteria")
        target_dir.mkdir(parents=True, exist_ok=True)
        path = target_dir / f"{_safe_filename(row['display_name'], row['strain_or_isolate'], accession)}.fasta"
        url = _efetch_url(accession)
        if not path.exists():
            fasta = _download_text(url)
            if not fasta.startswith(">"):
                raise ValueError(f"Downloaded record for {accession} is not FASTA")
            path.write_text(fasta, encoding="utf-8")
            time.sleep(args.sleep)
        stats = _fasta_stats(path)
        local_path = path.as_posix()
        linkage.at[index, "local_path"] = local_path
        if accession not in existing_accessions:
            new_records.append(
                {
                    "record_type": entity_type,
                    "name": _record_name(row),
                    "accession": accession,
                    "local_path": local_path,
                    "record_count": stats["record_count"],
                    "total_bp": stats["total_bp"],
                    "source_url": url,
                    "download_date": date.today().isoformat(),
                    "notes": row.get("notes", ""),
                }
            )
            existing_accessions.add(accession)

    linkage = linkage[ACCESSION_LINKAGE_COLUMNS]
    validate_accession_linkage_table(linkage)
    linkage.to_csv(args.linkage, sep="\t", index=False)

    if new_records:
        downloaded = pd.concat(
            [downloaded, pd.DataFrame(new_records, columns=DOWNLOADED_RECORD_COLUMNS)],
            ignore_index=True,
        )
        downloaded.to_csv(args.downloaded_records, sep="\t", index=False)
    print(f"downloaded_or_linked\t{len(new_records)}")


def _load_downloaded_records(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=DOWNLOADED_RECORD_COLUMNS)
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    for column in DOWNLOADED_RECORD_COLUMNS:
        if column not in table.columns:
            table[column] = ""
    return table[DOWNLOADED_RECORD_COLUMNS]


def _efetch_url(accession: str) -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }
    return f"{EUTILS_BASE}?{urllib.parse.urlencode(params)}"


def _download_text(url: str) -> str:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "SABR-curation/0.1 (genome-linkage-download)"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read().decode("utf-8")


def _fasta_stats(path: Path) -> dict[str, int]:
    record_count = 0
    total_bp = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            record_count += 1
        else:
            total_bp += len(line.strip())
    return {"record_count": record_count, "total_bp": total_bp}


def _record_name(row: pd.Series) -> str:
    if row["entity_type"] == "phage":
        return str(row["display_name"])
    strain = str(row["strain_or_isolate"])
    if strain == "*":
        return f"{row['display_name']} reference proxy"
    return f"{row['display_name']} {strain}"


def _safe_filename(display_name: object, strain_or_isolate: object, accession: str) -> str:
    base = f"{display_name}_{strain_or_isolate}_{accession}"
    return re.sub(r"[^A-Za-z0-9._-]+", "_", base).strip("_")


if __name__ == "__main__":
    main()
