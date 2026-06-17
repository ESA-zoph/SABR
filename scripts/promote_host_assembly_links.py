from __future__ import annotations

import argparse
from datetime import date
import gzip
from pathlib import Path
import re
import sys
import time
import urllib.parse
import urllib.request
import urllib.error
import json

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.accession_linkage import (
    ACCESSION_LINKAGE_COLUMNS,
    load_accession_linkage_table,
    validate_accession_linkage_table,
)
from scripts.download_accession_linkage_genomes import (
    DOWNLOADED_RECORD_COLUMNS,
    _fasta_stats,
    _load_downloaded_records,
)


ASSEMBLY_CANDIDATE_COLUMNS = [
    "source_key",
    "bacterium",
    "strain",
    "query",
    "candidate_rank",
    "assembly_uid",
    "assembly_accession",
    "assembly_name",
    "assembly_level",
    "organism",
    "ftp_path",
    "reported_strain",
    "review_status",
    "notes",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Promote exact host strain assembly links and download their genome FASTA files."
    )
    parser.add_argument(
        "--interactions",
        type=Path,
        default=Path("data/curation/phage_host_interactions.tsv"),
    )
    parser.add_argument(
        "--coverage",
        type=Path,
        default=Path("data/curation/accession_linkage_coverage.tsv"),
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
        "--candidates",
        type=Path,
        default=Path("data/curation/host_assembly_candidates.tsv"),
    )
    parser.add_argument(
        "--downloads-root",
        type=Path,
        default=Path("data/curation/downloads/bacteria"),
    )
    parser.add_argument("--max-strains", type=int, default=200)
    parser.add_argument("--retmax", type=int, default=5)
    parser.add_argument("--sleep", type=float, default=0.35)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write candidates but do not update linkage/downloaded-record tables.",
    )
    args = parser.parse_args()

    interactions = pd.read_csv(args.interactions, sep="\t", dtype=str).fillna("")
    coverage = pd.read_csv(args.coverage, sep="\t", dtype=str).fillna("")
    linkage = load_accession_linkage_table(args.linkage)
    downloaded = _load_downloaded_records(args.downloaded_records)

    unresolved_ids = set(
        coverage.loc[
            ~coverage["bacterium_genome_linked"].map(_truthy),
            "interaction_id",
        ].astype(str)
    )
    unresolved = interactions[interactions["interaction_id"].astype(str).isin(unresolved_ids)]
    strains = (
        unresolved[["source_key", "bacterium", "strain"]]
        .drop_duplicates()
        .head(args.max_strains)
    )

    candidate_rows = []
    promoted_rows = []
    existing_keys = {
        (
            row["entity_type"],
            row["source_key"],
            row["display_name"],
            row["strain_or_isolate"],
        )
        for _, row in linkage.iterrows()
    }
    existing_downloads = set(downloaded["accession"].astype(str))

    args.downloads_root.mkdir(parents=True, exist_ok=True)

    for _, strain_row in strains.iterrows():
        query = _query_for(strain_row["bacterium"], strain_row["strain"])
        ids = _esearch_assembly(query, retmax=args.retmax)
        if not ids:
            candidate_rows.append(_candidate_row(strain_row, query, rank=0))
            time.sleep(args.sleep)
            continue

        summaries = _esummary_assembly(ids)
        accepted_for_strain = False
        for rank, uid in enumerate(ids, start=1):
            summary = summaries.get(uid, {})
            ftp_path = _best_ftp_path(summary)
            report = _download_assembly_report(ftp_path) if ftp_path else ""
            reported_strain = _reported_strain(report)
            accepted, notes = _review_candidate(
                expected_bacterium=strain_row["bacterium"],
                expected_strain=strain_row["strain"],
                summary=summary,
                report=report,
                accepted_for_strain=accepted_for_strain,
            )
            candidate_rows.append(
                _candidate_row(
                    strain_row,
                    query=query,
                    rank=rank,
                    uid=uid,
                    summary=summary,
                    ftp_path=ftp_path,
                    reported_strain=reported_strain,
                    review_status="promoted_exact_assembly" if accepted else "rejected",
                    notes=notes,
                )
            )
            if not accepted:
                continue

            accepted_for_strain = True
            key = (
                "bacterium",
                strain_row["source_key"],
                strain_row["bacterium"],
                strain_row["strain"],
            )
            if key in existing_keys:
                continue

            accession = str(summary.get("assemblyaccession", "")).strip()
            local_path = ""
            if not args.dry_run:
                local_path = _download_assembly_fasta(
                    ftp_path=ftp_path,
                    display_name=strain_row["bacterium"],
                    strain=strain_row["strain"],
                    accession=accession,
                    downloads_root=args.downloads_root,
                )
                stats = _fasta_stats(Path(local_path))
                if accession not in existing_downloads:
                    downloaded = pd.concat(
                        [
                            downloaded,
                            pd.DataFrame(
                                [
                                    {
                                        "record_type": "bacterium",
                                        "name": f"{strain_row['bacterium']} {strain_row['strain']}",
                                        "accession": accession,
                                        "local_path": local_path,
                                        "record_count": stats["record_count"],
                                        "total_bp": stats["total_bp"],
                                        "source_url": _assembly_fasta_url(ftp_path),
                                        "download_date": date.today().isoformat(),
                                        "notes": notes,
                                    }
                                ],
                                columns=DOWNLOADED_RECORD_COLUMNS,
                            ),
                        ],
                        ignore_index=True,
                    )
                    existing_downloads.add(accession)

            promoted_rows.append(
                {
                    "linkage_id": _linkage_id(strain_row["source_key"], strain_row["strain"]),
                    "entity_type": "bacterium",
                    "source_key": strain_row["source_key"],
                    "display_name": strain_row["bacterium"],
                    "strain_or_isolate": strain_row["strain"],
                    "accession": accession,
                    "accession_database": "GenBank",
                    "assembly_level": _assembly_level(summary),
                    "sequence_status": "available",
                    "linkage_status": "exact",
                    "confidence": "high",
                    "local_path": local_path,
                    "notes": notes,
                }
            )
            existing_keys.add(key)
            time.sleep(args.sleep)

    candidates = pd.DataFrame(candidate_rows, columns=ASSEMBLY_CANDIDATE_COLUMNS)
    args.candidates.parent.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(args.candidates, sep="\t", index=False)

    if promoted_rows and not args.dry_run:
        linkage = pd.concat(
            [linkage, pd.DataFrame(promoted_rows, columns=ACCESSION_LINKAGE_COLUMNS)],
            ignore_index=True,
        )
        linkage = linkage[ACCESSION_LINKAGE_COLUMNS]
        validate_accession_linkage_table(linkage)
        linkage.to_csv(args.linkage, sep="\t", index=False)
        downloaded = downloaded[DOWNLOADED_RECORD_COLUMNS]
        downloaded.to_csv(args.downloaded_records, sep="\t", index=False)

    print(f"queried_strains\t{len(strains)}")
    print(f"candidate_rows\t{len(candidates)}")
    print(f"promoted_exact_assemblies\t{len(promoted_rows)}")
    print(f"dry_run\t{args.dry_run}")


def _query_for(bacterium: str, strain: str) -> str:
    clean_strain = str(strain).replace("_", " ")
    return f'"{bacterium}"[Organism] AND "{clean_strain}"[All Fields]'


def _esearch_assembly(query: str, retmax: int) -> list[str]:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "assembly", "term": query, "retmode": "json", "retmax": str(retmax)}
    )
    data = json.loads(_download_text(url))
    return data["esearchresult"].get("idlist", [])


def _esummary_assembly(ids: list[str]) -> dict[str, dict[str, object]]:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "assembly", "id": ",".join(ids), "retmode": "json"}
    )
    data = json.loads(_download_text(url))["result"]
    return {uid: data[uid] for uid in data.get("uids", [])}


def _download_assembly_report(ftp_path: str) -> str:
    if not ftp_path:
        return ""
    base = ftp_path.rstrip("/").rsplit("/", 1)[-1]
    url = ftp_path.replace("ftp://", "https://").rstrip("/") + f"/{base}_assembly_report.txt"
    return _download_text(url)


def _download_assembly_fasta(
    ftp_path: str,
    display_name: str,
    strain: str,
    accession: str,
    downloads_root: Path,
) -> str:
    url = _assembly_fasta_url(ftp_path)
    path = downloads_root / f"{_safe_filename(display_name, strain, accession)}.fasta"
    if path.exists():
        return path.as_posix()
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "SABR-curation/0.1 (host-assembly-download)"},
    )
    payload = b""
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read()
            break
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(2 * attempt)
    fasta = gzip.decompress(payload).decode("utf-8")
    if not fasta.startswith(">"):
        raise ValueError(f"Downloaded assembly {accession} is not FASTA")
    path.write_text(fasta, encoding="utf-8")
    return path.as_posix()


def _assembly_fasta_url(ftp_path: str) -> str:
    base = ftp_path.rstrip("/").rsplit("/", 1)[-1]
    return ftp_path.replace("ftp://", "https://").rstrip("/") + f"/{base}_genomic.fna.gz"


def _download_text(url: str) -> str:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "SABR-curation/0.1 (host-assembly-resolution)"},
    )
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(2 * attempt)
    raise RuntimeError(f"Unable to download {url}")


def _best_ftp_path(summary: dict[str, object]) -> str:
    return str(summary.get("ftppath_refseq") or summary.get("ftppath_genbank") or "").strip()


def _reported_strain(report: str) -> str:
    match = re.search(r"^# Infraspecific name:\s+strain=(.+)$", report, flags=re.MULTILINE)
    return match.group(1).strip() if match else ""


def _review_candidate(
    expected_bacterium: str,
    expected_strain: str,
    summary: dict[str, object],
    report: str,
    accepted_for_strain: bool,
) -> tuple[bool, str]:
    if accepted_for_strain:
        return False, "Rejected: another exact assembly was already promoted for this strain."
    accession = str(summary.get("assemblyaccession", "")).strip()
    ftp_path = _best_ftp_path(summary)
    if not accession or not ftp_path:
        return False, "Rejected: assembly accession or FTP path is missing."
    reported_strain = _reported_strain(report)
    if not _same_strain(expected_strain, reported_strain):
        return False, f"Rejected: assembly report strain '{reported_strain}' does not match."
    organism = str(summary.get("organism", "")).lower()
    if expected_bacterium.lower() not in organism:
        return False, f"Rejected: organism '{summary.get('organism', '')}' does not match."
    if "Genome representation: full" not in report:
        return False, "Rejected: assembly report does not state full genome representation."
    return True, (
        f"Exact host assembly promoted from NCBI Assembly on {date.today().isoformat()}; "
        f"assembly report gives strain={reported_strain} and full genome representation."
    )


def _same_strain(expected: object, reported: object) -> bool:
    return _normalize_strain(expected) == _normalize_strain(reported)


def _normalize_strain(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def _assembly_level(summary: dict[str, object]) -> str:
    raw = str(summary.get("assemblylevel") or summary.get("assemblystatus") or "").lower()
    if "complete" in raw:
        return "complete_genome"
    if "chromosome" in raw:
        return "chromosome"
    if "scaffold" in raw:
        return "scaffold"
    if "contig" in raw:
        return "contig"
    return "unknown"


def _linkage_id(source_key: object, strain: object) -> str:
    base = f"host_assembly_{source_key}_{strain}"
    return re.sub(r"[^a-z0-9]+", "_", base.lower()).strip("_")


def _safe_filename(display_name: object, strain: object, accession: str) -> str:
    base = f"{display_name}_{strain}_{accession}"
    return re.sub(r"[^A-Za-z0-9._-]+", "_", base).strip("_")


def _candidate_row(
    row: pd.Series,
    query: str,
    rank: int,
    uid: str = "",
    summary: dict[str, object] | None = None,
    ftp_path: str = "",
    reported_strain: str = "",
    review_status: str = "no_candidate",
    notes: str = "",
) -> dict[str, object]:
    summary = summary or {}
    return {
        "source_key": row["source_key"],
        "bacterium": row["bacterium"],
        "strain": row["strain"],
        "query": query,
        "candidate_rank": rank,
        "assembly_uid": uid,
        "assembly_accession": summary.get("assemblyaccession", ""),
        "assembly_name": summary.get("assemblyname", ""),
        "assembly_level": summary.get("assemblylevel") or summary.get("assemblystatus", ""),
        "organism": summary.get("organism", ""),
        "ftp_path": ftp_path,
        "reported_strain": reported_strain,
        "review_status": review_status,
        "notes": notes,
    }


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


if __name__ == "__main__":
    main()
