from __future__ import annotations

import argparse
import json
import time
import urllib.parse
import urllib.request
from http.client import IncompleteRead
from dataclasses import dataclass
from pathlib import Path
from urllib.error import URLError

import pandas as pd
from Bio import SeqIO

from crispr_phage_predictor.external.minced import detect_arrays_with_minced, minced_available
from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.ml.dataset import (
    REPEAT_CAS_DATASET_COLUMNS,
    cas_type_from_subtype,
    validate_repeat_cas_training_table,
)


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass(frozen=True)
class GenBankSubtypeProfile:
    subtype: str
    query: str
    required_any: tuple[tuple[str, ...], ...]


DEFAULT_PROFILES = [
    GenBankSubtypeProfile("I-E", "cse cas3 complete genome bacteria", (("cas3",), ("cse", "cas8e"))),
    GenBankSubtypeProfile("I-F", "csy cas3 complete genome bacteria", (("cas3",), ("csy", "cas8f"))),
    GenBankSubtypeProfile("II-A", "cas9 csn2 complete genome bacteria", (("cas9",), ("csn2",))),
    GenBankSubtypeProfile("III-A", "cas10 csm complete genome bacteria", (("cas10",), ("csm",))),
    GenBankSubtypeProfile("III-B", "cas10 cmr complete genome bacteria", (("cas10",), ("cmr",))),
    GenBankSubtypeProfile("V-A", "cas12a complete genome bacteria OR cpf1 complete genome bacteria", (("cas12a", "cpf1"),)),
    GenBankSubtypeProfile("VI-A", "cas13a complete genome bacteria OR c2c2 complete genome bacteria", (("cas13a", "c2c2"),)),
]


def collect_genbank_annotated_training_table(
    output_dir: str | Path,
    max_records_per_profile: int = 10,
    email: str | None = None,
    profiles: list[GenBankSubtypeProfile] | None = None,
) -> pd.DataFrame:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    seen_accessions: set[str] = set()

    for profile in profiles or DEFAULT_PROFILES:
        ids = _candidate_ids_for_profile(
            profile.query,
            max_records=max_records_per_profile,
            email=email,
        )
        for ncbi_id in ids:
            record = _fetch_genbank_record(ncbi_id, output_dir, email=email)
            accession = record.annotations.get("accessions", [record.id])[0]
            if accession in seen_accessions:
                continue
            seen_accessions.add(accession)

            matched_profiles = [
                candidate
                for candidate in DEFAULT_PROFILES
                if _record_matches_profile(record, candidate)
            ]
            if [candidate.subtype for candidate in matched_profiles] != [profile.subtype]:
                continue

            arrays = _detect_record_arrays(record, genome_id=accession)
            for array in arrays:
                rows.append(
                    {
                        "source": "genbank_signature_candidate",
                        "genome_id": accession,
                        "organism": record.annotations.get("organism", ""),
                        "taxonomy": ";".join(record.annotations.get("taxonomy", [])),
                        "assembly_level": "GenBank/RefSeq annotated record",
                        "contig_id": record.id,
                        "array_start": array.start,
                        "array_end": array.end,
                        "repeat_sequence": array.repeat_consensus,
                        "repeat_length": array.repeat_length,
                        "spacer_count": array.spacer_count,
                        "mean_spacer_length": round(array.mean_spacer_length, 6),
                        "cas_type": cas_type_from_subtype(profile.subtype),
                        "cas_subtype": profile.subtype,
                        "label_source": "GenBank_annotation_signature_single_profile_match",
                        "label_confidence": "genbank_signature_candidate",
                        "pam_rule": "",
                    }
                )
            time.sleep(0.12)

    table = pd.DataFrame(rows, columns=REPEAT_CAS_DATASET_COLUMNS)
    validate_repeat_cas_training_table(table)
    return table


def _candidate_ids_for_profile(query: str, max_records: int, email: str | None) -> list[str]:
    ids = _esearch(query, retmax=max(50, max_records * 10), email=email)
    summaries = _esummary(ids, email=email)
    candidates = [
        ncbi_id
        for ncbi_id in ids
        if _looks_like_complete_genomic_record(summaries.get(ncbi_id, {}))
    ]
    return candidates[:max_records]


def _esearch(query: str, retmax: int, email: str | None) -> list[str]:
    params = {
        "db": "nuccore",
        "term": query,
        "retmode": "json",
        "retmax": str(retmax),
    }
    if email:
        params["email"] = email
    url = f"{EUTILS_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    data = json.loads(_read_url(url, timeout=60).decode("utf-8"))
    return data["esearchresult"].get("idlist", [])


def _esummary(ids: list[str], email: str | None, batch_size: int = 100) -> dict[str, dict]:
    if not ids:
        return {}
    summaries: dict[str, dict] = {}
    for start in range(0, len(ids), batch_size):
        batch_ids = ids[start : start + batch_size]
        params = {
            "db": "nuccore",
            "id": ",".join(batch_ids),
            "retmode": "json",
        }
        if email:
            params["email"] = email
        url = f"{EUTILS_BASE}/esummary.fcgi?{urllib.parse.urlencode(params)}"
        data = json.loads(_read_url(url, timeout=60).decode("utf-8"))["result"]
        summaries.update(
            {ncbi_id: data[ncbi_id] for ncbi_id in batch_ids if ncbi_id in data}
        )
        time.sleep(0.12)
    return summaries


def _looks_like_complete_genomic_record(summary: dict) -> bool:
    title = str(summary.get("title", "")).lower()
    length = int(summary.get("slen", 0) or 0)
    if not (100_000 <= length <= 10_000_000):
        return False
    excluded = ["shotgun", "project", "mag:", "tpa_asm", "assembly", "scaffold"]
    if any(term in title for term in excluded):
        return False
    complete_terms = ["complete genome", "complete sequence", "chromosome, complete"]
    return any(term in title for term in complete_terms)


def _fetch_genbank_record(ncbi_id: str, output_dir: Path, email: str | None):
    path = output_dir / f"{ncbi_id}.gb"
    if not path.exists():
        params = {
            "db": "nuccore",
            "id": ncbi_id,
            "rettype": "gbwithparts",
            "retmode": "text",
        }
        if email:
            params["email"] = email
        url = f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
        temp_path = path.with_suffix(".tmp")
        temp_path.write_bytes(_read_url(url, timeout=180))
        temp_path.replace(path)
    return SeqIO.read(path, "genbank")


def _read_url(url: str, timeout: int, attempts: int = 3) -> bytes:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return response.read()
        except (IncompleteRead, TimeoutError, URLError) as error:
            last_error = error
            time.sleep(1.0 * attempt)
    assert last_error is not None
    raise last_error


def _record_matches_profile(record, profile: GenBankSubtypeProfile) -> bool:
    annotation_text = _record_annotation_text(record)
    return all(
        any(term.lower() in annotation_text for term in term_group)
        for term_group in profile.required_any
    )


def _record_annotation_text(record) -> str:
    parts = [record.description, record.annotations.get("organism", "")]
    for feature in record.features:
        if feature.type != "CDS":
            continue
        for key in ["gene", "product", "note", "locus_tag"]:
            parts.extend(feature.qualifiers.get(key, []))
    return " ".join(parts).lower()


def _detect_record_arrays(record, genome_id: str):
    if not minced_available():
        return []
    fasta_record = FastaRecord(
        source_file=genome_id,
        record_id=record.id,
        description=record.description,
        sequence=str(record.seq).upper(),
    )
    return detect_arrays_with_minced([fasta_record])


def _selected_profiles(profile_text: str) -> list[GenBankSubtypeProfile] | None:
    requested = {item.strip() for item in profile_text.split(",") if item.strip()}
    if not requested:
        return None
    selected = [profile for profile in DEFAULT_PROFILES if profile.subtype in requested]
    missing = sorted(requested - {profile.subtype for profile in selected})
    if missing:
        raise ValueError("Unknown profiles: " + ", ".join(missing))
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Collect candidate repeat/Cas rows from GenBank records with subtype signature annotations."
    )
    parser.add_argument(
        "--genbank-dir",
        default="data/training/genbank_sources",
        help="Directory for downloaded GenBank records.",
    )
    parser.add_argument(
        "--output",
        default="data/training/repeats_cas_types_genbank_candidate.csv",
        help="Output training CSV path.",
    )
    parser.add_argument(
        "--max-records-per-profile",
        type=int,
        default=10,
        help="Maximum NCBI records to inspect per subtype profile.",
    )
    parser.add_argument("--email", default=None, help="Optional email for NCBI E-utilities.")
    parser.add_argument(
        "--profiles",
        default="",
        help="Optional comma-separated subtype profiles to run, e.g. III-A,V-A,VI-A.",
    )
    args = parser.parse_args()

    table = collect_genbank_annotated_training_table(
        output_dir=args.genbank_dir,
        max_records_per_profile=args.max_records_per_profile,
        email=args.email,
        profiles=_selected_profiles(args.profiles),
    )
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    print(f"Wrote {len(table)} GenBank candidate rows to {output_path}")


if __name__ == "__main__":
    main()
