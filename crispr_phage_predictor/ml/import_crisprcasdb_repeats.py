from __future__ import annotations

import argparse
import zipfile
from io import TextIOWrapper
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from crispr_phage_predictor.io import extract_accession, sequence_hash


IUPAC_DNA_BASES = set("ACGTRYSWKMBDHVN")
SABR_REPEAT_FEATURE_BASES = set("ACGTN")

CRISPRCASDB_REPEAT_COLUMNS = [
    "source",
    "release",
    "fasta_member",
    "record_id",
    "description",
    "first_accession",
    "accession_count",
    "repeat_sequence",
    "repeat_length",
    "sequence_hash",
    "valid_iupac_dna",
    "usable_for_sabr_repeat_features",
]


def import_crisprcasdb_direct_repeats(
    zip_path: str | Path,
    member: str = "direct_repeat_seqName.fsa",
    release: str = "34",
    only_usable: bool = False,
    max_records: int | None = None,
) -> pd.DataFrame:
    """Import CRISPRCasdb direct-repeat FASTA records as an unlabeled inventory.

    The direct-repeat FASTA exports do not by themselves provide a reliable
    Cas type/subtype label, so this function intentionally does not emit the
    SABR training-table schema. Use the output as provenance/audit input or as
    an intermediate table for later SQL-derived Cas-cluster joins.
    """
    rows = []
    with zipfile.ZipFile(zip_path) as archive:
        if member not in archive.namelist():
            available = ", ".join(archive.namelist())
            raise ValueError(f"{member!r} not found in {zip_path}. Available members: {available}")
        with archive.open(member) as raw_handle:
            text_handle = TextIOWrapper(raw_handle, encoding="utf-8", errors="replace")
            for index, record in enumerate(SeqIO.parse(text_handle, "fasta"), start=1):
                repeat = str(record.seq).upper().strip()
                accession_values = _extract_accessions(record.id, record.description)
                usable = _is_usable_for_sabr_repeat_features(repeat)
                if only_usable and not usable:
                    continue
                rows.append(
                    {
                        "source": "crisprcasdb_direct_repeat_fasta",
                        "release": str(release),
                        "fasta_member": member,
                        "record_id": record.id,
                        "description": record.description,
                        "first_accession": accession_values[0] if accession_values else "",
                        "accession_count": len(accession_values),
                        "repeat_sequence": repeat,
                        "repeat_length": len(repeat),
                        "sequence_hash": sequence_hash(repeat),
                        "valid_iupac_dna": set(repeat).issubset(IUPAC_DNA_BASES),
                        "usable_for_sabr_repeat_features": usable,
                    }
                )
                if max_records is not None and index >= max_records:
                    break
    return pd.DataFrame(rows, columns=CRISPRCASDB_REPEAT_COLUMNS)


def _extract_accessions(record_id: str, description: str) -> list[str]:
    values: list[str] = []
    seen: set[str] = set()
    for text in [record_id, description]:
        for part in str(text or "").replace(";", "+").split("+"):
            accession = extract_accession(part)
            if accession and accession not in seen:
                seen.add(accession)
                values.append(accession)
    return values


def _is_usable_for_sabr_repeat_features(repeat: str) -> bool:
    return 23 <= len(repeat) <= 47 and set(repeat).issubset(SABR_REPEAT_FEATURE_BASES)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create an unlabeled CRISPRCasdb direct-repeat inventory from a FASTA ZIP export."
    )
    parser.add_argument("zip_path", help="Path to CRISPRCasdb dr_34.zip or equivalent direct-repeat FASTA ZIP.")
    parser.add_argument(
        "--member",
        default="direct_repeat_seqName.fsa",
        help="FASTA member inside the ZIP to import.",
    )
    parser.add_argument(
        "--release",
        default="34",
        help="CRISPRCasdb release/version label to write into the output.",
    )
    parser.add_argument(
        "--output",
        default="data/training/crisprcasdb_34_direct_repeats_inventory.csv",
        help="Output CSV path.",
    )
    parser.add_argument(
        "--only-usable",
        action="store_true",
        help="Keep only repeats compatible with current SABR repeat features.",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=None,
        help="Optional cap for quick inspections.",
    )
    args = parser.parse_args()

    table = import_crisprcasdb_direct_repeats(
        args.zip_path,
        member=args.member,
        release=args.release,
        only_usable=args.only_usable,
        max_records=args.max_records,
    )
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    print(f"Wrote {len(table)} direct-repeat inventory rows to {output_path}")


if __name__ == "__main__":
    main()
