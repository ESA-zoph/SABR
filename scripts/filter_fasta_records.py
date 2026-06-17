from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO


def filter_fasta_records(
    input_path: Path,
    output_path: Path,
    min_length: int,
    max_records: int | None = None,
) -> tuple[int, int]:
    records = []
    total = 0
    for record in SeqIO.parse(input_path, "fasta"):
        total += 1
        if len(record.seq) >= min_length:
            records.append(record)

    records.sort(key=lambda record: len(record.seq), reverse=True)
    if max_records is not None:
        records = records[:max_records]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, output_path, "fasta")
    return total, len(records)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter a fragmented FASTA to longer records for SABR testing."
    )
    parser.add_argument("input", type=Path, help="Input FASTA path.")
    parser.add_argument("output", type=Path, help="Output FASTA path.")
    parser.add_argument(
        "--min-length",
        type=int,
        default=10_000,
        help="Minimum record length to keep. Default: 10000.",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=500,
        help="Maximum number of longest records to keep. Use 0 for no cap. Default: 500.",
    )
    args = parser.parse_args()

    max_records = None if args.max_records == 0 else args.max_records
    total, kept = filter_fasta_records(
        input_path=args.input,
        output_path=args.output,
        min_length=args.min_length,
        max_records=max_records,
    )
    print(
        f"Read {total:,} FASTA records; kept {kept:,} records "
        f"with length >= {args.min_length:,} bp."
    )
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
