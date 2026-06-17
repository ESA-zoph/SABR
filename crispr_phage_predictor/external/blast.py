from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.matching import extract_protospacer_context
from crispr_phage_predictor.matching import SpacerHit


BLASTN_COMMAND = "blastn"
MAKEBLASTDB_COMMAND = "makeblastdb"


def find_spacer_hits_with_blast(
    crispr_arrays: list[CrisprArray],
    phage_records: list[FastaRecord],
    min_identity: float = 0.9,
    min_coverage: float = 0.95,
    require_full_query: bool = True,
) -> list[SpacerHit]:
    """Find spacer-protospacer hits with local BLAST+.

    BLAST is used as an optional approximate matcher. By default, alignments
    must span the full spacer so short local matches do not inflate evidence.
    """
    spacer_entries = _build_spacer_entries(crispr_arrays)
    if not spacer_entries or not phage_records:
        return []

    with tempfile.TemporaryDirectory(prefix="crispr_blast_") as temp_dir:
        work_dir = Path(temp_dir)
        spacer_fasta = work_dir / "spacers.fasta"
        phage_fasta = work_dir / "phages.fasta"
        db_prefix = work_dir / "phage_db"
        output_path = work_dir / "blast.tsv"

        spacer_map = _write_spacer_fasta(spacer_fasta, spacer_entries)
        phage_map = _write_phage_fasta(phage_fasta, phage_records)

        subprocess.run(
            [
                MAKEBLASTDB_COMMAND,
                "-in",
                str(phage_fasta),
                "-dbtype",
                "nucl",
                "-out",
                str(db_prefix),
                "-parse_seqids",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        subprocess.run(
            [
                BLASTN_COMMAND,
                "-task",
                "blastn-short",
                "-query",
                str(spacer_fasta),
                "-db",
                str(db_prefix),
                "-dust",
                "no",
                "-outfmt",
                (
                    "6 qseqid sseqid pident length mismatch gapopen qstart qend "
                    "sstart send evalue bitscore qseq sseq"
                ),
                "-out",
                str(output_path),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        return _parse_blast_hits(
            output_path=output_path,
            spacer_map=spacer_map,
            phage_map=phage_map,
            min_identity=min_identity,
            min_coverage=min_coverage,
            require_full_query=require_full_query,
        )


def _build_spacer_entries(crispr_arrays: list[CrisprArray]) -> list[tuple[CrisprArray, int, str]]:
    entries = []
    for array in crispr_arrays:
        for spacer_index, spacer in enumerate(array.spacers, start=1):
            entries.append((array, spacer_index, spacer))
    return entries


def _write_spacer_fasta(
    path: Path,
    spacer_entries: list[tuple[CrisprArray, int, str]],
) -> dict[str, tuple[CrisprArray, int, str]]:
    spacer_map = {}
    with path.open("w", encoding="utf-8") as handle:
        for index, (array, spacer_index, spacer) in enumerate(spacer_entries, start=1):
            query_id = f"spacer_{index}"
            spacer_map[query_id] = (array, spacer_index, spacer.upper())
            handle.write(f">{query_id}\n{spacer.upper()}\n")
    return spacer_map


def _write_phage_fasta(path: Path, phage_records: list[FastaRecord]) -> dict[str, FastaRecord]:
    phage_map = {}
    with path.open("w", encoding="utf-8") as handle:
        for index, record in enumerate(phage_records, start=1):
            subject_id = f"phage_{index}"
            phage_map[subject_id] = record
            handle.write(f">{subject_id}\n{record.sequence.upper()}\n")
    return phage_map


def _parse_blast_hits(
    output_path: Path,
    spacer_map: dict[str, tuple[CrisprArray, int, str]],
    phage_map: dict[str, FastaRecord],
    min_identity: float,
    min_coverage: float,
    require_full_query: bool,
) -> list[SpacerHit]:
    hits = []
    if not output_path.exists():
        return hits

    for line in output_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) == 14:
            (
                query_id,
                subject_id,
                pident,
                alignment_length,
                mismatch,
                _gapopen,
                _qstart,
                _qend,
                sstart,
                send,
                evalue,
                bitscore,
                query_sequence,
                subject_sequence,
            ) = fields
        elif len(fields) == 13:
            (
                query_id,
                subject_id,
                pident,
                alignment_length,
                mismatch,
                _gapopen,
                _qstart,
                _qend,
                sstart,
                send,
                evalue,
                bitscore,
                subject_sequence,
            ) = fields
            query_sequence = spacer_map[query_id][2]
        else:
            continue
        array, spacer_index, spacer_sequence = spacer_map[query_id]
        phage_record = phage_map[subject_id]
        identity = float(pident) / 100
        aligned_length = int(alignment_length)
        coverage = aligned_length / len(spacer_sequence) if spacer_sequence else 0.0
        if identity < min_identity:
            continue
        if require_full_query and aligned_length != len(spacer_sequence):
            continue
        if not require_full_query and coverage < min_coverage:
            continue

        start = int(sstart)
        end = int(send)
        strand = "+" if start <= end else "-"
        normalized_start = min(start, end)
        normalized_end = max(start, end)
        context = extract_protospacer_context(
            phage_sequence=phage_record.sequence,
            start=normalized_start,
            end=normalized_end,
            strand=strand,
        )
        hits.append(
            SpacerHit(
                bacterium_id=array.genome_id,
                array_id=array.array_id,
                phage_id=phage_record.source_file,
                spacer_id=f"{array.array_id}|spacer_{spacer_index}",
                phage_contig_id=phage_record.record_id,
                start=normalized_start,
                end=normalized_end,
                strand=strand,
                identity=identity,
                mismatches=int(mismatch),
                alignment_length=aligned_length,
                spacer_length=len(spacer_sequence),
                coverage=coverage,
                evalue=float(evalue),
                bitscore=float(bitscore),
                spacer_sequence=spacer_sequence,
                aligned_spacer_sequence=query_sequence.upper(),
                aligned_protospacer_sequence=subject_sequence.upper(),
                protospacer_sequence=subject_sequence.replace("-", "").upper(),
                protospacer_5p_flank=context.protospacer_5p_flank,
                protospacer_3p_flank=context.protospacer_3p_flank,
                genomic_upstream_flank=context.genomic_upstream_flank,
                genomic_downstream_flank=context.genomic_downstream_flank,
            )
        )
    return hits
