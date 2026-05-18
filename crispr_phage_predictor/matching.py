from __future__ import annotations

from dataclasses import dataclass

from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.io import FastaRecord


@dataclass(frozen=True)
class SpacerHit:
    bacterium_id: str
    array_id: str
    phage_id: str
    spacer_id: str
    phage_contig_id: str
    start: int
    end: int
    strand: str
    identity: float
    mismatches: int
    spacer_sequence: str
    protospacer_sequence: str


def find_spacer_hits(
    crispr_arrays: list[CrisprArray],
    phage_records: list[FastaRecord],
) -> list[SpacerHit]:
    hits: list[SpacerHit] = []
    for array in crispr_arrays:
        for spacer_index, spacer in enumerate(array.spacers, start=1):
            spacer_id = f"{array.array_id}|spacer_{spacer_index}"
            for phage_record in phage_records:
                hits.extend(
                    _find_exact_hits(
                        bacterium_id=array.genome_id,
                        array_id=array.array_id,
                        spacer_id=spacer_id,
                        spacer_sequence=spacer,
                        phage_id=phage_record.source_file,
                        phage_contig_id=phage_record.record_id,
                        phage_sequence=phage_record.sequence,
                    )
                )
    return hits


def _find_exact_hits(
    bacterium_id: str,
    array_id: str,
    spacer_id: str,
    spacer_sequence: str,
    phage_id: str,
    phage_contig_id: str,
    phage_sequence: str,
) -> list[SpacerHit]:
    hits: list[SpacerHit] = []
    query = spacer_sequence.upper()
    target = phage_sequence.upper()
    reverse_query = reverse_complement(query)

    hits.extend(
        _scan_one_strand(
            query=query,
            strand="+",
            bacterium_id=bacterium_id,
            array_id=array_id,
            spacer_id=spacer_id,
            spacer_sequence=spacer_sequence,
            phage_id=phage_id,
            phage_contig_id=phage_contig_id,
            phage_sequence=target,
        )
    )
    if reverse_query != query:
        hits.extend(
            _scan_one_strand(
                query=reverse_query,
                strand="-",
                bacterium_id=bacterium_id,
                array_id=array_id,
                spacer_id=spacer_id,
                spacer_sequence=spacer_sequence,
                phage_id=phage_id,
                phage_contig_id=phage_contig_id,
                phage_sequence=target,
            )
        )
    return hits


def _scan_one_strand(
    query: str,
    strand: str,
    bacterium_id: str,
    array_id: str,
    spacer_id: str,
    spacer_sequence: str,
    phage_id: str,
    phage_contig_id: str,
    phage_sequence: str,
) -> list[SpacerHit]:
    hits: list[SpacerHit] = []
    start_index = 0
    while True:
        match_index = phage_sequence.find(query, start_index)
        if match_index == -1:
            break
        hits.append(
            SpacerHit(
                bacterium_id=bacterium_id,
                array_id=array_id,
                phage_id=phage_id,
                spacer_id=spacer_id,
                phage_contig_id=phage_contig_id,
                start=match_index + 1,
                end=match_index + len(query),
                strand=strand,
                identity=1.0,
                mismatches=0,
                spacer_sequence=spacer_sequence,
                protospacer_sequence=phage_sequence[match_index : match_index + len(query)],
            )
        )
        start_index = match_index + 1
    return hits


def reverse_complement(sequence: str) -> str:
    translation = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(translation)[::-1].upper()
