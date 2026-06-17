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
    alignment_length: int
    spacer_length: int
    coverage: float
    evalue: float | None
    bitscore: float | None
    spacer_sequence: str
    protospacer_sequence: str
    aligned_spacer_sequence: str = ""
    aligned_protospacer_sequence: str = ""
    protospacer_5p_flank: str = ""
    protospacer_3p_flank: str = ""
    genomic_upstream_flank: str = ""
    genomic_downstream_flank: str = ""
    predicted_cas_subtype: str = ""
    cas_subtype_confidence: float | None = None
    cas_subtype_prediction_method: str = ""
    pam_rule: str = ""
    pam_rule_source: str = ""
    pam_sequence: str = ""
    pam_match: bool | None = None
    pam_support_level: str = "not_evaluated"
    pam_compatibility_score: float | None = None
    pam_offset_from_protospacer: int | None = None
    seed_region: str = ""
    seed_length: int | None = None
    seed_mismatches: int | None = None
    seed_mismatch_positions: str = ""


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
        start = match_index + 1
        end = match_index + len(query)
        context = extract_protospacer_context(
            phage_sequence=phage_sequence,
            start=start,
            end=end,
            strand=strand,
        )
        hits.append(
            SpacerHit(
                bacterium_id=bacterium_id,
                array_id=array_id,
                phage_id=phage_id,
                spacer_id=spacer_id,
                phage_contig_id=phage_contig_id,
                start=start,
                end=end,
                strand=strand,
                identity=1.0,
                mismatches=0,
                alignment_length=len(query),
                spacer_length=len(spacer_sequence),
                coverage=1.0,
                evalue=None,
                bitscore=None,
                spacer_sequence=spacer_sequence,
                protospacer_sequence=phage_sequence[match_index : match_index + len(query)],
                protospacer_5p_flank=context.protospacer_5p_flank,
                protospacer_3p_flank=context.protospacer_3p_flank,
                genomic_upstream_flank=context.genomic_upstream_flank,
                genomic_downstream_flank=context.genomic_downstream_flank,
            )
        )
        start_index = match_index + 1
    return hits


@dataclass(frozen=True)
class ProtospacerContext:
    protospacer_5p_flank: str
    protospacer_3p_flank: str
    genomic_upstream_flank: str
    genomic_downstream_flank: str


def extract_protospacer_context(
    phage_sequence: str,
    start: int,
    end: int,
    strand: str,
    flank_length: int = 10,
) -> ProtospacerContext:
    """Return flanks around a one-based inclusive protospacer interval.

    Genomic flanks are reported in phage-reference orientation. Protospacer
    flanks are oriented relative to the strand matched by the spacer, which is
    the orientation needed for PAM/PFS checks.
    """
    sequence = phage_sequence.upper()
    zero_based_start = max(start - 1, 0)
    zero_based_end = min(end, len(sequence))
    upstream_start = max(zero_based_start - flank_length, 0)
    downstream_end = min(zero_based_end + flank_length, len(sequence))

    genomic_upstream = sequence[upstream_start:zero_based_start]
    genomic_downstream = sequence[zero_based_end:downstream_end]

    if strand == "-":
        protospacer_5p = reverse_complement(genomic_downstream)
        protospacer_3p = reverse_complement(genomic_upstream)
    else:
        protospacer_5p = genomic_upstream
        protospacer_3p = genomic_downstream

    return ProtospacerContext(
        protospacer_5p_flank=protospacer_5p,
        protospacer_3p_flank=protospacer_3p,
        genomic_upstream_flank=genomic_upstream,
        genomic_downstream_flank=genomic_downstream,
    )


@dataclass(frozen=True)
class SeedMismatchSummary:
    seed_region: str
    seed_length: int
    seed_mismatches: int
    seed_mismatch_positions: str


def summarize_seed_mismatches(
    spacer_sequence: str,
    protospacer_sequence: str,
    pam_rule: str | None,
    seed_length: int = 8,
) -> SeedMismatchSummary | None:
    side = _pam_rule_side(pam_rule)
    if not side:
        return None
    spacer = spacer_sequence.upper().replace("-", "")
    protospacer = protospacer_sequence.upper().replace("-", "")
    comparable_length = min(len(spacer), len(protospacer))
    if comparable_length <= 0:
        return None
    effective_seed_length = min(seed_length, comparable_length)

    if side == "5prime":
        seed_start = 0
        seed_end = effective_seed_length
    else:
        seed_start = comparable_length - effective_seed_length
        seed_end = comparable_length

    spacer_seed = spacer[seed_start:seed_end]
    protospacer_seed = protospacer[seed_start:seed_end]
    mismatch_positions = [
        str(index)
        for index, (spacer_base, protospacer_base) in enumerate(
            zip(spacer_seed, protospacer_seed),
            start=1,
        )
        if spacer_base != protospacer_base
    ]
    return SeedMismatchSummary(
        seed_region=f"{side}:{seed_start + 1}-{seed_end}",
        seed_length=effective_seed_length,
        seed_mismatches=len(mismatch_positions),
        seed_mismatch_positions=",".join(mismatch_positions),
    )


def _pam_rule_side(pam_rule: str | None) -> str:
    if not pam_rule or ":" not in pam_rule:
        return ""
    side = pam_rule.split(":", maxsplit=1)[0].strip().upper()
    if side in {"5", "5P", "5PRIME", "UPSTREAM"}:
        return "5prime"
    if side in {"3", "3P", "3PRIME", "DOWNSTREAM"}:
        return "3prime"
    if side in {"GENOMIC_5", "GENOMIC5", "GENOMIC_5P", "GENOMIC_5PRIME", "GENOMIC_UPSTREAM"}:
        return "5prime"
    if side in {"GENOMIC_3", "GENOMIC3", "GENOMIC_3P", "GENOMIC_3PRIME", "GENOMIC_DOWNSTREAM"}:
        return "3prime"
    return ""


def reverse_complement(sequence: str) -> str:
    translation = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(translation)[::-1].upper()
