from pathlib import Path

from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.external.blast import _parse_blast_hits
from crispr_phage_predictor.io import FastaRecord


def _maps():
    array = CrisprArray(
        array_id="b.fna|contig|array_1",
        genome_id="b.fna",
        contig_id="contig",
        start=1,
        end=100,
        repeat_consensus="GTTCACTGCCGTACAGGCAGCTTAGAAA",
        spacers=["A" * 20],
    )
    phage = FastaRecord("p.fna", "p1", "p1", "C" * 100)
    return {"spacer_1": (array, 1, "A" * 20)}, {"phage_1": phage}


def test_blast_parser_filters_by_minimum_coverage_when_not_full_query(tmp_path: Path):
    output = tmp_path / "blast.tsv"
    output.write_text(
        "spacer_1\tphage_1\t100.0\t18\t0\t0\t1\t18\t5\t22\t1e-5\t50\t" + "A" * 18 + "\n",
        encoding="utf-8",
    )
    spacer_map, phage_map = _maps()

    hits = _parse_blast_hits(
        output_path=output,
        spacer_map=spacer_map,
        phage_map=phage_map,
        min_identity=0.9,
        min_coverage=0.95,
        require_full_query=False,
    )

    assert hits == []


def test_blast_parser_accepts_high_coverage_partial_hit(tmp_path: Path):
    output = tmp_path / "blast.tsv"
    output.write_text(
        "spacer_1\tphage_1\t95.0\t19\t1\t0\t1\t19\t5\t23\t1e-5\t50\t" + "A" * 19 + "\n",
        encoding="utf-8",
    )
    spacer_map, phage_map = _maps()

    hits = _parse_blast_hits(
        output_path=output,
        spacer_map=spacer_map,
        phage_map=phage_map,
        min_identity=0.9,
        min_coverage=0.95,
        require_full_query=False,
    )

    assert len(hits) == 1
    assert hits[0].coverage == 0.95


def test_blast_parser_preserves_aligned_query_and_subject_sequences(tmp_path: Path):
    output = tmp_path / "blast.tsv"
    output.write_text(
        "spacer_1\tphage_1\t93.103\t29\t2\t0\t2\t30\t5\t33\t1e-5\t50\t"
        "GCATCAAGCACGTTCGAGTTTACTGTTTC\t"
        "GCATCAAGCACGTTCGAGTTTACTGTTTC\n",
        encoding="utf-8",
    )
    spacer_map, phage_map = _maps()

    hits = _parse_blast_hits(
        output_path=output,
        spacer_map=spacer_map,
        phage_map=phage_map,
        min_identity=0.9,
        min_coverage=0.0,
        require_full_query=False,
    )

    assert len(hits) == 1
    assert hits[0].spacer_sequence == "A" * 20
    assert hits[0].aligned_spacer_sequence == "GCATCAAGCACGTTCGAGTTTACTGTTTC"
    assert hits[0].aligned_protospacer_sequence == "GCATCAAGCACGTTCGAGTTTACTGTTTC"
