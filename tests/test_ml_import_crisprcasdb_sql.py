from __future__ import annotations

from crispr_phage_predictor.ml.import_crisprcasdb_sql import (
    import_crisprcasdb_sql_candidate_labels,
)


def test_import_crisprcasdb_sql_candidate_labels_links_locus_to_nearest_cluster(tmp_path):
    sql_path = tmp_path / "crisprcasdb.sql"
    sql_path.write_text(
        "\n".join(
            [
                "COPY public.clustercas (id, sequence, start, length, class) FROM stdin;",
                "cluster-near\tseq-1\t2000\t3000\tCAS-TypeI-F",
                "cluster-far\tseq-1\t50000\t3000\tCAS-TypeIII-B",
                r"\.",
                "COPY public.crisprlocus (id, sequence, start, length, orientation, trusted, evidencelevel, drconsensus, drconservation, spacerconservation, potentialorientation, evidencelevelreeval, blastscore) FROM stdin;",
                "locus-1\tseq-1\t1000\t500\t1\tf\t4\trepeat-1\t99.0\t90.0\t1\t\\N\t\\N",
                r"\.",
                "COPY public.crisprlocus_region (crisprlocus, region, start, length) FROM stdin;",
                "locus-1\trepeat-1\t1000\t29",
                "locus-1\tspacer-1\t1030\t32",
                "locus-1\tspacer-2\t1092\t34",
                r"\.",
                "COPY public.region (id, sequence, category) FROM stdin;",
                "repeat-1\tGTTTCAATGCTGCTTCGCCTGCAATGGGTTTAGTAT\t1",
                "spacer-1\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t3",
                "spacer-2\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\t3",
                r"\.",
                "COPY public.sequence (id, strain, category, length, ncount, description, job) FROM stdin;",
                "seq-1\tstrain-1\t1\t100000\t0\tChromosome\t\\N",
                r"\.",
                "COPY public.strain (id, genbank, refseq, taxon, gb_release_date, release_level, assembly_status) FROM stdin;",
                "strain-1\tGCA_000001.1\tGCF_000001.1\t1\t2020-01-01 00:00:00\tmajor\tComplete Genome",
                r"\.",
            ]
        ),
        encoding="utf-8",
    )

    table = import_crisprcasdb_sql_candidate_labels(sql_path)

    assert len(table) == 1
    row = table.iloc[0]
    assert row["genome_id"] == "GCF_000001.1"
    assert row["repeat_sequence"] == "GTTTCAATGCTGCTTCGCCTGCAATGGGTTTAGTAT"
    assert row["spacer_count"] == 2
    assert row["mean_spacer_length"] == 33
    assert row["cas_type"] == "Type I"
    assert row["cas_subtype"] == "I-F"
    assert row["label_confidence"] == "computational_nearby_cas_cluster"


def test_import_crisprcasdb_sql_candidate_labels_filters_ambiguous_or_distant_clusters(tmp_path):
    sql_path = tmp_path / "crisprcasdb.sql"
    sql_path.write_text(
        "\n".join(
            [
                "COPY public.clustercas (id, sequence, start, length, class) FROM stdin;",
                "cluster-ambiguous\tseq-1\t2000\t3000\tCAS",
                "cluster-distant\tseq-2\t100000\t3000\tCAS-TypeII-A",
                r"\.",
                "COPY public.crisprlocus (id, sequence, start, length, orientation, trusted, evidencelevel, drconsensus, drconservation, spacerconservation, potentialorientation, evidencelevelreeval, blastscore) FROM stdin;",
                "locus-1\tseq-1\t1000\t500\t1\tf\t4\trepeat-1\t99.0\t90.0\t1\t\\N\t\\N",
                "locus-2\tseq-2\t1000\t500\t1\tf\t4\trepeat-2\t99.0\t90.0\t1\t\\N\t\\N",
                r"\.",
                "COPY public.crisprlocus_region (crisprlocus, region, start, length) FROM stdin;",
                "locus-1\trepeat-1\t1000\t29",
                "locus-2\trepeat-2\t1000\t29",
                r"\.",
                "COPY public.region (id, sequence, category) FROM stdin;",
                "repeat-1\tGTTTCAATGCTGCTTCGCCTGCAATGGGTTTAGTAT\t1",
                "repeat-2\tGTTTCAATGCTGCTTCGCCTGCAATGGGTTTAGTAT\t1",
                r"\.",
                "COPY public.sequence (id, strain, category, length, ncount, description, job) FROM stdin;",
                "seq-1\tstrain-1\t1\t100000\t0\tChromosome\t\\N",
                "seq-2\tstrain-1\t1\t100000\t0\tChromosome\t\\N",
                r"\.",
                "COPY public.strain (id, genbank, refseq, taxon, gb_release_date, release_level, assembly_status) FROM stdin;",
                "strain-1\tGCA_000001.1\tGCF_000001.1\t1\t2020-01-01 00:00:00\tmajor\tComplete Genome",
                r"\.",
            ]
        ),
        encoding="utf-8",
    )

    table = import_crisprcasdb_sql_candidate_labels(sql_path, max_cas_distance_bp=20_000)

    assert table.empty
