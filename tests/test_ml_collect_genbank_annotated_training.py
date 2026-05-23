from __future__ import annotations

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from crispr_phage_predictor.ml.collect_genbank_annotated_training import (
    DEFAULT_PROFILES,
    _esummary,
    _looks_like_complete_genomic_record,
    _record_matches_profile,
)


def test_record_matches_single_subtype_signature_profile():
    record = SeqRecord(Seq("A" * 1000), id="TEST.1", description="test record")
    record.annotations["organism"] = "Example bacterium"
    record.features = [
        SeqFeature(
            FeatureLocation(10, 100),
            type="CDS",
            qualifiers={"gene": ["cas9"], "product": ["CRISPR-associated protein Cas9"]},
        ),
        SeqFeature(
            FeatureLocation(120, 200),
            type="CDS",
            qualifiers={"gene": ["csn2"], "product": ["CRISPR-associated protein Csn2"]},
        ),
    ]

    matches = [
        profile.subtype
        for profile in DEFAULT_PROFILES
        if _record_matches_profile(record, profile)
    ]

    assert matches == ["II-A"]


def test_filters_complete_genomic_records_before_download():
    assert _looks_like_complete_genomic_record(
        {"title": "Example bacterium chromosome, complete genome", "slen": 2_000_000}
    )
    assert not _looks_like_complete_genomic_record(
        {"title": "Example bacterium whole genome shotgun sequencing project", "slen": 2_000_000}
    )
    assert not _looks_like_complete_genomic_record(
        {"title": "MAG: Example bacterium contig", "slen": 150_000}
    )
    assert not _looks_like_complete_genomic_record(
        {"title": "Example bacterium chromosome, complete genome", "slen": 50}
    )


def test_esummary_batches_large_id_lists(monkeypatch):
    requested_urls = []

    def fake_read_url(url: str, timeout: int) -> bytes:
        requested_urls.append(url)
        ids = url.split("id=", maxsplit=1)[1].split("&", maxsplit=1)[0].split("%2C")
        result = {"result": {ncbi_id: {"uid": ncbi_id} for ncbi_id in ids}}
        import json

        return json.dumps(result).encode("utf-8")

    monkeypatch.setattr(
        "crispr_phage_predictor.ml.collect_genbank_annotated_training._read_url",
        fake_read_url,
    )
    monkeypatch.setattr(
        "crispr_phage_predictor.ml.collect_genbank_annotated_training.time.sleep",
        lambda _: None,
    )

    summaries = _esummary([str(index) for index in range(5)], email=None, batch_size=2)

    assert len(requested_urls) == 3
    assert set(summaries) == {"0", "1", "2", "3", "4"}
