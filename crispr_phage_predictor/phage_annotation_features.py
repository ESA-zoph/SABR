from __future__ import annotations

from pathlib import Path
import re

import pandas as pd
from Bio import SeqIO


PHAGE_ANNOTATION_FEATURES = [
    "phage_cds_count",
    "phage_hypothetical_cds_fraction",
    "phage_integrase_count",
    "phage_recombinase_count",
    "phage_excisionase_count",
    "phage_repressor_count",
    "phage_terminase_count",
    "phage_portal_count",
    "phage_capsid_count",
    "phage_tail_count",
    "phage_tail_fiber_count",
    "phage_baseplate_count",
    "phage_lysin_count",
    "phage_holin_count",
    "phage_endolysin_count",
    "phage_depolymerase_count",
    "phage_dna_methyltransferase_count",
    "phage_restriction_evasion_count",
    "phage_anti_crispr_keyword_count",
    "phage_temperate_marker_count",
    "phage_structural_marker_count",
    "phage_lysis_marker_count",
]


KEYWORD_GROUPS = {
    "integrase": ("integrase",),
    "recombinase": ("recombinase", "resolvase"),
    "excisionase": ("excisionase", "xis"),
    "repressor": ("repressor", "cro", "ci protein", "ci-like"),
    "terminase": ("terminase",),
    "portal": ("portal",),
    "capsid": ("capsid", "head protein", "major head"),
    "tail": ("tail", "tape measure"),
    "tail_fiber": ("tail fiber", "tail fibre", "receptor binding", "receptor-binding", "rbp"),
    "baseplate": ("baseplate",),
    "lysin": ("lysin", "lysozyme", "muramidase"),
    "holin": ("holin",),
    "endolysin": ("endolysin",),
    "depolymerase": ("depolymerase", "polysaccharide lyase", "capsular polysaccharide"),
    "dna_methyltransferase": ("methyltransferase", "dna methylase", "methylase"),
    "restriction_evasion": ("anti-restriction", "antirestriction", "ocr protein", "dar"),
    "anti_crispr_keyword": ("anti-crispr", "anti crispr", "acrif", "acrie", "acr"),
}


def add_phage_annotation_features(
    feature_table: pd.DataFrame,
    genbank_dir: str | Path = "data/curation/downloads/phages_genbank",
) -> pd.DataFrame:
    cache: dict[str, dict[str, float]] = {}
    rows = []
    for _, row in feature_table.iterrows():
        accession = str(row.get("phage_accession", "")).strip()
        features = _features_for_accession(accession, Path(genbank_dir), cache)
        merged = row.to_dict()
        merged.update(features)
        rows.append(merged)
    return pd.DataFrame(rows)


def phage_annotation_features_from_genbank(path: str | Path) -> dict[str, float]:
    records = list(SeqIO.parse(path, "genbank"))
    product_texts = []
    cds_count = 0
    hypothetical_count = 0
    for record in records:
        for feature in record.features:
            if feature.type != "CDS":
                continue
            cds_count += 1
            text = _feature_text(feature)
            product_texts.append(text)
            if "hypothetical protein" in text or text.strip() in {"hypothetical", "unknown"}:
                hypothetical_count += 1
    counts = {
        f"phage_{name}_count": _count_keyword_hits(product_texts, keywords)
        for name, keywords in KEYWORD_GROUPS.items()
    }
    counts["phage_cds_count"] = cds_count
    counts["phage_hypothetical_cds_fraction"] = (
        round(hypothetical_count / cds_count, 6) if cds_count else 0.0
    )
    counts["phage_temperate_marker_count"] = (
        counts["phage_integrase_count"]
        + counts["phage_recombinase_count"]
        + counts["phage_excisionase_count"]
        + counts["phage_repressor_count"]
    )
    counts["phage_structural_marker_count"] = (
        counts["phage_terminase_count"]
        + counts["phage_portal_count"]
        + counts["phage_capsid_count"]
        + counts["phage_tail_count"]
        + counts["phage_tail_fiber_count"]
        + counts["phage_baseplate_count"]
    )
    counts["phage_lysis_marker_count"] = (
        counts["phage_lysin_count"]
        + counts["phage_holin_count"]
        + counts["phage_endolysin_count"]
    )
    return {column: counts.get(column, 0.0) for column in PHAGE_ANNOTATION_FEATURES}


def _features_for_accession(
    accession: str,
    genbank_dir: Path,
    cache: dict[str, dict[str, float]],
) -> dict[str, float]:
    if not accession:
        return _empty_features()
    if accession not in cache:
        path = _genbank_path(genbank_dir, accession)
        cache[accession] = (
            phage_annotation_features_from_genbank(path) if path.exists() else _empty_features()
        )
    return cache[accession]


def _genbank_path(genbank_dir: Path, accession: str) -> Path:
    return genbank_dir / f"{_safe_filename(accession)}.gb"


def _feature_text(feature) -> str:
    values = []
    for key in ["product", "gene", "note", "function"]:
        values.extend(str(value) for value in feature.qualifiers.get(key, []))
    return " ".join(values).lower()


def _count_keyword_hits(texts: list[str], keywords: tuple[str, ...]) -> int:
    return sum(1 for text in texts if any(keyword in text for keyword in keywords))


def _empty_features() -> dict[str, float]:
    return {column: 0.0 for column in PHAGE_ANNOTATION_FEATURES}


def _safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")
