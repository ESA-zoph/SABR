from __future__ import annotations

import base64
import json
from time import perf_counter
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import streamlit as st

from crispr_phage_predictor.external.blast import BLASTN_COMMAND, MAKEBLASTDB_COMMAND
from crispr_phage_predictor.external.minced import (
    MINCED_COMMAND,
    active_minced_backend,
    minced_available,
)
from crispr_phage_predictor.external.tools import missing_tool_message, tool_available
from crispr_phage_predictor.matching import reverse_complement
from crispr_phage_predictor.cas_prediction import predict_array_cas_subtypes
from crispr_phage_predictor.io import (
    deduplicate_records,
    parse_uploaded_fastas,
    summarize_accession_conflicts,
    summarize_duplicate_records,
    summarize_records,
    summarize_uploaded_files,
)
from crispr_phage_predictor.ml.model_artifact import DEFAULT_MODEL_PATH, model_artifact_metadata
from crispr_phage_predictor.output import save_analysis_run
from crispr_phage_predictor.pipeline import (
    annotate_spacer_hits_with_pam,
    build_crispr_targeting_evidence_matrix,
    build_exact_match_heatmap,
    build_initial_run_summary,
    detect_arrays_for_records,
    find_spacer_hits_for_records,
    summarize_crispr_arrays,
    summarize_pam_subtype_support,
    summarize_spacer_hits,
    summarize_spacers,
)


st.set_page_config(
    page_title="SABR",
    layout="wide",
)

LOGO_PATH = Path("assets") / "aub-logo.png"
FRAGMENTED_FASTA_WARNING_RECORDS = 1_000
FRAGMENTED_FASTA_INTERNAL_STOP_RECORDS = 5_000
DEMO_INPUT_DIR = Path("data") / "examples" / "real_demo"
REPORT_FILE_NAME = "sabr_report.md"


class DemoUploadedFile:
    def __init__(self, path: Path):
        self.path = path
        self.name = path.name
        self._bytes = path.read_bytes()

    def getvalue(self) -> bytes:
        return self._bytes


def image_data_uri(path: Path) -> str:
    if not path.exists():
        base64_path = path.with_name(f"{path.name}.b64")
        if not base64_path.exists():
            return ""
        encoded = "".join(base64_path.read_text(encoding="ascii").split())
        return f"data:image/png;base64,{encoded}" if encoded else ""
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def render_brand_header() -> None:
    logo_uri = image_data_uri(LOGO_PATH)
    logo_markup = (
        f"<img class='sabr-logo' src='{logo_uri}' alt='The Phage Lab logo'>"
        if logo_uri
        else "<div class='sabr-logo-fallback'>The Phage Lab</div>"
    )
    st.markdown(
        f"""
        <style>
            .block-container {{
                padding-top: 3rem;
            }}
            .sabr-header {{
                background:
                    linear-gradient(115deg, #0f4c81 0%, #8b2633 58%, #7a1f2b 100%);
                border-radius: 8px;
                color: #ffffff;
                display: grid;
                grid-template-columns: minmax(0, 1fr) minmax(360px, 460px);
                align-items: center;
                gap: 2rem;
                min-height: 166px;
                overflow: visible;
                padding: 2.05rem 2rem 1.45rem;
                margin-bottom: 1.2rem;
                box-shadow: 0 8px 24px rgba(15, 76, 129, 0.16);
            }}
            .sabr-copy {{
                min-width: 0;
            }}
            .sabr-title {{
                font-size: 2.45rem;
                line-height: 1.18;
                font-weight: 800;
                letter-spacing: 0;
                margin: 0;
            }}
            .sabr-subtitle {{
                font-size: 1rem;
                line-height: 1.35;
                margin-top: 0.35rem;
                color: rgba(255, 255, 255, 0.92);
            }}
            .sabr-lab {{
                font-size: 0.86rem;
                margin-top: 0.45rem;
                color: rgba(255, 255, 255, 0.78);
            }}
            .sabr-logo-wrap {{
                display: flex;
                align-items: center;
                justify-content: flex-end;
                min-width: 0;
            }}
            .sabr-logo {{
                display: block;
                max-height: 118px;
                max-width: 440px;
                object-fit: contain;
                width: auto;
                height: auto;
            }}
            .sabr-logo-fallback {{
                border: 1px solid rgba(255,255,255,0.45);
                border-radius: 6px;
                padding: 0.55rem 0.75rem;
                font-weight: 700;
                white-space: nowrap;
            }}
            @media (max-width: 720px) {{
                .sabr-header {{
                    align-items: flex-start;
                    grid-template-columns: 1fr;
                    gap: 1rem;
                    min-height: 0;
                    padding: 1.25rem;
                }}
                .sabr-title {{
                    font-size: 1.9rem;
                }}
                .sabr-logo-wrap {{
                    justify-content: flex-start;
                }}
                .sabr-logo {{
                    max-height: 82px;
                    max-width: min(100%, 360px);
                }}
            }}
        </style>
        <div class="sabr-header">
            <div class="sabr-copy">
                <div class="sabr-title">SABR</div>
                <div class="sabr-subtitle">Spacer Alignment-Based Recognition</div>
                <div class="sabr-lab">CRISPR-phage targeting evidence mapper | The Phage Lab, Faculty of Medicine, AUB</div>
            </div>
            <div class="sabr-logo-wrap">{logo_markup}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


render_brand_header()

if "page" not in st.session_state:
    st.session_state.page = "analysis"
if "latest_result" not in st.session_state:
    st.session_state.latest_result = None
if "latest_run_id" not in st.session_state:
    st.session_state.latest_run_id = None
if "selected_pair" not in st.session_state:
    st.session_state.selected_pair = None

if st.query_params.get("view") == "hit_details":
    st.session_state.page = "hit_details"
    query_bacterium = st.query_params.get("bacterium")
    query_phage = st.query_params.get("phage")
    query_run = st.query_params.get("run")
    if query_bacterium and query_phage:
        st.session_state.selected_pair = {
            "bacterium": query_bacterium,
            "phage": query_phage,
        }
    if query_run:
        st.session_state.latest_run_id = query_run


def _format_seconds(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    minutes, remainder = divmod(int(seconds), 60)
    if minutes < 60:
        return f"{minutes}m {remainder}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes}m"


def _clean_saved_value(value):
    if pd.isna(value):
        return ""
    return value


def _optional_float(value):
    value = _clean_saved_value(value)
    if value == "":
        return None
    return float(value)


def _optional_int(value):
    value = _clean_saved_value(value)
    if value == "":
        return None
    return int(float(value))


def _saved_hit_from_row(row) -> SimpleNamespace:
    return SimpleNamespace(
        bacterium_id=str(_clean_saved_value(row.get("bacterium", ""))),
        phage_id=str(_clean_saved_value(row.get("phage", ""))),
        array_id=str(_clean_saved_value(row.get("array_id", ""))),
        spacer_id=str(_clean_saved_value(row.get("spacer_id", ""))),
        phage_contig_id=str(_clean_saved_value(row.get("phage_contig", ""))),
        start=_optional_int(row.get("start", "")) or 0,
        end=_optional_int(row.get("end", "")) or 0,
        strand=str(_clean_saved_value(row.get("strand", ""))),
        identity=(_optional_float(row.get("identity_percent", "")) or 0.0) / 100,
        mismatches=_optional_int(row.get("mismatches", "")) or 0,
        alignment_length=_optional_int(row.get("alignment_length", "")) or 0,
        spacer_length=_optional_int(row.get("spacer_length", "")) or 0,
        coverage=(_optional_float(row.get("coverage_percent", "")) or 0.0) / 100,
        spacer_sequence=str(_clean_saved_value(row.get("spacer_sequence", ""))),
        aligned_spacer_sequence=str(_clean_saved_value(row.get("aligned_spacer_sequence", ""))),
        aligned_protospacer_sequence=str(
            _clean_saved_value(row.get("aligned_protospacer_sequence", ""))
        ),
        protospacer_sequence=str(_clean_saved_value(row.get("protospacer_sequence", ""))),
        protospacer_5p_flank=str(_clean_saved_value(row.get("protospacer_5p_flank", ""))),
        protospacer_3p_flank=str(_clean_saved_value(row.get("protospacer_3p_flank", ""))),
        genomic_upstream_flank=str(_clean_saved_value(row.get("genomic_upstream_flank", ""))),
        genomic_downstream_flank=str(_clean_saved_value(row.get("genomic_downstream_flank", ""))),
        predicted_cas_subtype=str(_clean_saved_value(row.get("predicted_cas_subtype", ""))),
        cas_subtype_confidence=_optional_float(row.get("cas_subtype_confidence", "")),
        pam_rule=str(_clean_saved_value(row.get("pam_rule", ""))),
        pam_sequence=str(_clean_saved_value(row.get("pam_sequence", ""))),
        pam_support_level=str(_clean_saved_value(row.get("pam_support_level", ""))),
        pam_offset_from_protospacer=_optional_int(row.get("pam_offset_from_protospacer", "")),
        seed_region=str(_clean_saved_value(row.get("seed_region", ""))),
        seed_mismatches=_optional_int(row.get("seed_mismatches", "")),
        seed_mismatch_positions=str(_clean_saved_value(row.get("seed_mismatch_positions", ""))),
    )


def _load_saved_result(run_id: str | None):
    if not run_id:
        return None
    run_dir = Path("outputs") / "runs" / run_id
    evidence_path = run_dir / "evidence_matrix.csv"
    hits_path = run_dir / "spacer_hits.csv"
    if not evidence_path.exists() or not hits_path.exists():
        return None
    evidence_matrix = _read_saved_table(evidence_path)
    hit_table = _read_saved_table(hits_path)
    metadata = {}
    metadata_path = run_dir / "run_metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    heatmap_path = run_dir / "heatmap.csv"
    heatmap = None
    if heatmap_path.exists():
        heatmap_table = pd.read_csv(heatmap_path)
        if not heatmap_table.empty:
            index_column = heatmap_table.columns[0]
            heatmap = heatmap_table.set_index(index_column)
    return {
        "evidence_matrix": evidence_matrix,
        "heatmap": heatmap,
        "spacer_hits": [_saved_hit_from_row(row) for _, row in hit_table.iterrows()],
        "output_dir": str(run_dir),
        "metadata": metadata,
    }


def _read_saved_table(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _clear_detail_query_params() -> None:
    if st.query_params.get("view") == "hit_details":
        st.query_params.clear()


def _bacterial_fragmentation_rows(records) -> list[dict[str, int | str]]:
    rows = []
    for source_file in sorted({record.source_file for record in records}):
        file_records = [record for record in records if record.source_file == source_file]
        record_count = len(file_records)
        if record_count < FRAGMENTED_FASTA_WARNING_RECORDS:
            continue
        rows.append(
            {
                "file_name": source_file,
                "parsed_records": record_count,
                "total_bp": sum(record.length for record in file_records),
                "longest_record_bp": max((record.length for record in file_records), default=0),
            }
        )
    return rows


def _demo_uploaded_files() -> tuple[list[DemoUploadedFile], list[DemoUploadedFile]]:
    bacteria_dir = DEMO_INPUT_DIR / "bacteria"
    phage_dir = DEMO_INPUT_DIR / "phages"
    if not bacteria_dir.exists() or not phage_dir.exists():
        return [], []
    bacteria_paths = sorted(bacteria_dir.glob("*.fasta"))
    phage_paths = sorted(phage_dir.glob("*.fasta"))
    if not bacteria_paths or not phage_paths:
        return [], []
    return (
        [DemoUploadedFile(path) for path in bacteria_paths],
        [DemoUploadedFile(path) for path in phage_paths],
    )


def _quality_warning_rows(records, role: str) -> list[dict[str, object]]:
    rows = []
    for source_file in sorted({record.source_file for record in records}):
        file_records = [record for record in records if record.source_file == source_file]
        total_bp = sum(record.length for record in file_records)
        if not file_records:
            rows.append(
                {
                    "file_name": source_file,
                    "input_type": role,
                    "issue": "no FASTA records parsed",
                    "detail": "Check that the file is FASTA, not GenBank or FASTQ.",
                }
            )
            continue
        ambiguous_bp = sum(
            sum(1 for base in record.sequence.upper() if base not in {"A", "C", "G", "T"})
            for record in file_records
        )
        ambiguous_percent = (ambiguous_bp / total_bp * 100) if total_bp else 0.0
        short_records = sum(record.length < 500 for record in file_records)
        if len(file_records) >= FRAGMENTED_FASTA_WARNING_RECORDS:
            rows.append(
                {
                    "file_name": source_file,
                    "input_type": role,
                    "issue": "highly fragmented upload",
                    "detail": f"{len(file_records):,} records; prefer complete/scaffolded assemblies.",
                }
            )
        if ambiguous_percent >= 5:
            rows.append(
                {
                    "file_name": source_file,
                    "input_type": role,
                    "issue": "high ambiguous-base content",
                    "detail": f"{ambiguous_percent:.1f}% non-ACGT bases.",
                }
            )
        if short_records and short_records / len(file_records) >= 0.5:
            rows.append(
                {
                    "file_name": source_file,
                    "input_type": role,
                    "issue": "many short records",
                    "detail": f"{short_records}/{len(file_records)} records are under 500 bp.",
                }
            )
    return rows


def _label_score(score: object, level: object) -> str:
    try:
        numeric = float(score)
    except (TypeError, ValueError):
        numeric = 0.0
    text = str(level or "")
    if text:
        return f"{text} ({numeric:.2f})"
    if numeric <= 0:
        return "no spacer-match evidence (0.00)"
    if numeric >= 75:
        return f"strong candidate CRISPR targeting evidence ({numeric:.2f})"
    if numeric >= 50:
        return f"moderate candidate CRISPR targeting evidence ({numeric:.2f})"
    return f"weak candidate CRISPR targeting evidence ({numeric:.2f})"


def _pair_score_note(pair_row) -> str:
    unique_spacers = int(pair_row.get("unique_matching_spacers", 0))
    score = float(pair_row.get("crispr_targeting_score", 0.0))
    pam_supported = int(pair_row.get("pam_supported_hits", 0))
    pam_evaluated = int(pair_row.get("pam_evaluated_hits", 0))
    identity = float(pair_row.get("best_identity_percent", 0.0))
    coverage = float(pair_row.get("best_coverage_percent", 0.0))
    if score <= 0:
        return "No spacer-protospacer match was detected for this bacteria-phage pair."
    if pam_evaluated == 0:
        return (
            f"This score is driven by {unique_spacers} matching spacer(s) with best "
            f"{identity:.2f}% identity and {coverage:.2f}% coverage. PAM/PFS and seed "
            "layers did not contribute because no PAM/PFS rule was evaluated for this run or hit."
        )
    if pam_supported > 0:
        return (
            f"This score includes spacer-match evidence plus PAM/PFS support in "
            f"{pam_supported}/{pam_evaluated} evaluated hit(s)."
        )
    return (
        f"This score includes spacer-match evidence, but none of the {pam_evaluated} "
        "PAM/PFS-evaluated hit(s) supported the selected rule."
    )


def _pam_evaluation_note(hit) -> str:
    pam_rule = str(getattr(hit, "pam_rule", "") or "").strip()
    pam_support = str(getattr(hit, "pam_support_level", "") or "").strip()
    predicted_subtype = str(getattr(hit, "predicted_cas_subtype", "") or "").strip()
    confidence = getattr(hit, "cas_subtype_confidence", None)
    if pam_rule:
        return ""
    if predicted_subtype:
        confidence_text = (
            f" with confidence {float(confidence):.3f}"
            if confidence is not None
            else ""
        )
        return (
            f"PAM/PFS was not evaluated because no curated rule was selected for the "
            f"predicted subtype {predicted_subtype}{confidence_text}."
        )
    if pam_support == "not_evaluated" or not pam_support:
        return (
            "PAM/PFS was not evaluated because this hit has no selected PAM/PFS rule. "
            "Use 'Auto from predicted subtype' or a manual expert override before running analysis."
        )
    return ""


def _generate_markdown_report(result: dict) -> str:
    evidence_matrix = result.get("evidence_matrix", pd.DataFrame())
    output_dir = result.get("output_dir", "")
    metadata = result.get("metadata") or {}
    lines = [
        "# SABR Analysis Report",
        "",
        f"- Run ID: `{Path(output_dir).name if output_dir else 'unsaved'}`",
        f"- Created at: `{metadata.get('created_at', '')}`",
        f"- Detection method: `{metadata.get('detection_method', '')}`",
        f"- Matching method: `{metadata.get('matching_method', '')}`",
        f"- PAM/PFS mode: `{metadata.get('pam_mode', '')}`",
        "",
        "SABR reports CRISPR spacer-targeting evidence. It does not prove phenotype-level phage resistance.",
        "",
        "## Summary",
        "",
        f"- Bacterial sequences: `{metadata.get('bacterial_sequence_count', '')}`",
        f"- Phage sequences: `{metadata.get('phage_sequence_count', '')}`",
        f"- Candidate arrays: `{metadata.get('candidate_array_count', '')}`",
        f"- Extracted spacers: `{metadata.get('extracted_spacer_count', '')}`",
        f"- Spacer hits: `{metadata.get('spacer_hit_count', '')}`",
        "",
        "## Highest Evidence Pairs",
        "",
    ]
    if evidence_matrix.empty:
        lines.append("No bacteria-phage evidence rows were available.")
    else:
        top_rows = evidence_matrix.sort_values(
            by=["crispr_targeting_score", "unique_matching_spacers", "spacer_hits"],
            ascending=False,
        ).head(10)
        lines.extend(
            [
                "| Bacterium | Phage | Evidence | Unique spacers | PAM/PFS |",
                "| --- | --- | --- | ---: | --- |",
            ]
        )
        for _, row in top_rows.iterrows():
            lines.append(
                "| "
                f"{row.get('bacterium', '')} | {row.get('phage', '')} | "
                f"{_label_score(row.get('crispr_targeting_score', 0), row.get('current_evidence_level', ''))} | "
                f"{row.get('unique_matching_spacers', '')} | {row.get('pam_support_level', '')} |"
            )
    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- Spacer/PAM/PFS evidence is computational and requires biological context.",
            "- The current PAM/PFS rule set is incomplete.",
            "- Cas subtype prediction is repeat/array-derived and does not confirm active Cas function.",
            "- Phage escape, anti-CRISPR genes, expression, and experimental conditions are not fully evaluated.",
        ]
    )
    return "\n".join(lines) + "\n"


def _write_report_if_possible(result: dict) -> None:
    output_dir = result.get("output_dir")
    if not output_dir:
        return
    Path(output_dir, REPORT_FILE_NAME).write_text(
        _generate_markdown_report(result),
        encoding="utf-8",
    )


def render_exact_match_heatmap(heatmap) -> None:
    max_value = int(heatmap.to_numpy().max()) if not heatmap.empty else 0
    scale_max = max(max_value, 1)

    html = [
        "<div style='overflow-x:auto;'>",
        "<table style='border-collapse:collapse;width:100%;font-size:0.9rem;'>",
        "<thead><tr>",
        "<th style='text-align:left;padding:0.5rem;border-bottom:1px solid #ddd;'>Bacterium</th>",
    ]
    for phage in heatmap.columns:
        html.append(
            "<th style='text-align:center;padding:0.5rem;border-bottom:1px solid #ddd;'>"
            f"{phage}</th>"
        )
    html.append("</tr></thead><tbody>")

    for bacterium, row in heatmap.iterrows():
        html.append("<tr>")
        html.append(
            "<td style='font-weight:600;padding:0.5rem;border-bottom:1px solid #eee;'>"
            f"{bacterium}</td>"
        )
        for value in row:
            intensity = int(255 - (int(value) / scale_max) * 165)
            background = f"rgb(255,{intensity},{max(80, intensity - 80)})" if value else "#f7f7f7"
            color = "#111" if value else "#666"
            html.append(
                "<td style='text-align:center;padding:0.5rem;border-bottom:1px solid #eee;"
                f"background:{background};color:{color};font-weight:600;'>"
                f"{int(value)}</td>"
            )
        html.append("</tr>")

    html.append("</tbody></table></div>")
    st.markdown("".join(html), unsafe_allow_html=True)


def _sequence_alignment_lines(spacer: str, protospacer: str) -> tuple[str, str, str]:
    left = (spacer or "").upper()
    right = (protospacer or "").upper()
    comparable = min(len(left), len(right))
    marker = "".join("|" if left[index] == right[index] else "." for index in range(comparable))
    if len(left) > comparable:
        marker += " " * (len(left) - comparable)
    elif len(right) > comparable:
        marker += " " * (len(right) - comparable)
    return left, marker, right


def _format_hit_context(hit) -> str:
    upstream = hit.genomic_upstream_flank or ""
    protospacer = hit.protospacer_sequence or ""
    downstream = hit.genomic_downstream_flank or ""
    return f"{upstream}[{protospacer}]{downstream}"


def _display_alignment_sequences(hit) -> tuple[str, str]:
    displayed_spacer = hit.aligned_spacer_sequence or hit.spacer_sequence
    displayed_protospacer = hit.aligned_protospacer_sequence or hit.protospacer_sequence
    if not hit.aligned_protospacer_sequence and hit.strand == "-":
        displayed_protospacer = reverse_complement(displayed_protospacer)
    return displayed_spacer, displayed_protospacer


def _positive_pair_options(evidence_matrix) -> list[tuple[str, str, str]]:
    if evidence_matrix.empty:
        return []

    positive_rows = evidence_matrix[
        evidence_matrix["unique_matching_spacers"].fillna(0).astype(int) > 0
    ].copy()
    if positive_rows.empty:
        return []

    positive_rows = positive_rows.sort_values(
        by=["crispr_targeting_score", "unique_matching_spacers", "spacer_hits"],
        ascending=False,
    )
    options = []
    for _, row in positive_rows.iterrows():
        label = (
            f"{row['bacterium']} -> {row['phage']} | "
            f"{int(row['unique_matching_spacers'])} spacers | "
            f"score {float(row['crispr_targeting_score']):.2f}"
        )
        options.append((label, row["bacterium"], row["phage"]))
    return options


def render_open_hit_details_control(evidence_matrix) -> None:
    options = _positive_pair_options(evidence_matrix)
    if not options:
        st.info("No nonzero heatmap cells are available for detailed hit inspection.")
        return

    selected_label = st.selectbox(
        "Select heatmap cell",
        [option[0] for option in options],
        help="Choose a bacteria-phage cell, then open a dedicated hit-details page.",
    )
    _, selected_bacterium, selected_phage = next(
        option for option in options if option[0] == selected_label
    )
    st.session_state.selected_pair = {
        "bacterium": selected_bacterium,
        "phage": selected_phage,
    }
    st.caption(f"Selected for details: {selected_bacterium} -> {selected_phage}")
    if st.button("Open hit details", type="primary"):
        st.session_state.page = "hit_details"
        st.rerun()


def render_latest_analysis_result(result, title: str = "Latest Analysis Results") -> None:
    evidence_matrix = result.get("evidence_matrix")
    heatmap = result.get("heatmap")
    if heatmap is None and evidence_matrix is not None:
        heatmap = build_exact_match_heatmap(evidence_matrix)

    st.subheader(title)
    output_dir = result.get("output_dir")
    if output_dir:
        st.caption(f"Analysis outputs saved to {output_dir}")
    st.caption(
        "Cell values are unique bacterial spacers with protospacer matches in each phage. "
        "This is evidence of candidate CRISPR targeting, not a final resistance call."
    )
    if evidence_matrix is not None and not evidence_matrix.empty:
        top_rows = evidence_matrix.sort_values(
            by=["crispr_targeting_score", "unique_matching_spacers", "spacer_hits"],
            ascending=False,
        ).head(8)
        top_view = top_rows[
            [
                "bacterium",
                "phage",
                "crispr_targeting_score",
                "current_evidence_level",
                "unique_matching_spacers",
                "pam_support_level",
                "evidence_summary",
            ]
        ].copy()
        top_view["evidence"] = top_view.apply(
            lambda row: _label_score(
                row["crispr_targeting_score"],
                row["current_evidence_level"],
            ),
            axis=1,
        )
        top_view = top_view[
            [
                "bacterium",
                "phage",
                "evidence",
                "unique_matching_spacers",
                "pam_support_level",
                "evidence_summary",
            ]
        ]
        with st.expander("Highest-evidence pairs", expanded=True):
            st.dataframe(top_view, use_container_width=True)
    report_text = _generate_markdown_report(result)
    st.download_button(
        "Download Markdown report",
        data=report_text,
        file_name=f"{Path(output_dir).name if output_dir else 'sabr'}_report.md",
        mime="text/markdown",
    )
    _write_report_if_possible(result)

    if heatmap is None or heatmap.empty:
        st.info("No nonzero heatmap is available for the latest analysis.")
        return
    render_exact_match_heatmap(heatmap)
    render_open_hit_details_control(evidence_matrix)


def render_pair_hit_details_page(evidence_matrix, spacer_hits) -> None:
    selected_pair = st.session_state.get("selected_pair") or {}
    selected_bacterium = selected_pair.get("bacterium")
    selected_phage = selected_pair.get("phage")
    if not selected_bacterium or not selected_phage:
        st.warning("No heatmap cell is selected.")
        if st.button("Back to results"):
            st.session_state.page = "analysis"
            _clear_detail_query_params()
            st.rerun()
        return

    if st.button("Back to results"):
        st.session_state.page = "analysis"
        _clear_detail_query_params()
        st.rerun()

    st.subheader("Spacer-Hit Details")
    st.caption(
        f"{selected_bacterium} -> {selected_phage}. "
        "Coordinates are one-based positions on the uploaded phage FASTA record."
    )

    positive_rows = evidence_matrix[
        evidence_matrix["unique_matching_spacers"].fillna(0).astype(int) > 0
    ].copy()
    pair_hits = [
        hit
        for hit in spacer_hits
        if hit.bacterium_id == selected_bacterium and hit.phage_id == selected_phage
    ]
    pair_row = positive_rows[
        (positive_rows["bacterium"] == selected_bacterium)
        & (positive_rows["phage"] == selected_phage)
    ]
    if pair_row.empty:
        st.warning("The selected heatmap cell is not available in the current run.")
        return
    pair_row = pair_row.iloc[0]

    col_a, col_b, col_c, col_d = st.columns(4)
    col_a.metric("Targeting score", f"{float(pair_row['crispr_targeting_score']):.2f}")
    col_b.metric("Unique spacers", int(pair_row["unique_matching_spacers"]))
    col_c.metric("Spacer hits", int(pair_row["spacer_hits"]))
    col_d.metric("PAM/PFS", str(pair_row["pam_support_level"]))

    st.markdown(f"**Evidence label:** {_label_score(pair_row['crispr_targeting_score'], pair_row.get('current_evidence_level', ''))}")
    st.caption(_pair_score_note(pair_row))
    st.info(str(pair_row.get("interpretation", "")))

    component_rows = [
        {
            "component": "best identity",
            "value": f"{float(pair_row.get('best_identity_percent', 0.0)):.2f}%",
            "meaning": "Best spacer/protospacer identity among hits.",
        },
        {
            "component": "best coverage",
            "value": f"{float(pair_row.get('best_coverage_percent', 0.0)):.2f}%",
            "meaning": "Best aligned spacer coverage among hits.",
        },
        {
            "component": "PAM/PFS-supported hits",
            "value": (
                f"{int(pair_row.get('pam_supported_hits', 0))}/"
                f"{int(pair_row.get('pam_evaluated_hits', 0))}"
            ),
            "meaning": "Hits compatible with the predicted or supplied PAM/PFS rule.",
        },
        {
            "component": "best PAM/PFS compatibility",
            "value": pair_row.get("best_pam_compatibility_score", ""),
            "meaning": "Best motif compatibility score when PAM/PFS was evaluated.",
        },
        {
            "component": "best seed mismatches",
            "value": pair_row.get("best_seed_mismatches", ""),
            "meaning": "Lowest mismatch count in the PAM/PFS-proximal seed region.",
        },
    ]
    with st.expander("Score components", expanded=True):
        st.dataframe(component_rows, use_container_width=True)

    hit_rows = []
    for index, hit in enumerate(pair_hits, start=1):
        hit_rows.append(
            {
                "hit": index,
                "array_id": hit.array_id,
                "spacer_id": hit.spacer_id.rsplit("|", 1)[-1],
                "phage_contig": hit.phage_contig_id,
                "start": hit.start,
                "end": hit.end,
                "strand": hit.strand,
                "identity_percent": round(hit.identity * 100, 2),
                "coverage_percent": round(hit.coverage * 100, 2),
                "predicted_subtype": hit.predicted_cas_subtype,
                "subtype_confidence": hit.cas_subtype_confidence,
                "pam_rule": hit.pam_rule,
                "pam_sequence": hit.pam_sequence,
                "pam_support": hit.pam_support_level,
                "pam_offset_from_protospacer": hit.pam_offset_from_protospacer,
                "seed_mismatches": hit.seed_mismatches,
            }
        )
    st.dataframe(hit_rows, use_container_width=True)

    for index, hit in enumerate(pair_hits, start=1):
        title = (
            f"Hit {index}: {hit.phage_contig_id}:{hit.start}-{hit.end} "
            f"({hit.strand}) | {round(hit.identity * 100, 2)}% identity"
        )
        with st.expander(title, expanded=index == 1):
            displayed_spacer, displayed_protospacer = _display_alignment_sequences(hit)
            left, marker, right = _sequence_alignment_lines(
                displayed_spacer,
                displayed_protospacer,
            )
            st.code(
                "\n".join(
                    [
                        f"Spacer      {left}",
                        f"            {marker}",
                        f"Protospacer {right}",
                    ]
                ),
                language="text",
            )
            if hit.aligned_spacer_sequence:
                st.caption(
                    "Alignment is the BLAST aligned segment. The full original spacer is "
                    f"{len(hit.spacer_sequence)} bp; this alignment spans "
                    f"{hit.alignment_length} bp."
                )
            elif hit.strand == "-":
                st.caption(
                    "This hit is on the reverse strand. The alignment shows the "
                    "strand-oriented protospacer so it can be compared directly with "
                    "the spacer; the context block below remains in phage-reference orientation."
                )
            st.markdown("**Phage-reference context**")
            st.code(_format_hit_context(hit), language="text")
            pam_note = _pam_evaluation_note(hit)
            if pam_note:
                st.warning(pam_note)
            flank_rows = [
                {"field": "protospacer_5p_flank", "value": hit.protospacer_5p_flank},
                {"field": "protospacer_3p_flank", "value": hit.protospacer_3p_flank},
                {"field": "genomic_upstream_flank", "value": hit.genomic_upstream_flank},
                {"field": "genomic_downstream_flank", "value": hit.genomic_downstream_flank},
                {"field": "PAM/PFS rule", "value": hit.pam_rule},
                {"field": "PAM/PFS sequence", "value": hit.pam_sequence},
                {"field": "PAM/PFS support", "value": hit.pam_support_level},
                {
                    "field": "PAM/PFS offset from protospacer",
                    "value": (
                        hit.pam_offset_from_protospacer
                        if hit.pam_offset_from_protospacer is not None
                        else ""
                    ),
                },
                {"field": "seed region", "value": hit.seed_region},
                {"field": "seed mismatch positions", "value": hit.seed_mismatch_positions},
            ]
            st.dataframe(flank_rows, use_container_width=True)

with st.sidebar:
    st.header("Inputs")
    use_demo_inputs = st.checkbox(
        "Use built-in real demo",
        value=False,
        help="Loads a small real PA14 CRISPR-locus panel with JBD18 positive and Lambda negative-control phages.",
    )
    bacterial_files = st.file_uploader(
        "Bacterial FASTA files",
        type=["fasta", "fa", "fna", "ffn", "txt", "gz"],
        accept_multiple_files=True,
    )
    phage_files = st.file_uploader(
        "Phage FASTA files",
        type=["fasta", "fa", "fna", "ffn", "txt", "gz"],
        accept_multiple_files=True,
    )
    st.header("Methods")
    detection_label = st.selectbox(
        "CRISPR detection",
        ["Auto recommended", "Internal exact-repeat MVP", "MinCED-compatible"],
        index=1,
    )
    matching_label = st.selectbox(
        "Spacer-phage matching",
        ["Auto recommended", "Internal exact match", "BLASTN"],
        index=1,
    )

    st.header("External tools")
    minced_backend_available = minced_available()
    minced_backend_name = active_minced_backend()
    blast_available = tool_available(BLASTN_COMMAND) and tool_available(MAKEBLASTDB_COMMAND)
    st.caption(
        "MinCED-compatible detection: "
        f"{minced_backend_name if minced_backend_available else 'not found'}"
    )
    st.caption(f"BLAST+: {'available' if blast_available else 'not found'}")

    if detection_label == "Auto recommended":
        detection_method = "minced" if minced_backend_available else "internal"
    else:
        detection_method = "minced" if detection_label == "MinCED-compatible" else "internal"

    if matching_label == "Auto recommended":
        matching_method = "blast" if blast_available else "internal"
    else:
        matching_method = "blast" if matching_label == "BLASTN" else "internal"

    st.caption(f"Selected detection backend: {detection_method}")
    st.caption(f"Selected matching backend: {matching_method}")
    st.header("Subtype model")
    if DEFAULT_MODEL_PATH.exists():
        model_size_mb = DEFAULT_MODEL_PATH.stat().st_size / (1024 * 1024)
        st.caption(f"ExtraTrees artifact: available ({model_size_mb:.1f} MB)")
    else:
        st.warning(
            "ExtraTrees subtype model artifact is missing. SABR will use a nearest-repeat "
            "fallback only if a local training table is available."
        )
    blast_min_identity_percent = 90
    blast_min_coverage_percent = 95
    blast_require_full_query = False
    if matching_method == "blast":
        blast_min_identity_percent = st.slider(
            "BLASTN minimum identity",
            min_value=80,
            max_value=100,
            value=90,
            step=1,
        )
        blast_min_coverage_percent = st.slider(
            "BLASTN minimum spacer coverage",
            min_value=50,
            max_value=100,
            value=95,
            step=1,
        )
        blast_require_full_query = st.checkbox(
            "Require full spacer alignment",
            value=False,
        )
    pam_mode = st.selectbox(
        "PAM/PFS evaluation",
        ["Auto from predicted subtype", "Do not evaluate", "Manual expert override"],
    )
    default_pam_rule = ""
    if pam_mode == "Manual expert override":
        default_pam_rule = st.text_input(
            "Manual PAM/PFS rule",
            value="",
            placeholder="Example: 3prime:NGG or 5prime:AWG",
            help="Applies one candidate rule to all hits in this run.",
        ).strip()
    seed_length = st.number_input(
        "Seed-region length",
        min_value=4,
        max_value=20,
        value=8,
        step=1,
        help="Used only when PAM/PFS side is known. Reports mismatches nearest the PAM/PFS.",
    )
    run_analysis = st.button("Run SABR analysis", type="primary")
    run_progress_slot = st.empty()
    run_status_slot = st.empty()
    matching_progress_slot = st.empty()

if st.session_state.page == "hit_details":
    latest_result = st.session_state.latest_result
    if latest_result is None:
        latest_result = _load_saved_result(st.session_state.latest_run_id)
        if latest_result is not None:
            st.session_state.latest_result = latest_result
            st.session_state.latest_run_id = Path(latest_result["output_dir"]).name
    if latest_result is None:
        st.warning("No analysis result is available yet.")
        if st.button("Back to analysis"):
            st.session_state.page = "analysis"
            st.rerun()
        st.stop()
    render_pair_hit_details_page(
        latest_result["evidence_matrix"],
        latest_result["spacer_hits"],
    )
    st.stop()

if detection_method == "internal" and matching_method == "internal":
    st.warning(
        "This run uses internal MVP methods only. Results are preliminary and should be "
        "benchmarked before scientific interpretation."
    )
else:
    st.warning(
        "External-tool results are still evidence of candidate CRISPR targeting, not final "
        "biological resistance. PAM, Cas function, phage escape, and anti-CRISPR genes are not fully evaluated."
    )

st.info(
    "Public demo privacy note: do not upload sensitive unpublished genomes to a hosted SABR demo. "
    "Use local Docker for private lab analyses."
)

if use_demo_inputs:
    demo_bacterial_files, demo_phage_files = _demo_uploaded_files()
    if not demo_bacterial_files or not demo_phage_files:
        st.error("The built-in demo FASTA files were not found under data/examples/real_demo/.")
        st.stop()
    bacterial_files = demo_bacterial_files
    phage_files = demo_phage_files
    st.caption("Built-in real demo FASTA files are loaded: PA14 CRISPR region, JBD18 positive control, and Lambda negative control.")
    st.info(
        "Demo expected result: PA14_CRISPR_region.fasta vs JBD18_positive_control.fasta "
        "should show spacer-targeting evidence; with the subtype model available, it "
        "should include PAM-supported strong evidence. PA14_CRISPR_region.fasta vs "
        "Lambda_negative_control.fasta should show no spacer-match evidence."
    )

if bacterial_files or phage_files:
    bacteria_records = parse_uploaded_fastas(bacterial_files)
    phage_records = parse_uploaded_fastas(phage_files)
    bacterial_duplicate_table = summarize_duplicate_records(bacteria_records)
    phage_duplicate_table = summarize_duplicate_records(phage_records)
    bacterial_accession_conflicts = summarize_accession_conflicts(bacteria_records)
    phage_accession_conflicts = summarize_accession_conflicts(phage_records)
    unique_bacteria_records = deduplicate_records(bacteria_records)
    unique_phage_records = deduplicate_records(phage_records)
    fragmented_bacterial_uploads = _bacterial_fragmentation_rows(unique_bacteria_records)
    quality_warning_table = pd.DataFrame(
        _quality_warning_rows(bacteria_records, "bacterial")
        + _quality_warning_rows(phage_records, "phage")
    )
    run_summary = build_initial_run_summary(bacteria_records, phage_records)

    st.subheader("Batch Summary")
    col_a, col_b, col_c, col_d = st.columns(4)
    col_a.metric("Bacterial files", run_summary.bacterial_file_count)
    col_b.metric("Bacterial sequences", run_summary.bacterial_sequence_count)
    col_c.metric("Phage files", run_summary.phage_file_count)
    col_d.metric("Phage sequences", run_summary.phage_sequence_count)

    st.subheader("Upload Diagnostics")
    diag_left, diag_right = st.columns(2)
    with diag_left:
        st.markdown("**Bacterial uploads**")
        bacterial_uploads = summarize_uploaded_files(bacterial_files, bacteria_records)
        if bacterial_uploads.empty:
            st.info("No bacterial files uploaded.")
        else:
            st.dataframe(bacterial_uploads, use_container_width=True)
    with diag_right:
        st.markdown("**Phage uploads**")
        phage_uploads = summarize_uploaded_files(phage_files, phage_records)
        if phage_uploads.empty:
            st.info("No phage files uploaded.")
        else:
            st.dataframe(phage_uploads, use_container_width=True)

    st.subheader("Bacterial FASTA Records")
    bacteria_table = summarize_records(bacteria_records)
    if bacteria_table.empty:
        st.info("No bacterial FASTA records loaded.")
    else:
        with st.expander("View bacterial FASTA records", expanded=False):
            st.dataframe(bacteria_table, use_container_width=True)
    if not bacterial_duplicate_table.empty:
        with st.expander("Duplicate bacterial records excluded from analysis", expanded=False):
            st.dataframe(bacterial_duplicate_table, use_container_width=True)
    if not bacterial_accession_conflicts.empty:
        with st.expander("Bacterial accession conflicts", expanded=False):
            st.warning(
                "One or more accessions appear with different sequence hashes. "
                "These records are not automatically merged."
            )
            st.dataframe(bacterial_accession_conflicts, use_container_width=True)

    st.subheader("Phage FASTA Records")
    phage_table = summarize_records(phage_records)
    if phage_table.empty:
        st.info("No phage FASTA records loaded.")
    else:
        with st.expander("View phage FASTA records", expanded=False):
            st.dataframe(phage_table, use_container_width=True)
    if not phage_duplicate_table.empty:
        with st.expander("Duplicate phage records excluded from analysis", expanded=False):
            st.dataframe(phage_duplicate_table, use_container_width=True)
    if not phage_accession_conflicts.empty:
        with st.expander("Phage accession conflicts", expanded=False):
            st.warning(
                "One or more accessions appear with different sequence hashes. "
                "These records are not automatically merged."
            )
            st.dataframe(phage_accession_conflicts, use_container_width=True)

    if fragmented_bacterial_uploads:
        st.warning(
            "One or more bacterial uploads contain many FASTA records. Highly fragmented "
            "assemblies can make internal CRISPR detection very slow. For public-beta runs, "
            "use a better assembled genome or filter to longer contigs before analysis."
        )
        with st.expander("Fragmented bacterial uploads", expanded=True):
            st.dataframe(fragmented_bacterial_uploads, use_container_width=True)

    if not quality_warning_table.empty:
        with st.expander("Input quality warnings", expanded=True):
            st.warning(
                "These warnings do not block analysis, but they can affect CRISPR detection, "
                "spacer matching, or interpretation."
            )
            st.dataframe(quality_warning_table, use_container_width=True)

    if not run_analysis:
        latest_result = st.session_state.latest_result
        if latest_result is not None:
            render_latest_analysis_result(latest_result)
        else:
            st.info(
                "Files are parsed. Press 'Run SABR analysis' in the sidebar to detect "
                "candidate arrays and match spacers against phages."
            )
        st.stop()

    if not unique_bacteria_records:
        st.error(
            "No bacterial FASTA records were parsed. Check that the bacterial file starts "
            "with FASTA headers like '>contig_1' and is not GenBank, FASTQ, or another format."
        )
        st.stop()

    if not unique_phage_records:
        st.error(
            "No phage FASTA records were parsed. Check that the phage file starts with FASTA headers."
        )
        st.stop()

    if detection_method == "minced" and not minced_backend_available:
        st.error(
            missing_tool_message(MINCED_COMMAND)
            + " The Python package 'diced' can also provide MinCED-compatible detection."
        )
        st.stop()

    if matching_method == "blast" and not blast_available:
        st.error(
            "BLAST+ was selected, but 'blastn' and/or 'makeblastdb' were not found on PATH. "
            "Install BLAST+ or select the internal exact matcher."
        )
        st.stop()

    internal_blocked_uploads = [
        row
        for row in fragmented_bacterial_uploads
        if row["parsed_records"] >= FRAGMENTED_FASTA_INTERNAL_STOP_RECORDS
    ]
    if detection_method == "internal" and internal_blocked_uploads:
        st.error(
            "Internal CRISPR detection was not started because at least one bacterial file "
            f"has {FRAGMENTED_FASTA_INTERNAL_STOP_RECORDS:,} or more FASTA records. "
            "This usually indicates a highly fragmented assembly and can take hours. "
            "Use a complete/scaffolded genome, filter short contigs, or add an external "
            "MinCED-compatible detector to the container before retrying."
        )
        st.dataframe(internal_blocked_uploads, use_container_width=True)
        st.stop()

    detection_progress = run_progress_slot.progress(
        0,
        text="Preparing CRISPR array detection...",
    )
    detection_status = run_status_slot
    analysis_start = perf_counter()

    def update_detection_progress(completed: int, total: int, record) -> None:
        elapsed = perf_counter() - analysis_start
        average = elapsed / completed if completed else 0
        remaining = max(total - completed, 0) * average
        detection_progress.progress(
            completed / total if total else 1.0,
            text=(
                f"CRISPR detection: {completed}/{total} scan steps "
                f"({record.source_file})"
            ),
        )
        detection_status.caption(
            f"Elapsed: {_format_seconds(elapsed)} | Estimated remaining: "
            f"{_format_seconds(remaining)}"
        )

    try:
        crispr_arrays = detect_arrays_for_records(
            unique_bacteria_records,
            method=detection_method,
            progress_callback=update_detection_progress,
        )
    except Exception as exc:
        st.error(f"Analysis failed during CRISPR detection: {exc}")
        st.stop()

    detection_elapsed = perf_counter() - analysis_start
    detection_progress.progress(
        1.0,
        text=f"CRISPR detection complete in {_format_seconds(detection_elapsed)}",
    )

    cas_predictions_by_array = {}
    if pam_mode == "Auto from predicted subtype":
        cas_predictions_by_array = predict_array_cas_subtypes(crispr_arrays)

    matching_progress = matching_progress_slot.progress(
        0,
        text="Preparing spacer-phage matching...",
    )
    matching_start = perf_counter()
    matching_unit_count = sum(array.spacer_count for array in crispr_arrays) * len(unique_phage_records)
    matching_backend_label = "BLASTN" if matching_method == "blast" else "internal exact matching"
    matching_progress.progress(
        0.2,
        text=(
            f"Running {matching_backend_label} for {sum(array.spacer_count for array in crispr_arrays)} "
            f"spacers across {len(unique_phage_records)} phage sequences"
        ),
    )
    try:
        spacer_hits = find_spacer_hits_for_records(
            crispr_arrays,
            unique_phage_records,
            method=matching_method,
            blast_min_identity=blast_min_identity_percent / 100,
            blast_min_coverage=blast_min_coverage_percent / 100,
            blast_require_full_query=blast_require_full_query,
        )
        spacer_hits = annotate_spacer_hits_with_pam(
            spacer_hits,
            cas_predictions_by_array=cas_predictions_by_array,
            default_pam_rule=default_pam_rule or None,
            seed_length=int(seed_length),
        )
    except Exception as exc:
        st.error(f"Analysis failed during spacer-phage matching: {exc}")
        st.stop()
    matching_elapsed = perf_counter() - matching_start
    matching_progress.progress(
        1.0,
        text=(
            f"Spacer-phage matching complete in {_format_seconds(matching_elapsed)} "
            f"({matching_unit_count} spacer-phage comparisons)"
        ),
    )

    evidence_matrix = build_crispr_targeting_evidence_matrix(
        bacteria_records=unique_bacteria_records,
        phage_records=unique_phage_records,
        hits=spacer_hits,
    )
    heatmap = build_exact_match_heatmap(evidence_matrix)
    output_dir = save_analysis_run(
        bacteria_records=unique_bacteria_records,
        phage_records=unique_phage_records,
        crispr_arrays=crispr_arrays,
        spacer_hits=spacer_hits,
        evidence_matrix=evidence_matrix,
        heatmap=heatmap,
        detection_method=detection_method,
        matching_method=matching_method,
        detection_backend_detail=active_minced_backend() if detection_method == "minced" else None,
        blast_min_identity=blast_min_identity_percent / 100 if matching_method == "blast" else None,
        blast_min_coverage=blast_min_coverage_percent / 100 if matching_method == "blast" else None,
        blast_require_full_query=blast_require_full_query if matching_method == "blast" else None,
        pam_rule=default_pam_rule or None,
        pam_mode=pam_mode,
        cas_prediction_count=len(cas_predictions_by_array),
        cas_model_artifact=model_artifact_metadata(),
        seed_length=int(seed_length),
        bacterial_duplicate_record_count=len(bacteria_records) - len(unique_bacteria_records),
        phage_duplicate_record_count=len(phage_records) - len(unique_phage_records),
        detection_elapsed_seconds=detection_elapsed,
        matching_elapsed_seconds=matching_elapsed,
        total_elapsed_seconds=perf_counter() - analysis_start,
    )
    st.session_state.latest_result = {
        "evidence_matrix": evidence_matrix,
        "heatmap": heatmap,
        "spacer_hits": spacer_hits,
        "output_dir": str(output_dir),
        "metadata": json.loads((Path(output_dir) / "run_metadata.json").read_text(encoding="utf-8")),
    }
    st.session_state.latest_run_id = Path(output_dir).name
    st.session_state.page = "analysis"

    st.success(f"Analysis outputs saved to {output_dir}")

    st.subheader("Bacteria-Phage Spacer-Hit Heatmap")
    st.caption(
        "Cell values are unique bacterial spacers with protospacer matches in each phage. "
        "This is evidence of candidate CRISPR targeting, not a final resistance call."
    )
    if heatmap.empty:
        st.info("Upload at least one bacterial FASTA and one phage FASTA to build the heatmap.")
    else:
        render_exact_match_heatmap(heatmap)
        render_open_hit_details_control(evidence_matrix)

    fresh_result = st.session_state.latest_result
    report_text = _generate_markdown_report(fresh_result)
    _write_report_if_possible(fresh_result)
    st.download_button(
        "Download Markdown report",
        data=report_text,
        file_name=f"{Path(output_dir).name}_report.md",
        mime="text/markdown",
    )

    st.subheader("Summary")
    col_a, col_b, col_c, col_d = st.columns(4)
    col_a.metric("Candidate arrays", len(crispr_arrays))
    col_b.metric("Extracted spacers", sum(array.spacer_count for array in crispr_arrays))
    col_c.metric("Spacer hits", len(spacer_hits))
    col_d.metric(
        "Duplicates excluded",
        (len(bacteria_records) - len(unique_bacteria_records))
        + (len(phage_records) - len(unique_phage_records)),
    )
    if pam_mode == "Auto from predicted subtype":
        st.caption(
            f"Automatic subtype prediction evaluated {len(cas_predictions_by_array)} arrays. "
            "PAM/PFS rules are applied only for high-confidence subtypes with curated rules."
        )

    with st.expander("Bacteria-phage evidence table", expanded=False):
        if evidence_matrix.empty:
            st.info("Upload at least one bacterial FASTA and one phage FASTA to build the matrix.")
        else:
            st.dataframe(evidence_matrix, use_container_width=True)

    with st.expander("Candidate CRISPR arrays", expanded=False):
        array_table = summarize_crispr_arrays(crispr_arrays)
        if array_table.empty:
            st.info("No candidate CRISPR arrays detected with the current exact-repeat MVP detector.")
        else:
            st.dataframe(array_table, use_container_width=True)

    with st.expander("Extracted spacers", expanded=False):
        spacer_table = summarize_spacers(crispr_arrays)
        if spacer_table.empty:
            st.info("No spacers extracted yet.")
        else:
            st.dataframe(spacer_table, use_container_width=True)

    with st.expander("Detailed spacer-phage hits", expanded=False):
        hit_table = summarize_spacer_hits(spacer_hits)
        if hit_table.empty:
            st.info("No exact spacer-protospacer hits detected in the uploaded phage genomes.")
        else:
            st.dataframe(hit_table, use_container_width=True)

    with st.expander("Exploratory PAM/PFS subtype support", expanded=False):
        pam_subtype_table = summarize_pam_subtype_support(spacer_hits)
        st.caption(
            "Diagnostic only: observed protospacer flanks are tested against curated subtype "
            "PAM/PFS rules and compared with the repeat-based subtype prediction. This does "
            "not change scoring or the primary Cas subtype call."
        )
        if pam_subtype_table.empty:
            st.info("No spacer hits available for PAM/PFS subtype support analysis.")
        else:
            st.dataframe(pam_subtype_table, use_container_width=True)

    with st.expander("Pipeline stages", expanded=False):
        st.dataframe(run_summary.stage_table(), use_container_width=True)
else:
    latest_result = st.session_state.latest_result
    if latest_result is not None:
        render_latest_analysis_result(latest_result)
    else:
        st.info("Upload bacterial and phage FASTA files to begin.")
