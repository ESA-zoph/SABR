from __future__ import annotations

import base64
from time import perf_counter
from pathlib import Path

import streamlit as st

from crispr_phage_predictor.external.blast import BLASTN_COMMAND, MAKEBLASTDB_COMMAND
from crispr_phage_predictor.external.minced import (
    MINCED_COMMAND,
    active_minced_backend,
    minced_available,
)
from crispr_phage_predictor.external.tools import missing_tool_message, tool_available
from crispr_phage_predictor.cas_prediction import predict_array_cas_subtypes
from crispr_phage_predictor.io import (
    deduplicate_records,
    parse_uploaded_fastas,
    summarize_accession_conflicts,
    summarize_duplicate_records,
    summarize_records,
    summarize_uploaded_files,
)
from crispr_phage_predictor.ml.model_artifact import model_artifact_metadata
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


def image_data_uri(path: Path) -> str:
    if not path.exists():
        return ""
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
                <div class="sabr-lab">CRISPR-phage targeting evidence mapper | The Phage Lab, Faculty of Medecine, AUB</div>
            </div>
            <div class="sabr-logo-wrap">{logo_markup}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


render_brand_header()


def _format_seconds(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    minutes, remainder = divmod(int(seconds), 60)
    if minutes < 60:
        return f"{minutes}m {remainder}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes}m"


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

with st.sidebar:
    st.header("Inputs")
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
    )
    matching_label = st.selectbox(
        "Spacer-phage matching",
        ["Auto recommended", "Internal exact match", "BLASTN"],
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
        detection_method = "internal"
    else:
        detection_method = "minced" if detection_label == "MinCED-compatible" else "internal"

    if matching_label == "Auto recommended":
        matching_method = "blast" if blast_available else "internal"
    else:
        matching_method = "blast" if matching_label == "BLASTN" else "internal"

    st.caption(f"Selected detection backend: {detection_method}")
    st.caption(f"Selected matching backend: {matching_method}")
    if detection_method == "minced":
        st.warning(
            "MinCED-compatible detection is available for benchmarking, but it may be slow on "
            "whole-genome uploads. Start with a small genome or subset before large batches."
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
    run_analysis = st.button("Run CRISPR-phage analysis", type="primary")

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

if bacterial_files or phage_files:
    bacteria_records = parse_uploaded_fastas(bacterial_files)
    phage_records = parse_uploaded_fastas(phage_files)
    bacterial_duplicate_table = summarize_duplicate_records(bacteria_records)
    phage_duplicate_table = summarize_duplicate_records(phage_records)
    bacterial_accession_conflicts = summarize_accession_conflicts(bacteria_records)
    phage_accession_conflicts = summarize_accession_conflicts(phage_records)
    unique_bacteria_records = deduplicate_records(bacteria_records)
    unique_phage_records = deduplicate_records(phage_records)
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

    if not run_analysis:
        st.info(
            "Files are parsed. Press 'Run CRISPR-phage analysis' in the sidebar to detect "
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

    detection_progress = st.progress(0, text="Preparing CRISPR array detection...")
    detection_status = st.empty()
    analysis_start = perf_counter()

    def update_detection_progress(completed: int, total: int, record) -> None:
        elapsed = perf_counter() - analysis_start
        average = elapsed / completed if completed else 0
        remaining = max(total - completed, 0) * average
        detection_progress.progress(
            completed / total if total else 1.0,
            text=(
                f"CRISPR detection: {completed}/{total} bacterial sequences "
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

    matching_progress = st.progress(0, text="Preparing spacer-phage matching...")
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
    st.info("Upload bacterial and phage FASTA files to begin.")
