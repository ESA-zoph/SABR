from __future__ import annotations

import streamlit as st

from crispr_phage_predictor.io import (
    parse_uploaded_fastas,
    summarize_records,
    summarize_uploaded_files,
)
from crispr_phage_predictor.matching import find_spacer_hits
from crispr_phage_predictor.pipeline import (
    build_exact_match_heatmap,
    build_resistance_evidence_matrix,
    build_initial_run_summary,
    detect_arrays_for_records,
    summarize_crispr_arrays,
    summarize_spacer_hits,
    summarize_spacers,
)


st.set_page_config(
    page_title="CRISPR-Phage Resistance Predictor",
    layout="wide",
)

st.title("CRISPR-Phage Resistance Predictor")
st.caption("Early scaffold for batch CRISPR-phage resistance likelihood analysis.")


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
    run_analysis = st.button("Run CRISPR-phage analysis", type="primary")

st.warning(
    "This version uses an exact-repeat MVP detector for candidate CRISPR arrays. "
    "Results are preliminary and should be benchmarked before scientific interpretation."
)

if bacterial_files or phage_files:
    bacteria_records = parse_uploaded_fastas(bacterial_files)
    phage_records = parse_uploaded_fastas(phage_files)
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

    st.subheader("Phage FASTA Records")
    phage_table = summarize_records(phage_records)
    if phage_table.empty:
        st.info("No phage FASTA records loaded.")
    else:
        with st.expander("View phage FASTA records", expanded=False):
            st.dataframe(phage_table, use_container_width=True)

    if not run_analysis:
        st.info(
            "Files are parsed. Press 'Run CRISPR-phage analysis' in the sidebar to detect "
            "candidate arrays and match spacers against phages."
        )
        st.stop()

    if not bacteria_records:
        st.error(
            "No bacterial FASTA records were parsed. Check that the bacterial file starts "
            "with FASTA headers like '>contig_1' and is not GenBank, FASTQ, or another format."
        )
        st.stop()

    if not phage_records:
        st.error(
            "No phage FASTA records were parsed. Check that the phage file starts with FASTA headers."
        )
        st.stop()

    with st.spinner("Detecting candidate CRISPR arrays and matching spacers..."):
        crispr_arrays = detect_arrays_for_records(bacteria_records)
        spacer_hits = find_spacer_hits(crispr_arrays, phage_records)

    evidence_matrix = build_resistance_evidence_matrix(
        bacteria_records=bacteria_records,
        phage_records=phage_records,
        hits=spacer_hits,
    )
    heatmap = build_exact_match_heatmap(evidence_matrix)

    st.subheader("Bacteria-Phage Exact-Match Heatmap")
    st.caption(
        "Cell values are unique bacterial spacers with exact protospacer matches in each phage. "
        "This is evidence of candidate CRISPR targeting, not a final resistance probability."
    )
    if heatmap.empty:
        st.info("Upload at least one bacterial FASTA and one phage FASTA to build the heatmap.")
    else:
        render_exact_match_heatmap(heatmap)

    st.subheader("Summary")
    col_a, col_b, col_c = st.columns(3)
    col_a.metric("Candidate arrays", len(crispr_arrays))
    col_b.metric("Extracted spacers", sum(array.spacer_count for array in crispr_arrays))
    col_c.metric("Exact spacer hits", len(spacer_hits))

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

    with st.expander("Pipeline stages", expanded=False):
        st.dataframe(run_summary.stage_table(), use_container_width=True)
else:
    st.info("Upload bacterial and phage FASTA files to begin.")
