from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.interactions import (
    INTERACTION_COLUMNS,
    eop_class_from_value,
    parse_eop,
    susceptibility_from_eop_class,
    validate_interaction_table,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Import small traceable EOP tables embedded from open literature."
    )
    parser.add_argument("output_tsv", type=Path)
    args = parser.parse_args()

    rows = []
    rows.extend(_kpp5_rows())
    rows.extend(_k23_klebsiella_rows())
    rows.extend(_kpkp_cocktail_rows())
    rows.extend(_stm2_salmonella_rows())
    rows.extend(_vabwu2101_rows())
    rows.extend(_vecpw8_rows())
    rows.extend(_kp1_kp12_rows())
    rows.extend(_f48_rows())
    rows.extend(_xp4_rows())
    table = pd.DataFrame(rows, columns=INTERACTION_COLUMNS)
    validate_interaction_table(table)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_tsv, sep="\t", index=False)


def _kpp5_rows() -> list[dict[str, str]]:
    rows = []
    entries = [
        ("Klebsiella pneumoniae", "K_pneumoniae_1_host", "food", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_2", "food", "M"),
        ("Klebsiella pneumoniae", "K_pneumoniae_3", "food", "M"),
        ("Klebsiella pneumoniae", "K_pneumoniae_4", "food", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_5", "food", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_6", "clinical", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_7", "clinical", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_8", "clinical", "M"),
        ("Klebsiella pneumoniae", "K_pneumoniae_9", "clinical", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_10", "clinical", "M"),
        ("Klebsiella pneumoniae", "K_pneumoniae_11", "clinical", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_12", "clinical", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_13", "clinical", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_14", "clinical", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_15", "clinical", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_16", "clinical", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_17", "clinical", "M"),
        ("Klebsiella pneumoniae", "K_pneumoniae_18", "clinical", "H"),
        ("Klebsiella pneumoniae", "K_pneumoniae_19", "clinical", "H"),
        ("Escherichia coli", "E_coli_1", "clinical", "N"),
        ("Escherichia coli", "E_coli_2", "clinical", "N"),
        ("Escherichia coli", "E_coli_3", "clinical", "N"),
        ("Escherichia coli", "E_coli_4", "clinical", "N"),
        ("Escherichia coli", "E_coli_5", "clinical", "N"),
        ("Escherichia coli", "E_coli_6", "clinical", "N"),
        ("Salmonella enterica", "S_Typhimurium_1", "food", "N"),
        ("Salmonella enterica", "S_Typhimurium_2", "food", "N"),
        ("Salmonella enterica", "S_Typhimurium_3", "food", "N"),
        ("Pseudomonas aeruginosa", "P_aeruginosa_1", "clinical", "N"),
        ("Pseudomonas aeruginosa", "P_aeruginosa_2", "clinical", "N"),
        ("Pseudomonas aeruginosa", "P_aeruginosa_3", "food", "N"),
    ]
    class_map = {
        "H": ("high", "H; 0.5-1.0"),
        "M": ("medium", "M; 0.2-0.4"),
        "L": ("low", "L; 0.001-0.1"),
        "N": ("none", "N; <0.001"),
    }
    for bacterium, strain, source, code in entries:
        eop_class, raw = class_map[code]
        rows.append(
            _row(
                interaction_id=f"sofy2021_{_safe(strain)}_kpp5_eop",
                source_key="Sofy2021_KPP5",
                source_type="paper",
                pmid="33805298",
                doi="10.3390/biomedicines9040342",
                assay_type="eop",
                bacterium=bacterium,
                strain=strain,
                phage="KPP-5",
                raw_eop=raw,
                eop_relation="range",
                eop_value="",
                eop_class=eop_class,
                plaque_result=_plaque_result(eop_class),
                experimental_conditions=(
                    "Spot test and EOP host-range table; H 0.5-1.0, M 0.2-0.4, "
                    "L 0.001-0.1, N <0.001."
                ),
                notes=f"Table 1 source category: {source}.",
            )
        )
    return rows


def _k23_klebsiella_rows() -> list[dict[str, str]]:
    strains = [
        ("Kp-9068", "K23", "ST11", ["1.0", "1.0", "1.0"]),
        ("KPi4275", "K23", "ST11", ["1.0", "1.0", "0.006"]),
        ("KPB2304-15", "K23", "ST11", ["1.0", "1.0", "0.006"]),
        ("KPB536-17-2", "K23", "ST1869", ["0.03", "0.01", "0.01"]),
        ("KPi1748", "K2", "ST65", ["0", "0", "0"]),
        ("KPi3014", "K2", "ST2174", ["0", "0", "0"]),
        ("KPi8289", "K57", "ST218", ["0", "0", "0"]),
        ("KPB2580", "K1", "ST23", ["0", "0", "0"]),
    ]
    phages = ["KpS8", "vB_KpnM_Seu621", "vB_KpnP_Dlv622"]
    rows = []
    for strain, cps_type, mlst, values in strains:
        for phage, raw_eop in zip(phages, values):
            relation, value = parse_eop(raw_eop)
            eop_class = eop_class_from_value(value, relation)
            rows.append(
                _row(
                    interaction_id=f"latka2021_{_safe(strain)}_{_safe(phage)}_eop",
                    source_key="Latka2021_K23Klebsiella",
                    source_type="paper",
                    pmid="34093432",
                    doi="10.3389/fmicb.2021.669618",
                    assay_type="eop",
                    bacterium="Klebsiella pneumoniae",
                    strain=strain,
                    phage=phage,
                    raw_eop=raw_eop,
                    eop_relation=relation,
                    eop_value="" if value is None else f"{value:.8g}",
                    eop_class=eop_class,
                    plaque_result=_plaque_result(eop_class),
                    experimental_conditions=(
                        "Table 2 EOP analysis for K23-specific Klebsiella phages."
                    ),
                    notes=f"cps-type {cps_type}; MLST {mlst}.",
                )
            )
    return rows


def _kpkp_cocktail_rows() -> list[dict[str, str]]:
    strains = {
        "KPCI1": {"KPKp": "0.64", "KSKp": "0.026", "KPKp_KSKp_cocktail": "0.85"},
        "KPCI2": {"KPKp": "0.05", "KSKp": "6.78", "KPKp_KSKp_cocktail": "0.90"},
        "KPCI3": {"KPKp": "0", "KSKp": "0.26", "KPKp_KSKp_cocktail": "0.33"},
        "KPCI4": {"KPKp": "0.55", "KSKp": "0.0002", "KPKp_KSKp_cocktail": "0.70"},
    }
    rows = []
    for strain, phage_values in strains.items():
        for phage, raw_eop in phage_values.items():
            relation, value = parse_eop(raw_eop)
            eop_class = eop_class_from_value(value, relation)
            rows.append(
                _row(
                    interaction_id=f"singh2025_{_safe(strain)}_{_safe(phage)}_eop",
                    source_key="Singh2025_KpCocktail",
                    source_type="paper",
                    pmid="",
                    doi="10.3389/fmicb.2025.1588472",
                    assay_type="eop",
                    bacterium="Klebsiella pneumoniae",
                    strain=strain,
                    phage=phage,
                    raw_eop=raw_eop,
                    eop_relation=relation,
                    eop_value="" if value is None else f"{value:.8g}",
                    eop_class=eop_class,
                    plaque_result=_plaque_result(eop_class),
                    experimental_conditions=(
                        "Table 3 EOP against MDR K. pneumoniae clinical isolates, "
                        "compared with K. pneumoniae ATCC 700603."
                    ),
                    notes="Individual phages KPKp/KSKp and combined cocktail reported.",
                )
            )
    return rows


def _stm2_salmonella_rows() -> list[dict[str, str]]:
    entries = [
        ("Salmonella enterica", "S_Typhimurium_ST_4_host", ["M", "H", "H"]),
        ("Salmonella enterica", "S_Typhimurium_ST_7", ["H", "H", "M"]),
        ("Salmonella enterica", "S_Typhimurium_ST_9", ["N", "M", "N"]),
        ("Salmonella enterica", "S_Typhimurium_ST_14", ["N", "H", "N"]),
        ("Salmonella enterica", "S_Typhimurium_ST_16", ["H", "H", "M"]),
        ("Salmonella enterica", "S_Typhimurium_ST_19", ["M", "H", "M"]),
        ("Salmonella enterica", "S_Typhimurium_ST_21", ["N", "M", "M"]),
        ("Salmonella enterica", "S_Typhimurium_ST_22", ["N", "H", "N"]),
        ("Salmonella enterica", "S_Typhimurium_ST_27", ["M", "H", "M"]),
        ("Salmonella enterica", "S_Typhimurium_ST_30", ["H", "H", "N"]),
        ("Salmonella enterica", "S_Typhimurium_ST_33", ["N", "M", "N"]),
        ("Salmonella enterica", "S_Typhimurium_ST_35", ["H", "H", "H"]),
        ("Salmonella enterica", "S_Typhimurium_ST_36", ["N", "M", "N"]),
        ("Salmonella enterica", "S_Typhimurium_ST_37", ["M", "H", "N"]),
        ("Salmonella enterica", "S_Typhimurium_ST_38", ["M", "H", "M"]),
        ("Salmonella enterica", "S_Typhimurium_ST_41", ["L", "H", "N"]),
        ("Salmonella enterica", "S_Typhimurium_ST_43", ["L", "H", "N"]),
        ("Salmonella enterica", "S_Typhimurium_ST_45", ["M", "H", "L"]),
        ("Salmonella enterica", "S_Typhimurium_ST_47", ["L", "H", "L"]),
        ("Salmonella enterica", "S_Typhimurium_ST_49", ["N", "M", "N"]),
        ("Salmonella enterica", "S_Typhimurium_ST_53", ["H", "H", "M"]),
        ("Salmonella enterica", "S_Typhimurium_ST_55", ["N", "M", "N"]),
        ("Salmonella enterica", "S_Typhimurium_ST_56", ["H", "H", "M"]),
        ("Pseudomonas aeruginosa", "PsaCI_1", ["N", "M", "N"]),
        ("Pseudomonas aeruginosa", "PsaCI_2", ["N", "M", "N"]),
        ("Pseudomonas aeruginosa", "PsaFI_1", ["N", "H", "N"]),
        ("Pseudomonas aeruginosa", "PsaFI_2", ["N", "H", "N"]),
        ("Pseudomonas aeruginosa", "PsaFI_1_table_cont", ["N", "H", "N"]),
        ("Staphylococcus aureus", "SaFI_1", ["N", "L", "N"]),
        ("Staphylococcus aureus", "SaFI_2", ["N", "M", "N"]),
        ("Staphylococcus aureus", "SaFI_3", ["N", "L", "N"]),
        ("Escherichia coli", "EcCI_1", ["N", "M", "N"]),
        ("Escherichia coli", "EcCI_2", ["N", "M", "N"]),
        ("Escherichia coli", "EcCI_3", ["N", "M", "N"]),
        ("Escherichia coli", "EcCI_4", ["N", "M", "N"]),
        ("Klebsiella pneumoniae", "KpCI_1", ["N", "N", "N"]),
    ]
    phages = ["vB_STS_1", "vB_STM_2", "vB_STS_3"]
    rows = []
    for bacterium, strain, classes in entries:
        for phage, code in zip(phages, classes):
            eop_class, raw = _category_to_class_and_raw(code)
            rows.append(
                _row(
                    interaction_id=f"abdelhadi2021_{_safe(strain)}_{_safe(phage)}_eop",
                    source_key="Abdelhadi2021_STM2",
                    source_type="paper",
                    pmid="",
                    doi="10.3390/su132111602",
                    assay_type="eop",
                    bacterium=bacterium,
                    strain=strain,
                    phage=phage,
                    raw_eop=raw,
                    eop_relation="range",
                    eop_value="",
                    eop_class=eop_class,
                    plaque_result=_plaque_result(eop_class),
                    experimental_conditions=(
                        "Table 2 spot test and EOP categories; H 0.5-1.0, "
                        "M 0.2-0.4, L 0.001-0.1, N inefficient."
                    ),
                    notes=(
                        "The source table repeats one P. aeruginosa label as PsaFI-1; "
                        "the second occurrence is kept as PsaFI_1_table_cont pending "
                        "manual author-table verification."
                    )
                    if strain == "PsaFI_1_table_cont"
                    else "Imported from Table 2 category matrix.",
                )
            )
    return rows


def _vabwu2101_rows() -> list[dict[str, str]]:
    entries = [
        ("Acinetobacter baumannii", "ABPW0181", "0.64"),
        ("Acinetobacter baumannii", "ABPW0182", "0.55"),
        ("Acinetobacter baumannii", "ABPW0183", "0"),
        ("Acinetobacter baumannii", "ABPW0184", "0"),
        ("Acinetobacter baumannii", "ABPW0185_host", "1"),
        ("Acinetobacter baumannii", "ABPW0186", "0.17"),
        ("Acinetobacter baumannii", "ABPW0187", "0.21"),
        ("Acinetobacter baumannii", "ABPW0188", "0"),
        ("Acinetobacter baumannii", "ABPW0189", "0.08"),
        ("Acinetobacter baumannii", "ABPW0190", "0.91"),
        ("Acinetobacter baumannii", "ABPW0191", "0.09"),
        ("Acinetobacter baumannii", "ABPW0192", "0.32"),
        ("Acinetobacter baumannii", "ABPW0193", "0.48"),
        ("Acinetobacter baumannii", "ABPW0194", "0"),
        ("Acinetobacter baumannii", "ABPW0195", "0.72"),
        ("Acinetobacter baumannii", "ABPW0196", "0.27"),
        ("Acinetobacter baumannii", "ABPW0197", "0.93"),
        ("Acinetobacter baumannii", "ABPW0198", "0"),
        ("Acinetobacter baumannii", "ABPW0199", "0"),
        ("Acinetobacter baumannii", "ABPW0200", "0.06"),
        ("Klebsiella pneumoniae", "K_pneumoniae_control", "0"),
        ("Staphylococcus aureus", "MRSA_control", "0"),
    ]
    rows = []
    for bacterium, strain, raw_eop in entries:
        relation, value = parse_eop(raw_eop)
        eop_class = eop_class_from_value(value, relation)
        rows.append(
            _row(
                interaction_id=f"wintachai2022_{_safe(strain)}_vabwu2101_eop",
                source_key="Wintachai2022_vABWU2101",
                source_type="paper",
                pmid="35215915",
                doi="10.3390/v14020194",
                assay_type="eop",
                bacterium=bacterium,
                strain=strain,
                phage="vABWU2101",
                raw_eop=raw_eop,
                eop_relation=relation,
                eop_value="" if value is None else f"{value:.8g}",
                eop_class=eop_class,
                plaque_result=_plaque_result(eop_class),
                experimental_conditions=(
                    "Table 1 host range and EOP; duplicate independent experiments "
                    "with duplicate assay."
                ),
                notes="MDR A. baumannii SSTI isolate panel plus K. pneumoniae/MRSA controls.",
            )
        )
    return rows


def _vecpw8_rows() -> list[dict[str, str]]:
    entries = [
        ("PW001", "0"),
        ("PW002", "0.8"),
        ("PW003", "0.65"),
        ("PW004", "0"),
        ("PW005_host", "1"),
        ("PW006", "0.49"),
        ("PW007", "0"),
        ("PW008", "0.57"),
        ("PW009", "0"),
        ("PW010", "0.32"),
        ("PW011", "0"),
        ("PW012", "0.44"),
        ("PW013", "0.02"),
        ("PW014", "0"),
        ("PW015", "0.41"),
        ("PW016", "0.79"),
        ("PW017", "0"),
        ("PW018", "0.03"),
        ("PW019", "0"),
        ("PW020", "0"),
    ]
    rows = []
    for strain, raw_eop in entries:
        relation, value = parse_eop(raw_eop)
        eop_class = eop_class_from_value(value, relation)
        rows.append(
            _row(
                interaction_id=f"wintachai2024_{_safe(strain)}_vecpw8_eop",
                source_key="Wintachai2024_vECPW8",
                source_type="paper",
                pmid="39595283",
                doi="10.3390/antibiotics13111083",
                assay_type="eop",
                bacterium="Escherichia coli",
                strain=f"MDR_APEC_{strain}",
                phage="vECPW8",
                raw_eop=raw_eop,
                eop_relation=relation,
                eop_value="" if value is None else f"{value:.8g}",
                eop_class=eop_class,
                plaque_result=_plaque_result(eop_class),
                experimental_conditions=(
                    "Table 1 bacterial lysis efficacy and EOP; triplicate "
                    "experiments with duplicate plaque assays, repeated twice."
                ),
                notes="MDR avian pathogenic E. coli isolate panel.",
            )
        )
    return rows


def _kp1_kp12_rows() -> list[dict[str, str]]:
    entries = [
        ("Klebsiella pneumoniae", "K16-KPN-13-022", "100", "100"),
        ("Klebsiella pneumoniae", "K01-KPN-13-134", "40.44", "49.08"),
        ("Klebsiella pneumoniae", "K01-KPN-13-149", "21.98", "0"),
        ("Klebsiella pneumoniae", "K07-KPN-13-002", "97.65", "81.06"),
        ("Klebsiella pneumoniae", "K14-KPN-13-016", "83.39", "75.15"),
        ("Klebsiella pneumoniae", "K16-KPN-13-008", "0", "51.73"),
        ("Klebsiella pneumoniae", "K20-KPN-12-057", "14.77", "0"),
        ("Klebsiella pneumoniae", "K20-KPN-12-067", "44.30", "53.77"),
        ("Klebsiella pneumoniae", "K21-KPN-12-013", "62.42", "96.95"),
        ("Klebsiella pneumoniae", "K22-KPN-13-007", "0", "0"),
        ("Klebsiella pneumoniae", "K22-KPN-13-013", "14.60", "26.68"),
        ("Acinetobacter baumannii", "ATCC17978", "0", "0"),
        ("Citrobacter freundii", "clinical_isolate_15-0628", "0", "0"),
        ("Cronobacter sakazakii", "KCTC_2949", "0", "0"),
        ("Escherichia coli", "KCTC_1039", "0", "0"),
        ("Escherichia coli", "K01-ECO12-052", "0", "0"),
        ("Proteus mirabilis", "KCTC_2566", "0", "0"),
        ("Pseudomonas aeruginosa", "KCTC2004", "0", "0"),
        ("Salmonella enterica", "Typhimurium_ATCC14028", "0", "0"),
        ("Salmonella enterica", "Enteritidis_KCCM12021", "0", "0"),
    ]
    rows = []
    for bacterium, strain, kp1_percent, kp12_percent in entries:
        for phage, percent in [("KP1", kp1_percent), ("KP12", kp12_percent)]:
            value = float(percent) / 100.0
            eop_class = eop_class_from_value(value, "=")
            rows.append(
                _row(
                    interaction_id=f"kim2022_{_safe(strain)}_{_safe(phage)}_eop",
                    source_key="Kim2022_KP1_KP12",
                    source_type="paper",
                    pmid="36704203",
                    doi="10.3389/fmicb.2022.990910",
                    assay_type="eop",
                    bacterium=bacterium,
                    strain=strain,
                    phage=phage,
                    raw_eop=f"{percent}%",
                    eop_relation="=",
                    eop_value=f"{value:.8g}",
                    eop_class=eop_class,
                    plaque_result=_plaque_result(eop_class),
                    experimental_conditions=(
                        "Table 1 host specificity; EOP reported as percent with "
                        "K16-KPN-13-022 considered 100%; triplicate experiment."
                    ),
                    notes="Percent EOP converted to ratio for SABR schema.",
                )
            )
    return rows


def _f48_rows() -> list[dict[str, str]]:
    entries = [
        ("12C47", "ST101", "1"),
        ("12C73", "ST101", "0.9 ± 0.2"),
        ("KPC174", "ST1633", "0.44 ± 0.06"),
        ("5559", "ST101", "0.5 ± 0.2"),
        ("5583", "ST2502", "0"),
        ("C002", "ST101", "1.9 ± 0.7"),
        ("K13", "ST101", "2.62 ± 0.01"),
        ("K18", "ST101", "<0.001"),
        ("6071", "ST2502", "0"),
        ("12C29", "ST101", "2.2 ± 0.7"),
        ("5546", "ST101", "3.9 ± 0.3"),
        ("12C72", "ST101", "3.1 ± 0.3"),
        ("KPC220", "ST101", "0"),
        ("494647", "ST101", "0.9 ± 0.2"),
    ]
    rows = []
    for strain, sequence_type, raw_eop in entries:
        relation, value = parse_eop(raw_eop)
        eop_class = eop_class_from_value(value, relation)
        rows.append(
            _row(
                interaction_id=f"ciacci2018_{_safe(strain)}_vb_kpn_f48_eop",
                source_key="Ciacci2018_vB_Kpn_F48",
                source_type="paper",
                pmid="30205589",
                doi="10.3390/v10090482",
                assay_type="eop",
                bacterium="Klebsiella pneumoniae",
                strain=strain,
                phage="vB_Kpn_F48",
                reference_host="12C47",
                raw_eop=raw_eop,
                eop_relation=relation,
                eop_value="" if value is None else f"{value:.8g}",
                eop_class=eop_class,
                plaque_result=_plaque_result(eop_class),
                experimental_conditions=(
                    "Table 2 EOP for spot-test-sensitive K. pneumoniae isolates; "
                    "mean of three observations with production category."
                ),
                notes=f"Sequence type {sequence_type}.",
            )
        )
    return rows


def _xp4_rows() -> list[dict[str, str]]:
    entries = [
        ("P4", "K1", False, False, True, "100"),
        ("1", "K64", True, False, False, "0"),
        ("2", "K64", True, False, False, "0"),
        ("3", "K64", True, False, False, "0"),
        ("4", "K64", True, False, False, "0"),
        ("5", "K64", True, False, False, "0"),
        ("6", "K19", True, False, False, "0"),
        ("7", "K64", True, False, False, "0"),
        ("8", "K64", True, False, False, "0"),
        ("9", "K64", True, False, False, "0"),
        ("10", "K19", True, False, True, "51"),
        ("11", "K64", True, False, False, "0"),
        ("12", "K149", True, False, False, "0"),
        ("13", "K19", True, False, False, "0"),
        ("14", "K64", True, False, False, "0"),
        ("15", "K19", False, False, False, "0"),
        ("16", "K57", False, False, False, "0"),
        ("17", "K125", False, False, False, "0"),
        ("18", "K64", True, False, False, "0"),
        ("19", "K102", False, True, False, "0"),
        ("20", "K64", True, False, False, "0"),
    ]
    rows = []
    for strain, k_type, kpc, ndm, sensitive, percent in entries:
        value = float(percent) / 100.0
        eop_class = eop_class_from_value(value, "=")
        rows.append(
            _row(
                interaction_id=f"peng2025_{_safe(strain)}_vb_kp_xp4_eop",
                source_key="Peng2025_vB_Kp_XP4",
                source_type="paper",
                pmid="",
                doi="10.3389/fmicb.2025.1491961",
                assay_type="eop",
                bacterium="Klebsiella pneumoniae",
                strain=f"XP4_panel_{strain}",
                phage="vB_Kp_XP4",
                reference_host="P4",
                raw_eop=f"{percent}%",
                eop_relation="=",
                eop_value=f"{value:.8g}",
                eop_class=eop_class,
                plaque_result=_plaque_result(eop_class),
                experimental_conditions=(
                    "Table 2 host-range details; EOP reported as percent relative "
                    "to strain P4."
                ),
                notes=(
                    f"K-type {k_type}; KPC {'+' if kpc else '-'}; "
                    f"NDM {'+' if ndm else '-'}; sensitivity {'+' if sensitive else '-'}."
                ),
            )
        )
    return rows


def _row(**kwargs: str) -> dict[str, str]:
    row = {column: "" for column in INTERACTION_COLUMNS}
    row.update(
        {
            "reference_host": "reported source reference host",
            "susceptibility_label": susceptibility_from_eop_class(kwargs["eop_class"]),
            "anti_crispr_status": "not_evaluated",
            "crispr_interference_evidence": "not_evaluated",
            "curation_status": "curated",
            "curation_confidence": "high",
        }
    )
    row.update(kwargs)
    return row


def _plaque_result(eop_class: str) -> str:
    if eop_class == "none":
        return "no_plaques"
    if eop_class in {"trace", "low"}:
        return "pinpoint_plaques"
    if eop_class in {"medium", "high"}:
        return "clear_plaques"
    return "not_reported"


def _safe(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def _category_to_class_and_raw(code: str) -> tuple[str, str]:
    mapping = {
        "H": ("high", "H; 0.5-1.0"),
        "M": ("medium", "M; 0.2-0.4"),
        "L": ("low", "L; 0.001-0.1"),
        "N": ("none", "N; inefficient"),
    }
    return mapping[code]


if __name__ == "__main__":
    main()
