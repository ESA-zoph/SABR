from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.ml.classifier import NearestRepeatClassifier
from crispr_phage_predictor.ml.dataset import load_repeat_cas_training_table
from crispr_phage_predictor.ml.features import build_repeat_feature_table
from crispr_phage_predictor.ml.model_artifact import (
    DEFAULT_MODEL_PATH,
    load_artifact,
    prediction_table_for_arrays,
)


DEFAULT_RUNTIME_TRAINING_TABLE = (
    Path("data") / "training" / "repeats_cas_types_augmented_vink_genbank_targeted.csv"
)


CURATED_PAM_RULES_BY_SUBTYPE = {
    "I-E": "5prime:AWG",
    "I-F": "genomic_3prime:GG",
    "I-A": "5prime:CCN",
    "II-A": "3prime:NGG",
    "V-A": "5prime:TTTN",
}


@dataclass(frozen=True)
class ArrayCasPrediction:
    array_id: str
    cas_subtype: str
    cas_subtype_confidence: float
    prediction_method: str
    pam_rule: str
    pam_rule_source: str


def predict_array_cas_subtypes(
    arrays: list[CrisprArray],
    training_table_path: Path = DEFAULT_RUNTIME_TRAINING_TABLE,
    model_artifact_path: Path = DEFAULT_MODEL_PATH,
    min_confidence_for_pam: float = 0.9,
) -> dict[str, ArrayCasPrediction]:
    if not arrays:
        return {}
    artifact_predictions = _predict_with_artifact(
        arrays=arrays,
        model_artifact_path=model_artifact_path,
        min_confidence_for_pam=min_confidence_for_pam,
    )
    if artifact_predictions is not None:
        return artifact_predictions
    return _predict_with_nearest_repeat(
        arrays=arrays,
        training_table_path=training_table_path,
        min_confidence_for_pam=min_confidence_for_pam,
    )


def _predict_with_artifact(
    arrays: list[CrisprArray],
    model_artifact_path: Path,
    min_confidence_for_pam: float,
) -> dict[str, ArrayCasPrediction] | None:
    if not model_artifact_path.exists():
        return None
    artifact = load_artifact(model_artifact_path)
    query_table = prediction_table_for_arrays(arrays)
    feature_table = build_repeat_feature_table(query_table)
    aligned_features = feature_table.reindex(columns=artifact.feature_names, fill_value=0.0)
    probabilities = artifact.model.predict_proba(aligned_features)
    predictions: dict[str, ArrayCasPrediction] = {}
    for index, array in enumerate(arrays):
        probability_map = {
            str(label): float(probability)
            for label, probability in zip(artifact.model.classes_, probabilities[index])
        }
        cas_subtype = max(probability_map, key=probability_map.get)
        confidence = probability_map[cas_subtype]
        predictions[array.array_id] = _array_prediction(
            array=array,
            cas_subtype=cas_subtype,
            confidence=confidence,
            method=str(artifact.metadata.get("method", "extra_trees")),
            min_confidence_for_pam=min_confidence_for_pam,
        )
    return predictions


def _predict_with_nearest_repeat(
    arrays: list[CrisprArray],
    training_table_path: Path,
    min_confidence_for_pam: float,
) -> dict[str, ArrayCasPrediction]:
    if not arrays or not training_table_path.exists():
        return {}

    training_table = load_repeat_cas_training_table(training_table_path)
    classifier = NearestRepeatClassifier().fit(training_table)
    predictions: dict[str, ArrayCasPrediction] = {}
    for array in arrays:
        prediction = classifier.predict_one(array.repeat_consensus)
        predictions[array.array_id] = _array_prediction(
            array=array,
            cas_subtype=prediction.cas_subtype,
            confidence=prediction.confidence,
            method="nearest_repeat",
            min_confidence_for_pam=min_confidence_for_pam,
        )
    return predictions


def _array_prediction(
    array: CrisprArray,
    cas_subtype: str,
    confidence: float,
    method: str,
    min_confidence_for_pam: float,
) -> ArrayCasPrediction:
    pam_rule = ""
    pam_rule_source = "not_available"
    if confidence >= min_confidence_for_pam:
        pam_rule = CURATED_PAM_RULES_BY_SUBTYPE.get(cas_subtype, "")
        pam_rule_source = "curated_subtype_catalog" if pam_rule else "no_curated_rule"
    else:
        pam_rule_source = "subtype_confidence_below_threshold"

    return ArrayCasPrediction(
        array_id=array.array_id,
        cas_subtype=cas_subtype,
        cas_subtype_confidence=confidence,
        prediction_method=method,
        pam_rule=pam_rule,
        pam_rule_source=pam_rule_source,
    )
