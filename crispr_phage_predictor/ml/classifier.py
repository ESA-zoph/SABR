from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from sklearn.ensemble import RandomForestClassifier

from crispr_phage_predictor.ml.dataset import load_repeat_cas_training_table
from crispr_phage_predictor.ml.features import (
    DEFAULT_KMER_SIZES,
    build_repeat_feature_table,
    feature_columns,
)


@dataclass(frozen=True)
class CasSubtypePrediction:
    cas_subtype: str
    confidence: float
    probabilities: dict[str, float]


@dataclass(frozen=True)
class SimilarityPrediction:
    cas_subtype: str
    confidence: float
    best_identity: float
    matched_repeat: str


class RepeatCasSubtypeClassifier:
    """Baseline repeat/array feature classifier for CRISPR-Cas subtype."""

    def __init__(
        self,
        kmer_sizes: tuple[int, ...] = DEFAULT_KMER_SIZES,
        random_state: int = 42,
        n_estimators: int = 200,
    ) -> None:
        self.kmer_sizes = kmer_sizes
        self.model = RandomForestClassifier(
            n_estimators=n_estimators,
            random_state=random_state,
            class_weight="balanced",
        )
        self._feature_names: list[str] = []
        self._classes: list[str] = []

    @property
    def is_fitted(self) -> bool:
        return bool(self._feature_names and self._classes)

    def fit(self, training_table: pd.DataFrame) -> "RepeatCasSubtypeClassifier":
        feature_table = build_repeat_feature_table(training_table, kmer_sizes=self.kmer_sizes)
        self._feature_names = feature_columns(feature_table)
        self._classes = sorted(str(label) for label in feature_table["cas_subtype"].unique())
        if len(self._classes) < 2:
            raise ValueError("At least two Cas subtypes are required to train the classifier")
        self.model.fit(feature_table[self._feature_names], feature_table["cas_subtype"])
        return self

    def fit_csv(self, training_table_path: str | Path) -> "RepeatCasSubtypeClassifier":
        return self.fit(load_repeat_cas_training_table(training_table_path))

    def predict_table(self, training_like_table: pd.DataFrame) -> list[CasSubtypePrediction]:
        if not self.is_fitted:
            raise ValueError("Classifier must be fitted before prediction")
        feature_table = build_repeat_feature_table(training_like_table, kmer_sizes=self.kmer_sizes)
        aligned_features = feature_table.reindex(columns=self._feature_names, fill_value=0.0)
        probabilities = self.model.predict_proba(aligned_features)
        predictions: list[CasSubtypePrediction] = []
        for row_probabilities in probabilities:
            probability_map = {
                str(label): float(probability)
                for label, probability in zip(self.model.classes_, row_probabilities)
            }
            cas_subtype = max(probability_map, key=probability_map.get)
            predictions.append(
                CasSubtypePrediction(
                    cas_subtype=cas_subtype,
                    confidence=probability_map[cas_subtype],
                    probabilities=probability_map,
                )
            )
        return predictions

    def predict_one(
        self,
        repeat_sequence: str,
        spacer_count: int = 0,
        mean_spacer_length: float = 0.0,
    ) -> CasSubtypePrediction:
        repeat = repeat_sequence.upper()
        table = pd.DataFrame(
            [
                {
                    "genome_id": "query",
                    "contig_id": "query",
                    "repeat_sequence": repeat,
                    "repeat_length": len(repeat),
                    "spacer_count": spacer_count,
                    "mean_spacer_length": mean_spacer_length,
                    "cas_type": "",
                    "cas_subtype": "",
                }
            ]
        )
        return self.predict_table(table)[0]


class NearestRepeatClassifier:
    """Interpretable nearest-repeat baseline for Cas subtype prediction."""

    def __init__(self) -> None:
        self._examples: list[tuple[str, str]] = []

    @property
    def is_fitted(self) -> bool:
        return bool(self._examples)

    def fit(self, training_table: pd.DataFrame) -> "NearestRepeatClassifier":
        self._examples = [
            (str(row["repeat_sequence"]).upper(), str(row["cas_subtype"]))
            for _, row in training_table.iterrows()
            if str(row.get("repeat_sequence", "")).strip()
            and str(row.get("cas_subtype", "")).strip()
        ]
        if not self._examples:
            raise ValueError("At least one labeled repeat is required")
        return self

    def predict_table(self, training_like_table: pd.DataFrame) -> list[SimilarityPrediction]:
        if not self.is_fitted:
            raise ValueError("Classifier must be fitted before prediction")
        return [
            self.predict_one(str(row["repeat_sequence"]).upper())
            for _, row in training_like_table.iterrows()
        ]

    def predict_one(self, repeat_sequence: str) -> SimilarityPrediction:
        if not self.is_fitted:
            raise ValueError("Classifier must be fitted before prediction")
        query = repeat_sequence.upper()
        scored = [
            (_global_identity(query, repeat), repeat, subtype)
            for repeat, subtype in self._examples
        ]
        scored.sort(key=lambda item: item[0], reverse=True)
        best_identity, matched_repeat, cas_subtype = scored[0]
        return SimilarityPrediction(
            cas_subtype=cas_subtype,
            confidence=best_identity,
            best_identity=best_identity,
            matched_repeat=matched_repeat,
        )


def _global_identity(left: str, right: str) -> float:
    if not left and not right:
        return 1.0
    max_length = max(len(left), len(right))
    if max_length == 0:
        return 0.0
    matches = sum(1 for left_base, right_base in zip(left, right) if left_base == right_base)
    return matches / max_length


CasTypeClassifier = RepeatCasSubtypeClassifier
