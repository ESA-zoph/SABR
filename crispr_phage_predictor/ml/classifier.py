from __future__ import annotations


class CasTypeClassifier:
    """Interface for future CRISPR-Cas type/subtype classifiers."""

    def fit(self, training_table_path: str) -> None:
        raise NotImplementedError("Classifier training is not implemented yet.")

    def predict(self, repeat_sequence: str) -> dict[str, object]:
        raise NotImplementedError("Classifier prediction is not implemented yet.")
