import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, ConfusionMatrixDisplay
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder

from rules_engine import classify_metabolizer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/sample_genotypes.csv")
    ap.add_argument("--outdir", default="reports/figures")
    args = ap.parse_args()

    df = pd.read_csv(args.data)

    # Labels still come from the 3 core pharmacogenomic SNPs
    df["label"] = df.apply(
        lambda r: classify_metabolizer(
            int(r["rs4244285"]),
            int(r["rs12248560"]),
            int(r["rs4986893"])
        ),
        axis=1
    )

    le_sex = LabelEncoder()
    le_pop = LabelEncoder()

    df["sex_encoded"] = le_sex.fit_transform(df["sex"])
    df["pop_encoded"] = le_pop.fit_transform(df["super_pop"])

    print(df.groupby(["super_pop", "label"]).size())

    # Use all 5 SNPs + sex/pop
    X = df[
        [
            "rs4244285",
            "rs4986893",
            "rs12248560",
            "rs28399504",
            "rs41291556",
            "sex_encoded",
            "pop_encoded",
        ]
    ].to_numpy(dtype=np.int64)

    y = df["label"].to_numpy()

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )

    models = {
        "rf": RandomForestClassifier(n_estimators=200, random_state=42),
    }

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    Path("models").mkdir(parents=True, exist_ok=True)

    for name, model in models.items():
        model.fit(X_train, y_train)
        joblib.dump(model, f"models/{name}_model.joblib")
        preds = model.predict(X_test)
        probs = model.predict_proba(X_test[:5])

        print("\n=== Example prediction probabilities ===")
        for i, row in enumerate(probs):
            print(f"\nSample {i+1}")
            for cls, p in zip(model.classes_, row):
                print(f"{cls:15s} {p:.4f}")

        print(f"\n=== {name} ===")
        print(classification_report(y_test, preds, zero_division=0))

        disp = ConfusionMatrixDisplay.from_predictions(y_test, preds, xticks_rotation=45)
        plt.title(f"Confusion Matrix - {name}")
        plt.savefig(outdir / f"confusion_{name}.png", dpi=200, bbox_inches="tight")
        plt.close()

        feature_names = [
            "rs4244285",
            "rs4986893",
            "rs12248560",
            "rs28399504",
            "rs41291556",
            "sex_encoded",
            "pop_encoded",
            ]

        print("\n=== Random Forest Feature Importance ===")
        for name, imp in sorted(zip(feature_names, model.feature_importances_), key=lambda x: x[1], reverse=True):
            print(f"{name:15s} {imp:.4f}")        


if __name__ == "__main__":
    main()