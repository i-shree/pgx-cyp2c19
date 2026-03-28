import argparse
from pathlib import Path

import joblib
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.preprocessing import LabelEncoder


def build_activity_score(df: pd.DataFrame) -> pd.Series:
    return (
        -2.0 * df["rs4244285"]
        -1.5 * df["rs4986893"]
        -1.2 * df["rs28399504"]
        -1.2 * df["rs41291556"]
        +1.0 * df["rs12248560"]
    )


def score_to_label(score: float) -> str:
    if score <= -2:
        return "Poor"
    elif score < 0:
        return "Intermediate"
    elif score == 0:
        return "Normal"
    else:
        return "Rapid"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/cyp2c19_1000g_5snp.csv")
    ap.add_argument("--outdir", default="reports/figures")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    Path("models").mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.data)

    le_sex = LabelEncoder()
    le_pop = LabelEncoder()

    df["sex_encoded"] = le_sex.fit_transform(df["sex"])
    df["pop_encoded"] = le_pop.fit_transform(df["super_pop"])

    df["activity_score"] = build_activity_score(df)

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
    ]

    y = df["activity_score"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42
    )

    model = RandomForestRegressor(n_estimators=300, random_state=42)
    model.fit(X_train, y_train)

    preds = model.predict(X_test)

    print("MAE:", mean_absolute_error(y_test, preds))
    print("R2 :", r2_score(y_test, preds))

    print("\n=== Feature Importance ===")
    for name, imp in sorted(zip(X.columns, model.feature_importances_), key=lambda x: x[1], reverse=True):
        print(f"{name:15s} {imp:.4f}")

    # Optional: convert predicted score back to classes
    pred_labels = [score_to_label(s) for s in preds[:10]]
    true_labels = [score_to_label(s) for s in y_test.iloc[:10]]

    print("\nExample predictions:")
    for t, p, ts, ps in zip(true_labels, pred_labels, y_test.iloc[:10], preds[:10]):
        print(f"true_score={ts:5.2f} pred_score={ps:5.2f} true={t:12s} pred={p}")

    joblib.dump(model, "models/rf_score_model.joblib")

    plt.figure(figsize=(6, 6))
    plt.scatter(y_test, preds, alpha=0.5)
    plt.xlabel("True activity score")
    plt.ylabel("Predicted activity score")
    plt.title("Predicted vs True Activity Score")
    plt.tight_layout()
    plt.savefig(outdir / "rf_score_scatter.png", dpi=200)
    plt.close()


if __name__ == "__main__":
    main()