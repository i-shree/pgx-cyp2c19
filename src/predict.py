import argparse
import joblib
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True)
    parser.add_argument("--rs4244285", type=int, required=True)
    parser.add_argument("--rs4986893", type=int, required=True)
    parser.add_argument("--rs12248560", type=int, required=True)
    parser.add_argument("--rs28399504", type=int, default=0)
    parser.add_argument("--rs41291556", type=int, default=0)
    parser.add_argument("--sex_encoded", type=int, default=0)
    parser.add_argument("--pop_encoded", type=int, default=0)
    args = parser.parse_args()

    model = joblib.load(args.model)

    X_new = np.array([[
        args.rs4244285,
        args.rs4986893,
        args.rs12248560,
        args.rs28399504,
        args.rs41291556,
        args.sex_encoded,
        args.pop_encoded,
    ]])

    pred = model.predict(X_new)[0]
    probs = model.predict_proba(X_new)[0]

    print(f"\nPredicted metabolizer class: {pred}")
    print("\nClass probabilities:")
    for cls, p in zip(model.classes_, probs):
        print(f"{cls:15s} {p:.4f}")


if __name__ == "__main__":
    main()