import argparse
import joblib
import numpy as np
from rules_engine import recommendation

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True, help="Path to saved model")
    parser.add_argument("--rs4244285", type=int, required=True)
    parser.add_argument("--rs12248560", type=int, required=True)
    parser.add_argument("--rs4986893", type=int, required=True)
    args = parser.parse_args()

    model = joblib.load(args.model)

    X_new = np.array([[args.rs4244285, args.rs12248560, args.rs4986893]])
    prediction = model.predict(X_new)[0]

    print(f"\nPredicted metabolizer class: {prediction}")
    print(f"Recommendation: {recommendation(prediction)}")

if __name__ == "__main__":
    main()
