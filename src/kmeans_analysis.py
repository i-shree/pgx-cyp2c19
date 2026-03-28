import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.cluster import KMeans
from sklearn.preprocessing import LabelEncoder

from rules_engine import classify_metabolizer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/cyp2c19_1000g.csv")
    ap.add_argument("--outdir", default="reports/figures")
    ap.add_argument("--n_clusters", type=int, default=4,
                    help="Number of clusters (default: 4 matching Poor/Intermediate/Normal/Rapid)")
    args = ap.parse_args()

    # ── Load data ────────────────────────────────────────────────────────────
    df = pd.read_csv(args.data)

    # ── Apply rules engine to get known labels ───────────────────────────────
    df["label"] = df.apply(
        lambda r: classify_metabolizer(
            int(r["rs4244285"]), int(r["rs12248560"]), int(r["rs4986893"])
        ),
        axis=1,
    )

    # ── Encode categorical columns ───────────────────────────────────────────
    le_sex = LabelEncoder()
    le_pop = LabelEncoder()
    df["sex_encoded"] = le_sex.fit_transform(df["sex"])
    df["pop_encoded"] = le_pop.fit_transform(df["super_pop"])

    print("Sex encoding:      ", dict(zip(le_sex.classes_, le_sex.transform(le_sex.classes_))))
    print("Population encoding:", dict(zip(le_pop.classes_, le_pop.transform(le_pop.classes_))))

    # ── Feature matrix — SNPs only (no sex/population) ──────────────────────
    feature_cols = ["rs4244285", "rs4986893", "rs12248560",
                    "rs28399504", "rs41291556"]
    X = df[feature_cols].to_numpy()

    # ── KMeans clustering ────────────────────────────────────────────────────
    kmeans = KMeans(n_clusters=args.n_clusters, random_state=42, n_init=10)
    df["cluster"] = kmeans.fit_predict(X)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── OUTPUT 1: Cluster vs label comparison table ──────────────────────────
    print("\n=== Cluster vs Known Label (Rules Engine) ===")
    cluster_label = df.groupby(["cluster", "label"]).size().unstack(fill_value=0)
    print(cluster_label)

    # Map each cluster to its dominant label
    cluster_map = cluster_label.idxmax(axis=1).to_dict()
    print("\nDominant label per cluster:")
    for cluster, label in cluster_map.items():
        print(f"  Cluster {cluster} → {label}")

    # ── OUTPUT 2: Confusion-style heatmap of clusters vs labels ─────────────
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(
        cluster_label,
        annot=True,
        fmt="d",
        cmap="YlOrRd",
        ax=ax,
        linewidths=0.5,
    )
    ax.set_title("KMeans Clusters vs Known Metabolizer Labels", fontsize=13)
    ax.set_xlabel("Known Label (Rules Engine)")
    ax.set_ylabel("KMeans Cluster")
    plt.tight_layout()
    plt.savefig(outdir / "kmeans_cluster_vs_label.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outdir}/kmeans_cluster_vs_label.png")

    # ── OUTPUT 3: Population breakdown per cluster ───────────────────────────
    print("\n=== Population Breakdown per Cluster ===")
    pop_cluster = df.groupby(["cluster", "super_pop"]).size().unstack(fill_value=0)
    print(pop_cluster)

    fig, ax = plt.subplots(figsize=(9, 5))
    pop_cluster.plot(kind="bar", ax=ax, colormap="tab10", edgecolor="white")
    ax.set_title("Population Distribution per KMeans Cluster", fontsize=13)
    ax.set_xlabel("KMeans Cluster")
    ax.set_ylabel("Number of Individuals")
    ax.set_xticklabels([f"Cluster {i}\n({cluster_map.get(i, '?')})"
                        for i in pop_cluster.index], rotation=0)
    ax.legend(title="Population", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(outdir / "kmeans_population_per_cluster.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outdir}/kmeans_population_per_cluster.png")

    # ── OUTPUT 4: Sex breakdown per cluster ─────────────────────────────────
    print("\n=== Sex Breakdown per Cluster ===")
    sex_cluster = df.groupby(["cluster", "sex"]).size().unstack(fill_value=0)
    print(sex_cluster)

    fig, ax = plt.subplots(figsize=(8, 5))
    sex_cluster.plot(kind="bar", ax=ax, colormap="Set2", edgecolor="white")
    ax.set_title("Sex Distribution per KMeans Cluster", fontsize=13)
    ax.set_xlabel("KMeans Cluster")
    ax.set_ylabel("Number of Individuals")
    ax.set_xticklabels([f"Cluster {i}\n({cluster_map.get(i, '?')})"
                        for i in sex_cluster.index], rotation=0)
    ax.legend(title="Sex", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(outdir / "kmeans_sex_per_cluster.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outdir}/kmeans_sex_per_cluster.png")

    # ── SUMMARY ──────────────────────────────────────────────────────────────
    print("\n=== Summary ===")
    print(f"Total samples:   {len(df)}")
    print(f"Number of clusters: {args.n_clusters}")
    print(f"Features used:   {feature_cols}")
    print("\nCluster sizes:")
    print(df["cluster"].value_counts().sort_index())


if __name__ == "__main__":
    main()