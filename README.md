# pgx-cyp2c19


This project demonstrates a complete pharmacogenomics workflow that predicts clopidogrel metabolizer status from CYP2C19 genetic variants using real 1000 Genomes population data.

It integrates whole-genome VCF processing, targeted SNP extraction, genotype encoding, rule-based phenotype classification, machine learning validation, unsupervised clustering analysis, and CLI-based patient prediction — all in one reproducible pipeline.

---

## Complete Start-to-Finish Setup

### 1) Install Python

Download Python 3 from:
https://www.python.org/downloads/

Verify installation:
    python3 --version

---

### 2) Create Project and Environment

    mkdir pgx-cyp2c19
    cd pgx-cyp2c19

    python3 -m venv .venv
    source .venv/bin/activate

    pip install --upgrade pip
    pip install pandas numpy scikit-learn matplotlib seaborn joblib cyvcf2

---

### 3) Create Project Structure

    mkdir -p data/raw
    mkdir -p src
    mkdir -p models
    mkdir -p reports/figures

---

### 4) Download Real Genomic Data (Chromosome 10 – 1000 Genomes Phase 3)

    curl -L -o data/raw/chr10.vcf.gz \
    https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

    curl -L -o data/raw/chr10.vcf.gz.tbi \
    https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi

    curl -L -o data/raw/panel.txt \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

File sizes:
- chr10.vcf.gz ≈ 778 MB compressed
- chr10.vcf.gz.tbi ≈ few MB
- panel.txt ≈ few KB

---

### 5) Extract CYP2C19 SNPs

    python3 src/vcf_to_csv.py \
      --vcf data/raw/chr10.vcf.gz \
      --out data/cyp2c19_1000g_5snp.csv \
      --panel data/raw/panel.txt

This generates:
- 2504 individuals
- 5 SNP features:
    rs4244285  (*2,  loss-of-function)
    rs4986893  (*3,  loss-of-function)
    rs12248560 (*17, gain-of-function)
    rs28399504 (*4,  loss-of-function — rare, 4 samples)
    rs41291556 (*8,  loss-of-function — rare, 5 samples)
- sex (male/female)
- pop (26 specific populations e.g. GBR, YRI)
- super_pop (5 continental groups: AFR, AMR, EAS, EUR, SAS)

Genotypes encoded as:
0 = no alternate allele
1 = heterozygous
2 = homozygous alternate

---

### 6) K-Means Clustering Analysis (Dataset Exploration)

    python3 src/kmeans_analysis.py \
      --data data/cyp2c19_1000g_5snp.csv \
      --outdir reports/figures

K-Means is an unsupervised algorithm — no labels are provided.
It groups individuals based only on SNP pattern similarity.
After clustering, results are compared against known metabolizer labels.

Key findings:

Cluster vs Known Label:

| Cluster | Dominant Label | Purity  |
|---------|---------------|---------|
| 0       | Intermediate  | 97%     |
| 1       | Normal        | 96%     |
| 2       | Rapid         | 99%     |
| 3       | Poor          | 100%    |

Without being told any labels, K-Means independently rediscovered
the same four metabolizer classes established by clinical research.

Population finding:
- SAS (South Asian): 16.4% Poor metabolizers — highest of all populations
- EAS (East Asian):  11.7% Poor metabolizers
- EUR (European):     1.2% Poor metabolizers
- AFR (African):   highest Rapid metabolizer rate

Sex finding:
- Male/female distribution is nearly identical across all clusters
- Sex has no influence on CYP2C19 metabolizer status
- Confirmed by chromosome 10 being autosomal (not a sex chromosome)

Outputs saved to reports/figures/:
- kmeans_cluster_vs_label.png     — heatmap of clusters vs known labels
- kmeans_population_per_cluster.png — population distribution per cluster
- kmeans_sex_per_cluster.png      — sex distribution per cluster

---

### 7) Train Machine Learning Models

    python3 src/train_and_evaluate.py \
      --data data/cyp2c19_1000g_5snp.csv

Training details:
- Rules engine (CPIC guidelines) assigns metabolizer labels first
- Stratified 70/30 train-test split
- Logistic Regression
- Random Forest (200 trees)

Dataset size: 2504 samples
Test set: 752 samples

Example test distribution:

| Class        | Count |
|--------------|-------|
| Intermediate | 259   |
| Normal       | 278   |
| Poor         | 46    |
| Rapid        | 169   |

Accuracy ≈ 100%

High accuracy occurs because metabolizer labels are generated
deterministically from the same SNP features used for training.
The model learned clinically established CPIC rules — not new biology.

---

### 8) Interpret Confusion Matrix

- Rows = True metabolizer class
- Columns = Predicted class
- Diagonal cells = Correct predictions
- Off-diagonal cells = Model confusion
- Color bar = Number of samples per cell

Since genotype → phenotype mapping is rule-based and noise-free,
predictions are nearly perfect.

---

### 9) Predict New Patient (CLI)

    python3 src/predict.py \
      --model models/rf_model.joblib \
      --rs4244285 1 \
      --rs12248560 0 \
      --rs4986893 0

Example output:

    Predicted metabolizer class: Intermediate
    Recommendation: Consider alternative antiplatelet therapy.

---

## Scientific Background

Clopidogrel is a prodrug that requires activation by the CYP2C19 enzyme.

Key variants:
- rs4244285  (*2)  → reduced enzyme function
- rs4986893  (*3)  → reduced enzyme function
- rs12248560 (*17) → increased enzyme function
- rs28399504 (*4)  → reduced enzyme function (rare — 0.16% frequency)
- rs41291556 (*8)  → reduced enzyme function (rare — 0.20% frequency)

Individuals are classified as:
- Poor metabolizers        — two loss-of-function alleles
- Intermediate metabolizers — one loss-of-function allele
- Normal metabolizers      — no loss-of-function, no gain-of-function
- Rapid metabolizers       — one or more gain-of-function alleles

This project demonstrates how genomic variation can be translated
into personalized drug response predictions.

---

## Two complementary ML approaches

| Approach       | Algorithm     | Labels needed | Purpose                        |
|----------------|---------------|---------------|--------------------------------|
| Supervised     | Random Forest | Yes (CPIC)    | Predict new patient's class    |
| Unsupervised   | K-Means       | No labels     | Explore natural patterns in data |

K-Means confirmed that SNP patterns alone are strong enough for an
algorithm with no prior knowledge to group patients the same way
doctors do — independently validating the clinical classification system.

---

## Conclusion

Whole-genome sequencing data
→ Pharmacogenetic variant extraction (5 SNPs + sex + population)
→ Dataset exploration with K-Means clustering
→ Metabolizer classification via CPIC rules engine
→ Machine learning validation with Random Forest
→ Clinical decision support via command-line predictor

This repository demonstrates a complete precision medicine workflow
suitable for educational research purposes.