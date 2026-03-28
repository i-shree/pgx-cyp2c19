
# pgx-cyp2c19

> Educational use only. Not for clinical decision making.

This project demonstrates a complete pharmacogenomics workflow that predicts clopidogrel metabolizer status from CYP2C19 genetic variants using real 1000 Genomes population data.

It integrates whole-genome VCF processing, targeted SNP extraction, genotype encoding, rule-based phenotype classification, machine learning validation, and CLI-based patient prediction — all in one reproducible pipeline.

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
    pip install pandas numpy scikit-learn matplotlib joblib cyvcf2

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
- chr10.vcf.gz ≈ 1.5 GB
- chr10.vcf.gz.tbi ≈ few MB
- panel.txt
---

### 5) Extract CYP2C19 SNPs

    python3 src/vcf_to_csv.py \
      --vcf data/raw/chr10.vcf.gz \
      --out data/cyp2c19_1000g.csv

This generates:
- 2504 individuals
- 3 SNP features:
    rs4244285 (*2, loss-of-function)
    rs4986893 (*3, loss-of-function)
    rs12248560 (*17, gain-of-function)

Genotypes encoded as:
0 = no alternate allele  
1 = heterozygous  
2 = homozygous alternate  

---

### 6) Train Machine Learning Models

    python3 src/train_and_evaluate.py \
      --data data/cyp2c19_1000g.csv

Training details:
- Stratified 70/30 train-test split
- Logistic Regression
- Random Forest

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

High accuracy occurs because metabolizer labels are generated deterministically from the same SNP features used for training.

---

### 7) Interpret Confusion Matrix

- Rows = True metabolizer class  
- Columns = Predicted class  
- Diagonal cells = Correct predictions  
- Off-diagonal cells = Model confusion  
- Color bar = Number of samples per cell  

Since genotype → phenotype mapping is rule-based and noise-free, predictions are nearly perfect.

---

### 8) Predict New Patient (CLI)

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
- rs4244285 (*2) → reduced enzyme function
- rs4986893 (*3) → reduced enzyme function
- rs12248560 (*17) → increased enzyme function

Individuals are classified as:
- Poor metabolizers
- Intermediate metabolizers
- Normal metabolizers
- Rapid metabolizers

This project demonstrates how genomic variation can be translated into personalized drug response predictions.

---

## Conclusion

Whole-genome sequencing data  
→ Pharmacogenetic variant extraction  
→ Metabolizer classification  
→ Machine learning validation  
→ Clinical decision support  

This repository demonstrates a complete precision medicine workflow suitable for educational research purposes.
