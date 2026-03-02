import argparse
import csv
from cyvcf2 import VCF

RSIDS = ["rs4244285", "rs12248560", "rs4986893"]  # *2, *17, *3

def gt_to_alt_count(gt):
    # cyvcf2 returns (a1, a2, phased) where a1/a2 in {0,1,2,...} or -1 for missing
    a1, a2, _phased = gt
    if a1 < 0 or a2 < 0:
        return ""  # missing
    # count how many alleles are non-reference (>=1)
    return int(a1 >= 1) + int(a2 >= 1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="Path to chr10 .vcf.gz")
    ap.add_argument("--out", required=True, help="Output CSV path")
    args = ap.parse_args()

    v = VCF(args.vcf)
    samples = v.samples  # list of sample IDs in header

    # We'll store per-RSID alt counts in a dict of dicts: data[sample][rsid] = 0/1/2
    data = {s: {rs: "" for rs in RSIDS} for s in samples}

    found = set()
    for variant in v:
        if variant.ID in RSIDS:
            found.add(variant.ID)
            gts = variant.genotypes  # list aligned with samples
            for s, gt in zip(samples, gts):
                data[s][variant.ID] = gt_to_alt_count(gt)

        if len(found) == len(RSIDS):
            break

    missing = [r for r in RSIDS if r not in found]
    if missing:
        raise RuntimeError(f"Did not find these RSIDs in VCF: {missing}")

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id"] + RSIDS)
        for s in samples:
            w.writerow([s] + [data[s][rs] for rs in RSIDS])

if __name__ == "__main__":
    main()
