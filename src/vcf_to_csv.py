import argparse
import csv
from cyvcf2 import VCF

RSIDS = [
    "rs4244285",   # *2  LOF
    "rs4986893",   # *3  LOF
    "rs12248560",  # *17 GOF
    "rs28399504",  # *4  LOF
    "rs41291556",  # *8  LOF
]

def gt_to_alt_count(gt):
    a1, a2, _phased = gt
    if a1 < 0 or a2 < 0:
        return ""  # missing
    return int(a1 >= 1) + int(a2 >= 1)

def load_panel(panel_path):
    """
    Returns a dict: sample_id -> {"sex": ..., "pop": ..., "super_pop": ...}
    Expected columns: sample  pop  super_pop  gender
    """
    panel = {}
    with open(panel_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            panel[row["sample"]] = {
                "sex":       row.get("gender", ""),
                "pop":       row.get("pop", ""),
                "super_pop": row.get("super_pop", ""),
            }
    return panel

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf",   required=True, help="Path to chr10 .vcf.gz")
    ap.add_argument("--out",   required=True, help="Output CSV path")
    ap.add_argument("--panel", required=False, default=None,
                    help="Path to 1000G panel file (TSV with sample/pop/super_pop/gender). "
                         "Download from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
                         "integrated_call_samples_v3.20130502.ALL.panel")
    args = ap.parse_args()

    # Load panel metadata if provided
    panel = load_panel(args.panel) if args.panel else {}

    v = VCF(args.vcf)
    samples = v.samples

    data = {s: {rs: "" for rs in RSIDS} for s in samples}

    found = set()
    for variant in v:
        if variant.ID in RSIDS:
            found.add(variant.ID)
            gts = variant.genotypes
            for s, gt in zip(samples, gts):
                data[s][variant.ID] = gt_to_alt_count(gt)
        if len(found) == len(RSIDS):
            break

    missing = [r for r in RSIDS if r not in found]
    if missing:
        print(f"Warning: these RSIDs were not found in the VCF: {missing}")

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)

        # Header — include metadata columns only if panel was provided
        if panel:
            w.writerow(["sample_id"] + RSIDS + ["sex", "pop", "super_pop"])
        else:
            w.writerow(["sample_id"] + RSIDS)

        for s in samples:
            row = [s] + [data[s][rs] for rs in RSIDS]
            if panel:
                meta = panel.get(s, {"sex": "", "pop": "", "super_pop": ""})
                row += [meta["sex"], meta["pop"], meta["super_pop"]]
            w.writerow(row)

    print(f"Done. Wrote {len(samples)} samples to {args.out}")
    if panel:
        matched = sum(1 for s in samples if s in panel)
        print(f"Panel metadata matched: {matched}/{len(samples)} samples")

if __name__ == "__main__":
    main()