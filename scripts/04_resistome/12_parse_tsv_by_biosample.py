import pandas as pd
from pathlib import Path

# Load master metadata (has GCF and BioSample)
master = pd.read_csv("data/metadata/master_metadata.tsv", sep="\t")
# Ensure BioSample column exists and is string
if "biosample" not in master.columns:
    raise KeyError("master_metadata.tsv missing 'biosample' column")
master = master[["accession", "biosample"]].dropna()
master["biosample"] = master["biosample"].astype(str)
print(f"Master: {len(master)} entries")

# Load TSV files
gh = pd.read_csv("data/processed/gh_complete_isolates.tsv", sep="\t", dtype=str, low_memory=False)
ng = pd.read_csv("data/processed/nig_complete_isolates.tsv", sep="\t", dtype=str, low_memory=False)
df = pd.concat([gh, ng], ignore_index=True)
print(f"TSV combined: {len(df)} rows")

# Identify BioSample column in TSV
biosample_col = None
for col in df.columns:
    if "biosample" in col.lower():
        biosample_col = col
        break
if biosample_col is None:
    raise KeyError("No BioSample column in TSV files")

# Identify AMR genotypes column
amr_col = None
for col in df.columns:
    if "amr genotype" in col.lower():
        amr_col = col
        break
if amr_col is None:
    raise KeyError("No AMR genotypes column in TSV files")

print(f"BioSample column: '{biosample_col}'")
print(f"AMR column: '{amr_col}'")

# Merge on BioSample
df = df.merge(master, left_on=biosample_col, right_on="biosample", how="inner")
print(f"After merge: {len(df)} rows for {df['accession'].nunique()} unique GCF accessions")

# Gene to class mapping
gene_to_class = {
    "blaTEM": "Beta-lactam", "blaCTX-M": "Beta-lactam", "blaOXA": "Beta-lactam",
    "blaSHV": "Beta-lactam", "blaCMY": "Beta-lactam", "blaNDM": "Carbapenem",
    "sul1": "Sulfonamide", "sul2": "Sulfonamide", "sul3": "Sulfonamide",
    "qnrA": "Fluoroquinolone", "qnrB": "Fluoroquinolone", "qnrS": "Fluoroquinolone",
    "aac(3)": "Aminoglycoside", "aac(6')": "Aminoglycoside", "aph(3')": "Aminoglycoside",
    "ant(2'')": "Aminoglycoside", "tetA": "Tetracycline", "tetB": "Tetracycline",
    "catA": "Chloramphenicol", "cmlA": "Chloramphenicol", "floR": "Chloramphenicol",
    "erm": "Macrolide", "mph": "Macrolide", "mef": "Macrolide",
    "mcr": "Colistin", "armA": "Aminoglycoside", "rmt": "Aminoglycoside",
    "dfrA": "Trimethoprim", "arr": "Rifampicin"
}

def extract_classes(amr_str):
    if pd.isna(amr_str):
        return []
    classes = set()
    for item in amr_str.split(","):
        if "=" not in item:
            continue
        gene, status = item.split("=", 1)
        gene = gene.strip()
        status = status.strip().upper()
        if status != "COMPLETE":
            continue
        if "_" in gene and any(x in gene for x in ["gyrA", "parC", "parE", "nfsA", "nfsB", "omp"]):
            continue
        for key, cls in gene_to_class.items():
            if gene.startswith(key):
                classes.add(cls)
                break
    return list(classes)

df["drug_classes"] = df[amr_col].apply(extract_classes)

# Group by accession (GCF)
grouped = df.groupby("accession")["drug_classes"].agg(lambda x: list(set(cls for sublist in x for cls in sublist)))
mdr_labels = {acc: 1 if len(classes) >= 3 else 0 for acc, classes in grouped.items()}

# Ensure all master accessions are present
all_acc = master["accession"].unique()
result = []
for acc in all_acc:
    result.append({"accession": acc, "MDR_amrfinder": mdr_labels.get(acc, 0)})

out_df = pd.DataFrame(result)
out_path = Path("results/tables/amrfinder_mdr_labels_independent.tsv")
out_path.parent.mkdir(parents=True, exist_ok=True)
out_df.to_csv(out_path, sep="\t", index=False)

print(f"\n✅ Saved to {out_path}")
print(f"MDR positive: {out_df['MDR_amrfinder'].sum()} / {len(out_df)} ({out_df['MDR_amrfinder'].mean()*100:.1f}%)")