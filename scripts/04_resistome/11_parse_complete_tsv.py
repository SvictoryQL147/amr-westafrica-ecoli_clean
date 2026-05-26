import pandas as pd
import json
from pathlib import Path

# Paths
gh_file = Path("data/processed/gh_complete_isolates.tsv")
ng_file = Path("data/processed/nig_complete_isolates.tsv")
meta_path = Path("data/metadata/master_metadata.tsv")
gca_map_path = Path("data/metadata/gca_to_gcf.json")
out_path = Path("results/tables/amrfinder_mdr_labels_independent.tsv")

# Load target GCFs
meta = pd.read_csv(meta_path, sep="\t")
target_gcf = set(meta["accession"].astype(str))
print(f"Target GCF accessions: {len(target_gcf)}")

# Load GCA->GCF mapping
with open(gca_map_path) as f:
    gca_to_gcf = json.load(f)
print(f"Loaded {len(gca_to_gcf)} GCA→GCF mappings")

# Gene to class mapping (same as before)
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
        parts = item.split("=", 1)
        if len(parts) != 2:
            continue
        gene, status = parts[0].strip(), parts[1].strip().upper()
        if status != "COMPLETE":
            continue
        # Skip chromosomal mutations
        if "_" in gene and any(x in gene for x in ["gyrA", "parC", "parE", "nfsA", "nfsB", "omp"]):
            continue
        for key, cls in gene_to_class.items():
            if gene.startswith(key):
                classes.add(cls)
                break
    return list(classes)

def process_tsv(file_path, country_name):
    print(f"\nProcessing {country_name} file: {file_path}")
    df = pd.read_csv(file_path, sep="\t", dtype=str, low_memory=False)
    # Find assembly column
    acc_col = None
    for col in df.columns:
        if "assembly" in col.lower():
            acc_col = col
            break
    if acc_col is None:
        raise KeyError(f"No assembly column in {file_path}. Columns: {df.columns.tolist()}")
    # Find AMR genotypes column
    amr_col = None
    for col in df.columns:
        if "amr genotype" in col.lower() or col == "AMR genotypes":
            amr_col = col
            break
    if amr_col is None:
        raise KeyError(f"No AMR column in {file_path}")
    print(f"  Using assembly column: '{acc_col}'")
    print(f"  Using AMR column: '{amr_col}'")
    # Map GCA to GCF
    df["assembly_gcf"] = df[acc_col].map(gca_to_gcf)
    # Keep only target GCFs
    df = df[df["assembly_gcf"].isin(target_gcf)]
    print(f"  Kept {len(df)} rows for {df['assembly_gcf'].nunique()} assemblies")
    # Extract classes
    df["drug_classes"] = df[amr_col].apply(extract_classes)
    return df

# Process both files
gh_df = process_tsv(gh_file, "Ghana")
ng_df = process_tsv(ng_file, "Nigeria")
combined = pd.concat([gh_df, ng_df], ignore_index=True)

# Group by assembly and count unique drug classes
grouped = combined.groupby("assembly_gcf")["drug_classes"].agg(
    lambda x: list(set(cls for sublist in x for cls in sublist))
)
mdr_labels = {acc: 1 if len(classes) >= 3 else 0 for acc, classes in grouped.items()}

# Ensure all targets are present
result = []
for acc in target_gcf:
    result.append({"accession": acc, "MDR_amrfinder": mdr_labels.get(acc, 0)})

out_df = pd.DataFrame(result)
out_df.to_csv(out_path, sep="\t", index=False)
print(f"\n✅ Saved MDR labels to {out_path}")
print(f"MDR positive: {out_df['MDR_amrfinder'].sum()} / {len(out_df)} ({out_df['MDR_amrfinder'].mean()*100:.1f}%)")