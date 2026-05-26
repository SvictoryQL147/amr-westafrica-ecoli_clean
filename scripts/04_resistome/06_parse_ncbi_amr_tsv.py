import pandas as pd
from pathlib import Path

# ---------- PATHS ----------
master_meta = Path("data/metadata/master_metadata.tsv")
gh_file = Path("data/processed/gh_isolates.tsv")
ng_file = Path("data/processed/isolates.tsv")
out_file = Path("results/tables/amrfinder_mdr_labels_independent.tsv")

# ---------- LOAD YOUR 235 ACCESSIONS ----------
meta = pd.read_csv(master_meta, sep="\t")
# Ensure the column name is correct – might be 'accession' or 'Assembly Accession'? Check.
# Assuming column is 'accession' as in your master_metadata example.
target_accessions = set(meta["accession"].astype(str))
print(f"Loaded {len(target_accessions)} target accessions from master_metadata.tsv")

# ---------- MAPPING: gene prefix -> drug class ----------
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

def extract_drug_classes(amr_str):
    """Parse AMR genotypes string and return list of unique drug classes."""
    if pd.isna(amr_str):
        return []
    classes = []
    for item in amr_str.split(","):
        if "=" not in item:
            continue
        gene = item.split("=")[0].strip()
        # Skip chromosomal point mutations (e.g., gyrA_S83L)
        if "_" in gene and any(x in gene for x in ["gyrA", "parC", "parE", "nfsA", "nfsB", "omp"]):
            continue
        for key, cls in gene_to_class.items():
            if gene.startswith(key):
                classes.append(cls)
                break
    return list(set(classes))

def parse_amr_tsv(file_path):
    df = pd.read_csv(file_path, sep="\t", dtype=str, low_memory=False)
    # Identify assembly column
    acc_col = None
    for col in df.columns:
        if "assembly" in col.lower():
            acc_col = col
            break
    if acc_col is None:
        raise KeyError(f"No assembly column in {file_path}. Columns: {df.columns.tolist()}")
    # Identify AMR genotypes column
    amr_col = None
    for col in df.columns:
        if "amr genotype" in col.lower() or "resistance" in col.lower():
            amr_col = col
            break
    if amr_col is None:
        raise KeyError(f"No AMR genotypes column in {file_path}.")
    
    # Keep only rows matching our target accessions
    df = df[df[acc_col].isin(target_accessions)]
    print(f"  {file_path.name}: kept {len(df)} rows for {df[acc_col].nunique()} unique assemblies.")
    
    # Extract drug classes per row
    df["drug_classes"] = df[amr_col].apply(extract_drug_classes)
    # Group by assembly and combine lists properly
    grouped = df.groupby(acc_col)["drug_classes"].apply(lambda x: list(set([item for sublist in x for item in sublist])))
    mdr = {acc: 1 if len(classes) >= 3 else 0 for acc, classes in grouped.items()}
    return mdr

# ---------- RUN PARSER ON BOTH FILES ----------
print("\nParsing Ghana isolates...")
gh_labels = parse_amr_tsv(gh_file)
print("\nParsing Nigeria isolates...")
ng_labels = parse_amr_tsv(ng_file)

all_labels = {**gh_labels, **ng_labels}
print(f"\nTotal assemblies with labels from TSVs: {len(all_labels)}")

# Fill any missing accessions (should be none, but just in case)
missing = target_accessions - set(all_labels.keys())
if missing:
    print(f"Warning: {len(missing)} accessions not found in TSVs. Setting to 0 (non-MDR).")
    for acc in missing:
        all_labels[acc] = 0

# Save result
df_out = pd.DataFrame([{"accession": acc, "MDR_amrfinder": label} for acc, label in all_labels.items()])
df_out.to_csv(out_file, sep="\t", index=False)
print(f"\nSaved MDR labels to {out_file}")
print(f"MDR positive: {df_out['MDR_amrfinder'].sum()} / {len(df_out)} ({df_out['MDR_amrfinder'].mean()*100:.1f}%)")