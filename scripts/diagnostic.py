import pandas as pd

# Load master metadata
meta = pd.read_csv("data/metadata/master_metadata.tsv", sep="\t")
print(f"Master metadata columns: {meta.columns.tolist()}")
print(f"Number of rows: {len(meta)}")
print(f"First 3 accessions:\n{meta['accession'].head(3).tolist()}")
print(f"Accession type: {type(meta['accession'].iloc[0])}")

# Load Ghana TSV
print("\n--- Ghana TSV ---")
df_gh = pd.read_csv("data/processed/gh_isolates.tsv", sep="\t", nrows=3, dtype=str)
print(f"Columns: {df_gh.columns.tolist()}")
# Find assembly column
acc_col = None
for col in df_gh.columns:
    if "assembly" in col.lower():
        acc_col = col
        break
if acc_col:
    print(f"Assembly column name: '{acc_col}'")
    print(f"First 3 values: {df_gh[acc_col].head(3).tolist()}")
else:
    print("No assembly column found!")

# Load Nigeria TSV
print("\n--- Nigeria TSV ---")
df_ng = pd.read_csv("data/processed/isolates.tsv", sep="\t", nrows=3, dtype=str)
print(f"Columns: {df_ng.columns.tolist()}")
for col in df_ng.columns:
    if "assembly" in col.lower():
        acc_col = col
        break
if acc_col:
    print(f"Assembly column name: '{acc_col}'")
    print(f"First 3 values: {df_ng[acc_col].head(3).tolist()}")
else:
    print("No assembly column found!")