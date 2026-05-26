import json
import time
import requests
import pandas as pd
from pathlib import Path

# Paths
gh_tsv = Path("data/processed/gh_complete_isolates.tsv")
ng_tsv = Path("data/processed/nig_complete_isolates.tsv")
out_json = Path("data/metadata/gca_to_gcf.json")

# Read both TSVs to get all GCA accessions
def get_gca_list(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    # Find assembly column
    for col in df.columns:
        if "assembly" in col.lower():
            return df[col].dropna().unique().tolist()
    return []

gca_list = set()
gca_list.update(get_gca_list(gh_tsv))
gca_list.update(get_gca_list(ng_tsv))
gca_list = list(gca_list)
print(f"Found {len(gca_list)} unique GCA accessions")

# Function to get GCF from GCA
def gca_to_gcf(gca):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{gca}/dataset_report"
    try:
        r = requests.get(url, timeout=20)
        if r.status_code == 200:
            data = r.json()
            for report in data.get("reports", []):
                gcf = report.get("assembly_info", {}).get("assembly_accession")
                if gcf:
                    return gcf
    except:
        pass
    return None

# Build mapping
mapping = {}
for i, gca in enumerate(gca_list):
    if i % 20 == 0:
        print(f"Processing {i}/{len(gca_list)}...")
    gcf = gca_to_gcf(gca)
    if gcf:
        mapping[gca] = gcf
    time.sleep(0.34)  # rate limit

print(f"Obtained {len(mapping)} mappings")
# Save
with open(out_json, "w") as f:
    json.dump(mapping, f, indent=2)
print(f"Saved to {out_json}")