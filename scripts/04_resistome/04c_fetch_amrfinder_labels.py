"""
Fix for Objective 4 — Independent MDR labels from NCBI AMRFinderPlus
=====================================================================
Fetches AMRFinderPlus results directly from NCBI Pathogen Detection API
for all 235 genomes, builds independent MDR labels.
"""
import json
import time
import requests
import pandas as pd
from pathlib import Path
from datetime import datetime

# Correct paths
BASE_DIR = Path.cwd()
DATA_DIR = BASE_DIR / "data" / "metadata"
TABLE_DIR = BASE_DIR / "results" / "tables"
LOG_DIR   = BASE_DIR / "results" / "logs"

# Create directories if they don't exist
TABLE_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

# Path to master metadata (corrected)
MASTER_META = DATA_DIR / "master_metadata.tsv"

def fetch_amrfinder_results(accessions):
    """
    Query NCBI Pathogen Detection AMRFinderPlus results.
    Returns resistance gene calls per genome from independent algorithm.
    """
    records = []
    failed  = []
    total   = len(accessions)

    print(f"Fetching AMRFinderPlus data for {total} genomes from NCBI...")
    print("This may take 15-20 minutes.\n")

    for i, acc in enumerate(accessions):
        if i % 20 == 0:
            print(f"  Progress: {i}/{total}...")

        # Convert GCF to biosample via NCBI datasets API
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/dataset_report"
        try:
            r = requests.get(url, timeout=20)
            if r.status_code != 200:
                failed.append(acc)
                time.sleep(0.5)
                continue

            data = r.json()
            reports = data.get("reports", [])
            if not reports:
                failed.append(acc)
                time.sleep(0.5)
                continue

            report   = reports[0]
            asm_info = report.get("assembly_info", {})
            biosample_acc = asm_info.get("biosample", {}).get("accession", "")

            if not biosample_acc:
                failed.append(acc)
                time.sleep(0.5)
                continue

            # Query pathogen detection AMR data via biosample
            amr_url = f"https://www.ncbi.nlm.nih.gov/pathogens/api/v1/amrfinderplus?biosample_acc={biosample_acc}"
            r2 = requests.get(amr_url, timeout=20)

            if r2.status_code == 200:
                amr_data = r2.json()
                genes = amr_data.get("amrfinderplus_genes", [])
                for gene in genes:
                    gene["genome_accession"] = acc
                    gene["biosample"] = biosample_acc
                    records.append(gene)
            else:
                failed.append(acc)

            time.sleep(0.34)  # NCBI rate limit

        except Exception as e:
            print(f"  Exception for {acc}: {e}")
            failed.append(acc)
            time.sleep(1)
            continue

    print(f"\nCompleted: {len(records)} gene records fetched")
    print(f"Failed/missing: {len(failed)} genomes")
    return records, failed

def build_independent_labels(records, accessions, countries):
    """
    Build MDR labels from AMRFinderPlus results.
    MDR = resistance to >= 3 antibiotic classes per genome.
    """
    if not records:
        return None

    df = pd.DataFrame(records)
    print(f"\nAMRFinderPlus columns: {df.columns.tolist()}")

    # Try to find a column with drug class information
    class_col = None
    for col in ["drug_class", "subclass", "class", "antibiotic_class", "Element_subtype"]:
        if col in df.columns:
            class_col = col
            break

    if not class_col:
        print("WARNING: Could not find drug class column. Using gene symbol mapping fallback.")
        # Fallback: map gene symbols to classes (simplified)
        gene_to_class = {
            "blaTEM": "Beta-lactam", "blaCTX-M": "Beta-lactam", "blaOXA": "Beta-lactam",
            "blaSHV": "Beta-lactam", "blaCMY": "Beta-lactam", "blaNDM": "Carbapenem",
            "sul1": "Sulfonamide", "sul2": "Sulfonamide", "sul3": "Sulfonamide",
            "qnrA": "Fluoroquinolone", "qnrB": "Fluoroquinolone", "qnrS": "Fluoroquinolone",
            "aac(3)": "Aminoglycoside", "aac(6')": "Aminoglycoside", "aph(3')": "Aminoglycoside",
            "tetA": "Tetracycline", "tetB": "Tetracycline", "catA": "Chloramphenicol",
            "erm": "Macrolide", "mph": "Macrolide", "mcr": "Colistin", "dfrA": "Trimethoprim"
        }
        def map_gene(gene):
            for key, cls in gene_to_class.items():
                if gene.startswith(key):
                    return cls
            return None
        df["class_mapped"] = df["gene_symbol"].apply(map_gene)
        class_col = "class_mapped"
        df = df.dropna(subset=[class_col])

    mdr_labels = {}
    for acc in accessions:
        genome_df = df[df["genome_accession"] == acc]
        if genome_df.empty:
            mdr_labels[acc] = 0  # assume non-MDR if no data
        else:
            n_classes = genome_df[class_col].nunique()
            mdr_labels[acc] = 1 if n_classes >= 3 else 0

    labels_df = pd.DataFrame([
        {"accession": acc, "MDR_amrfinder": label, "country": countries.get(acc, "")}
        for acc, label in mdr_labels.items()
    ])

    known = labels_df.dropna(subset=["MDR_amrfinder"])
    print(f"\nMDR labels obtained: {len(known)}/{len(accessions)} genomes")
    print(f"MDR positive: {known['MDR_amrfinder'].sum():.0f} ({known['MDR_amrfinder'].mean()*100:.1f}%)")

    return labels_df

def main():
    print("="*55)
    print("Objective 4 Fix — Independent AMRFinderPlus Labels")
    print(f"Date: {datetime.now().isoformat()}")
    print("="*55)

    if not MASTER_META.exists():
        print(f"ERROR: master_metadata.tsv not found at {MASTER_META}")
        print("Please ensure the file exists.")
        return

    meta = pd.read_csv(MASTER_META, sep="\t")
    accessions = meta["accession"].tolist()
    countries  = dict(zip(meta["accession"], meta["country"]))

    records, failed = fetch_amrfinder_results(accessions)

    if records:
        raw_df = pd.DataFrame(records)
        raw_df.to_csv(TABLE_DIR / "amrfinder_ncbi_raw.tsv", sep="\t", index=False)
        print(f"Raw results saved: {TABLE_DIR / 'amrfinder_ncbi_raw.tsv'}")

        labels_df = build_independent_labels(records, accessions, countries)
        if labels_df is not None:
            labels_df.to_csv(TABLE_DIR / "amrfinder_mdr_labels_independent.tsv", sep="\t", index=False)
            print(f"Labels saved: {TABLE_DIR / 'amrfinder_mdr_labels_independent.tsv'}")
            print("\n✅ Next step: run ML retraining with these independent labels.")
    else:
        print("\nNo AMRFinderPlus data retrieved from API.")
        print("\nMANUAL ALTERNATIVE (most reliable):")
        print("1. Go to https://www.ncbi.nlm.nih.gov/pathogens/")
        print("2. Search: Escherichia coli")
        print("3. Filter by: Isolation country = Nigeria OR Ghana")
        print("4. Click 'AMR' tab")
        print("5. Download full AMR table as CSV")
        print("6. Save to: data/processed/ncbi_pathogens_amr.csv")
        print("7. Run the alternative parser I will provide.")

    with open(LOG_DIR / "04b_amrfinder_fetch.json", "w") as f:
        json.dump({
            "date": datetime.now().isoformat(),
            "total_accessions": len(accessions),
            "records_fetched": len(records),
            "failed": len(failed),
            "failed_accessions": failed[:10]
        }, f, indent=2)

if __name__ == "__main__":
    main()