import json
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE_DIR  = Path(".")
META_DIR  = BASE_DIR / "data" / "metadata"
TABLE_DIR = BASE_DIR / "results" / "tables"
LOG_DIR   = BASE_DIR / "results" / "logs"

# Shorter 20-mer markers — more tolerant of strain variation
resistance_markers = {
    "blaTEM":     ["ATGAGTATTCAACATTTCCG", "TTACCAATGCTTAATCAGTG", "CAATGCTTAATCAGTGAGGC"],
    "blaCTX-M":   ["ATGTGCAGCACCAGTAAAGT", "TTTGCCGTCAGGTCGCGAAA", "GCACGATGACATTCACGCCG"],
    "blaSHV":     ["ATGCGTTATATTCGCCTGTG", "TTCCCGATGCCGCCGCCAGT", "TTATCGGAGGTATAATGAGC"],
    "blaOXA-1":   ["ATGAAAATAGAGCAATTTTT", "TTAAAATCATTTGCGGAGCC", "AAAGCTTGTCATACTTGTGG"],
    "tetA":       ["ATGAAATCTAACAATGCAGA", "CCCAATTTGTCAGTGGATAA", "GTATCATTGGGTTTGTCACG"],
    "tetB":       ["ATGAAATGGAATAAAATTCT", "CCGCAGAAGCGCTATGAAGT", "GATCTTGCCATTGGTGTAGT"],
    "sul1":       ["ATGGTGACGGTGTTCGGCAT", "CACCGGAAACATCGCTGCAC", "TTACCGAGCGTTTTGCGGAC"],
    "sul2":       ["ATGAATTTCTGGTTTTTGAT", "TTACCGTTTGCGGTCGTGGC", "CGGCATCGTCAACATAACC"],
    "dfrA1":      ["ATGAAATTAACCTTAGTTGT", "TTAATCGCTAATTTCGGCTG", "GCGAAAGCGGAGATCCTGAA"],
    "aac(3)-II":  ["ATGACTCTTAAAAAATCAGG", "TCAGAAAAACTCACCGAGGC", "GAAATCGTCTCCGACGAGAT"],
    "aac(6')-Ib": ["ATGACTAACGATCAATTGCT", "TTACCAGCACTGCGGAATCG", "AAGACCTCGAATCGGCATTG"],
    "qnrB":       ["ATGGATATTATTGATAAAAG", "TTAATGGCGAAATTTGCCCC", "CTGCGCAACATGGACGATGT"],
    "qnrS":       ["ATGGATATTATTGATGAAAG", "TTAGTCAGGATAAACAACAA", "CTTCATGGAGCAAGATTTCC"],
    "mcr-1":      ["ATGCGGTTGATGACGGTGAC", "TTAGACATCATGCCGGTGCG", "GCAGTGGCCGATGAAATCGG"],
    "aph(3')-Ia": ["ATGAGCCATATTCAACGGGA", "TCAGAAGAACTCGTCAAGAA", "GCGAATTTCTTATGATTTAT"],
    "rmtB":       ["ATGCACAAGCAATTTATTGA", "TCAGTCTTCGAATTTTTCAG", "AAGCGGCTTCGGCTTTTGAT"],
    "armA":       ["ATGGATAAGAAATTTCAAAG", "TCAGAACATCGGCAAGACGT", "TTCAAAGCTTAACAGAGCAA"],
    "blaCMY-2":   ["ATGATGAAAAAATCGTTATG", "TTAGCCGTTGCCGTCGATGT", "GAAGTGCTTGCCGTTGGAGT"],
    "blaKPC-2":   ["ATGTCACTGTATCGCCGTCT", "TTACTGCCCGTTGACGCCCA", "GAAGATCAAGTCAGTTTGTC"],
    "blaNDM-1":   ["ATGCGCTTCATTCACGCACT", "TCAGCGCTTGGCTGGGTCAT", "CACCACCACCACCACCACCA"],
}

class_map = {
    "beta_lactam":    ["blaTEM","blaCTX-M","blaSHV","blaOXA-1","blaCMY-2","blaKPC-2","blaNDM-1"],
    "tetracycline":   ["tetA","tetB"],
    "sulfonamide":    ["sul1","sul2"],
    "trimethoprim":   ["dfrA1"],
    "aminoglycoside": ["aac(3)-II","aac(6')-Ib","aph(3')-Ia","rmtB","armA"],
    "fluoroquinolone":["qnrB","qnrS"],
    "colistin":       ["mcr-1"],
}

def rc(seq):
    return seq.translate(str.maketrans("ATCG","TAGC"))[::-1]

def load_accessions():
    df = pd.read_csv(META_DIR / "master_metadata.tsv", sep="\t")
    return dict(zip(df["accession"], df["country"]))

def scan_genomes(countries):
    fna_files = (
        list(Path("data/raw/nigeria/ncbi_dataset/data").rglob("*.fna")) +
        list(Path("data/raw/ghana/ncbi_dataset/data").rglob("*.fna"))
    )
    print(f"Scanning {len(fna_files)} genomes for {len(resistance_markers)} resistance genes...")
    results = []
    for i, fna_path in enumerate(fna_files):
        if i % 25 == 0:
            print(f"  {i}/{len(fna_files)} genomes scanned...")
        acc = fna_path.parent.name
        country = countries.get(acc, countries.get(
            acc.replace("GCA_","GCF_").replace("GCF_","GCA_"), "Unknown"))
        with open(fna_path) as f:
            lines = f.readlines()
        seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
        row = {"accession": acc, "country": country}
        for gene, markers in resistance_markers.items():
            found = any(m in seq or rc(m) in seq for m in markers)
            row[gene] = 1 if found else 0
        results.append(row)
    return pd.DataFrame(results)

def main():
    print("="*55)
    print("STEP 04 - Resistome Profiling (improved markers)")
    print(f"Date: {datetime.now().isoformat()}")
    print("="*55)
    countries = load_accessions()
    df = scan_genomes(countries)
    gene_cols = list(resistance_markers.keys())

    def count_classes(row):
        return sum(1 for genes in class_map.values() if any(row.get(g,0)==1 for g in genes))

    df["resistance_classes"] = df.apply(count_classes, axis=1)
    df["MDR"] = (df["resistance_classes"] >= 3).astype(int)
    prevalence = df[gene_cols].mean().sort_values(ascending=False)

    print("\n--- Resistance gene prevalence ---")
    for gene, prev in prevalence.items():
        if prev > 0:
            bar = "+" * int(prev * 40)
            print(f"  {gene:<15} {prev*100:5.1f}%  {bar}")
        else:
            print(f"  {gene:<15}  0.0%")

    print("\n--- Resistance classes per genome ---")
    print(df["resistance_classes"].value_counts().sort_index().to_string())

    print("\n--- MDR prevalence by country ---")
    mdr = df.groupby("country")["MDR"].agg(["sum","mean"])
    mdr.columns = ["mdr_count","mdr_rate"]
    mdr["mdr_rate"] = (mdr["mdr_rate"]*100).round(1)
    print(mdr.to_string())
    print(f"\nTotal MDR genomes: {df['MDR'].sum()}/{len(df)} ({df['MDR'].mean()*100:.1f}%)")

    df.to_csv(TABLE_DIR / "amrfinder_raw.tsv", sep="\t", index=False)
    matrix = df.set_index("accession")[gene_cols + ["MDR","country"]]
    matrix.to_csv(TABLE_DIR / "resistome_matrix.tsv", sep="\t")
    prevalence.reset_index().rename(
        columns={"index":"gene", 0:"prevalence"}).to_csv(
        TABLE_DIR / "table2_gene_prevalence.tsv", sep="\t", index=False)

    summary = {
        "date": datetime.now().isoformat(),
        "genomes": len(df),
        "genes_scanned": len(gene_cols),
        "mdr_total": int(df["MDR"].sum()),
        "mdr_pct": round(float(df["MDR"].mean()*100),1),
        "prevalence": {g: round(float(v*100),1) for g,v in prevalence.items()},
    }
    with open(LOG_DIR / "04_resistome_summary.json","w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nSaved: results/tables/amrfinder_raw.tsv")
    print(f"Saved: results/tables/resistome_matrix.tsv")
    print(f"Saved: results/logs/04_resistome_summary.json")
    print("\nDone.")

main()
