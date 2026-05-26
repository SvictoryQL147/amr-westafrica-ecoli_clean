"""
STEP 02 — Metadata Parser
=========================
Parses NCBI metadata JSON files for Nigeria and Ghana E. coli genomes.
Produces:
  - data/metadata/master_metadata.tsv
  - results/tables/table1_sample_summary.tsv
  - results/logs/02_metadata_summary.txt
"""

import json
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE_DIR  = Path(".")
META_DIR  = BASE_DIR / "data" / "metadata"
TABLE_DIR = BASE_DIR / "results" / "tables"
LOG_DIR   = BASE_DIR / "results" / "logs"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

def parse_metadata(json_file, country):
    with open(json_file, encoding="utf-8-sig") as f:
        data = json.load(f)

    records = []
    for item in data.get("reports", []):
        asm  = item.get("assembly_info", {})
        org  = item.get("organism", {})
        seq  = item.get("assembly_stats", {})
        chk  = item.get("checkm_info", {})
        bio  = asm.get("biosample", {})

        records.append({
            "accession":         item.get("accession", ""),
            "current_accession": item.get("current_accession", ""),
            "country":           country,
            "organism":          org.get("organism_name", ""),
            "strain":            org.get("infraspecific_names", {}).get("strain", ""),
            "assembly_level":    asm.get("assembly_level", ""),
            "assembly_status":   asm.get("assembly_status", ""),
            "assembly_type":     asm.get("assembly_type", ""),
            "assembly_name":     asm.get("assembly_name", ""),
            "release_date":      asm.get("release_date", ""),
            "sequencing_tech":   asm.get("sequencing_tech", ""),
            "submitter":         asm.get("submitter", ""),
            "bioproject":        asm.get("bioproject_accession", ""),
            "biosample":         bio.get("accession", ""),
            "total_length":      seq.get("total_sequence_length", ""),
            "contig_n50":        seq.get("contig_n50", ""),
            "scaffold_n50":      seq.get("scaffold_n50", ""),
            "contig_count":      seq.get("number_of_contigs", ""),
            "gc_percent":        seq.get("gc_percent", ""),
            "checkm_completeness":   chk.get("completeness", ""),
            "checkm_contamination":  chk.get("contamination", ""),
        })
    return records

# ── Parse both countries ───────────────────────────────────────────────────────
all_records = []
for country, fname in [("Nigeria", "nigeria_ecoli_metadata.json"),
                        ("Ghana",   "ghana_ecoli_metadata.json")]:
    fpath = META_DIR / fname
    if fpath.exists():
        recs = parse_metadata(fpath, country)
        all_records.extend(recs)
        print(f"{country}: {len(recs)} records parsed")
    else:
        print(f"WARNING: {fpath} not found — skipping")

df = pd.DataFrame(all_records)

# ── Save master metadata ───────────────────────────────────────────────────────
master_path = META_DIR / "master_metadata.tsv"
df.to_csv(master_path, sep="\t", index=False)
print(f"\nMaster metadata saved: {master_path} ({len(df)} rows)")

# ── Assembly level breakdown ───────────────────────────────────────────────────
print("\n--- Assembly level breakdown ---")
lvl = df.groupby(["country", "assembly_level"]).size().reset_index(name="count")
print(lvl.to_string(index=False))

# ── Sequencing technology breakdown ───────────────────────────────────────────
print("\n--- Sequencing technology ---")
tech = df.groupby(["country", "sequencing_tech"]).size().reset_index(name="count")
print(tech.to_string(index=False))

# ── Release year breakdown ─────────────────────────────────────────────────────
df["release_year"] = pd.to_datetime(df["release_date"], errors="coerce").dt.year
print("\n--- Genomes released by year ---")
yr = df.groupby(["country", "release_year"]).size().reset_index(name="count")
print(yr.to_string(index=False))

# ── Numeric stats ──────────────────────────────────────────────────────────────
for col in ["total_length", "contig_n50", "contig_count", "gc_percent",
            "checkm_completeness", "checkm_contamination"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

print("\n--- Genome stats by country ---")
stats = df.groupby("country").agg(
    total_genomes   = ("accession",      "count"),
    mean_length_mb  = ("total_length",   lambda x: round(x.mean()/1e6, 2)),
    mean_n50_kb     = ("contig_n50",     lambda x: round(x.mean()/1e3, 1)),
    mean_contigs    = ("contig_count",   lambda x: round(x.mean(), 1)),
    mean_gc         = ("gc_percent",     lambda x: round(x.mean(), 2)),
    mean_completeness  = ("checkm_completeness",  lambda x: round(x.mean(), 2)),
    mean_contamination = ("checkm_contamination", lambda x: round(x.mean(), 2)),
).reset_index()
print(stats.to_string(index=False))

# ── Table 1 for manuscript ─────────────────────────────────────────────────────
table1 = df.groupby("country").agg(
    total_genomes      = ("accession",       "count"),
    scaffold           = ("assembly_level",  lambda x: (x=="Scaffold").sum()),
    contig             = ("assembly_level",  lambda x: (x=="Contig").sum()),
    chromosome         = ("assembly_level",  lambda x: (x=="Chromosome").sum()),
    mean_length_mb     = ("total_length",    lambda x: round(x.mean()/1e6, 2)),
    mean_n50_kb        = ("contig_n50",      lambda x: round(x.mean()/1e3, 1)),
    mean_gc_percent    = ("gc_percent",      lambda x: round(x.mean(), 2)),
).reset_index()

table1_path = TABLE_DIR / "table1_sample_summary.tsv"
table1.to_csv(table1_path, sep="\t", index=False)

print(f"\nTable 1 saved: {table1_path}")
print("\n=== TABLE 1 — MANUSCRIPT READY ===")
print(table1.to_string(index=False))

# ── Save full log ──────────────────────────────────────────────────────────────
log_path = LOG_DIR / "02_metadata_summary.txt"
with open(log_path, "w") as f:
    f.write(f"Metadata parsing — {datetime.now().isoformat()}\n")
    f.write(f"Total genomes: {len(df)}\n\n")
    f.write("Assembly levels:\n")
    f.write(lvl.to_string(index=False))
    f.write("\n\nSequencing tech:\n")
    f.write(tech.to_string(index=False))
    f.write("\n\nRelease years:\n")
    f.write(yr.to_string(index=False))
    f.write("\n\nGenome stats:\n")
    f.write(stats.to_string(index=False))
    f.write("\n\nTable 1:\n")
    f.write(table1.to_string(index=False))

print(f"\nLog saved: {log_path}")
print("\nDone. Copy Table 1 numbers into MANUSCRIPT_DRAFT.md and LAB_NOTEBOOK.md")
