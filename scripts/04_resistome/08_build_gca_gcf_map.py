import json
from pathlib import Path

mapping = {}

# Nigeria
ng_json = Path("data/raw/nigeria/ncbi_dataset/data/assembly_data_report.json")
if ng_json.exists():
    with open(ng_json) as f:
        data = json.load(f)
        for report in data.get("reports", []):
            asm = report.get("assembly_info", {})
            gcf = asm.get("assembly_accession")
            gca = asm.get("genbank_accession")
            if gcf and gca:
                mapping[gca] = gcf
                print(f"  {gca} → {gcf}")
    print(f"Nigeria: added {len([k for k in mapping if k.startswith('GCA')])} mappings")
else:
    print(f"⚠️ Nigeria JSON not found at {ng_json}")

# Ghana
gh_json = Path("data/raw/ghana/ncbi_dataset/data/assembly_data_report.json")
if gh_json.exists():
    with open(gh_json) as f:
        data = json.load(f)
        for report in data.get("reports", []):
            asm = report.get("assembly_info", {})
            gcf = asm.get("assembly_accession")
            gca = asm.get("genbank_accession")
            if gcf and gca:
                mapping[gca] = gcf
    print(f"Ghana: added {len([k for k in mapping if k.startswith('GCA')])} mappings")
else:
    print(f"⚠️ Ghana JSON not found at {gh_json}")

# Save mapping
out_path = Path("data/metadata/gca_to_gcf.json")
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, "w") as f:
    json.dump(mapping, f, indent=2)

print(f"\n✅ Saved {len(mapping)} GCA→GCF mappings to {out_path}")