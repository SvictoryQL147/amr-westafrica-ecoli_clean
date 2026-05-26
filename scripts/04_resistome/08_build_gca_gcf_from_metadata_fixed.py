import json
from pathlib import Path

mapping = {}

# Nigeria
ng_json = Path("data/metadata/nigeria_ecoli_metadata.json")
if ng_json.exists():
    with open(ng_json) as f:
        data = json.load(f)
        for report in data.get("reports", []):
            gcf = report.get("accession")
            asm_info = report.get("assembly_info", {})
            gca = asm_info.get("genbank_accession")
            if gcf and gca:
                mapping[gca] = gcf
    print(f"Nigeria: added {len([k for k in mapping if k.startswith('GCA')])} mappings")

# Ghana
gh_json = Path("data/metadata/ghana_ecoli_metadata.json")
if gh_json.exists():
    with open(gh_json) as f:
        data = json.load(f)
        for report in data.get("reports", []):
            gcf = report.get("accession")
            asm_info = report.get("assembly_info", {})
            gca = asm_info.get("genbank_accession")
            if gcf and gca:
                mapping[gca] = gcf
    print(f"Ghana: added {len([k for k in mapping if k.startswith('GCA')])} mappings")

# Save mapping
out_path = Path("data/metadata/gca_to_gcf.json")
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, "w") as f:
    json.dump(mapping, f, indent=2)

print(f"\n✅ Saved {len(mapping)} GCA→GCF mappings to {out_path}")