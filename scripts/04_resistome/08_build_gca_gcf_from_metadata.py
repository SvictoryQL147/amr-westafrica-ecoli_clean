import json
from pathlib import Path

mapping = {}

# Nigeria metadata
ng_json = Path("data/metadata/nigeria_ecoli_metadata.json")
if ng_json.exists():
    with open(ng_json) as f:
        data = json.load(f)
        # The JSON structure from datasets: it's a dictionary with key "assemblies"?
        # Let's handle both list and dict
        if isinstance(data, dict) and "assemblies" in data:
            assemblies = data["assemblies"]
        elif isinstance(data, list):
            assemblies = data
        else:
            assemblies = [data]
        for item in assemblies:
            asm_info = item.get("assembly_info", {})
            gcf = asm_info.get("assembly_accession")
            gca = asm_info.get("genbank_accession")
            if gcf and gca:
                mapping[gca] = gcf
    print(f"Nigeria: added {len([k for k in mapping if k.startswith('GCA')])} mappings")
else:
    print(f"⚠️ Nigeria metadata not found at {ng_json}")

# Ghana metadata
gh_json = Path("data/metadata/ghana_ecoli_metadata.json")
if gh_json.exists():
    with open(gh_json) as f:
        data = json.load(f)
        if isinstance(data, dict) and "assemblies" in data:
            assemblies = data["assemblies"]
        elif isinstance(data, list):
            assemblies = data
        else:
            assemblies = [data]
        for item in assemblies:
            asm_info = item.get("assembly_info", {})
            gcf = asm_info.get("assembly_accession")
            gca = asm_info.get("genbank_accession")
            if gcf and gca:
                mapping[gca] = gcf
    print(f"Ghana: added {len([k for k in mapping if k.startswith('GCA')])} mappings")
else:
    print(f"⚠️ Ghana metadata not found at {gh_json}")

# Save mapping
out_path = Path("data/metadata/gca_to_gcf.json")
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, "w") as f:
    json.dump(mapping, f, indent=2)

print(f"\n✅ Saved {len(mapping)} GCA→GCF mappings to {out_path}")