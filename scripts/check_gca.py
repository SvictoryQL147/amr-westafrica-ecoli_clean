import pandas as pd
meta = pd.read_csv("data/metadata/master_metadata.tsv", sep="\t")
print(meta[["accession", "current_accession"]].head(10))