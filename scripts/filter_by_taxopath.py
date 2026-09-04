import pandas as pd

df = pd.read_csv(snakemake.input.tsv, sep='\t')
filtered = df[df['taxopath'].str.contains(snakemake.params.search_term, na=False)]
filtered.to_csv(snakemake.output.tsv, sep='\t', index=False)
