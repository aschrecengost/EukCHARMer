import pandas as pd

perquery = pd.read_csv(snakemake.input.per_query, sep='\t')
lwr      = pd.read_csv(snakemake.input.lwr_list).rename(columns={'PqueryName': 'name'})
edpl     = pd.read_csv(snakemake.input.edpl_list).rename(columns={'Pquery': 'name'})

merged = perquery.merge(lwr, on='name').merge(edpl, on='name')
filtered = merged[merged['EDPL'] < snakemake.params.edpl_threshold]

drop_cols = [c for c in filtered.columns
             if c in {'Sample_x', 'Sample_y', 'Multiplicity_x', 'Multiplicity_y'}]
filtered = filtered.drop(columns=drop_cols)
filtered = filtered[['name'] + [c for c in filtered.columns if c != 'name']]

filtered.to_csv(snakemake.output.tsv, sep='\t', index=False)
