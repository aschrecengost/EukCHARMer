import pandas as pd

df = pd.read_csv(snakemake.input.tsv, sep='\t')
names = set(df['name'].astype(str))

with open(snakemake.input.fasta) as fin, open(snakemake.output.fasta, 'w') as fout:
    write = False
    for line in fin:
        if line.startswith('>'):
            seq_name = line[1:].strip().split()[0]
            write = seq_name in names
        if write:
            fout.write(line)
