import pandas as pd

df = pd.read_csv(snakemake.input.tsv, sep="\t", dtype=str)
names = set(df["name"].dropna().astype(str).str.strip())
names.discard("")

if len(df) and not names:
    raise ValueError(
        "The placement table contains rows but no usable sequence names"
    )

found = set()

with open(snakemake.input.fasta) as fin, \
     open(snakemake.output.fasta, "w") as fout:

    write = False
    seq_name = None

    for line in fin:
        if line.startswith(">"):
            seq_name = line[1:].strip().split()[0]
            write = seq_name in names
        elif write and line.strip():
            found.add(seq_name)

        if write:
            fout.write(line)

missing = names - found

if missing:
    examples = ", ".join(sorted(missing)[:5])
    raise ValueError(
        f"{len(missing)} requested sequence IDs were not found "
        f"with sequence data in {snakemake.input.fasta}: {examples}"
    )

print(
    f"Extracted {len(found)} sequences into "
    f"{snakemake.output.fasta}"
)
