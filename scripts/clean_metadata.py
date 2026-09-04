import pandas as pd

df = pd.read_csv(snakemake.input.metadata)

# Rename all aliases to their target column name
flat_rename = {alias: target
               for target, aliases in snakemake.params.rename.items()
               for alias in aliases
               if alias in df.columns and alias != target}
df.rename(columns=flat_rename, inplace=True)

# Merge duplicate columns that arose from multiple aliases mapping to the same target
for target in snakemake.params.rename:
    if df.columns.tolist().count(target) > 1:
        df[target] = df[[target]].bfill(axis=1).iloc[:, 0]
        df = df.loc[:, ~df.columns.duplicated()]

# Build lat_lon from raw lat/lon columns if present and not already created
if "lat_lon" not in df.columns:
    if ("geographic location (latitude)" in df.columns
            and "geographic location (longitude)" in df.columns):
        df["lat_lon"] = (df["geographic location (latitude)"].astype(str)
                         + " " + df["geographic location (longitude)"].astype(str))

keep = [c for c in snakemake.params.keep if c in df.columns]
df[keep].to_csv(snakemake.output.cleaned, index=False)
