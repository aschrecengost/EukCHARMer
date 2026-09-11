from pathlib import Path
import pandas as pd

input_files = list(snakemake.input.cleaned)
output_file = Path(str(snakemake.output.merged))

output_file.parent.mkdir(parents=True, exist_ok=True)

dataframes = [pd.read_csv(file) for file in input_files]
merged = pd.concat(dataframes, ignore_index=True, sort=False)
merged.to_csv(output_file, index=False)

print(f"Merged {len(input_files)} projects and {len(merged)} total rows")
