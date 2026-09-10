import pandas as pd
import os
is.makedirs("documents/merged", exist_ok=True)
df = pd.concat([pd.read_csv(f) for f in input.cleaned], ignore_index=True)
df.to_csv(output.merged, index=False)
print(f"Merged {len(input.cleaned)} projects, {len(df)} total rows")
