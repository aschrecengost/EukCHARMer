#!/usr/bin/env python
"""Write an empty QIIME2 artifact, used as a placeholder when a taxonomic
filter (qiime taxa filter-table/filter-seqs) matches zero features and
refuses to write an output. Downstream `qiime feature-table merge`/
`merge-seqs` accept an empty artifact and simply contribute nothing from it.

Usage: make_empty_qza.py <semantic-type> <output-path>
"""
import sys

import biom
import numpy as np
import pandas as pd
import qiime2

SEMANTIC_TYPE, OUTPUT_PATH = sys.argv[1], sys.argv[2]

if SEMANTIC_TYPE == "FeatureTable[Frequency]":
    view = biom.Table(np.empty((0, 0)), [], [])
elif SEMANTIC_TYPE == "FeatureData[Sequence]":
    view = pd.Series([], dtype=object, name="Sequence")
else:
    raise ValueError(f"Unsupported semantic type: {SEMANTIC_TYPE}")

qiime2.Artifact.import_data(SEMANTIC_TYPE, view).save(OUTPUT_PATH)
