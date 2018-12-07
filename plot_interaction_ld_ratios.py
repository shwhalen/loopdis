#!/usr/bin/env python

import feather
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from common import *

max_interaction_ld = feather.read_dataframe(f'{output_dir}/max_interaction_ld.feather')
max_interaction_ld['cell_type'] = max_interaction_ld['cell_type'].str.replace('IMR90', 'IMR-90')
print(max_interaction_ld['cell_type'].value_counts())

max_interaction_ld = (
    max_interaction_ld
    .groupby(['cell_type', 'super_pop', 'is_significant'])
    ['max_ld']
    .apply(np.mean)
)

(
    max_interaction_ld
    .rename(r'mean(max(\RSQUARED))')
    .rename_axis(['Cell Type', 'Super-population', 'Significant'])
    .reindex(encode_cell_types + blood_cell_types, level = 0)
    .to_frame()
    .to_latex(
        'output/interaction_ld_raw-table.tex',
        escape = False,
        float_format = '%.3f',
        longtable = True
    )
)

ld_ratios = (
    max_interaction_ld
    .groupby(level = ['cell_type', 'super_pop'])
    .apply(lambda x: np.log2(x.iloc[1] / x.iloc[0]))
    .rename('ld_ratio')
)

plt.figure(figsize = (11, 2))
cg = sns.heatmap(
    (
        ld_ratios
        .reset_index()
        .pivot('super_pop', 'cell_type', 'ld_ratio')
        .reindex(columns = encode_cell_types + blood_cell_types)
    ),
    cmap = 'RdBu_r',
    center = 0,
    cbar_kws= {'label': 'log$_2$(Interacting:Non-Interacting LD Ratio)'}
)
plt.setp(cg.get_xticklabels(), rotation = 90)
plt.setp(cg.get_yticklabels(), rotation = 0)
plt.ylabel('Super-population')
plt.xlabel('Cell Type')
plt.savefig('output/interaction_ld_ratios.pdf', bbox_inches = 'tight')
