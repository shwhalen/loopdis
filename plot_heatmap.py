#!/usr/bin/env python

import json
import matplotlib.gridspec as grid
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

from common import *
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle, Rectangle
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from pathlib import Path

# GR requires same font size within figure
general_fontsize = 10
plt.rc('legend', fontsize = general_fontsize)
plt.rc('xtick', labelsize = general_fontsize)
plt.rc('ytick', labelsize = general_fontsize)

config_fn = Path(sys.argv[1])
config = json.load(open(config_fn))
super_pop = sys.argv[2]
cell_type = sys.argv[3]
observed_over_expected = int(sys.argv[4])
show_legend = int(sys.argv[5])
show_colorbar = int(sys.argv[6])
print(json.dumps(config, indent = 4))

chrom = str(config['locus_chrom'])
locus_start = config['locus_start']
locus_end = config['locus_end']
ld_bin_size = config['ld_bin_size'] if 'ld_bin_size' in config else 5000
annotations = config['annotations'] if 'annotations' in config else None
oe_flag = 'oe' if observed_over_expected else 'o'
ld_window_r2 = 0.1 # r2 threshold for plot *only* (named after plink param)
lw = 1.3
ld_lw = lw * 0.5 if locus_end - locus_start > 1e6 else lw
alpha = 1
loop_color = 'orange'
interpolation = 'nearest' # matplotlib pdf renderer has issues: choose 'nearest' for non-blurry pixels, 'none' for smooth diagonals

min_bin = locus_start // ld_bin_size
max_bin = locus_end // ld_bin_size
bin_and_shift = lambda x: x // ld_bin_size - min_bin

print('\nreading significant hi-c loops')
loops_df = (
    get_rao_loops(cell_type)
    .query(
        'f1_chr == @chrom and \
        @locus_start <= f1_start and \
        f1_end <= @locus_end and \
        @locus_start <= f2_start and \
        f2_end <= @locus_end'
    )
)

loops_df['f1_start'] = bin_and_shift(loops_df['f1_start'])
loops_df['f1_end'] = bin_and_shift(loops_df['f1_end'])
loops_df['f2_start'] = bin_and_shift(loops_df['f2_start'])
loops_df['f2_end'] = bin_and_shift(loops_df['f2_end'])

print('\nreading hi-c contact domains')
domains_df = (
    get_rao_domains(cell_type)
    .query('domain_chr == @chrom')
)
domains_df['domain_start'] = bin_and_shift(domains_df['domain_start'])
domains_df['domain_end'] = bin_and_shift(domains_df['domain_end'])

print('\nreading hi-c contacts')
contacts_df = (
    get_rao_contacts(chrom, cell_type, observed_over_expected)
    .reset_index()
)
contacts_df['f1_start'] = bin_and_shift(contacts_df['f1_start'])
contacts_df['f2_start'] = bin_and_shift(contacts_df['f2_start'])

# swap coords to plot lower triangle
contacts_df.rename(
    columns = {
        'f1_start': 'f2_start',
        'f2_start': 'f1_start'
    },
    inplace = True
)
aggregate_contacts = (
    contacts_df
    .groupby(['f1_start', 'f2_start'], as_index = False)
    .mean()
)
del contacts_df

print('\nreading pairwise ld')
ld_df = (
    get_plink_ld(chrom, super_pop)
    .query(
        '@locus_start <= snp1_start <= @locus_end and \
        @locus_start <= snp2_start <= @locus_end'
    )
    .reset_index()
)
ld_df['snp1_start'] = bin_and_shift(ld_df['snp1_start'])
ld_df['snp2_start'] = bin_and_shift(ld_df['snp2_start'])
aggregate_ld = (
    ld_df
    .groupby(['snp1_start', 'snp2_start'], as_index = False)
    .mean()
)
del ld_df

print('\nreading ld blocks')
ld_blocks_df = (
    get_plink_ld_blocks(chrom, super_pop)
    .query(
        '@locus_start <= ld_block_start and \
        ld_block_end <= @locus_end'
    )
)
ld_blocks_df['ld_block_start'] = bin_and_shift(ld_blocks_df['ld_block_start'])
ld_blocks_df['ld_block_end'] = bin_and_shift(ld_blocks_df['ld_block_end'])

print('\npivoting')
index = range(max_bin - min_bin)
ld_heatmap_df = (
    aggregate_ld
    .query('r2 >= @ld_window_r2')
    .pivot('snp1_start', 'snp2_start', 'r2')
    .reindex(index = index, columns = index)
)
contacts_heatmap_df = (
    aggregate_contacts
    .pivot('f1_start', 'f2_start', 'contact_value')
    .reindex(index = index, columns = index)
)
# if not observed_over_expected:
#     contacts_heatmap_df = np.log2(contacts_heatmap_df)

if contact_bin_size != ld_bin_size:
    print('\ninterpolating hi-c bins')
    stride = contact_bin_size // ld_bin_size
    for i, j in zip(*np.where(contacts_heatmap_df.notnull())):
        contacts_heatmap_df.iloc[i:(i + stride), j:(j + stride)] = contacts_heatmap_df.iloc[i, j]
    upper_triangle_mask = np.triu(np.ones_like(contacts_heatmap_df)).astype(bool)
    contacts_heatmap_df.values[upper_triangle_mask] = np.nan

print('\nplotting')
fig = plt.figure(figsize = (9, 6))
gs = grid.GridSpec(
    1,
    3,
    width_ratios = [1, 1, 20],
    wspace = 0.55
)

main_ax = plt.subplot(gs[2])
ld_plot = main_ax.imshow(
    ld_heatmap_df,
    cmap = 'Greens',
    interpolation = interpolation
)
contacts_plot = main_ax.imshow(
    contacts_heatmap_df,
    cmap = 'Purples',
    interpolation = interpolation,
    norm = LogNorm()
)

locus_size_kb = (locus_end - locus_start) // 1000
print(locus_size_kb)
scalebar_size = 100 if locus_size_kb < 5e3 else 1000
scalebar = AnchoredSizeBar(
    transform = main_ax.transData,
    size = scalebar_size,
    label = f'Scale: {scalebar_size} kb',
    loc = 'upper right',
    borderpad = 1,
    pad = 0.1,
    frameon = False,
    label_top = False
)
main_ax.add_artist(scalebar)

if show_colorbar:
    ld_colorbar_ax = plt.subplot(gs[0])
    ld_colorbar_ax.tick_params(labelsize = general_fontsize)
    ld_colorbar_label = f'Mean LD ({ld_bin_size // 1000}kb bins, MAF {maf:.2f})'
    ld_colorbar = fig.colorbar(
        ld_plot,
        cax = ld_colorbar_ax,
        orientation = 'vertical'
    )
    ld_colorbar.set_label(ld_colorbar_label, size = general_fontsize)

    contacts_colorbar_ax = plt.subplot(gs[1])
    contacts_colorbar_ax.tick_params(labelsize = general_fontsize)
    contacts_colorbar_label = f'Normalized Hi-C Contact Count ({contact_bin_size // 1000}kb bins)' #if observed_over_expected \
    #    else f'log$_2$ Normalized Hi-C Contact Count ({contact_bin_size // 1000}kb bins)'
    contacts_colorbar = fig.colorbar(
        contacts_plot,
        cax = contacts_colorbar_ax,
        orientation = 'vertical'
    )
    contacts_colorbar.set_label(contacts_colorbar_label, size = general_fontsize)

contact_patch = mpatches.Patch(
    color = plt.cm.Purples(127),
    label = f'Hi-C Contact Count ({contact_bin_size // 1000}kb bins)'
)
domain_patch = mlines.Line2D(
    [],
    [],
    color = plt.cm.Purples(255),
    linestyle = '--',
    label = 'Contact Domain'
)
ld_patch = mpatches.Patch(
    color = plt.cm.Greens(127),
    label = f'Mean LD ({ld_bin_size // 1000}kb bins, MAF {maf:.2f})'
)
ld_block_patch = mlines.Line2D(
    [],
    [],
    color = plt.cm.Greens(255),
    label = 'LD Block'
)
loop_patch = mlines.Line2D(
    [],
    [],
    color = loop_color,
    label = 'Significant Loop'
)
if show_legend:
    handles = [domain_patch, ld_block_patch, loop_patch] if show_colorbar \
        else [contact_patch, domain_patch, ld_patch, ld_block_patch, loop_patch]
    main_ax.legend(
        handles = handles,
        loc = 'center left',
        bbox_to_anchor = (1, 0.5),
        ncol = 1
    )

main_ax.spines['top'].set_visible(False)
main_ax.spines['right'].set_visible(False)
main_ax.spines['bottom'].set_visible(False)
main_ax.spines['left'].set_visible(False)
main_ax.xaxis.set_label_position('top')

if annotations:
    print('\nadding annotations')

    annotations_df = pd.DataFrame(annotations, columns = ['start', 'label'])
    annotations_df['start'] = bin_and_shift(annotations_df['start'])
    annotations_df = (
        annotations_df
        .groupby('start')
        ['label']
        .apply(
            lambda x: ' +\n'.join(x)
        )
    )
    annotations = annotations_df.to_dict()

    tick_locations = list(annotations.keys())
    tick_labels = list(annotations.values())

    main_ax.set_xticks(tick_locations)
    main_ax.set_yticks(tick_locations)

    main_ax.set_xticklabels(tick_labels, rotation = 90)
    main_ax.set_yticklabels(tick_labels)

    main_ax.tick_params(
        axis = 'both',
        which = 'both',
        length = 5
    )
    main_ax.xaxis.tick_top()
    main_ax.yaxis.tick_left()
else:
    main_ax.xaxis.set_visible(False)
    main_ax.yaxis.set_visible(False)

print('\nadding hi-c contact domains to heatmap')
kwargs = {
    'fill': False,
    'edgecolor': plt.cm.Purples(255),
    'lw': lw,
    'linestyle': ':',
    'alpha': alpha
}
for _, x1, x2 in domains_df.itertuples(index = False):
    w = x2 - x1
    main_ax.add_patch(
        Rectangle(
            (x1, x1),
            w,
            w,
            **kwargs
        )
    )

print('\nadding ld blocks to heatmap')
kwargs = {
    'fill': False,
    'edgecolor': plt.cm.Greens(255),
    'lw': ld_lw,
    'alpha': alpha
}
for _, x1, x2, _, _ in ld_blocks_df.itertuples(index = False):
    w = x2 - x1
    main_ax.add_patch(
        Rectangle(
            (x1, x1),
            w,
            w,
            **kwargs
        )
    )

print('\nadding hi-c loops to heatmap')
kwargs = {
    'fill': True,
    'edgecolor': loop_color,
    'facecolor': loop_color,
    'lw': lw,
    'alpha': alpha
}
for _, x1, x2, _, y1, y2 in loops_df.itertuples(index = False):
    main_ax.add_patch(
        Rectangle(
            (x1, y1),
            x2 - x1,
            y2 - y1,
            **kwargs
        )
    )

print('\nsaving')
output_fn = f'output/hic_ld-{config_fn.stem}-{super_pop}-{cell_type}-{oe_flag}.pdf'
plt.savefig(output_fn, bbox_inches = 'tight')
