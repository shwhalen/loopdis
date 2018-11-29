#!/bin/bash

chroms=`seq 1 22`
super_pops="AFR AMR EAS EUR SAS"
encode_cell_types="K562 GM12878 IMR90 NHEK HUVEC"
blood_cell_types="Mon Mac0 Mac1 Mac2 Neu MK EP Ery FoeT nCD4 tCD4 aCD4 naCD4 nCD8 tCD8 nB tB"
cell_types="${encode_cell_types} ${blood_cell_types}"

# generate plink data common to all scripts (~7.5hrs)
# the majority of variants should start with 'rs', a few are submitted snps ('ss'), a handful have no id ('.')
parallel --bar ./generate_plink_bed.sh ::: ${chroms} ::: ${super_pops}
parallel --bar ./generate_ld_blocks.sh ::: ${chroms} ::: ${super_pops}
parallel --bar ./generate_ld.sh ::: ${chroms} ::: ${super_pops}

# convert plink ld to binary format (~12 hrs, mostly due to sorting and dropping duplicates)
parallel --bar --jobs 2 'python -c "from common import get_plink_ld; get_plink_ld({1}, \"{2}\");"' ::: ${chroms} ::: ${super_pops}

# convert rao contacts to binary format (~15m)
parallel --bar 'python -c "from common import get_rao_contacts; get_rao_contacts({1}, \"{2}\", True);"' ::: ${chroms} :::  ${encode_cell_types}
parallel --bar 'python -c "from common import get_rao_contacts; get_rao_contacts({1}, \"{2}\", False);"' ::: ${chroms} ::: ${encode_cell_types}

# significant generation runtimes:
# generate_interaction_ld.py: ~1 week, limited to single core due to ld ram requirements
# generate_domain_snps.py: ~1.5h
# generate_go_stats: ~12h
# generate_hic_ld_concordance: days
parallel --bar --jobs 1 ./generate_interaction_ld.py ::: ${chroms} ::: ${super_pops} ::: ${cell_types}
./generate_max_interaction_ld.py
./generate_dataset_stats.py
parallel --bar --jobs 1 ./generate_domain_snps.py ::: AFR AMR ::: K562 GM12878 ::: 0 1
./generate_eqtl_stats.py
./generate_go_stats.py
parallel --bar --jobs 1 ./generate_hic_ld_concordance.py ::: ${encode_cell_types}
./generate_interaction_crossings.py
./generate_target_stats.py
./preprocess_gencode_expression.py

# significant plot runtimes:
# plot_architectures.R: ~30m
# plot_hic_scaling-rao.R: ~30m
./plot_architectures.R
parallel --bar --jobs 4 ./plot_domain_crossing_ld.R ::: AFR AMR ::: K562 GM12878
./plot_eqtl_enrichment.R
./plot_eqtl_support.R
./plot_go_stats.R
./plot_heatmap.py config/chr14.json EAS GM12878 0 1 1
./plot_heatmap.py config/SPATA7.json EAS GM12878 1 0 0
./plot_heatmap.py config/SPATA7.json EAS NHEK 1 0 0
./plot_heatmap.py config/chr4.json EAS NHEK 1 1 1
./plot_heatmap.py config/SLC9B2.json EAS NHEK 1 0 0
./plot_heatmap.py config/SLC9B2.json EAS HUVEC 1 0 0
./plot_hic_ld_concordance.R
parallel --bar --jobs 4 ./plot_hic_scaling-javierre.R ::: ${blood_cell_types}
parallel --bar --jobs 2 ./plot_hic_scaling-rao.R ::: ${encode_cell_types} ::: 0 1
./plot_interaction_crossings.R
./plot_interaction_ld_ratios.py
./plot_interaction_ld_scaling.R
./plot_raw_ld_scaling.R
