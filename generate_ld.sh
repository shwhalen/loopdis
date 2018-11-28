#!/bin/bash

chrom=$1
super_pop=$2
bedfile_prefix=${LD_OUTPUT_DIR}/plink/${chrom}-${super_pop}

if [ -f ${bedfile_prefix}.ld.gz ]; then
    exit
fi

plink \
    --threads 1 \
    --memory 4096 \
    --bfile ${bedfile_prefix} \
    --r2 gz dprime \
    --ld-window-r2 0.01 \
    --ld-window 10000 \
    --ld-window-kb 2000 \
    --out ${bedfile_prefix}
