#!/bin/bash

chrom=$1
super_pop=$2
bedfile_prefix=${LD_OUTPUT_DIR}/plink/${chrom}-${super_pop}

if [ -f ${bedfile_prefix}.blocks.det ]; then
    exit
fi

plink \
    --threads 1 \
    --memory 4096 \
    --bfile ${bedfile_prefix} \
    --blocks no-pheno-req no-small-max-span \
    --blocks-max-kb 2000 \
    --out ${bedfile_prefix}
