#!/bin/bash

chrom=$1
super_pop=$2
bcf_fn=${DB_DIR}/1kg/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v*.20130502.genotypes.bcf
samples_fn=${DB_DIR}/1kg/integrated_call_samples_v3.20130502.ALL.panel
bedfile_prefix=${LD_OUTPUT_DIR}/plink/${chrom}-${super_pop}

if [ -f ${bedfile_prefix}.bed ]; then
    exit
fi

plink \
    --threads 1 \
    --memory 4096 \
    --bcf ${bcf_fn} \
    --allow-extra-chr \
    --snps-only \
    --biallelic-only strict \
    --maf 0.05 \
    --make-bed \
    --out ${bedfile_prefix} \
    --filter <( \
        grep ${super_pop} ${samples_fn} \
        | cut -f1 \
        | awk '{ printf("%s\t%s\t0\n", $1, $1) }' \
    ) 0
