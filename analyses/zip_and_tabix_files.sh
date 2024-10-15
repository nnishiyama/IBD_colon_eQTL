#!/usr/bin/bash

# Zips and tabix indexes bed & vcf files for QTL input

# Argument 1 = BED file with gene expression data
# Argument 2 = VCF file with genotype data

if test ! -f "${1}.gz"
then 
    sbatch -t 1- --mem=4g --wrap="bgzip -c ${1} > ${1}.gz;
    tabix -f -p bed ${1}.gz"
else
    echo "bed file is zipped"
fi

if test ! -f "${2}.tbi"
then 
    sbatch -t 1- --mem=10g --wrap="tabix -f -p vcf ${2}"
else
    echo "vcf file index exists"
fi