#!/usr/bin/bash

# Runs cis-QTL analysis on all covariance matrices to determine optimal number of RUV factors

# Argument 1 = covariate file directory
# Argument 2 = VCF file with genotype data
# Argument 3 = BED file with gene expression data
# Argument 4 = cis-QTL results output directory

if test -d "${4%/*}";
then 
    echo "Directory ${4%/*} exists.";
else
    echo "Directory ${4%/*} does not exists. Making directory...";
    mkdir ${4%/*};
fi

module load qtltools/1.3.1;

for FILE in ${1}*;
do file=${FILE%.*};
sbatch -t 2- --mem=10g --wrap="QTLtools cis --vcf ${2} --bed ${3}.gz --cov $FILE --permute 1000 --out ${4}.RUV.${file##*.}.txt --std-err";
done;