#!/usr/bin/bash

# Renames bed file chr annotation to match vcf format

# Argument 1 = bed file to convert

awk -F'\t' -vOFS='\t' 'NR>1 {$1="chr" $1} 1' ${1}>temp.txt && mv temp.txt ${1}
