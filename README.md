IBD colon eQTL analysis
=======================

eQTL analysis pipeline


1. Starting with processed RNA-seq data, BAM files are converted into individual QTLtools BED formatted files:
```
QTLtools quan --bam sample.bam --gtf gencode.v39.annotation.gtf.gz --sample sample --out-prefix sample --no-hash
```

2. Create a raw count matrix:
```
analyses/concat_RNA_bed.pl <path_to_BED_files> <path_to_raw_count_matrix> <path_to_analyses_dir>
```

3. Calculate RUV factors:
```
```

4. Create covariate matrices:
```
```

5. Normalize counts:
```
```

6. Zip and index files:
```
```

7. Perform cis-eQTL mapping:
```
```

8. Determine number of RUV factors to include in the final model:
```
```

9. Run multiple hypothesis testing correct and calculate p-value thresholds
```
Rscript analyses/qtltools_runFDR_cis.R IBD_cis-eQTL_RUV_21.txt 0.05 IBD_cis-eQTL_lead
```

10. Conditional cis-eQTL mapping:
```
QTLtools cis --vcf vcf.gz  --bed normalized.bed.gz  --cov covariates.txt  --mapping IBD_cis-eQTL_lead.thresholds.txt  --out IBD_cis-eQTL_conditional.txt
```

11. Full cis-eQTL summary statistics:
```
QTLtools cis --vcf vcf.gz --bed normalized.bed.gz --cov covariates.txt --std-err --nominal 1.0 --out IBD_cis-eQTL_full_summary_stats.txt
```
