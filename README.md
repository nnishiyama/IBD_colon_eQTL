IBD colon eQTL analysis
=======================

eQTL analysis pipeline


1. Starting with processed RNA-seq data, BAM files are converted into individual QTLtools BED formatted files:
```
QTLtools quan --bam <sample.bam> --gtf gencode.v39.annotation.gtf.gz --sample <sample_ID> --out-prefix <sample_ID> --no-hash
```

2. Create a raw count matrix:
```
analyses/concat_RNA_bed.pl <path_to_BED_files> <path_to_raw_count_matrix> <path_to_analyses_dir>
```

3. Calculate RUV factors:
```
Rscript analyses/RUVseq.R <path_to_raw_count_matrix> IBD_coldata_UIDs.txt
```

4. Create covariate matrices:
```
python3 analyses/make_covariance_matrices.py --RNA_covariates IBD_coldata_UIDs.txt --ruv RUV_factors.txt --out_prefix covariate_matrices/covariates_RUV
```

5. Normalize counts:
```
Rscript analyses/normalize_counts.R <path_to_raw_count_matrix> IBD_coldata_UIDs.txt > normalized.bed
```

6. Zip and index files:
```
analyses/zip_and_tabix_files.sh normalized.bed vcf.gz
```

7. Perform cis-eQTL mapping:
```
analyses/QTLtools_batch_run.sh covariate_matrices/covariates_RUV <path_to_vcf.gz> normalized.bed <path_to_results_directory/results_file_prefix>
```

8. Determine number of RUV factors to include in the final model:
```
Rscript analyses/kneedle.R <path_to_results_directory>
```

9. Run multiple hypothesis testing correct and calculate p-value thresholds
```
Rscript analyses/qtltools_runFDR_cis.R IBD_cis-eQTL_RUV_21.txt 0.05 IBD_cis-eQTL_lead
```

10. Conditional cis-eQTL mapping:
```
QTLtools cis --vcf <path_to_vcf.gz>  --bed normalized.bed.gz  --cov <path_to_covariate_file.txt>  --mapping IBD_cis-eQTL_lead.thresholds.txt  --out IBD_cis-eQTL_conditional.txt
```

11. Full cis-eQTL summary statistics:
```
QTLtools cis --vcf <path_to_vcf.gz> --bed normalized.bed.gz --cov <path_to_covariate_file.txt> --std-err --nominal 1.0 --out IBD_cis-eQTL_full_summary_stats.txt
```
