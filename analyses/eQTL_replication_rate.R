# IBD eQTL replication rate in GTEx and BarcUVa
library(data.table)
library(tidyverse)
library(qvalue)

# load lead eQTL summary statistics
ibd <- fread('data/UNC/Supplementary_Table_1_IBD_cis-eQTL_lead.significant.txt')
# load full GTEx summary statistics
gtex <- fread('data/GTEx/GTEx_Analysis_v8_eQTL_all_associations/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Colon_Transverse.allpairs.txt.gz')
# load full BarcUVa summary statistics
barcuva <- fread('data/BarcUVa/barcuvaseq.eqtls.allpairs.sumstats.hg38.txt')
# reformat gene IDs
ibd <- ibd %>% separate_wider_delim(phe_id, delim = '.', names = c('gene', 'ver'))
gtex <- gtex %>% separate_wider_delim(gene_id, delim = '.', names = c('gene', 'ver'))
barcuva <- barcuva %>% separate_wider_delim(gene_id, delim = '.', names = c('gene', 'ver'))
# reformat GTEx variants
gtex$variant_id <- gsub('_b38', '', gtex$variant_id)
gtex$variant_id <- gsub('_', ':', gtex$variant_id)
# create eSNP-eGene pairs
ibd <- ibd %>% unite('eqtl', c(gene, var_id), sep = '-', remove = FALSE)
gtex <- gtex %>% unite('eqtl', c(gene, variant_id), sep = '-', remove = FALSE)
barcuva <- barcuva %>% unite('eqtl', c(gene, variant_id), sep = '-', remove = FALSE)
# find shared eSNP-eGene pairs
pairs.gtex <- intersect(ibd$eqtl, gtex$eqtl)
pairs.barcuva <- intersect(ibd$eqtl, barcuva$eqtl)
# subset shared pairs in GTEx and BarcUVa
g <- subset(gtex, eqtl %in% pairs.gtex)
b <- subset(barcuva, eqtl %in% pairs.barcuva)
# calculate pi0, estimate of the proportion of null p-values
q.gtex <- qvalue(g$pval_nominal)
q.barcuva <- qvalue(b$pval_nominal)
# calculate pi1, estimate of the proportion of significant p-values aka replication rate
pi1.gtex <- 1 - q.gtex$pi0
pi1.barcuva <- 1 - q.barcuva$pi0
