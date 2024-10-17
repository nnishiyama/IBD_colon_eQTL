# Annotating full summary statistics with rsIDs

library(data.table)
library(tidyverse)

# load dbsnp key
dbsnp <- fread('data/other/dbsnp_b156_rsIDs.txt.gz', header = F)
# reformat
dbsnp <- dbsnp %>% separate_longer_delim(V4, delim = ',')
dbsnp <- dbsnp %>% unite('snp', V1:V4, sep = ':')
colnames(dbsnp) <- c('snp', 'rsID')
# load GWAS full summary stats
cd <- fread('data/GWAS/ibd_EAS_EUR_SiKJEF_meta_CD.TBL.txt.gz')
ibd <- fread('data/GWAS/ibd_EAS_EUR_SiKJEF_meta_IBD.TBL.txt.gz')
uc <- fread('data/GWAS/ibd_EAS_EUR_SiKJEF_meta_UC.TBL.txt.gz')
# annotate GWAS full summary stats
cd <- merge(dbsnp, cd, by.x = 'snp', by.y = 'MarkerName', all.y = TRUE)
ibd <- merge(dbsnp, ibd, by.x = 'snp', by.y = 'MarkerName', all.y = TRUE)
uc <- merge(dbsnp, uc, by.x = 'snp', by.y = 'MarkerName', all.y = TRUE)
# load eQTL full summary stats
eqtl.barcuva <- fread('/work/users/n/n/nnishi/BarcUVa/barcuvaseq.eqtls.allpairs.sumstats.hg38.txt')
eqtl.gtex <- fread('/work/users/n/n/nnishi/GTEx/full_summary_stats/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Colon_Transverse.allpairs.txt.gz')
eqtl.ibd <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/cis-eqtl.RUV.IBD.geno.corrected.MAF002.RUV.21.20230724.nominal.primary.txt')
# reformat GTEx eQTL
eqtl.gtex <- eqtl.gtex %>% na.omit()
eqtl.gtex$variant_id <- gsub('_b38', '', eqtl.gtex$variant_id)
eqtl.gtex$variant_id <- gsub('_', ':', eqtl.gtex$variant_id)
# annotate eQTL full summary stats
eqtl.barcuva <- merge(dbsnp, eqtl.barcuva, by.x = 'snp', by.y = 'variant_id', all.y = TRUE)
eqtl.gtex <- merge(dbsnp, eqtl.gtex, by.x = 'snp', by.y = 'variant_id', all.y = TRUE)
eqtl.ibd <- merge(dbsnp, eqtl.ibd, by.x = 'snp', by.y = 'var_id', all.y = TRUE)
# save
remove(dbsnp)
save.image('/work/users/n/n/nnishi/eqtl/freeze/final/RData/full_summary_stats_rsIDs.RData')