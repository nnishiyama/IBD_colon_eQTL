# apply mash fit to colocalizing eQTL signals
library(data.table)
library(tidyverse)
library(mashr)

# load mash fit
load('data/UNC/RData/mashr_colon_eqtl_random_null_corr.RData')
# load summary statistics
load('/work/users/n/n/nnishi/eqtl/freeze/final/RData/full_summary_stats_rsIDs.RData')
# drop GWAS data
remove(cd, ibd, uc)
# load coloc results
coloc.barcuva <- fread('data/BarcUVa/Supplementary_Table_9_coloc_H4_BarcUVa_eQTL_GWAS.txt')
coloc.gtex <- fread('data/GTEx/Supplementary_Table_8_coloc_H4_GTEx_eQTL_colon_transverse_GWAS.txt')
coloc.ibd <- fread('data/UNC/Supplementary_Table_7_coloc_h4_IBD_eQTL_GWAS.txt')
# reformat
coloc.barcuva <- coloc.barcuva %>% separate(gene, c('gene', 'ver'))
coloc.gtex <- coloc.gtex %>% separate(gene, c('gene', 'ver'))
coloc.ibd <- coloc.ibd %>% separate(gene, c('gene', 'ver'))
# create a union list of colocalizing eGenes
genes <- unique(c(coloc.gtex$gene, coloc.ibd$gene))
# reformat ENSEMBL IDs
eqtl.barcuva <- eqtl.barcuva %>% separate(gene_id, c('gene', 'ver'))
eqtl.gtex <- eqtl.gtex %>% separate(gene_id, c('gene', 'ver'))
eqtl.ibd <- eqtl.ibd %>% separate(phe_id, c('gene', 'ver'))
# annotate eqtl
eqtl.barcuva$eqtl <- paste(eqtl.barcuva$gene, eqtl.barcuva$snp, sep = '-')
eqtl.gtex$eqtl <- paste(eqtl.gtex$gene, eqtl.gtex$snp, sep = '-')
eqtl.ibd$eqtl <- paste(eqtl.ibd$gene, eqtl.ibd$snp, sep = '-')
# subset data to gene list
eb <- subset(eqtl.barcuva, gene %in% genes) %>% select(eqtl, slope, slope_se, rsID)
eg <- subset(eqtl.gtex, gene %in% genes) %>% select(eqtl, slope, slope_se, rsID)
ei <- subset(eqtl.ibd, gene %in% genes) %>% select(eqtl, slope, slope_se, rsID)
colnames(ei) <- c('eqtl', 'slope_IBD', 'slope_se_IBD', 'rsID')
# drop unnecessary data
remove(eqtl.barcuva, eqtl.gtex, eqtl.ibd)
# merge
df.coloc <- merge(eb, eg, by = c('eqtl', 'rsID'), all = TRUE, suffixes = c('_BarcUVa', '_GTEx'))
df.coloc <- merge(df.coloc, ei, by = c('eqtl', 'rsID'), all = TRUE)
# reformat
bhat <- cbind(df.coloc$slope_BarcUVa, df.coloc$slope_GTEx, df.coloc$slope_IBD)
colnames(bhat) <- c('barcuva', 'gtex', 'ibd')
row.names(bhat) <- df.coloc$eqtl
sehat <- cbind(df.coloc$slope_se_BarcUVa, df.coloc$slope_se_GTEx, df.coloc$slope_se_IBD)
colnames(sehat) <- c('barcuva', 'gtex', 'ibd')
row.names(sehat) <- df.coloc$eqtl
dof <- cbind(rep(443, nrow(df.coloc)), rep(366, nrow(df.coloc)), rep(250, nrow(df.coloc)))
colnames(dof) <- c('barcuva', 'gtex', 'ibd')
row.names(dof) <- df.coloc$eqtl
# create a mash object
data.coloc <- mash_set_data(Bhat = bhat, Shat = sehat, df = dof, V = data.strong$V)
# apply fit
m.coloc <- mash(data.coloc, g = get_fitted_g(m), fixg = TRUE)
# total number of effects that mash considers significant
ic <- get_significant_results(m.coloc)
# significant results by condition
length(get_significant_results(m.coloc, conditions = 1))
length(get_significant_results(m.coloc, conditions = 2))
length(get_significant_results(m.coloc, conditions = 3))
# calculate the proportion of pairwise sharing between studies
get_pairwise_sharing(m.coloc)
# extract posterior matrices & subset to mash significant eQTL
coloc.beta <- as.data.frame(m.coloc$result$PosteriorMean)
coloc.sd <- as.data.frame(m.coloc$result$PosteriorSD)
coloc.lfsr <- as.data.frame(m.coloc$result$lfsr)
# reformat
coloc.beta$rsID <- df.coloc$rsID
coloc.beta$eqtl <- df.coloc$eqtl
coloc.beta <- coloc.beta %>% separate(eqtl, c('gene', 'snp'), sep = '-', remove = FALSE)
# pull out significant results for each condition
mc.b <- coloc.beta[get_significant_results(m.coloc, conditions =1),]
mc.g <- coloc.beta[get_significant_results(m.coloc, conditions =2),]
mc.i <- coloc.beta[get_significant_results(m.coloc, conditions =3),]
################################################################################
# save workspace
save.image(file = 'data/UNC/RData/mashr_coloc_fit_full_signal.RData')
################################################################################
# format for writing output
colnames(coloc.lfsr) <- paste0(colnames(coloc.lfsr), '_lfsr')
coloc.lfsr$eqtl <- coloc.beta$eqtl
coloc.sd$eqtl <- coloc.beta$eqtl
out <- merge(coloc.beta, coloc.sd, by = 'eqtl', suffixes = c('_beta', '_sd')) %>% select(-c(gene:snp))
out <- merge(out, coloc.lfsr, by = 'eqtl') %>% distinct()
################################################################################
# write out
write.table(out, sep = '\t', quote = FALSE, row.names = FALSE, 
            file = 'data/UNC/Supplementary_Table_13_mashr_coloc_fit_full_signal_colocalizing_eQTL.txt')
################################################################################