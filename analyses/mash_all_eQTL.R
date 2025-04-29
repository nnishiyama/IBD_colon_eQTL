# apply mash to all eQTL
library(data.table)
library(tidyverse)
library(mashr)
library(qvalue)
library(ggplot2)
library(ggsignif)

# load mash fit data
load('data/UNC/RData/mashr_colon_eqtl_random_null_corr.RData')
# load summary statistics
load('data/UNC/RData/full_summary_stats_rsIDs.RData')
# drop GWAS data
remove(cd, ibd, uc)
# remove duplicates
eqtl.barcuva <- eqtl.barcuva %>% distinct(gene_id, snp, .keep_all = TRUE)
eqtl.gtex <- eqtl.gtex %>% distinct(gene_id, snp, .keep_all = TRUE)
eqtl.ibd <- eqtl.ibd %>% distinct(phe_id, snp, .keep_all = TRUE)
# reformat ENSEMBL IDs
eqtl.barcuva <- eqtl.barcuva %>% separate(gene_id, c('gene', 'ver'))
eqtl.gtex <- eqtl.gtex %>% separate(gene_id, c('gene', 'ver'))
eqtl.ibd <- eqtl.ibd %>% separate(phe_id, c('gene', 'ver'))
# annotate eqtl
eqtl.barcuva$eqtl <- paste(eqtl.barcuva$gene, eqtl.barcuva$snp, sep = '-')
eqtl.gtex$eqtl <- paste(eqtl.gtex$gene, eqtl.gtex$snp, sep = '-')
eqtl.ibd$eqtl <- paste(eqtl.ibd$gene, eqtl.ibd$snp, sep = '-')
# subset data
eb <- eqtl.barcuva %>% select(eqtl, slope, slope_se, rsID)
eg <- eqtl.gtex %>% select(eqtl, slope, slope_se, rsID)
ei <- eqtl.ibd %>% select(eqtl, slope, slope_se, rsID)
colnames(ei) <- c('eqtl', 'slope_IBD', 'slope_se_IBD', 'rsID')
# merge
df.all <- merge(eb, eg, by = c('eqtl', 'rsID'), all = TRUE, suffixes = c('_BarcUVa', '_GTEx'))
df.all <- merge(df.all, ei, by = c('eqtl', 'rsID'), all = TRUE) %>% distinct(eqtl, .keep_all = TRUE)
# reformat
bhat <- cbind(df.all$slope_BarcUVa, df.all$slope_GTEx, df.all$slope_IBD)
colnames(bhat) <- c('barcuva', 'gtex', 'ibd')
row.names(bhat) <- df.all$eqtl
sehat <- cbind(df.all$slope_se_BarcUVa, df.all$slope_se_GTEx, df.all$slope_se_IBD)
colnames(sehat) <- c('barcuva', 'gtex', 'ibd')
row.names(sehat) <- df.all$eqtl
dof <- cbind(rep(443, nrow(df.all)), rep(366, nrow(df.all)), rep(250, nrow(df.all)))
colnames(dof) <- c('barcuva', 'gtex', 'ibd')
row.names(dof) <- df.all$eqtl

# create a mash object
data.all <- mash_set_data(Bhat = bhat, Shat = sehat, df = dof, V = data.strong$V)
# apply fit
m.all <- mash(data.all, g = get_fitted_g(m), fixg = TRUE)
# total number of effects that mash considers significant
ic.all <- get_significant_results(m.all)
# significant results by condition
length(get_significant_results(m.all, conditions = 1))
length(get_significant_results(m.all, conditions = 2))
length(get_significant_results(m.all, conditions = 3))
# calculate the proportion of pairwise sharing between studies
get_pairwise_sharing(m.all)
# extract posterior matrices & subset to mash significant eQTL
all.beta <- as.data.frame(m.all$result$PosteriorMean)
all.sd <- as.data.frame(m.all$result$PosteriorSD)
all.lfsr <- as.data.frame(m.all$result$lfsr)
# reformat
all.beta$rsID <- df.all$rsID
all.beta$eqtl <- df.all$eqtl
all.beta <- all.beta %>% separate(eqtl, c('gene', 'snp'), sep = '-', remove = FALSE)
# pull out significant results for each condition
ma.b <- all.beta[get_significant_results(m.all, conditions =1),]
ma.g <- all.beta[get_significant_results(m.all, conditions =2),]
ma.i <- all.beta[get_significant_results(m.all, conditions =3),]
# save workspace
save(data.all, m.all, all.beta, all.sd, all.lfsr, ma.b, ma.g, ma.i, 
     file = 'data/UNC/RData/mashr_all_fit_full_signal.RData')
# format for writing output
colnames(all.lfsr) <- paste0(colnames(all.lfsr), '_lfsr')
all.lfsr$eqtl <- all.beta$eqtl
all.sd$eqtl <- all.beta$eqtl
out.all <- merge(all.beta, all.sd, by = 'eqtl', suffixes = c('_beta', '_sd')) %>% select(-c(gene:snp))
out.all <- merge(out.all, all.lfsr, by = 'eqtl') %>% distinct()
# write out
write.table(out.all, sep = '\t', quote = FALSE, row.names = FALSE, 
            file = 'data/UNC/mashr_all_fit_full_signal_eQTL.txt')