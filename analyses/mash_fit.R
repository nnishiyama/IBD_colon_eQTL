# mash analysis across disease states in colon tissue
library(data.table)
library(tidyverse)
library(mashr)

set.seed(12345)

# load full sumary statistics
load('/work/users/n/n/nnishi/eqtl/freeze/final/coloc/coloc_all_data.RData')
# load top eQTL 
barcuva <- fread('data/BarcUVa/BarcUVA-seq_colon_eQTLs.csv')
gtex <- fread('data/GTEx/GTEx_Analysis_v8_eQTL/Colon_Transverse.v8.egenes.txt.gz')
ibd <- fread('data/UNC/Supplementary_Table_1_IBD_cis-eQTL_lead.significant.txt')
# drop GWAS data
remove(cd, ibd, uc)
# reformat GTEx
eqtl.gtex <- eqtl.gtex %>% na.omit()
eqtl.gtex$variant_id <- gsub('_b38', '', eqtl.gtex$variant_id)
eqtl.gtex$variant_id <- gsub('_', ':', eqtl.gtex$variant_id)
# reformat ENSEMBL IDs
eqtl.barcuva <- eqtl.barcuva %>% separate(gene_id, c('gene', 'ver'))
eqtl.gtex <- eqtl.gtex %>% separate(gene_id, c('gene', 'ver'))
eqtl.ibd <- eqtl.ibd %>% separate(phe_id, c('gene', 'ver'))
barcuva <- barcuva %>% separate(`Ensembl gene id (GENCODE v19)`, c('gene', 'ver'))
gtex <- gtex %>% separate(gene_id, c('gene','ver'))
ibd <- ibd %>% separate(phe_id, c('gene', 'ver'))
# reformat SNP
barcuva <- barcuva %>% unite('snp', BED_hg38:ALT, sep = ':', remove = FALSE)
gtex <- gtex %>% unite('snp', chr:alt, sep = ':', remove = FALSE)
ibd$snp <- ibd$var_id
# annotate eQTL
barcuva$eqtl <- paste(barcuva$gene, barcuva$snp, sep = '-')
gtex$eqtl <- paste(gtex$gene, gtex$snp, sep = '-')
ibd$eqtl <- paste(ibd$gene, ibd$snp, sep = '-')
# separate signficant eQTL
barcuva_sig <- barcuva[barcuva$qval < 0.05,]
barcuva_non <- barcuva[barcuva$qval >= 0.05,]
gtex_sig <- gtex[gtex$qval < 0.05,]
gtex_non <- gtex[gtex$qval >= 0.05,]
# list of significant eQTL
pairs.sig <- Reduce(union, list(barcuva_sig$eqtl, gtex_sig$eqtl, ibd$eqtl))
# find a common set of genes
genes <- Reduce(union, list(barcuva_sig$gene, gtex_sig$gene, ibd$gene))
# subset full summary stats
eb <- subset(eqtl.barcuva, gene %in% genes)
eg <- subset(eqtl.gtex, gene %in% genes)
ei <- subset(eqtl.ibd, gene %in% genes)
# create eSNP-eGene pairs
eb$eqtl <- paste(eb$gene, eb$variant_id, sep = '-')
eg$eqtl <- paste(eg$gene, eg$variant_id, sep = '-')
ei$eqtl <- paste(ei$gene, ei$var_id, sep = '-')
# find a common set of eSNP-eGene pairs
pairs <- Reduce(intersect, list(eb$eqtl, eg$eqtl, ei$eqtl))
# subset
eb <- subset(eb, eqtl %in% pairs) %>% arrange(match(eqtl, pairs))
eg <- subset(eg, eqtl %in% pairs) %>% arrange(match(eqtl, pairs))
ei <- subset(ei, eqtl %in% pairs) %>% arrange(match(eqtl, pairs))
# drop duplicate rows
eb <- eb[!duplicated(eb$eqtl), ]
# create a matrix of effect size estimates
bhat <- cbind(eb$slope, eg$slope, ei$slope)
colnames(bhat) <- c('barcuva', 'gtex', 'ibd')
row.names(bhat) <- pairs
# create a matrix of standard error estimates
sehat <- cbind(eb$slope_se, eg$slope_se, ei$slope_se)
colnames(sehat) <- c('barcuva', 'gtex', 'ibd')
row.names(sehat) <- pairs
# create a degrees of freedom matrix
dof <- cbind(rep(443, length(pairs)), rep(366, length(pairs)), rep(250, length(pairs)))
colnames(dof) <- c('barcuva', 'gtex', 'ibd')
row.names(dof) <- pairs
# delete unnecessary data
remove(eb, eg, ei)
# sample data to create a testing set
random.subset <- sample(1:length(pairs), 5e6)
# learn correlation structure
data.temp <- mash_set_data(Bhat = bhat[random.subset,], Shat = sehat[random.subset,], 
                           df = dof[random.subset,])
Vhat <- estimate_null_correlation_simple(data.temp)
rm(data.temp)
# create mash matrix using random data
data.random <- mash_set_data(Bhat = bhat[random.subset,], Shat = sehat[random.subset,], 
                             df = dof[random.subset,], V = Vhat)
# create canonical covariance matrix
U.c <- cov_canonical(data.random)
# find significant pairs tested in all data sets
pairs.strong <- intersect(pairs.sig, pairs)
# find indices of "strong" signals
strong <- match(pairs.strong, pairs)
# create mash matrix using strong signals
data.strong <- mash_set_data(Bhat = bhat[strong,], Shat = sehat[strong,], 
                             df = dof[strong,], V = Vhat)
# create data-drive covariance matrices
U.pca <- cov_pca(data.strong, npc = 2)
U.ed <- cov_ed(data.strong, U.pca)
# run mash on random data to create fit
m <- mash(data.random, c(U.c, U.ed), outputlevel = 1)
# apply fit to strong set
m.strong <- mash(data.strong, g = get_fitted_g(m), fixg = TRUE)
################################################################################
# save mash object
save(m, m.strong, data.random, data.strong, file = 'data/UNC/RData/mashr_colon_eqtl_random_null_corr.RData')
################################################################################
# total number of effects that mash considers significant
ind <- get_significant_results(m.strong)
# extract posterior matrices & subset to mash significant eQTL
mash.beta <- as.data.frame(m.strong$result$PosteriorMean[ind,])
mash.sd <- as.data.frame(m.strong$result$PosteriorSD[ind,])
mash.lfsr <- as.data.frame(m.strong$result$lfsr[ind,])
# reformat
mash.beta$eqtl <- row.names(mash.beta)
mash.beta <- mash.beta %>% separate(eqtl, c('gene', 'snp'), sep = '-', remove = FALSE)
# calculate the proportion of pairwise sharing between studies
get_pairwise_sharing(m.strong)
# format for writing out
colnames(mash.lfsr) <- paste0(colnames(mash.lfsr), '_lfsr')
mash.lfsr$eqtl <- row.names(mash.lfsr)
mash.sd$eqtl <- row.names(mash.sd)
out <- merge(mash.beta, mash.sd, by = 'eqtl', suffixes = c('_beta', '_sd')) %>% select(-c(gene:snp))
out <- merge(out, mash.lfsr, by = 'eqtl')
################################################################################
write.table(out, sep = '\t', quote = FALSE, row.names = FALSE, 
            file = 'data/UNC/Supplementary_Table_12_mash_top_eQTL.txt')

