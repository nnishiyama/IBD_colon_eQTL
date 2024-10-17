# colocalize eQTL with Liu et al. GWAS
library(data.table)
library(tidyverse)
library(coloc)

# load pairs info
pairs.barcuva <- fread('data/BarcUVa/BarcUVa_eQTL_GWAS_pairs_for_coloc.txt')
pairs.gtex <- fread('data/GTEx/GTEx_eQTL_colon_transverse_GWAS_pairs_for_coloc.txt')
pairs.ibd <- fread('data/UNC/IBD_eQTL_GWAS_pairs_for_coloc.txt')
# load full summary stats
eqtl.barcuva <- fread('data/BarcUVa/barcuvaseq.eqtls.allpairs.sumstats.hg38.txt')
eqtl.gtex <- fread('data/GTEx/GTEx_Analysis_v8_eQTL_all_associations/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Colon_Transverse.allpairs.txt.gz')
eqtl.ibd <- fread('data/UNC/IBD_cis-eQTL_full_summary_stats_primary.txt.gz')
# load GWAS full summary stats
cd <- fread('data/GWAS/ibd_EAS_EUR_SiKJEF_meta_CD.TBL.txt.gz')
ibd <- fread('data/GWAS/ibd_EAS_EUR_SiKJEF_meta_IBD.TBL.txt.gz')
uc <- fread('data/GWAS/ibd_EAS_EUR_SiKJEF_meta_UC.TBL.txt.gz')
# load lead GWAS summary stats
lead <- fread('data/GWAS/Liu_2023_lead_summary_stats.txt')
# reformat pairs
pairs.barcuva <- pairs.barcuva %>% unite('gene', gene:ver, sep = '.')
pairs.gtex <- pairs.gtex %>% unite('gene', gene:ver, sep = '.')
pairs.ibd <- pairs.ibd %>% unite('gene', gene:ver, sep = '.')
# reformat lead
lead$Chr_index <- gsub('^', 'chr', lead$Chr_index)
# reformat GTEx eQTL
eqtl.gtex <- eqtl.gtex %>% na.omit()
eqtl.gtex$variant_id <- gsub('_b38', '', eqtl.gtex$variant_id)
eqtl.gtex$variant_id <- gsub('_', ':', eqtl.gtex$variant_id)
# reformat GWAS
cd <- cd %>% separate(snp, c('chr', 'pos', 'ref', 'alt'), sep = ':', remove = FALSE)
ibd <- ibd %>% separate(snp, c('chr', 'pos', 'ref', 'alt'), sep = ':', remove = FALSE)
uc <- uc %>% separate(snp, c('chr', 'pos', 'ref', 'alt'), sep = ':', remove = FALSE)
################################################################################
# BarcUVa
# assign variables to be filled
pairs.barcuva$GWAS_index_snp_slope <- NA
pairs.baxrcuva$best_trait <- NA
pairs.barcuva$BarcUVa_snps <- NA
pairs.barcuva$GWAS_snps <- NA
pairs.barcuva$shared_snps <- NA
pairs.barcuva$eSNP_in_GWAS <- NA
pairs.barcuva$GWAS_index_SNP_in_eQTL <- NA
pairs.barcuva$H1_hypothesis_eQTL <- NA
pairs.barcuva$H2_hypothesis_GWAS <- NA
pairs.barcuva$H3_hypothesis_different <- NA
pairs.barcuva$H4_hypothesis_shared <- NA
pairs.barcuva$number_snps_credible_set <- NA
pairs.barcuva$H4_snp <- NA
pairs.barcuva$coloc_likely_hypothesis <- NA
# loop through BarcUVa pairs
for (x in 1:nrow(pairs.barcuva)) {
  print(x)
  # subset eQTL data
  e <- eqtl.barcuva[eqtl.barcuva$gene_id == pairs.barcuva$gene[x],]
  e <- e %>% select(slope, slope_se, snp, pval_nominal)
  colnames(e) <- c('beta', 'varbeta', 'snp', 'pvalue')
  e$varbeta <- e$varbeta ^2
  e$sdY <- 1
  e <- e %>% distinct()
  e <- as.list(e)
  e$type <- 'quant'
  # subset GWAS lead data
  snp <- pairs.barcuva$GWAS_index_snp[x]
  l <- lead[lead$SNP_index == snp,]
  # select GWAS trait
  trait <- l$Phenotype_loci
  if (trait == 'CD'){
    g <- cd
  } else if (trait == 'IBD'){
    g <- ibd
  } else {
    g <- uc
  }
  # add trait
  pairs.barcuva$best_trait[x] <- trait
  # subset GWAS data
  g <- g[g$chr == l$Chr_index,]
  g$pos <- as.integer(g$pos)
  g <- g[g$pos >= l$BP_left_loci & g$pos <= l$BP_right_loci,]
  g <- g %>% select(Effect, StdErr, snp, `P-value`)
  colnames(g) <- c('beta', 'varbeta', 'snp', 'pvalue')
  g$varbeta <- g$varbeta ^2
  g$sdY <- 1
  g <- g %>% distinct()
  g <- as.list(g)
  g$type <- 'cc'
  # add SNP counts
  pairs.barcuva$GWAS_index_snp_slope[x] <- g$beta[which(g$snp == snp)]
  pairs.barcuva$BarcUVa_snps[x] <- length(e$snp)
  pairs.barcuva$GWAS_snps[x] <- length(g$snp)
  pairs.barcuva$shared_snps[x] <- length(intersect(e$snp, g$snp))
  if (pairs.barcuva$shared_snps[x] == 0)
    next
  # check if lead eSNPs exist in the other data set
  e.snp <- pairs.barcuva$snp[x]
  g.snp <- pairs.barcuva$GWAS_index_snp[x]
  pairs.barcuva$eSNP_in_GWAS[x] <- e.snp %in% g$snp
  pairs.barcuva$GWAS_index_SNP_in_eQTL[x] <- g.snp %in% e$snp
  # run coloc
  res <- coloc.abf(dataset1 = e, dataset2 = g)
  # pull out PPs
  pairs.barcuva$H1_hypothesis_eQTL[x] <- res$summary[[3]]
  pairs.barcuva$H2_hypothesis_GWAS[x] <- res$summary[[4]]
  pairs.barcuva$H3_hypothesis_different[x] <- res$summary[[5]]
  pairs.barcuva$H4_hypothesis_shared[x] <- res$summary[[6]]
  # order SNPs in order of H4
  o <- order(res$results$SNP.PP.H4,decreasing=TRUE)
  # pull out top H4 SNP
  pairs.barcuva$H4_snp[x] <- res$results$snp[o[1]]
  # define credible set of SNPs
  cs <- cumsum(res$results$SNP.PP.H4[o])
  pairs.barcuva$number_snps_credible_set[x] <- which(cs > 0.95)[1]
  # add most likely hypothesis
  pairs.barcuva$coloc_likely_hypothesis[x] <- names(which.max(res$summary[-1]))
}
h4.barcuva <- pairs.barcuva[pairs.barcuva$H4_hypothesis_shared > 0.5,]
length(unique(h4.barcuva$GWAS_index_snp))
################################################################################
# GTEx
# assign variables to be filled
pairs.gtex$GWAS_index_snp_slope <- NA
pairs.gtex$best_trait <- NA
pairs.gtex$GTEx_snps <- NA
pairs.gtex$GWAS_snps <- NA
pairs.gtex$shared_snps <- NA
pairs.gtex$eSNP_in_GWAS <- NA
pairs.gtex$GWAS_index_SNP_in_eQTL <- NA
pairs.gtex$H1_hypothesis_eQTL <- NA
pairs.gtex$H2_hypothesis_GWAS <- NA
pairs.gtex$H3_hypothesis_different <- NA
pairs.gtex$H4_hypothesis_shared <- NA
pairs.gtex$number_snps_credible_set <- NA
pairs.gtex$H4_snp <- NA
pairs.gtex$coloc_likely_hypothesis <- NA
# loop through GTEx pairs
for (x in 1:nrow(pairs.gtex)) {
  print(x)
  # subset IBD data
  e <- eqtl.gtex[eqtl.gtex$gene_id == pairs.gtex$gene[x],]
  e <- e %>% select(slope, slope_se, snp, pval_nominal)
  colnames(e) <- c('beta', 'varbeta', 'snp', 'pvalue')
  e$varbeta <- e$varbeta ^2
  e$sdY <- 1
  e <- e %>% distinct()
  e <- as.list(e)
  e$type <- 'quant'
  # subset GWAS lead data
  snp <- pairs.gtex$GWAS_index_snp[x]
  l <- lead[lead$SNP_index == snp,]
  # select GWAS trait
  trait <- l$Phenotype_loci
  if (trait == 'CD'){
    g <- cd
  } else if (trait == 'IBD'){
    g <- ibd
  } else {
    g <- uc
  }
  # add trait
  pairs.gtex$best_trait[x] <- trait
  # subset GWAS data
  g <- g[g$chr == l$Chr_index,]
  g$pos <- as.integer(g$pos)
  g <- g[g$pos >= l$BP_left_loci & g$pos <= l$BP_right_loci,]
  g <- g %>% select(Effect, StdErr, snp, `P-value`)
  colnames(g) <- c('beta', 'varbeta', 'snp', 'pvalue')
  g$varbeta <- g$varbeta ^2
  g$sdY <- 1
  g <- g %>% distinct()
  g <- as.list(g)
  g$type <- 'cc'
  # add SNP counts
  pairs.gtex$GWAS_index_snp_slope[x] <- g$beta[which(g$snp == snp)]
  pairs.gtex$GTEx_snps[x] <- length(e$snp)
  pairs.gtex$GWAS_snps[x] <- length(g$snp)
  pairs.gtex$shared_snps[x] <- length(intersect(e$snp, g$snp))
  if (pairs.gtex$shared_snps[x] == 0)
    next
  # check if lead eSNPs exist in the other data set
  e.snp <- pairs.gtex$snp[x]
  g.snp <- pairs.gtex$GWAS_index_snp[x]
  pairs.gtex$eSNP_in_GWAS[x] <- e.snp %in% g$snp
  pairs.gtex$GWAS_index_SNP_in_eQTL[x] <- g.snp %in% e$snp
  # run coloc
  res <- coloc.abf(dataset1 = e, dataset2 = g)
  # pull out PPs
  pairs.gtex$H1_hypothesis_eQTL[x] <- res$summary[[3]]
  pairs.gtex$H2_hypothesis_GWAS[x] <- res$summary[[4]]
  pairs.gtex$H3_hypothesis_different[x] <- res$summary[[5]]
  pairs.gtex$H4_hypothesis_shared[x] <- res$summary[[6]]
  # order SNPs in order of H4
  o <- order(res$results$SNP.PP.H4,decreasing=TRUE)
  # pull out top H4 SNP
  pairs.gtex$H4_snp[x] <- res$results$snp[o[1]]
  # define credible set of SNPs
  cs <- cumsum(res$results$SNP.PP.H4[o])
  pairs.gtex$number_snps_credible_set[x] <- which(cs > 0.95)[1]
  # add most likely hypothesis
  pairs.gtex$coloc_likely_hypothesis[x] <- names(which.max(res$summary[-1]))
}
h4.gtex <- pairs.gtex[pairs.gtex$H4_hypothesis_shared > 0.5,]
length(unique(h4.gtex$GWAS_index_snp))
################################################################################
# IBD
# assign variables to be filled
pairs.ibd$GWAS_index_snp_slope <- NA
pairs.ibd$best_trait <- NA
pairs.ibd$IBD_snps <- NA
pairs.ibd$GWAS_snps <- NA
pairs.ibd$shared_snps <- NA
pairs.ibd$eSNP_in_GWAS <- NA
pairs.ibd$GWAS_index_SNP_in_eQTL <- NA
pairs.ibd$H1_hypothesis_eQTL <- NA
pairs.ibd$H2_hypothesis_GWAS <- NA
pairs.ibd$H3_hypothesis_different <- NA
pairs.ibd$H4_hypothesis_shared <- NA
pairs.ibd$number_snps_credible_set <- NA
pairs.ibd$H4_snp <- NA
pairs.ibd$coloc_likely_hypothesis <- NA
# loop through IBD pairs
for (x in 1:nrow(pairs.ibd)) {
  print(x)
  # subset IBD data
  e <- eqtl.ibd[eqtl.ibd$phe_id == pairs.ibd$gene[x],]
  e <- e %>% select(slope, slope_se, snp, nom_pval)
  colnames(e) <- c('beta', 'varbeta', 'snp', 'pvalue')
  e$varbeta <- e$varbeta ^2
  e$sdY <- 1
  e <- e %>% distinct()
  e <- as.list(e)
  e$type <- 'quant'
  # subset GWAS lead data
  snp <- pairs.ibd$GWAS_index_snp[x]
  l <- lead[lead$SNP_index == snp,]
  # select GWAS trait
  trait <- l$Phenotype_loci
  if (trait == 'CD'){
    g <- cd
  } else if (trait == 'IBD'){
    g <- ibd
  } else {
    g <- uc
  }
  # add trait
  pairs.ibd$best_trait[x] <- trait
  # subset GWAS data
  g <- g[g$chr == l$Chr_index,]
  g$pos <- as.integer(g$pos)
  g <- g[g$pos >= l$BP_left_loci & g$pos <= l$BP_right_loci,]
  g <- g %>% select(Effect, StdErr, snp, `P-value`)
  colnames(g) <- c('beta', 'varbeta', 'snp', 'pvalue')
  g$varbeta <- g$varbeta ^2
  g$sdY <- 1
  g <- g %>% distinct()
  g <- as.list(g)
  g$type <- 'cc'
  # add SNP counts
  pairs.ibd$GWAS_index_snp_slope[x] <- g$beta[which(g$snp == snp)]
  pairs.ibd$IBD_snps[x] <- length(e$snp)
  pairs.ibd$GWAS_snps[x] <- length(g$snp)
  pairs.ibd$shared_snps[x] <- length(intersect(e$snp, g$snp))
  if (pairs.ibd$shared_snps[x] == 0)
    next
  # check if lead eSNPs exist in the other data set
  e.snp <- pairs.ibd$snp[x]
  g.snp <- pairs.ibd$GWAS_index_snp[x]
  pairs.ibd$eSNP_in_GWAS[x] <- e.snp %in% g$snp
  pairs.ibd$GWAS_index_SNP_in_eQTL[x] <- g.snp %in% e$snp
  # run coloc
  res <- coloc.abf(dataset1 = e, dataset2 = g)
  # pull out PPs
  pairs.ibd$H1_hypothesis_eQTL[x] <- res$summary[[3]]
  pairs.ibd$H2_hypothesis_GWAS[x] <- res$summary[[4]]
  pairs.ibd$H3_hypothesis_different[x] <- res$summary[[5]]
  pairs.ibd$H4_hypothesis_shared[x] <- res$summary[[6]]
  # order SNPs in order of H4
  o <- order(res$results$SNP.PP.H4,decreasing=TRUE)
  # pull out top H4 SNP
  pairs.ibd$H4_snp[x] <- res$results$snp[o[1]]
  # define credible set of SNPs
  cs <- cumsum(res$results$SNP.PP.H4[o])
  pairs.ibd$number_snps_credible_set[x] <- which(cs > 0.95)[1]
  # add most likely hypothesis
  pairs.ibd$coloc_likely_hypothesis[x] <- names(which.max(res$summary[-1]))
}
h4.ibd <- pairs.ibd[pairs.ibd$H4_hypothesis_shared > 0.5,]
length(unique(h4.ibd$GWAS_index_snp))
################################################################################
# check whether IBD H4 eQTL in GTEx/BarcUVa also colocalize with GWAS (whether or not in LD)

# GTEx
h4.ibd$shared_snps_GTEx <- NA
h4.ibd$GTEx_likely_hypothesis <- NA
h4.ibd$GTEx_PP <- NA

for (x in 1:nrow(h4.ibd)) {
  print(x)
  # subset GTEx data
  e <- eqtl.gtex[eqtl.gtex$gene_id %like% strsplit(h4.ibd$gene[x], split = '[.]')[[1]][1],]
  if (nrow(e) == 0)
    next
  e <- e %>% select(slope, slope_se, snp, pval_nominal)
  colnames(e) <- c('beta', 'varbeta', 'snp', 'pvalue')
  e$varbeta <- e$varbeta ^2
  e$sdY <- 1
  e <- e %>% distinct()
  e <- as.list(e)
  e$type <- 'quant'
  # subset GWAS lead data
  l <- lead[lead$SNP_index == h4.ibd$GWAS_index_snp[x],]
  # select GWAS trait
  trait <- l$Phenotype_loci
  if (trait == 'CD'){
    g <- cd
  } else if (trait == 'IBD'){
    g <- ibd
  } else {
    g <- uc
  }
  # subset GWAS data
  g <- g[g$chr == l$Chr_index,]
  g$pos <- as.integer(g$pos)
  g <- g[g$pos >= l$BP_left_loci & g$pos <= l$BP_right_loci,]
  g <- g %>% select(Effect, StdErr, snp, `P-value`)
  colnames(g) <- c('beta', 'varbeta', 'snp', 'pvalue')
  g$varbeta <- g$varbeta ^2
  g$sdY <- 1
  g <- g %>% distinct()
  g <- as.list(g)
  g$type <- 'cc'
  # skip if gene not annotated
  h4.ibd$shared_snps_GTEx[x] <- length(intersect(e$snp, g$snp))
  if (length(intersect(e$snp, g$snp)) == 0)
    next
  # run coloc
  res <- coloc.abf(dataset1 = e, dataset2 = g)
  # add most likely hypothesis
  h4.ibd$GTEx_likely_hypothesis[x] <- names(which.max(res$summary[-1]))
  # add PP
  h4.ibd$GTEx_PP[x] <- res$summary[[which.max(res$summary[-1])+1]]
}

# BarcUVa
h4.ibd$shared_snps_BarcUVa <- NA
h4.ibd$BarcUVa_likely_hypothesis <- NA
h4.ibd$BarcUVa_PP <- NA

for (x in 1:nrow(h4.ibd)) {
  print(x)
  # subset BarcUVa data
  e <- eqtl.barcuva[eqtl.barcuva$gene_id %like% strsplit(h4.ibd$gene[x], split = '[.]')[[1]][1],]
  if (nrow(e) == 0)
    next
  e <- e %>% select(slope, slope_se, snp, pval_nominal)
  colnames(e) <- c('beta', 'varbeta', 'snp', 'pvalue')
  e$varbeta <- e$varbeta ^2
  e$sdY <- 1
  e <- e %>% distinct()
  e <- as.list(e)
  e$type <- 'quant'
  # subset GWAS lead data
  l <- lead[lead$SNP_index == h4.ibd$GWAS_index_snp[x],]
  # select GWAS trait
  trait <- l$Phenotype_loci
  if (trait == 'CD'){
    g <- cd
  } else if (trait == 'IBD'){
    g <- ibd
  } else {
    g <- uc
  }
  # subset GWAS data
  g <- g[g$chr == l$Chr_index,]
  g$pos <- as.integer(g$pos)
  g <- g[g$pos >= l$BP_left_loci & g$pos <= l$BP_right_loci,]
  g <- g %>% select(Effect, StdErr, snp, `P-value`)
  colnames(g) <- c('beta', 'varbeta', 'snp', 'pvalue')
  g$varbeta <- g$varbeta ^2
  g$sdY <- 1
  g <- g %>% distinct()
  g <- as.list(g)
  g$type <- 'cc'
  # skip if gene not annotated
  h4.ibd$shared_snps_BarcUVa[x] <- length(intersect(e$snp, g$snp))
  if (length(intersect(e$snp, g$snp)) == 0)
    next
  # run coloc
  res <- coloc.abf(dataset1 = e, dataset2 = g)
  # add most likely hypothesis
  h4.ibd$BarcUVa_likely_hypothesis[x] <- names(which.max(res$summary[-1]))
  # add PP
  h4.ibd$BarcUVa_PP[x] <- res$summary[[which.max(res$summary[-1])+1]]
}
################################################################################
write.table(h4.barcuva, sep = '\t', quote = FALSE, row.names = FALSE,
            file = 'data/BarcUVa/Supplementary_Table_9_coloc_H4_BarcUVa_eQTL_colon_GWAS.txt')
write.table(h4.gtex, sep = '\t', quote = FALSE, row.names = FALSE,
            file = 'data/GTEx/Supplementary_Table_8_coloc_H4_GTEx_eQTL_colon_transverse_GWAS.txt')
write.table(h4.ibd, sep = '\t', quote = FALSE, row.names = FALSE,
            file = 'data/UNC/Supplementary_Table_7_coloc_h4_IBD_eQTL_GWAS.txt')
################################################################################
