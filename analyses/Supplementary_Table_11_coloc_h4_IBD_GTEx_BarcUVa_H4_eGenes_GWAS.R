# test GTEx and BarcUVa H4 eGenes in IBD eQTL data set
library(data.table)
library(tidyverse)
library(coloc)

# load coloc results
coloc.barcuva <- fread('data/BarcUVa/Supplementary_Table_9_coloc_H4_BarcUVa_eQTL_colon_GWAS.txt')
coloc.gtex <- fread('data/GTEx/Supplementary_Table_8_coloc_H4_GTEx_eQTL_colon_transverse_GWAS.txt')
coloc.ibd <- fread('data/UNC/Supplementary_Table_7_coloc_h4_IBD_eQTL_GWAS.txt')
# load summary stats
ss <- fread('data/UNC/IBD_cis-eQTL_RUV_21.txt')
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
# reformat 
coloc.ibd <- coloc.ibd %>% separate(gene, c('gene', 'ver'))
coloc.barcuva <- coloc.barcuva %>% separate(gene, c('gene', 'ver'))
coloc.gtex <- coloc.gtex %>% separate(gene, c('gene', 'ver'))
ss <- ss %>% separate(gene, c('gene', 'ver'))
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
# pull out GTEx and BarcUVa H4 eGenes to test for colocalization in IBD eQTL
# pull out non-IBD H4 eGenes
ni.egenes <- setdiff(c(coloc.barcuva$gene, coloc.gtex$gene), coloc.ibd$gene)
# subset IBD data
ibd.ni <- subset(ss, gene %in% ni.egenes)
# subset non-IBD results
b <- subset(coloc.barcuva, gene %in% ni.egenes)
b <- b %>% select(gene, GWAS_index_snp, shared_snps, H4_hypothesis_shared)
b <- b %>% unite('pairs', gene:GWAS_index_snp, sep = '-')
g <- subset(coloc.gtex, gene %in% ni.egenes)
g <- g %>% select(gene, GWAS_index_snp, shared_snps, H4_hypothesis_shared)
g <- g %>% unite('pairs', gene:GWAS_index_snp, sep = '-')
# merge
ni <- merge(b, g, by = 'pairs', all = TRUE)
colnames(ni) <- c('pairs', 'BarcUVa_shared_snps', 'BarcUVa_PP', 'GTEx_shared_snps', 'GTEx_PP')
ni <- ni %>% separate(pairs, c('gene', 'GWAS_index_snp'), sep = '-')
ni <- ni %>% mutate(Barcuva_likely_hypothesis = ifelse(is.na(BarcUVa_PP), NA, 'PP.H4.abf'))
ni <- ni %>% mutate(GTEx_likely_hypothesis = ifelse(is.na(GTEx_PP), NA, 'PP.H4.abf'))
# merge with IBD summary stats
pairs.ibd.ni <- merge(ibd.ni, ni, by = 'gene', all.x = TRUE)
pairs.ibd.ni <- pairs.ibd.ni %>% unite('gene', c(gene, ver), sep = '.')
# loop through pairs for coloc
for (x in 1:nrow(pairs.ibd.ni)) {
  print(x)
  # subset IBD data
  e <- eqtl.ibd[eqtl.ibd$phe_id == pairs.ibd.ni$gene[x],]
  e <- e %>% select(slope, slope_se, var_id, nom_pval)
  colnames(e) <- c('beta', 'varbeta', 'snp', 'pvalue')
  e$varbeta <- e$varbeta ^2
  e$sdY <- 1
  e <- as.list(e)
  e$type <- 'quant'
  # subset GWAS lead data
  snp <- pairs.ibd.ni$GWAS_index_snp[x]
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
  pairs.ibd.ni$best_trait[x] <- trait
  # subset GWAS data
  g <- g[g$chr == l$Chr_index,]
  g$pos <- as.integer(g$pos)
  g <- g[g$pos >= l$BP_left_loci & g$pos <= l$BP_right_loci,]
  g <- g %>% select(Effect, StdErr, MarkerName, `P-value`)
  colnames(g) <- c('beta', 'varbeta', 'snp', 'pvalue')
  g$varbeta <- g$varbeta ^2
  g$sdY <- 1
  g <- as.list(g)
  g$type <- 'cc'
  # add SNP counts
  pairs.ibd.ni$IBD_snps[x] <- length(e$snp)
  pairs.ibd.ni$GWAS_snps[x] <- length(g$snp)
  pairs.ibd.ni$shared_snps[x] <- length(intersect(e$snp, g$snp))
  if (pairs.ibd.ni$shared_snps[x] == 0)
    next
  # check if lead eSNPs exist in the other data set
  e.snp <- pairs.ibd.ni$snp[x]
  g.snp <- pairs.ibd.ni$GWAS_index_snp[x]
  pairs.ibd.ni$eSNP_in_GWAS[x] <- e.snp %in% g$snp
  pairs.ibd.ni$GWAS_index_SNP_in_eQTL[x] <- g.snp %in% e$snp
  # run coloc
  res <- coloc.abf(dataset1 = e, dataset2 = g)
  # pull out PPs
  pairs.ibd.ni$H1_hypothesis_eQTL[x] <- res$summary[[3]]
  pairs.ibd.ni$H2_hypothesis_GWAS[x] <- res$summary[[4]]
  pairs.ibd.ni$H3_hypothesis_different[x] <- res$summary[[5]]
  pairs.ibd.ni$H4_hypothesis_shared[x] <- res$summary[[6]]
  # order SNPs in order of H4
  o <- order(res$results$SNP.PP.H4,decreasing=TRUE)
  # pull out top H4 SNP
  pairs.ibd.ni$H4_snp[x] <- res$results$snp[o[1]]
  # define credible set of SNPs
  cs <- cumsum(res$results$SNP.PP.H4[o])
  pairs.ibd.ni$number_snps_credible_set[x] <- which(cs > 0.95)[1]
  # add most likely hypothesis
  pairs.ibd.ni$coloc_likely_hypothesis[x] <- names(which.max(res$summary[-1]))
}
table(pairs.ibd.ni$coloc_likely_hypothesis)
h4 <- pairs.ibd.ni[pairs.ibd.ni$coloc_likely_hypothesis == 'PP.H4.abf',]
################################################################################
write.table(h4, sep ='\t', quote = FALSE, row.names = FALSE,
            file = 'data/UNC/Supplementary_Table_11_coloc_h4_IBD_GTEx_BarcUVa_H4_eGenes_GWAS.txt')
################################################################################
