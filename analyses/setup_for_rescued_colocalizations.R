

################################################################################
# pull out GTEx and BarcUVa H4 eGenes to test for colocalization in IBD eQTL

ibd <- fread('~/Downloads/IBD.geno.corrected.cis-eqtl.MAF002.RUV.21.txt')
raw <- fread('~/Desktop/UNC/Rotations_2020-2021/Furey_S2021/eQTL/freeze/RNAseq_counts_raw_20220830.bed')
# annotate genes
raw <- raw %>% separate(info, c('length', 'type', 'position', 'name'), sep = ';',
                        remove = TRUE, convert = TRUE, extra = 'drop', fill = 'right')
raw <- raw %>% mutate(type = gsub('T=', '', type))
raw <- raw %>% mutate(name = gsub('N=', '', name))
raw <- raw %>% separate(gene, c('gene', 'ver'))
raw <- raw %>% select(gene, ver, type, name)
ibd <- ibd %>% separate(phe_id, c('gene', 'ver'))
ibd$snp <- ibd$var_id
ibd <- right_join(raw, ibd)
# drop genes not tested (aka no variants in cis)
ibd <- ibd[ibd$n_var_in_cis > 0,]
# calculate q-values for IBD eQTL
q <- qvalue(ibd$adj_beta_pval)
ibd$qval <- q$qvalues
ibd_sig <- ibd[ibd$qval < 0.05,]
ibd_non <- ibd[ibd$qval >= 0.05,]
# pull out non-IBD eGenes
ni.egenes <- c(extract_comb(mat.egenes, '011'),extract_comb(mat.egenes, '010'),extract_comb(mat.egenes, '001'))
# subset IBD data
ibd.ni <- subset(ibd, gene %in% ni.egenes)
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
# merge with IBD data
out <- merge(ibd.ni, ni, by = 'gene', all.x = TRUE)
################################################################################
write.table(out, sep = '\t', quote = FALSE, row.names = FALSE,
            file = '~/Desktop/UNC/Rotations_2020-2021/Furey_S2021/eQTL/freeze/final/IBD_eQTL_GTEx_BarcUVa_H4_eGenes_for_coloc_20230929.txt')
################################################################################