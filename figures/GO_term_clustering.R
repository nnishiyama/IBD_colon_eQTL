# pull colocalizing genes for a high level GO analysis
library(data.table)
library(tidyverse)

# load data
coloc.barcuva <- fread('/work/users/n/n/nnishi/gwas/liu/coloc_H4_BarcUVa_eQTL_colon_GWAS_Liu_all_pop_20240111.txt')
coloc.gtex <- fread('/work/users/n/n/nnishi/gwas/liu/coloc_H4_GTEx_eQTL_colon_transverse_GWAS_Liu_all_pop_20240111.txt')
coloc.ibd <- fread('/work/users/n/n/nnishi/gwas/liu/coloc_h4_IBD_MAF002_GWAS_Liu_all_pop_20240111.txt')
# drop NAs
coloc.barcuva <- coloc.barcuva %>% drop_na(GWAS_index_snp)
# reformat
coloc.barcuva <- coloc.barcuva %>% separate(gene, c('gene', 'ver'))
coloc.gtex <- coloc.gtex %>% separate(gene, c('gene', 'ver'))
coloc.ibd <- coloc.ibd %>% separate(gene, c('gene', 'ver'))
# create gene list
genes <- Reduce(union, list(coloc.barcuva$gene, coloc.gtex$gene, coloc.ibd$gene))
# create a gene table
barcuva.genes <- coloc.barcuva %>% select(GWAS_index_snp, gene, gene_name)
gtex.genes <- coloc.gtex %>% select(GWAS_index_snp, gene, gene_name)
ibd.genes <- coloc.ibd %>% select(GWAS_index_snp, gene, name)
colnames(barcuva.genes) <- c('GWAS_index_snp', 'ENSEMBL_ID', 'gene_symbol')
colnames(gtex.genes) <- c('GWAS_index_snp', 'ENSEMBL_ID', 'gene_symbol')
colnames(ibd.genes) <- c('GWAS_index_snp', 'ENSEMBL_ID', 'gene_symbol')
barcuva.genes$group <- 'BarcUVa'
gtex.genes$group <- 'GTEx'
ibd.genes$group <- 'IBD'
df <- full_join(barcuva.genes, gtex.genes)
df <- full_join(df, ibd.genes)
# test pivot
test <- df %>% pivot_wider(names_from = group, values_from = gene_symbol)
hist(table(test$GWAS_index_snp), breaks = seq(0,12), main = 'Number of colocalizing eGenes per loci')
table(test$GWAS_index_snp)[order(table(test$GWAS_index_snp))]
# pull out unique eGenes
barcuva.unique <- setdiff(coloc.barcuva$gene, c(coloc.gtex$gene, coloc.ibd$gene))
gtex.unique <- setdiff(coloc.gtex$gene, c(coloc.barcuva$gene, coloc.ibd$gene))
ibd.unique <- setdiff(coloc.ibd$gene, c(coloc.barcuva$gene, coloc.gtex$gene))
# write out unique eGene symbols
write.table(unique(df$gene_symbol), row.names = FALSE, quote = FALSE, col.names = FALSE,
            file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/All_colocalizing_eGenes_for_GO.txt')
write.table(coloc.ibd$name, row.names = FALSE, quote = FALSE, col.names = FALSE,
            file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/IBD_colocalizing_eGenes_all_for_GO.txt')
write.table(unique(df$ENSEMBL_ID), row.names = FALSE, quote = FALSE, col.names = FALSE,
            file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/All_colocalizing_eGenes_ENSEMBL_ID_for_GO.txt')
write.table(coloc.ibd$gene, row.names = FALSE, quote = FALSE, col.names = FALSE, 
            file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/IBD_colocalizing_eGenes_all_ENSEMBL_ID_for_GO.txt')
write.table(coloc.gtex$gene, row.names = FALSE, quote = FALSE, col.names = FALSE, 
            file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/GTEx_colocalizing_eGenes_all_ENSEMBL_ID_for_GO.txt')
write.table(coloc.barcuva$gene, row.names = FALSE, quote = FALSE, col.names = FALSE, 
            file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/BarcUVa_colocalizing_eGenes_all_ENSEMBL_ID_for_GO.txt')
write.table(barcuva.unique, row.names = FALSE, quote = FALSE, col.names = FALSE, 
            file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/BarcUVa_colocalizing_eGenes_unique_ENSEMBL_ID_for_GO.txt')
write.table(gtex.unique, row.names = FALSE, quote = FALSE, col.names = FALSE, 
            file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/GTEx_colocalizing_eGenes_unique_ENSEMBL_ID_for_GO.txt')
write.table(ibd.unique, row.names = FALSE, quote = FALSE, col.names = FALSE, 
            file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/IBD_colocalizing_eGenes_unique_ENSEMBL_ID_for_GO.txt')
################################################################################
# # GO term enrichment analysis
# library(gprofiler2)
# 
# # all eGenes
# go.res <- gost(list(all_eGenes = unique(df$ENSEMBL_ID)), evcodes = TRUE)
# publish_gostplot(gostplot(go.res, interactive = FALSE), highlight_terms = go.res$result$term_id)
# # IBD eGenes
# go.ibd <- gost(list(IBD_eGenes = coloc.ibd$gene), evcodes = TRUE)
# publish_gostplot(gostplot(go.ibd, interactive = FALSE), highlight_terms = go.ibd$result$term_id)
# View(go.ibd$result)
# # GTEx eGenes
# go.gtex <- gost(list(GTEx_eGenes = coloc.gtex$gene), evcodes = TRUE)
# publish_gostplot(gostplot(go.gtex, interactive = FALSE), highlight_terms = go.gtex$result$term_id)
# View(go.gtex$result)
# # BarcUVa eGenes
# go.barcuva <- gost(list(BarcUVa_eGenes = coloc.barcuva$gene), evcodes = TRUE)
# publish_gostplot(gostplot(go.barcuva, interactive = FALSE), highlight_terms = go.barcuva$result$term_id)
# View(go.barcuva$result)
################################################################################
# # another package
# library(clusterProfiler)
# ego.res <- enrichGO(unique(df$ENSEMBL_ID), ont = 'BP', 'org.Hs.eg.db', keyType = 'ENSEMBL', qvalueCutoff = 0.1, pvalueCutoff = 0.1)
# View(ego.res@result)
# ego.ibd <- enrichGO(coloc.ibd$gene, ont = 'BP', 'org.Hs.eg.db', keyType = 'ENSEMBL', qvalueCutoff = 0.1, pvalueCutoff = 0.1)
# View(ego.ibd@result)
# goplot(ego.ibd)
# 
# ego.gtex <- enrichGO(coloc.gtex$gene, ont = 'BP', 'org.Hs.eg.db', keyType = 'ENSEMBL', qvalueCutoff = 0.1, pvalueCutoff = 0.1)
# View(ego.gtex@result)
# 
# ego.barcuva <- enrichGO(coloc.barcuva$gene, ont = 'BP', 'org.Hs.eg.db', keyType = 'ENSEMBL', qvalueCutoff = 0.1, pvalueCutoff = 0.1)
# View(ego.barcuva@result)
################################################################################
# pull out GO term IDs
################################################################################
################################################################################
# # cluster GO terms into higher order categories
# library(data.table)
# library(tidyverse)
# library(rrvgo)
# # load GO terms from DAVID
# david.all <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/All_colocalizing_eGenes_ENSEMBL_ID_GO_terms_DAVID.txt')
# david.barcuva <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/BarcUVa_colocalizing_eGenes_all_ENSEMBL_ID_GO_terms_DAVID.txt')
# david.gtex <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/GTEx_colocalizing_eGenes_all_ENSEMBL_ID_GO_terms_DAVID.txt')
# david.ibd <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/IBD_colocalizing_eGenes_all_ENSEMBL_ID_GO_terms_DAVID.txt')
# david.barcuva.unique <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/BarcUVa_colocalizing_eGenes_unique_ENSEMBL_ID_GO_terms_DAVID.txt')
# david.gtex.unique <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/GTEx_colocalizing_eGenes_unique_ENSEMBL_ID_GO_terms_DAVID.txt')
# david.ibd.unique <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/IBD_colocalizing_eGenes_unique_ENSEMBL_ID_GO_terms_DAVID.txt')
# # reformat
# david.all <- david.all %>% mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ',GO:')) %>% 
#   unnest(GOTERM_BP_DIRECT) %>% separate(GOTERM_BP_DIRECT, c('GO_ID', 'GO_term'), sep = '~') %>% 
#   mutate(GO_ID = gsub('GO:', '', GO_ID)) %>% mutate(GO_ID = gsub('^', 'GO:', GO_ID))
# david.barcuva <- david.barcuva %>% mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ',GO:')) %>% 
#   unnest(GOTERM_BP_DIRECT) %>% separate(GOTERM_BP_DIRECT, c('GO_ID', 'GO_term'), sep = '~') %>% 
#   mutate(GO_ID = gsub('GO:', '', GO_ID)) %>% mutate(GO_ID = gsub('^', 'GO:', GO_ID))
# david.gtex <- david.gtex %>% mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ',GO:')) %>% 
#   unnest(GOTERM_BP_DIRECT) %>% separate(GOTERM_BP_DIRECT, c('GO_ID', 'GO_term'), sep = '~') %>% 
#   mutate(GO_ID = gsub('GO:', '', GO_ID)) %>% mutate(GO_ID = gsub('^', 'GO:', GO_ID))
# david.ibd <- david.ibd %>% mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ',GO:')) %>% 
#   unnest(GOTERM_BP_DIRECT) %>% separate(GOTERM_BP_DIRECT, c('GO_ID', 'GO_term'), sep = '~') %>% 
#   mutate(GO_ID = gsub('GO:', '', GO_ID)) %>% mutate(GO_ID = gsub('^', 'GO:', GO_ID))
# david.barcuva.unique <- david.barcuva.unique %>% mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ',GO:')) %>% 
#   unnest(GOTERM_BP_DIRECT) %>% separate(GOTERM_BP_DIRECT, c('GO_ID', 'GO_term'), sep = '~') %>% 
#   mutate(GO_ID = gsub('GO:', '', GO_ID)) %>% mutate(GO_ID = gsub('^', 'GO:', GO_ID))
# david.gtex.unique <- david.gtex.unique %>% mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ',GO:')) %>% 
#   unnest(GOTERM_BP_DIRECT) %>% separate(GOTERM_BP_DIRECT, c('GO_ID', 'GO_term'), sep = '~') %>% 
#   mutate(GO_ID = gsub('GO:', '', GO_ID)) %>% mutate(GO_ID = gsub('^', 'GO:', GO_ID))
# david.ibd.unique <- david.ibd.unique %>% mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ',GO:')) %>% 
#   unnest(GOTERM_BP_DIRECT) %>% separate(GOTERM_BP_DIRECT, c('GO_ID', 'GO_term'), sep = '~') %>% 
#   mutate(GO_ID = gsub('GO:', '', GO_ID)) %>% mutate(GO_ID = gsub('^', 'GO:', GO_ID))
# # calculate similarity scores
# simmat.all <- calculateSimMatrix(david.all$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
# simmat.b <- calculateSimMatrix(david.barcuva$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
# simmat.g <- calculateSimMatrix(david.gtex$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
# simmat.i <- calculateSimMatrix(david.ibd$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
# simmat.b.unique <- calculateSimMatrix(david.barcuva.unique$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
# simmat.g.unique <- calculateSimMatrix(david.gtex.unique$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
# simmat.i.unique <- calculateSimMatrix(david.ibd.unique$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
# # test out different IC methods
# # default method is Resnik
# simmat.wang <- calculateSimMatrix(david.all$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
# simmat.resnick <- calculateSimMatrix(david.all$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Resnik')
# simmat.lin <- calculateSimMatrix(david.all$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Lin')
# simmat.rel <- calculateSimMatrix(david.all$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Rel')
# simmat.jiang <- calculateSimMatrix(david.all$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Jiang')
# red.wang <- reduceSimMatrix(simmat.wang, threshold=0.99, orgdb="org.Hs.eg.db")
# red.resnick <- reduceSimMatrix(simmat.resnick, threshold=0.99, orgdb="org.Hs.eg.db")
# red.lin <- reduceSimMatrix(simmat.lin, threshold=0.99, orgdb="org.Hs.eg.db")
# red.rel <- reduceSimMatrix(simmat.rel, threshold=0.99, orgdb="org.Hs.eg.db")
# red.jiang <- reduceSimMatrix(simmat.jiang, threshold=0.99, orgdb="org.Hs.eg.db")
# max(red.wang$cluster) # lowest
# max(red.resnick$cluster)
# max(red.lin$cluster)
# max(red.rel$cluster)
# max(red.jiang$cluster) # way too many
# table(red.wang$parentTerm)
# table(red.resnick$parentTerm)
# table(red.lin$parentTerm)
# table(red.rel$parentTerm)
# table(red.jiang$parentTerm)
# # group terms based on similarity
# red.all <- reduceSimMatrix(simmat.all, threshold=0.99, orgdb="org.Hs.eg.db")
# red.b <- reduceSimMatrix(simmat.b, threshold=0.98, orgdb="org.Hs.eg.db")
# red.g <- reduceSimMatrix(simmat.g, threshold=0.98, orgdb="org.Hs.eg.db")
# red.i <- reduceSimMatrix(simmat.i, threshold=0.98, orgdb="org.Hs.eg.db")
# red.b.unique <- reduceSimMatrix(simmat.b.unique, threshold=0.98, orgdb="org.Hs.eg.db")
# red.g.unique <- reduceSimMatrix(simmat.g.unique, threshold=0.98, orgdb="org.Hs.eg.db")
# red.i.unique <- reduceSimMatrix(simmat.i.unique, threshold=0.98, orgdb="org.Hs.eg.db")
# # how many clusters
# range(red.all$cluster)
# range(red.b$cluster)
# range(red.g$cluster)
# range(red.i$cluster)
# range(red.b.unique$cluster)
# range(red.g.unique$cluster)
# range(red.i.unique$cluster)
# # parent cluster terms
# table(red.all$parentTerm)
# table(red.b$parentTerm)
# table(red.g$parentTerm)
# table(red.i$parentTerm)
# table(red.b.unique$parentTerm)
# table(red.g.unique$parentTerm)
# table(red.i.unique$parentTerm)
# # visualize grouping
# scatterPlot(simmat.all, red.all)
# #file = paste('/work/users/n/n/nnishi/eqtl/freeze/final/coloc/h3_eqtl/H3_', gene, '_eQTL.png', sep = '')
# png(filename = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/All_colocalizing_eGenes_GO_clustering.png', res = 300, units = 'in', height = 6, width = 6)
# treemapPlot(red.all)
# dev.off()
# png(filename = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/BarcUVa_colocalizing_eGenes_all_GO_clustering.png', res = 300, units = 'in', height = 6, width = 6)
# treemapPlot(red.b)
# dev.off()
# png(filename = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/GTEx_colocalizing_eGenes_all_GO_clustering.png', res = 300, units = 'in', height = 6, width = 6)
# treemapPlot(red.g)
# dev.off()
# png(filename = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/IBD_colocalizing_eGenes_all_GO_clustering.png', res = 300, units = 'in', height = 6, width = 6)
# treemapPlot(red.i)
# dev.off()
# png(filename = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/BarcUVa_colocalizing_eGenes_unique_GO_clustering.png', res = 300, units = 'in', height = 6, width = 6)
# treemapPlot(red.b.unique)
# dev.off()
# png(filename = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/GTEx_colocalizing_eGenes_unique_GO_clustering.png', res = 300, units = 'in', height = 6, width = 6)
# treemapPlot(red.g.unique)
# dev.off()
# png(filename = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/IBD_colocalizing_eGenes_unique_GO_clustering.png', res = 300, units = 'in', height = 6, width = 6)
# treemapPlot(red.i.unique)
# dev.off()
# # reformat
# red.all <- red.all %>% select(-term) %>% relocate(parentTerm, .after = parent)
# red.b <- red.b %>% select(-term) %>% relocate(parentTerm, .after = parent)
# red.g <- red.g %>% select(-term) %>% relocate(parentTerm, .after = parent)
# red.i <- red.i %>% select(-term) %>% relocate(parentTerm, .after = parent)
# red.b.unique <- red.b.unique %>% select(-term) %>% relocate(parentTerm, .after = parent)
# red.g.unique <- red.g.unique %>% select(-term) %>% relocate(parentTerm, .after = parent)
# red.i.unique <- red.i.unique %>% select(-term) %>% relocate(parentTerm, .after = parent)
# # map parent terms back to genes
# david.all <- merge(david.all, red.all, by.x = 'GO_ID', by.y = 'go', all.x = TRUE) %>% relocate(GO_ID, .before = GO_term) %>% arrange(cluster, ID)
# david.barcuva <- merge(david.barcuva, red.b, by.x = 'GO_ID', by.y = 'go', all.x = TRUE) %>% relocate(GO_ID, .before = GO_term) %>% arrange(cluster, ID)
# david.gtex <- merge(david.gtex, red.g, by.x = 'GO_ID', by.y = 'go', all.x = TRUE) %>% relocate(GO_ID, .before = GO_term) %>% arrange(cluster, ID)
# david.ibd <- merge(david.ibd, red.i, by.x = 'GO_ID', by.y = 'go', all.x = TRUE) %>% relocate(GO_ID, .before = GO_term) %>% arrange(cluster, ID)
# david.barcuva.unique <- merge(david.barcuva.unique, red.b.unique, by.x = 'GO_ID', by.y = 'go', all.x = TRUE) %>% relocate(GO_ID, .before = GO_term) %>% arrange(cluster, ID)
# david.gtex.unique <- merge(david.gtex.unique, red.g.unique, by.x = 'GO_ID', by.y = 'go', all.x = TRUE) %>% relocate(GO_ID, .before = GO_term) %>% arrange(cluster, ID)
# david.ibd.unique <- merge(david.ibd.unique, red.i.unique, by.x = 'GO_ID', by.y = 'go', all.x = TRUE) %>% relocate(GO_ID, .before = GO_term) %>% arrange(cluster, ID)
# # write out results
# write.table(david.all, sep = '\t', quote = FALSE, row.names = FALSE,
#             file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/All_colocalizing_eGenes_GO_clustering.txt')
# write.table(david.barcuva, sep = '\t', quote = FALSE, row.names = FALSE,
#             file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/BarcUVa_colocalizing_eGenes_all_GO_clustering.txt')
# write.table(david.gtex, sep = '\t', quote = FALSE, row.names = FALSE,
#             file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/GTEx_colocalizing_eGenes_all_GO_clustering.png.txt')
# write.table(david.ibd, sep = '\t', quote = FALSE, row.names = FALSE,
#             file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/IBD_colocalizing_eGenes_all_GO_clustering.txt')
# write.table(david.barcuva.unique, sep = '\t', quote = FALSE, row.names = FALSE,
#             file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/BarcUVa_colocalizing_eGenes_unique_GO_clustering.txt')
# write.table(david.gtex.unique, sep = '\t', quote = FALSE, row.names = FALSE,
#             file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/GTEx_colocalizing_eGenes_unique_GO_clustering.txt')
# write.table(david.ibd.unique, sep = '\t', quote = FALSE, row.names = FALSE,
#             file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/IBD_colocalizing_eGenes_unique_GO_clustering.txt')
################################################################################
# # other potential libraries
# library(clusterProfiler)
# library(simplifyEnrichment)
# 
# # simplifyEnrichment
# set.seed(888)
# go_id = random_GO(500)
# mat = GO_similarity(go_id, ont = 'BP')
# mat <- GO_similarity(david.all$GO_ID, ont = 'BP')
# df = simplifyGO(mat)
# 
# compare_clustering_methods(mat)
# cluster_terms(mat, method = "binary_cut")
# se <- simplifyEnrichment(mat)
################################################################################
# recursively cluster GO terms to create own hierarchy
library(data.table)
library(tidyverse)
library(rrvgo)
# load GO terms from DAVID
david.all <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/All_colocalizing_eGenes_ENSEMBL_ID_GO_terms_DAVID.txt')
# reformat into long format
david.all <- david.all %>% mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ',GO:')) %>% 
  unnest(GOTERM_BP_DIRECT) %>% separate(GOTERM_BP_DIRECT, c('GO_ID', 'GO_term'), sep = '~') %>% 
  mutate(GO_ID = gsub('GO:', '', GO_ID)) %>% mutate(GO_ID = gsub('^', 'GO:', GO_ID))
# summary stats
range(table(david.all$ID))
median(table(david.all$ID))
# plot number of GO terms per eGene
png(res = 300, units = 'in', height = 6, width = 6,
    filename = '/work/users/n/n/nnishi/eqtl/freeze/final/figures/Supplementary_Figure7_histogram_GO_terms_per_eGene.png')
hist(table(david.all$ID), xlab = 'number of GO terms associated with eGene', 
     main = '')
dev.off()
# calculate similarity scores
simmat.all <- calculateSimMatrix(david.all$GO_ID, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
# group terms based on similarity
red.all <- reduceSimMatrix(simmat.all, threshold=0.99, orgdb="org.Hs.eg.db")
table(red.all$parentTerm)
# drop primary GO terms
red.filtered <- subset(red.all, !(term %in% unique(red.all$parentTerm)))
# pull out clusters
clu.1 <- subset(red.filtered, cluster == 1)
clu.2 <- subset(red.filtered, cluster == 2)
clu.3 <- subset(red.filtered, cluster == 3)
clu.4 <- subset(red.filtered, cluster == 4)
clu.5 <- subset(red.filtered, cluster == 5)
clu.6 <- subset(red.filtered, cluster == 6)

# calculate similarity scores within each cluster
simmat.1 <- calculateSimMatrix(clu.1$go, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
simmat.2 <- calculateSimMatrix(clu.2$go, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
simmat.3 <- calculateSimMatrix(clu.3$go, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
simmat.4 <- calculateSimMatrix(clu.4$go, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
simmat.5 <- calculateSimMatrix(clu.5$go, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
simmat.6 <- calculateSimMatrix(clu.6$go, 'org.Hs.eg.db', ont = 'BP', method = 'Wang')
# group secondary terms
t <- 0.90
red.1 <- reduceSimMatrix(simmat.1, threshold=t, orgdb="org.Hs.eg.db")
red.2 <- reduceSimMatrix(simmat.2, threshold=t, orgdb="org.Hs.eg.db")
red.3 <- reduceSimMatrix(simmat.3, threshold=t, orgdb="org.Hs.eg.db")
red.4 <- reduceSimMatrix(simmat.4, threshold=t, orgdb="org.Hs.eg.db")
red.5 <- reduceSimMatrix(simmat.5, threshold=t, orgdb="org.Hs.eg.db")
red.6 <- reduceSimMatrix(simmat.6, threshold=t, orgdb="org.Hs.eg.db")

max(red.1$cluster)
max(red.2$cluster)
max(red.3$cluster)
max(red.4$cluster)
max(red.5$cluster)
max(red.6$cluster)

unique(red.1$parentTerm)
unique(red.2$parentTerm)
unique(red.3$parentTerm)
unique(red.4$parentTerm)
unique(red.5$parentTerm)
unique(red.6$parentTerm)
# merge results back together
red.merge <- do.call('rbind', list(red.1, red.2, red.3, red.4, red.5, red.6))
# rename columns to merge back with primary terms
colnames(red.merge) <- gsub('parent', 'secondary', colnames(red.merge))
colnames(red.merge) <- gsub('cluster', 'secondaryCluster', colnames(red.merge))
colnames(red.all) <- gsub('cluster', 'primaryCluster', colnames(red.all))
colnames(red.all) <- gsub('parent', 'primary', colnames(red.all))
# merge
red.out <- full_join(red.merge, red.all)
# rearrange columns & rows
red.out <- red.out %>% relocate(term, .after = go) %>% 
  relocate(c(secondarySimScore,secondaryCluster, secondary, secondaryTerm, primarySimScore), .after = size) %>%
  arrange(primaryCluster, secondaryCluster)
# map back to gene IDs
red.out <- merge(david.all, red.out, by.x = 'GO_ID', by.y = 'go', all.x = TRUE) %>% 
  select(-c(term)) %>% relocate(GO_ID, .before = GO_term) %>% arrange(primaryCluster, secondaryCluster, ID)
# write out results
write.table(red.out, sep = '\t', quote = FALSE, row.names = FALSE,
            file = '/work/users/n/n/nnishi/eqtl/freeze/final/GO_analysis/All_colocalizing_eGenes_GO_clustering_hierarchy_threshold_0.9.txt')
################################################################################
# create pie/donut charts based off of GO clustering
library(ggplot2)
library(webr)
# load data
coloc.barcuva <- fread('/work/users/n/n/nnishi/gwas/liu/coloc_H4_BarcUVa_eQTL_colon_GWAS_Liu_all_pop_20240111.txt')
coloc.gtex <- fread('/work/users/n/n/nnishi/gwas/liu/coloc_H4_GTEx_eQTL_colon_transverse_GWAS_Liu_all_pop_20240111.txt')
coloc.ibd <- fread('/work/users/n/n/nnishi/gwas/liu/coloc_h4_IBD_MAF002_GWAS_Liu_all_pop_20240111.txt')
# drop NAs
coloc.barcuva <- coloc.barcuva %>% drop_na(GWAS_index_snp)
# reformat
coloc.barcuva <- coloc.barcuva %>% separate(gene, c('gene', 'ver'))
coloc.gtex <- coloc.gtex %>% separate(gene, c('gene', 'ver'))
coloc.ibd <- coloc.ibd %>% separate(gene, c('gene', 'ver'))
# unnest ENSEMBL IDs
red.out <- red.out %>% mutate(ID = strsplit(as.character(ID), ',')) %>% unnest(ID)
# map GO clusters back to eQTL data sets
coloc.barcuva <- merge(coloc.barcuva, red.out, by.x = 'gene', by.y = 'ID', all.x = TRUE)
coloc.gtex <- merge(coloc.gtex, red.out, by.x = 'gene', by.y = 'ID', all.x = TRUE)
coloc.ibd <- merge(coloc.ibd, red.out, by.x = 'gene', by.y = 'ID', all.x = TRUE)
# relabel
coloc.barcuva <- coloc.barcuva %>% mutate(secondaryTerm = ifelse(GO_ID == primary, primaryTerm, secondaryTerm))
coloc.gtex <- coloc.gtex %>% mutate(secondaryTerm = ifelse(GO_ID == primary, primaryTerm, secondaryTerm))
coloc.ibd <- coloc.ibd %>% mutate(secondaryTerm = ifelse(GO_ID == primary, primaryTerm, secondaryTerm))
# label unassigned genes
coloc.barcuva <- coloc.barcuva %>% replace_na(list(primaryCluster = 0, primaryTerm = 'unassigned',
                                 secondaryCluster = 0, secondaryTerm = 'unassigned'))
coloc.gtex <- coloc.gtex %>% replace_na(list(primaryCluster = 0, primaryTerm = 'unassigned',
                                                   secondaryCluster = 0, secondaryTerm = 'unassigned'))
coloc.ibd <- coloc.ibd %>% replace_na(list(primaryCluster = 0, primaryTerm = 'unassigned',
                                                   secondaryCluster = 0, secondaryTerm = 'unassigned'))
# create a summary table quantifying the number of genes mapping to each category
# BarcUVa
b <- coloc.barcuva %>% group_by(primaryTerm, secondaryTerm) %>% summarize(n = n_distinct(gene)) %>% arrange(primaryTerm, secondaryTerm)
colnames(b)[1] <- 'BarcUVa'
b$BarcUVa <- gsub('positive regulation of cell population proliferation', 
                  'positive regulation of cell\npopulation proliferation', b$BarcUVa)
b$BarcUVa <- gsub('regulation of transcription by RNA polymerase II', 
                  'regulation of transcription\nby RNA polymerase II', b$BarcUVa)
b$secondaryTerm <- gsub('cell adhesion', 'cell adhesion\n', b$secondaryTerm)
b$secondaryTerm <- gsub('nervous system development', 
                        'nervous\nsystem\ndevelopment', b$secondaryTerm)
b$secondaryTerm <- gsub('cholesterol metabolic process', '\n\ncholesterol metabolic process', b$secondaryTerm)
b$secondaryTerm <- gsub('electron transport chain', '\nelectron transport chain', b$secondaryTerm)
b$secondaryTerm <- gsub('oligosaccharide biosynthetic process',
                        '\n\n\noligosaccharide biosynthetic process', b$secondaryTerm)
b$secondaryTerm <- gsub('positive regulation of gene expression', 
                        'positive regulation of\ngene expression', b$secondaryTerm)
b$secondaryTerm <- gsub('positive regulation of transcription by RNA polymerase II', 
                        'positive regulation of transcription\nby RNA polymerase II', b$secondaryTerm)
b$secondaryTerm <- gsub('proteasome-mediated ubiquitin-dependent protein catabolic process', 
                        'proteasome-mediated ubiquitin-\ndependent protein catabolic process', b$secondaryTerm)
b$secondaryTerm <- gsub('regulation of transcription by RNA polymerase II', 
                        'regulation of transcription\nby RNA polymerase II', b$secondaryTerm)
b$secondaryTerm <- gsub('complement activation, classical pathway', 
                        'complement activation,\nclassical pathway', b$secondaryTerm)
b$secondaryTerm <- gsub('cytokine-mediated signaling pathway', 
                        'cytokine-mediated\nsignaling pathway', b$secondaryTerm)
b$secondaryTerm <- gsub('G protein-coupled receptor signaling pathway', 
                        'G protein-coupled receptor\nsignaling pathway', b$secondaryTerm)
# pie-donut chart
png(res = 300, units = 'in', height = 12, width = 12,
    filename = '/work/users/n/n/nnishi/eqtl/freeze/final/figures/pie-donut_GO_clustering_BarcUVa_all_eGenes.png')
PieDonut(b, aes(BarcUVa, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.03, 
         titlesize = 10)
dev.off()
# GTEx
g <- coloc.gtex %>% group_by(primaryTerm, secondaryTerm) %>% summarize(n = n_distinct(gene)) %>% arrange(primaryTerm, secondaryTerm)
colnames(g)[1] <- 'GTEx'
g$GTEx <- gsub('positive regulation of cell population proliferation', 
               'positive regulation of cell\npopulation proliferation', g$GTEx)
g$GTEx <- gsub('regulation of transcription by RNA polymerase II', 
               'regulation of transcription\nby RNA polymerase II', g$GTEx)
g$secondaryTerm <- gsub('cell adhesion', 'cell adhesion\n', g$secondaryTerm)
g$secondaryTerm <- gsub('\\(SSU-rRNA, 5.8S rRNA, LSU-rRNA\\)', '', g$secondaryTerm)
g$secondaryTerm <- gsub('nervous system development', 
                        'nervous\nsystem\ndevelopment', g$secondaryTerm)
g$secondaryTerm <- gsub('electron transport chain', '\nelectron transport chain', g$secondaryTerm)
g$secondaryTerm <- gsub('cholesterol metabolic process', '\n\ncholesterol metabolic process', g$secondaryTerm)
g$secondaryTerm <- gsub('oligosaccharide biosynthetic process',
                        '\n\n\noligosaccharide biosynthetic process', g$secondaryTerm)
g$secondaryTerm <- gsub('positive regulation of transcription by RNA polymerase II', 
                        'positive regulation of transcription\nby RNA polymerase II', g$secondaryTerm)
g$secondaryTerm <- gsub('cytokine-mediated signaling pathway', 
                        'cytokine-mediated\nsignaling pathway', g$secondaryTerm)
g$secondaryTerm <- gsub('G protein-coupled receptor signaling pathway', 
                        'G protein-coupled receptor\nsignaling pathway', g$secondaryTerm)
g$secondaryTerm <- gsub('inflammatory response', 'inflammatory\nresponse', g$secondaryTerm)
# pie-donut chart
png(res = 300, units = 'in', height = 12, width = 12,
    filename = '/work/users/n/n/nnishi/eqtl/freeze/final/figures/pie-donut_GO_clustering_GTEx_all_eGenes.png')
PieDonut(g, aes(GTEx, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.03, 
         titlesize = 10)
dev.off()
# UNC
i <- coloc.ibd %>% group_by(primaryTerm, secondaryTerm) %>% summarize(n = n_distinct(gene)) %>% arrange(primaryTerm, secondaryTerm)
colnames(i)[1] <- 'UNC'
i$UNC <- gsub('positive regulation of cell population proliferation', 
              'positive regulation of cell\npopulation proliferation', i$UNC)
i$UNC <- gsub('regulation of transcription by RNA polymerase II', 
              'regulation of transcription\nby RNA polymerase II', i$UNC)
i$secondaryTerm <- gsub('cell adhesion', 'cell adhesion\n', i$secondaryTerm)
i$secondaryTerm <- gsub('\\(SSU-rRNA, 5.8S rRNA, LSU-rRNA\\)', '', i$secondaryTerm)
i$secondaryTerm <- gsub('nervous system development', 
                        'nervous\nsystem\ndevelopment', i$secondaryTerm)
i$secondaryTerm <- gsub('immune response', 'immune\nresponse', i$secondaryTerm)
i$secondaryTerm <- gsub('positive regulation of gene expression',
                        'positive regulation of\ngene expression', i$secondaryTerm)
i$secondaryTerm <- gsub('positive regulation of transcription by RNA polymerase II', 
                        'positive regulation of transcription\nby RNA polymerase II', i$secondaryTerm)
i$secondaryTerm <- gsub('cytokine-mediated signaling pathway', 
                        'cytokine-mediated\nsignaling pathway', i$secondaryTerm)
i$secondaryTerm <- gsub('G protein-coupled receptor signaling pathway', 
                        'G protein-coupled receptor\nsignaling pathway', i$secondaryTerm)




# pie-donut chart
png(res = 300, units = 'in', height = 12, width = 12,
    filename = '/work/users/n/n/nnishi/eqtl/freeze/final/figures/pie-donut_GO_clustering_IBD_all_eGenes.png')
PieDonut(i, aes(UNC, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.03, 
         titlesize = 10)
dev.off()
################################################################################
# create pie-donut charts for unique eGenes only

# label eGenes as shared or unique
coloc.barcuva <- coloc.barcuva %>% mutate(shared = ifelse(gene %in% c(coloc.gtex$gene, coloc.ibd$gene), 'shared', 'unique'))
coloc.gtex <- coloc.gtex %>% mutate(shared = ifelse(gene %in% c(coloc.barcuva$gene, coloc.ibd$gene), 'shared', 'unique'))
coloc.ibd <- coloc.ibd %>% mutate(shared = ifelse(gene %in% c(coloc.gtex$gene, coloc.barcuva$gene), 'shared', 'unique'))
# BarcUVa eGenes
# subset & summarize
b <- coloc.barcuva %>% filter(shared == 'unique') %>% group_by(primaryTerm, secondaryTerm) %>% summarize(n = n_distinct(gene)) %>% arrange(primaryTerm, secondaryTerm)
colnames(b)[1] <- 'BarcUVa'
# reformat for plotting
b$BarcUVa <- gsub('positive regulation of cell population proliferation', 
                  'positive regulation of cell\npopulation proliferation', b$BarcUVa)
b$BarcUVa <- gsub('regulation of transcription by RNA polymerase II', 
                  'regulation of transcription\nby RNA polymerase II', b$BarcUVa)
b$secondaryTerm <- gsub('nervous system development', 
                        'nervous\nsystem\ndevelopment', b$secondaryTerm)
b$secondaryTerm <- gsub('negative regulation of cell population proliferation', 
                        'negative regulation of\ncell population proliferation', b$secondaryTerm)
b$secondaryTerm <- gsub('protein stabilization', 'protein\nstabilization', b$secondaryTerm)
b$secondaryTerm <- gsub('cholesterol metabolic process', '\n\ncholesterol metabolic process', b$secondaryTerm)
b$secondaryTerm <- gsub('electron transport chain', '\nelectron transport chain', b$secondaryTerm)
b$secondaryTerm <- gsub('positive regulation of gene expression', 
                        'positive regulation of\ngene expression', b$secondaryTerm)
b$secondaryTerm <- gsub('positive regulation of transcription by RNA polymerase II', 
                        'positive regulation of transcription\nby RNA polymerase II', b$secondaryTerm)
b$secondaryTerm <- gsub('proteasome-mediated ubiquitin-dependent protein catabolic process', 
                        'proteasome-mediated ubiquitin-\ndependent protein catabolic process', b$secondaryTerm)
b$secondaryTerm <- gsub('regulation of transcription by RNA polymerase II', 
                        'regulation of transcription\nby RNA polymerase II', b$secondaryTerm)
b$secondaryTerm <- gsub('complement activation, classical pathway', 
                        'complement activation,\nclassical pathway', b$secondaryTerm)
b$secondaryTerm <- gsub('cytokine-mediated signaling pathway', 
                        'cytokine-mediated\nsignaling pathway', b$secondaryTerm)
b$secondaryTerm <- gsub('G protein-coupled receptor signaling pathway', 
                        'G protein-coupled receptor\nsignaling pathway', b$secondaryTerm)
#pie-donut chart
png(res = 300, units = 'in', height = 12, width = 12,
    filename = '/work/users/n/n/nnishi/eqtl/freeze/final/figures/pie-donut_GO_clustering_BarcUVa_unique_eGenes.png')
PieDonut(b, aes(BarcUVa, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.02, 
         titlesize = 10)
dev.off()
# GTEx eGenes
g <- coloc.gtex %>% filter(shared == 'unique') %>% group_by(primaryTerm, secondaryTerm) %>% summarize(n = n_distinct(gene)) %>% arrange(primaryTerm, secondaryTerm)
colnames(g)[1] <- 'GTEx'
g$GTEx <- gsub('positive regulation of cell population proliferation', 
                  'positive regulation of cell\npopulation proliferation', g$GTEx)
g$GTEx <- gsub('regulation of transcription by RNA polymerase II', 
                  'regulation of transcription\nby RNA polymerase II', g$GTEx)
g$secondaryTerm <- gsub('cell adhesion', 'cell adhesion\n', g$secondaryTerm)
g$secondaryTerm <- gsub('\\(SSU-rRNA, 5.8S rRNA, LSU-rRNA\\)', '', g$secondaryTerm)
g$secondaryTerm <- gsub('nervous system development', 
                        'nervous\nsystem\ndevelopment', g$secondaryTerm)
g$secondaryTerm <- gsub('electron transport chain', '\nelectron transport chain', g$secondaryTerm)
g$secondaryTerm <- gsub('positive regulation of transcription by RNA polymerase II', 
                        'positive regulation of transcription\nby RNA polymerase II', g$secondaryTerm)
g$secondaryTerm <- gsub('proteasome-mediated ubiquitin-dependent protein catabolic process', 
                        'proteasome-mediated ubiquitin-\ndependent protein catabolic process', g$secondaryTerm)
g$secondaryTerm <- gsub('regulation of transcription by RNA polymerase II', 
                        'regulation of transcription\nby RNA polymerase II', g$secondaryTerm)
g$secondaryTerm <- gsub('complement activation, classical pathway', 
                        'complement activation,\nclassical pathway', g$secondaryTerm)
g$secondaryTerm <- gsub('cytokine-mediated signaling pathway', 
                        'cytokine-mediated\nsignaling pathway', g$secondaryTerm)
g$secondaryTerm <- gsub('G protein-coupled receptor signaling pathway', 
                        'G protein-coupled receptor\nsignaling pathway', g$secondaryTerm)
g$secondaryTerm <- gsub('negative regulation of canonical Wnt signaling pathway', 
                        'negative regulation of canonical\nWnt signaling pathway', g$secondaryTerm)
g$secondaryTerm <- gsub('response to xenobiotic stimulus', 
                        'response to\nxenobiotic stimulus', g$secondaryTerm)
#pie-donut chart
png(res = 300, units = 'in', height = 12, width = 12,
    filename = '/work/users/n/n/nnishi/eqtl/freeze/final/figures/pie-donut_GO_clustering_GTEx_unique_eGenes.png')
PieDonut(g, aes(GTEx, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.02, 
         titlesize = 10)
dev.off()
# UNC eGenes
i <- coloc.ibd %>% filter(shared == 'unique') %>% group_by(primaryTerm, secondaryTerm) %>% summarize(n = n_distinct(gene)) %>% arrange(primaryTerm, secondaryTerm)
colnames(i)[1] <- 'UNC'
i$UNC <- gsub('positive regulation of cell population proliferation', 
               'positive regulation of cell\npopulation proliferation', i$UNC)
i$UNC <- gsub('regulation of transcription by RNA polymerase II', 
               'regulation of transcription\nby RNA polymerase II', i$UNC)
i$secondaryTerm <- gsub('cell adhesion', 'cell adhesion\n', i$secondaryTerm)
i$secondaryTerm <- gsub('\\(SSU-rRNA, 5.8S rRNA, LSU-rRNA\\)', '', i$secondaryTerm)
i$secondaryTerm <- gsub('positive regulation of inflammatory response', 
                        'positive regulation of\ninflammatory response', i$secondaryTerm)
i$secondaryTerm <- gsub('negative regulation of canonical Wnt signaling pathway', 
                        'negative regulation of canonical\nWnt signaling pathway', i$secondaryTerm)
png(res = 300, units = 'in', height = 12, width = 12,
    filename = '/work/users/n/n/nnishi/eqtl/freeze/final/figures/pie-donut_GO_clustering_IBD_unique_eGenes.png')
PieDonut(i, aes(UNC, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.05,
         titlesize = 10)
dev.off()
################################################################################
# run chi-square test to compare between studies

# pull out primary terms
pt <- sort(unique(coloc.barcuva$primaryTerm))
group <- c('barcuva', 'gtex', 'ibd')
# make contingency tables how many genes fall into a cluster vs not
contingency.tables <- list()
for (p in pt) {
  x <- list()
  for (g in group) {
    if (g == group[1]) {
      df <- coloc.barcuva
    } else if (g == group[2]) {
      df <- coloc.gtex
    } else df <- coloc.ibd
    df <- subset(df, shared == 'unique')
    y <- length(unique(subset(df, primaryTerm == p)$gene))
    z <- c(y, length(unique(df$gene))-y)
    x <- c(x, list(z))
  }
  names(x) <- group
  contingency.tables[[p]] <- x
  remove(x,y,z)
}
# secondary terms
st <- sort(unique(red.out$secondaryTerm))
contingency.tables <- list()
for (p in st) {
  x <- list()
  for (g in group) {
    if (g == group[1]) {
      df <- coloc.barcuva
    } else if (g == group[2]) {
      df <- coloc.gtex
    } else df <- coloc.ibd
    df <- subset(df, shared == 'unique')
    y <- length(unique(subset(df, secondaryTerm == p)$gene))
    z <- c(y, length(unique(df$gene))-y)
    x <- c(x, list(z))
  }
  names(x) <- group
  contingency.tables[[p]] <- x
  remove(x,y,z)
}
# run chi-squared test
pairs <- c('ibd vs barcuva', 'ibd vs gtex', 'gtex vs barcuva')
res <- list()
for (p in pt) {
  ft <- list()
  x <- as.data.frame(contingency.tables[[p]])
  y <- chisq.test(t(x))
  print(c(p, y$p.value))
  # add post-hoc pairwise fisher's exact comparisons
  for (q in pairs) {
    if (q == pairs[1]) {
      df <- as.data.frame(contingency.tables[[p]][c(3,1)])
    } else if (q == pairs[2]) {
      df <- as.data.frame(contingency.tables[[p]][c(3,2)])
    } else df <- as.data.frame(contingency.tables[[p]][c(2,1)])
    f <- fisher.test(t(df))
    ft <- c(ft, list(f))
    #print(c(q, f$p.value))
  }
  res[[p]] <- c(y, ft)
  names(res[[p]])[10:12] <- pairs
  remove(x,y)
}




# z-test for proportions
contingency.tables <- list()
for (p in pt) {
  x <- list()
  for (g in group) {
    if (g == group[1]) {
      df <- coloc.barcuva
    } else if (g == group[2]) {
      df <- coloc.gtex
    } else df <- coloc.ibd
    df <- subset(df, shared == 'unique')
    y <- length(unique(subset(df, primaryTerm == p)$gene))
    z <- c(y, length(unique(df$gene)))
    x <- c(x, list(z))
  }
  names(x) <- group
  contingency.tables[[p]] <- x
  remove(x,y,z)
}

res <- list()
for (p in pt) {
  zt <- list()
  x <- as.data.frame(contingency.tables[[p]])
  y <- chisq.test(t(x))
  print(c(p, y$p.value))
  # add post-hoc pairwise fisher's exact comparisons
  for (q in pairs) {
    if (q == pairs[1]) {
      z <- prop.test(c(contingency.tables[[p]][3][[1]][1], contingency.tables[[p]][1][[1]][1]),
                     c(contingency.tables[[p]][3][[1]][2], contingency.tables[[p]][1][[1]][2]))
    } else if (q == pairs[2]) {
      z <- prop.test(c(contingency.tables[[p]][3][[1]][1], contingency.tables[[p]][2][[1]][1]),
                     c(contingency.tables[[p]][3][[1]][2], contingency.tables[[p]][2][[1]][2]))
    } else z <- prop.test(c(contingency.tables[[p]][2][[1]][1], contingency.tables[[p]][1][[1]][1]),
                          c(contingency.tables[[p]][2][[1]][2], contingency.tables[[p]][1][[1]][2]))
    zt <- c(zt, list(z))
    print(c(q, z$p.value))
  }
  res[[p]] <- c(y, zt)
  names(res[[p]])[10:12] <- pairs
  remove(x,y,z)
}
################################################################################
# # make a list of gene set for each cluster
# gene.sets <- list()
# for (x in group) {
#   y <- list()
#     for (p in pt) {
#       if (x == 'barcuva') {
#         df <- coloc.barcuva
#       } else if (x == 'gtex') {
#         df <- coloc.gtex
#       } else df <- coloc.ibd
#       z <- unique(subset(df, primaryTerm == p)$gene)
#       y <- c(y, list(z))
#     }
#   names(y) <- pt
#   gene.sets[[x]] <- y
#   remove(y,z)
# }
# 
# gene.sets <- list()
# for (x in group) {
#   y <- list()
#   for (p in pt) {
#     if (x == 'barcuva') {
#       df <- coloc.barcuva
#     } else if (x == 'gtex') {
#       df <- coloc.gtex
#     } else df <- coloc.ibd
#     z <- length(unique(subset(df, primaryTerm == p)$gene))
#     zz <- c(z, length(unique(df$gene))-z)
#     y <- c(y, list(zz))
#   }
#   names(y) <- pt
#   gene.sets[[x]] <- y
#   remove(y,z, zz)
# }
# # number of unique genes that fall into cluster vs not (not study-unique genes)
# contingency.tables <- list()
# for (x in group) {
#   y <- list()
#   for (p in pt) {
#     n <- setdiff(pt, p)
#     z <- length(setdiff(gene.sets[[x]][[p]], Reduce(union, gene.sets[[x]][n])))
#     y <- c(y, list(z))
#     if (x == 'barcuva') {
#       df <- coloc.barcuva
#     } else if (x == 'gtex') {
#       df <- coloc.gtex
#     } else df <- coloc.ibd
#     z <- unique(subset(df, primaryTerm == p)$gene)
#   }
#   names(y) <- pt
#   contingency.tables[[x]] <- y
#   remove(y,z)
# }
################################################################################
# run fisher's exact test comparing all vs unique genes for each study


group <- c('all','unique')
# make contingency tables how many genes fall into a cluster vs not
contingency.tables <- list()
for (p in pt) {
  x <- list()
  df <- coloc.ibd
  for (g in group) {
    if (g == group[1]) {
      df <- df
    } else if (g == group[2]) {
      df <- subset(df, shared == 'unique')
    }
    y <- length(unique(subset(df, primaryTerm == p)$gene))
    z <- c(y, length(unique(df$gene))-y)
    x <- c(x, list(z))
  }
  names(x) <- group
  contingency.tables[[p]] <- x
  remove(x,y,z)
}
# run fisher's exact test
res <- list()
for (p in pt) {
  ft <- list()
  x <- as.data.frame(contingency.tables[[p]])
  y <- fisher.test(t(x))
  print(c(p, y$p.value))
  # # add post-hoc pairwise fisher's exact comparisons
  # for (q in pairs) {
  #   if (q == pairs[1]) {
  #     df <- as.data.frame(contingency.tables[[p]][c(3,1)])
  #   } else if (q == pairs[2]) {
  #     df <- as.data.frame(contingency.tables[[p]][c(3,2)])
  #   } else df <- as.data.frame(contingency.tables[[p]][c(2,1)])
  #   f <- fisher.test(t(df))
  #   ft <- c(ft, list(f))
  #   print(c(q, f$p.value))
  # }
  # res[[p]] <- c(y, ft)
  # names(res[[p]])[10:12] <- pairs
  res[[p]] <- y
  remove(x,y)
}




