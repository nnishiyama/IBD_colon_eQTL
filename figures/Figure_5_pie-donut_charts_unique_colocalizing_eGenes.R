# Figure 5. pie-donut charts unique colocalizing eGenes
library(data.table)
library(tidyverse)
library(ggplot2)
library(webr)

# load GO clusters
red.out <- fread('data/UNC/Supplementary_Table_10_GO_term_clustering.txt')
# load coloc data
coloc.barcuva <- fread('data/BarcUVa/Supplementary_Table_9_coloc_H4_BarcUVa_eQTL_GWAS.txt')
coloc.gtex <- fread('data/GTEx/Supplementary_Table_8_coloc_H4_GTEx_eQTL_colon_transverse_GWAS.txt')
coloc.ibd <- fread('data/UNC/Supplementary_Table_7_coloc_h4_IBD_eQTL_GWAS.txt')
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
# label eGenes as shared or unique
coloc.barcuva <- coloc.barcuva %>% mutate(shared = ifelse(gene %in% c(coloc.gtex$gene, coloc.ibd$gene), 'shared', 'unique'))
coloc.gtex <- coloc.gtex %>% mutate(shared = ifelse(gene %in% c(coloc.barcuva$gene, coloc.ibd$gene), 'shared', 'unique'))
coloc.ibd <- coloc.ibd %>% mutate(shared = ifelse(gene %in% c(coloc.gtex$gene, coloc.barcuva$gene), 'shared', 'unique'))
################################################################################
# create pie-donut charts for unique eGenes only
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
    filename = 'plots/Figure_5a_pie-donut_GO_clustering_BarcUVa_unique_eGenes.png')
PieDonut(b, aes(BarcUVa, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.02, 
         titlesize = 10)
dev.off()
################################################################################
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
    filename = 'plots/Figure_5a_pie-donut_GO_clustering_GTEx_unique_eGenes.png')
PieDonut(g, aes(GTEx, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.02, 
         titlesize = 10)
dev.off()
################################################################################
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
    filename = 'plots/Figure_5a_pie-donut_GO_clustering_IBD_unique_eGenes.png')
PieDonut(i, aes(UNC, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.05,
         titlesize = 10)
dev.off()
