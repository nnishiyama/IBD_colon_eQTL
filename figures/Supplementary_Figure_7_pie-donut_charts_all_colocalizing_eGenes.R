# Supplementary Figure 7. pie-donut charts all colocalizing eGenes
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
################################################################################
# create pie/donut charts based off of GO clustering
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
    filename = 'plots/Supplementary_Figure_7a_pie-donut_GO_clustering_BarcUVa_all_eGenes.png')
PieDonut(b, aes(BarcUVa, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.03, 
         titlesize = 10)
dev.off()
################################################################################
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
    filename = 'plots/Supplementary_Figure_7b_pie-donut_GO_clustering_GTEx_all_eGenes.png')
PieDonut(g, aes(GTEx, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.03, 
         titlesize = 10)
dev.off()
################################################################################
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
    filename = 'plots/Supplementary_Figure_7c_pie-donut_GO_clustering_IBD_all_eGenes.png')
PieDonut(i, aes(UNC, secondaryTerm, count = n), pieLabelSize = 4, donutLabelSize = 2, 
         showRatioThreshold = F, showRatioDonut = F, labelpositionThreshold = 0.03, 
         titlesize = 10)
dev.off()
