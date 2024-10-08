# Figure 1b. eQTL-eQTL colocalization hypotheses results
library(data.table)
library(tidyverse)
library(ggplot2)

# load coloc results
coloc.barcuva <- fread('data/UNC/Supplementary_Table_5_coloc_IBD_eQTL_BarcUVa_eQTL.txt')
coloc.gtex <- fread('data/UNC/Supplementary_Table_4_coloc_IBD_eQTL_GTEx_eQTL.txt')
# separate out ENSEMBL IDs
coloc.barcuva <- coloc.barcuva %>% separate(gene_id_IBD, c('gene', 'ver'), remove = FALSE)
coloc.gtex <- coloc.gtex %>% separate(gene_id_IBD, c('gene', 'ver'), remove = FALSE)
# number of genes tested
length(unique(c(coloc.barcuva$gene, coloc.gtex$gene)))
# limit results to genes with a significant eQTL in IBD tissue
coloc.barcuva <- coloc.barcuva[coloc.barcuva$qval_IBD < 0.05,]
coloc.gtex <- coloc.gtex[coloc.gtex$qval_IBD < 0.05,]
# number of genes tested, with a significant eQTL in IBD
length(unique(c(coloc.barcuva$gene, coloc.gtex$gene)))
# set PP thresholds
h1.b <- coloc.barcuva[coloc.barcuva$H1_hypothesis_IBD > 0.5,]
h2.b <- coloc.barcuva[coloc.barcuva$H2_hypothesis_BarcUVa > 0.5,]
h3.b <- coloc.barcuva[coloc.barcuva$H3_hypothesis_different > 0.5,]
h4.b <- coloc.barcuva[coloc.barcuva$H4_hypothesis_shared > 0.5,]
h1.g <- coloc.gtex[coloc.gtex$H1_hypothesis_IBD > 0.5,]
h2.g <- coloc.gtex[coloc.gtex$H2_hypothesis_GTEx > 0.5,]
h3.g <- coloc.gtex[coloc.gtex$H3_hypothesis_different > 0.5,]
h4.g <- coloc.gtex[coloc.gtex$H4_hypothesis_shared > 0.5,]
# count unique genes
length(unique(c(h4.b$gene, h4.g$gene)))
length(intersect(h3.b$gene, h3.g$gene))
length(unique(c(setdiff(h3.b$gene, h4.g$gene), setdiff(h3.g$gene, h4.b$gene))))
# merge
b <- rbind(h1.b, h2.b, h3.b, h4.b)
g <- rbind(h1.g, h2.g, h3.g, h4.g)
# combine and reformat results
cb <- b %>% select(gene, coloc_likely_hypothesis)
cg <- g %>% select(gene, coloc_likely_hypothesis)
cb$Comparison <- 'BarcUVa'
cg$Comparison <- 'GTEx'
df <- rbind(cb, cg)
df <- df %>% drop_na()
df$coloc_likely_hypothesis <- gsub('PP.', '', df$coloc_likely_hypothesis)
df$coloc_likely_hypothesis <- gsub('.abf', '', df$coloc_likely_hypothesis)
# plot
png(filename = 'plots/Figure_1b_eQTL-eQTL_coloc_hypotheses.png',
    res = 300, units = 'in', height = 6, width = 6)
ggplot(df, aes(x = coloc_likely_hypothesis, fill = Comparison)) + geom_bar(position = 'dodge') + 
  theme_minimal() + xlab('coloc Hypothesis (PP > 0.5)') + ylab('eQTL count') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = ..count..), stat = 'count', position = position_dodge(0.9), vjust = -0.3)
dev.off()
