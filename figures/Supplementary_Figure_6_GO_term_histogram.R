# Supplementary Figure 6. GO terms/eGene histogram
library(data.table)
library(tidyverse)
library(rrvgo)

# load GO terms from DAVID
david.all <- fread('data/DAVID_GO_terms_colocalizing_eGenes.txt')
# reformat into long format
david.all <- david.all %>% mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ',GO:')) %>% 
  unnest(GOTERM_BP_DIRECT) %>% separate(GOTERM_BP_DIRECT, c('GO_ID', 'GO_term'), sep = '~') %>% 
  mutate(GO_ID = gsub('GO:', '', GO_ID)) %>% mutate(GO_ID = gsub('^', 'GO:', GO_ID))
# summary stats
range(table(david.all$ID))
median(table(david.all$ID))
# plot number of GO terms per eGene
png(res = 300, units = 'in', height = 6, width = 6,
    filename = 'plots/Supplementary_Figure_6_histogram_GO_terms_per_eGene.png')
hist(table(david.all$ID), xlab = 'number of GO terms associated with eGene', 
     main = '')
dev.off()
