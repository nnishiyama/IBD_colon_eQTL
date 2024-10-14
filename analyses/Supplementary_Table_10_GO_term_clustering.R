# recursively cluster GO terms to create own hierarchy
library(data.table)
library(tidyverse)
library(rrvgo)

# load GO terms from DAVID
david.all <- fread('data/UNC/DAVID_GO_terms_colocalizing_eGenes.txt')
# reformat into long format
david.all <- david.all %>% mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ',GO:')) %>% 
  unnest(GOTERM_BP_DIRECT) %>% separate(GOTERM_BP_DIRECT, c('GO_ID', 'GO_term'), sep = '~') %>% 
  mutate(GO_ID = gsub('GO:', '', GO_ID)) %>% mutate(GO_ID = gsub('^', 'GO:', GO_ID))
# summary stats
range(table(david.all$ID))
median(table(david.all$ID))
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
            file = 'data/UNC/Supplementary_Table_10_GO_term_clustering.txt')
