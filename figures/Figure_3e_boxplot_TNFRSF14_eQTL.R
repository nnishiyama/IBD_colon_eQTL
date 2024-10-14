# Figure 3e. TNFRSF14 eQTL boxplot
library(data.table)
library(tidyverse)
library(ggplot2)

# load expression residuals
df <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/extract/IBD_extract_TNFRSF14_residuals.content.txt.gz')
# round risk allelic dosage
df <- df %>% mutate(ra_dosage = round(`chr1:2557169:T:C`))
df <- df %>% select(sample, ENSG00000157873.18, ra_dosage)
colnames(df) <- c('sample', 'TNFRSF14', 'ra_dosage')
# flip risk allele
df <- df %>% mutate(ra_dosage = ifelse(ra_dosage == 0, 2, 
                                       ifelse(ra_dosage == 2, 0, ra_dosage)))
df$ra_dosage <- factor(df$ra_dosage)
# load phenotype data
coldata <- fread('data/IBD_coldata_UIDs.txt')
# merge
df <- merge(df, coldata, by = 'sample', all.x = TRUE)
# plot
png(filename = 'plots/Figure_3e_boxplot_TNFRSF14_eQTL.png',
    res = 300, units = 'in', height = 7, width = 7)
ggplot(df, aes(ra_dosage, TNFRSF14)) + geom_boxplot(outlier.shape=NA, aes(group=ra_dosage), lwd=1) + 
  geom_jitter(height=0, width=0.3, size=4, aes(color=disease)) + 
  labs(x='Risk allelic dosage', y='TNFRSF14 rank-normalized residual expression', title='TNFRSF14-chr1:2557169:T:C eQTL',
       color='Disease') +
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20, face='bold'), 
        legend.text = element_text(size=20), legend.title = element_text(size=20), 
        legend.position = c(0.1,0.09), plot.title = element_text(hjust = 0.5, size = 20, face='bold'))
dev.off()
