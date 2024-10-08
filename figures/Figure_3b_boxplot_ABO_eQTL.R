# Figure 4b. ABO eQTL
library(data.table)
library(tidyverse)
library(ggplot2)

# load expression residuals
df <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/extract/IBD_extract_ABO_residuals.content.txt.gz')
# round risk allelic dosage
df <- df %>% mutate(ra_dosage = round(`chr9:133257521:T:TC`))
df <- df %>% select(sample, ENSG00000175164.16, ra_dosage)
colnames(df) <- c('sample', 'ABO', 'ra_dosage')
df$ra_dosage <- factor(df$ra_dosage)
# load phenotype data
coldata <- fread('/proj/fureylab/nnishi/eqtl/RNA-seq_All_genotyped_colon_2024-01-12.txt')
# merge
df <- merge(df, coldata, by = 'sample', all.x = TRUE)
# plot
png(filename = '/work/users/n/n/nnishi/eqtl/freeze/final/figures/Figure4b_boxplot_ABO_eQTL.png',
    res = 300, units = 'in', height = 6, width = 6)
ggplot(df, aes(ra_dosage, ABO)) + geom_boxplot(outlier.shape=NA, aes(group=ra_dosage), lwd=1) + 
  geom_jitter(height=0, width=0.3, size=4, aes(color=disease)) + 
  labs(x='Risk allelic dosage', y='ABO rank-normalized residual expression', title='ABO-chr9:133257521:T:TC eQTL',
       color='Disease') +
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20, face='bold'), 
        legend.text = element_text(size=20), legend.title = element_text(size=20), 
        legend.position = c(0.92,0.1), plot.title = element_text(hjust = 0.5, size = 20, face='bold'))
dev.off()
