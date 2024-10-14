# Figure 3c. ABO DE results: CL vs NIBD

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggsignif)

# load coldata
coldata <- fread('/nas/longleaf/home/nnishi/fureylab/projects/RNA-seq_processing/human/20230928_CL_nonIBD_sample_analysis/CL_vs_NIBD_coldata.txt')
# load counts
de <- fread('/nas/longleaf/home/nnishi/fureylab/projects/RNA-seq_processing/human/20230928_CL_nonIBD_sample_analysis/CL_vs_NIBD_colon_uninflamed_normalized_corrected_counts_wgeneid.txt')
# load results
res <- fread('/nas/longleaf/home/nnishi/fureylab/projects/RNA-seq_processing/human/20230928_CL_nonIBD_sample_analysis/DESeq2_results_all.txt')
# load genotypes
geno <- fread('/work/users/n/n/nnishi/eqtl/freeze/final/vcf/genotypes_for_DE_IBD_NIBD.txt', header = T)
# reformat genotypes
geno <- geno %>% separate(snp, c('chr', 'pos', 'ref', 'alt'), sep = ':', remove = F)
geno <- geno %>% mutate(across(colnames(geno)[-c(1:5)], gsub, pattern='\\|', replacement=''))
geno <- geno %>% mutate(across(colnames(geno)[-c(1:5)], str_replace_all, pattern='10', replacement='1'))
geno <- geno %>% mutate(across(colnames(geno)[-c(1:5)], str_replace_all, pattern='01', replacement='1'))
geno <- geno %>% mutate(across(colnames(geno)[-c(1:5)], str_replace_all, pattern='00', replacement='0'))
geno <- geno %>% mutate(across(colnames(geno)[-c(1:5)], str_replace_all, pattern='11', replacement='2'))
geno <- geno %>% select(-c(chr:alt))
geno <- as.data.frame(t(geno))
colnames(geno) <- geno[1,]
geno <- geno[-1,]
# subset genotype
g <- geno %>% select(`chr9:133257521:T:TC`)
colnames(g) <- c('Risk allelic dosage')
# subset to ABO
de <- de[de$gene_id == 'ABO',] %>% select(-ensembl)
# pull out stats for annotation
i <- which(res$gene_id == 'ABO')
stats <- c(round(res$log2FoldChange[i], 3), round(res$padj[i], 3))
annot <- paste('LFC = ', stats[1], ', padj = ', stats[2], sep ='')
# transpose
x <- as.data.frame(t(de))
colnames(x) <- x$V1[1]
x$sample <- row.names(x)
x <- x[-1,]
x$ABO <- as.numeric(x$ABO)
# merge
x <- merge(x, coldata, by = 'sample')
# add genotypes
x <- merge(x, g, by.x = 'sample', by.y = 'row.names')
# plot
png(filename = 'plots/Figure_3c_boxplot_ABO_DE.png',
    res = 300, units = 'in', height = 6, width = 5)
ggplot(x, aes(disease, ABO)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width=0.3, size = 2.5, aes(color = `Risk allelic dosage`)) + 
  theme_minimal() +
  theme(plot.title=element_text(face="bold", size=24, hjust = 0.5), 
        axis.title=element_text(face="bold", size=20),
        axis.text=element_text(size=14, color="black"), 
        panel.background = element_rect(fill="white",color="black"),
        legend.title=element_text(), legend.key=element_blank()) +
  labs(x="Disease phenotype", y="Normalized/Transformed Expression", title='ABO') +
  geom_signif(y_position=c(11.4), xmin=c(1), xmax=c(2),
              annotation=annot, tip_length=0) 
dev.off()
