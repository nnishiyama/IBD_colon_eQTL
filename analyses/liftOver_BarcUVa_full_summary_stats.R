# liftover BarcUVa eQTL to hg38
library(data.table)
library(tidyverse)
library(rtracklayer)

# load full summary stats with hg19 coordinates
barcuva <- fread('data/BarcUVa/barcuvaseq.eqtls.allpairs.sumstats.txt')
# load chain file
chain <- import.chain('/work/users/n/n/nnishi/liftover/hg19ToHg38.over.chain')
# reformat
barcuva$chr <- gsub('^', 'chr', barcuva$chr)
barcuva <- barcuva %>% unite('var_id_hg19', c(chr:ALT), sep = ':', remove = FALSE)
# create a GRanges object
b.gr <- GRanges(seqnames = Rle(barcuva$chr),
                ranges = IRanges(start = barcuva$pos, end = barcuva$pos),
                var_id_hg19 = barcuva$var_id_hg19, ref = barcuva$REF, alt = barcuva$ALT,
                rs_id = barcuva$variant_id, gene_id = barcuva$gene_id, 
                tss_distance = barcuva$tss_distance, ma_samples = barcuva$ma_samples, 
                ma_count = barcuva$ma_count, maf = barcuva$maf, 
                pval_nominal = barcuva$pval_nominal, slope = barcuva$slope,
                slope_se = barcuva$slope_se)
# liftover
b.hg38 <- unlist(liftOver(b.gr, chain))
# convert from Granges to df
b.hg38 <- as.data.frame(b.hg38)
# separate out coordinates that didn't liftover
snps <- setdiff(mcols(b.gr)$var_id_hg19, b.hg38$var_id_hg19)
b.hg19 <- as.data.frame(mcols(b.gr))
b.hg19 <- subset(b.hg19, var_id_hg19 %in% snps)
# check overlap with lead eSNPs
lead <- fread('/work/users/n/n/nnishi/BarcUVa/BarcUVA-seq_colon_eQTLs.csv')
lead <- lead[lead$qval < 0.05,]
lead <- lead %>% unite(var_id_hg19, c(chr, pos_hg19, REF, ALT), sep = ':', remove = FALSE)
length(intersect(lead$var_id_hg19, b.hg19$var_id_hg19))
# reformat
b.hg38 <- b.hg38 %>% unite('variant_id', c(seqnames, start, ref, alt), sep = ':')
b.hg38 <- b.hg38 %>% select(-c(end:var_id_hg19))
################################################################################
# save hg38 coordinates
write.table(b.hg38, sep = '\t', quote = FALSE, row.names = FALSE, 
            file = 'data/BarcUVa/barcuvaseq.eqtls.allpairs.sumstats.hg38.txt')
# save coordinates that didn't liftover
write.table(b.hg19, sep = '\t', quote = FALSE, row.names = FALSE,
            file = 'data/BarcUVa/barcuvaseq.eqtls.allpairs.sumstats.no_liftOver.txt')
