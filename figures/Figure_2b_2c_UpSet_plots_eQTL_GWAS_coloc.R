# Figure 2b, 2c. UpSet plots of eQTL-GWAS colocalization results
library(data.table)
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)

coloc.ibd <- fread('data/UNC/Supplementary_Table_7_coloc_h4_IBD_eQTL_GWAS.txt')
coloc.barcuva <- fread('data/GTEx/Supplementary_Table_8_coloc_H4_GTEx_eQTL_colon_transverse_GWAS.txt')
coloc.gtex <- fread('data/BarcUVa/Supplementary_Table_9_coloc_H4_BarcUVa_eQTL_GWAS.txt')
# reformat
coloc.ibd <- coloc.ibd %>% separate(gene, c('gene', 'ver'))
coloc.barcuva <- coloc.barcuva %>% separate(gene, c('gene', 'ver'))
coloc.gtex <- coloc.gtex %>% separate(gene, c('gene', 'ver'))
################################################################################
# UpSet plot GWAS loci
set.seed(12345)
# create a set of GWAS loci
set.loci <- list(UNC = coloc.ibd$GWAS_index_snp, 
                 GTEx = coloc.gtex$GWAS_index_snp, 
                 BarcUVa = coloc.barcuva$GWAS_index_snp)
# create combination matrix
mat.loci <- make_comb_mat(set.loci, mode = 'distinct')
# UpSet plot
UpSet(mat.loci, set_order = c('UNC', 'GTEx', 'BarcUVa'))
png(filename = 'plots/Figure_2b_UpSet_GWAS_loci.png',
    res = 300, units = 'in', height = 4, width = 6)
ht <- draw(UpSet(mat.loci, set_order = c('UNC', 'GTEx', 'BarcUVa'),
                 pt_size = unit(5, 'mm'), lwd = 3,
                 comb_col =c('#13294B','#13294B','#13294B','#13294B', '#7BAFD4','#13294B','#13294B'),
                 top_annotation = upset_top_annotation(mat.loci, extend = 0.2, annotation_name_gp = gpar(fontsize=15),axis_param = list(gp=gpar(fontsize=14))),
                 right_annotation = upset_right_annotation(mat.loci, extend = 1, 
                                                           annotation_name_gp = gpar(fontsize=15), 
                                                           axis_param = list(gp=gpar(fontsize=14)),
                                                           gp = gpar(col = c('#7BAFD4','#13294B','#13294B'),
                                                                     fill = c('#7BAFD4','#13294B','#13294B'))),
                 row_names_gp = gpar(fontsize=20)))
od <- column_order(ht)
cs <- comb_size(mat.loci)
rs <- set_size(mat.loci)
ro <- row_order(ht)
decorate_annotation("intersection_size", {grid.text(cs[od], x = seq_along(cs), 
                                                    y = unit(cs[od], "native") + unit(2, "pt"), default.units = "native", 
                                                    just = "bottom", gp = gpar(fontsize = 20, col = c('#13294B', '#13294B', '#13294B', '#13294B', '#7BAFD4', '#13294B', '#13294B')))})
decorate_annotation('set_size', {grid.text(rs[ro], x = unit(rs[ro], 'native') + unit(16, 'pt'),
                                           y = rev(seq_len(length(rs))), default.units = 'native', just = 'bottom', 
                                           gp = gpar(fontsize = 20, col = c('#7BAFD4','#13294B','#13294B')))})
dev.off()
################################################################################
# UpSet plot for colocalizing eGenes
set.seed(12345)
# create a list of colocalizing eGenes
set.egenes <- list(UNC = coloc.ibd$gene, 
                   GTEx = coloc.gtex$gene,
                   BarcUVa = coloc.barcuva$gene)
# create combination matrix
mat.egenes <- make_comb_mat(set.egenes, mode = 'distinct')
# UpSet plot
UpSet(mat.egenes, set_order = c('UNC', 'GTEx', 'BarcUVa'))
png(filename = 'plots/Figure_2c_UpSet_GWAS_eGenes.png',
    res = 300, units = 'in', height = 4, width = 6)
ht <- draw(UpSet(mat.egenes, set_order = c('UNC', 'GTEx', 'BarcUVa'),
                 pt_size = unit(5, 'mm'), lwd = 3,
                 comb_col =c('#13294B','#13294B','#13294B','#13294B', '#7BAFD4','#13294B','#13294B'),
                 top_annotation = upset_top_annotation(mat.egenes, extend = 0.2, annotation_name_gp = gpar(fontsize=15),axis_param = list(gp=gpar(fontsize=14))),
                 right_annotation = upset_right_annotation(mat.egenes, extend = 1, 
                                    annotation_name_gp = gpar(fontsize=15), 
                                    axis_param = list(gp=gpar(fontsize=14)),
                                    gp = gpar(col = c('#7BAFD4','#13294B','#13294B'),
                                              fill = c('#7BAFD4','#13294B','#13294B'))),
                 row_names_gp = gpar(fontsize=20)))
od <- column_order(ht)
cs <- comb_size(mat.egenes)
rs <- set_size(mat.egenes)
ro <- row_order(ht)
decorate_annotation("intersection_size", {grid.text(cs[od], x = seq_along(cs), 
                                            y = unit(cs[od], "native") + unit(2, "pt"), default.units = "native", 
                                            just = "bottom", gp = gpar(fontsize = 20, col = c('#13294B', '#13294B', '#13294B', '#13294B', '#7BAFD4', '#13294B', '#13294B')))})
decorate_annotation('set_size', {grid.text(rs[ro], x = unit(rs[ro], 'native') + unit(16, 'pt'),
                                           y = rev(seq_len(length(rs))), default.units = 'native', just = 'bottom', 
                                           gp = gpar(fontsize = 20, col = c('#7BAFD4','#13294B','#13294B')))})
dev.off()
