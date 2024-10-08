# Supplementary Figure 5. UpSet plots of eQTL-GWAS colocalization results for new GWAS loci
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

# focus on newly reported IBD GWAS loci
new <- c('chr11:91270280:A:G','chr15:79913012:A:G','chr17:47259971:T:G','chr2:105777237:G:A','chr17:62449865:G:A','chr9:88081422:G:A','chr8:141164479:C:T','chr12:116408331:C:T','chr3:112334718:G:A','chr12:37910758:G:A','chr2:172559679:C:G','chr5:112513193:C:T','chr5:172151212:C:T','chr11:118273990:A:G','chr1:56454342:C:T','chr3:30608835:T:C','chr12:69348988:AT:A','chr14:105686797:A:G','chr11:17311651:A:G','chr1:33337547:C:T','chr17:7074375:C:T','chr22:26641724:G:C','chr20:40169106:T:G','chr16:75463237:C:T','chr16:78768798:A:G','chr16:27384341:C:CT','chr2:127092119:T:C','chr5:80229545:A:G','chr17:1453285:C:A','chr1:205516162:T:A','chr5:156956411:A:G','chr8:115614948:T:C','chr20:41301389:T:G','chr17:81183249:G:T','chr10:101436625:AATAGATAGATAGATAG:AATAGATAGATAGATAGATAG','chr18:2734945:C:G','chr18:50977537:G:A','chr9:271455:T:C','chr1:226819948:C:T','chr2:111175424:G:T','chr11:36410895:C:T','chr3:52013158:TGGA:T','chr20:14618543:C:A','chr1:28131939:C:T','chr5:178217126:CAG:C','chr10:48182406:G:T','chr19:18522945:T:C','chr12:56076060:G:A','chr19:36991870:G:A','chr7:46834111:A:T','chr13:52062509:C:T','chr6:117522909:A:G','chr16:2117121:C:T','chr14:65142694:C:G','chr15:80675204:T:C','chr9:133257521:T:TC','chr15:60777789:C:T','chr7:41935319:C:T','chr1:108112400:A:C','chr19:40712706:G:A','chr17:48268286:C:T','chr1:93349933:T:C','chr1:179134699:T:C','chr5:142866576:C:T','chr12:110259525:G:A','chr11:56679337:A:T','chr19:54203491:G:A','chr1:234599210:C:T','chr4:36061370:A:G','chr12:96137691:A:G','chr8:144042819:T:C','chr7:74711703:C:T','chr10:6607973:C:T','chr10:89146767:C:G','chr9:92674850:A:G','chr7:37417880:A:G','chr1:24968082:G:A','chr18:44806588:T:C','chr7:944282:T:A','chr17:40619224:C:T','chr1:150706557:G:A')

barcuva.new <- subset(coloc.barcuva, GWAS_index_snp %in% new)
gtex.new <- subset(coloc.gtex, GWAS_index_snp %in% new)
ibd.new <- subset(coloc.ibd, GWAS_index_snp %in% new)
length(unique(c(barcuva.new$GWAS_index_snp, gtex.new$GWAS_index_snp, ibd.new$GWAS_index_snp)))
length(unique(c(barcuva.new$gene, gtex.new$gene, ibd.new$gene)))
# UpSet plot for new GWAS loci
set.seed(12345)
set.new.loci = list(UNC = ibd.new$GWAS_index_snp,
                    GTEx = gtex.new$GWAS_index_snp,
                    BarcUVa = barcuva.new$GWAS_index_snp)
# create combination matrix
mat.new.loci <- make_comb_mat(set.new.loci, mode = 'distinct')
# UpSet plot
UpSet(mat.new.loci, set_order = c('UNC', 'GTEx', 'BarcUVa'))
png(filename = 'plots/Supplementary_Figure_5a_UpSet_GWAS_new_loci.png',
    res = 300, units = 'in', height = 4, width = 6)
ht <- draw(UpSet(mat.new.loci, set_order = c('UNC', 'GTEx', 'BarcUVa'),
                 pt_size = unit(5, 'mm'), lwd = 3,
                 comb_col =c('#7BAFD4','#13294B','#13294B','#13294B', '#13294B','#13294B','#13294B'),
                 top_annotation = upset_top_annotation(mat.new.loci, extend = 0.2, annotation_name_gp = gpar(fontsize=15),axis_param = list(gp=gpar(fontsize=14))),
                 right_annotation = upset_right_annotation(mat.new.loci, extend = 1, 
                                                           annotation_name_gp = gpar(fontsize=15), 
                                                           axis_param = list(gp=gpar(fontsize=14)),
                                                           gp = gpar(col = c('#13294B','#13294B','#13294B'),
                                                                     fill = c('#13294B','#13294B','#13294B'))),
                 row_names_gp = gpar(fontsize=20)))
od <- column_order(ht)
cs <- comb_size(mat.new.loci)
rs <- set_size(mat.new.loci)
ro <- row_order(ht)
decorate_annotation("intersection_size", {grid.text(cs[od], x = seq_along(cs), 
                                                    y = unit(cs[od], "native") + unit(2, "pt"), default.units = "native", 
                                                    just = "bottom", gp = gpar(fontsize = 20, col = c('#7BAFD4', '#13294B', '#13294B', '#13294B', '#13294B', '#13294B', '#13294B')))})
decorate_annotation('set_size', {grid.text(rs[ro], x = unit(rs[ro], 'native') + unit(16, 'pt'),
                                           y = rev(seq_len(length(rs))), default.units = 'native', just = 'bottom', 
                                           gp = gpar(fontsize = 20, col = c('#13294B','#13294B','#13294B')))})
dev.off()

# UpSet plot for new GWAS loci: eGenes
set.seed(12345)
set.new.egenes = list(UNC = ibd.new$gene,
                    GTEx = gtex.new$gene,
                    BarcUVa = barcuva.new$gene)
# create combination matrix
mat.new.egenes <- make_comb_mat(set.new.egenes, mode = 'distinct')
# UpSet plot
UpSet(mat.new.egenes, set_order = c('UNC', 'GTEx', 'BarcUVa'))
png(filename = 'plots/Supplementary_Figure_5b_UpSet_GWAS_new_eGenes.png',
    res = 300, units = 'in', height = 4, width = 6)
ht <- draw(UpSet(mat.new.egenes, set_order = c('UNC', 'GTEx', 'BarcUVa'),
                 pt_size = unit(5, 'mm'), lwd = 3,
                 comb_col =c('#7BAFD4','#13294B','#13294B','#13294B', '#13294B','#13294B','#13294B'),
                 top_annotation = upset_top_annotation(mat.new.egenes, extend = 0.2, annotation_name_gp = gpar(fontsize=15),axis_param = list(gp=gpar(fontsize=14))),
                 right_annotation = upset_right_annotation(mat.new.egenes, extend = 1, 
                                                           annotation_name_gp = gpar(fontsize=15), 
                                                           axis_param = list(gp=gpar(fontsize=14)),
                                                           gp = gpar(col = c('#13294B','#13294B','#13294B'),
                                                                     fill = c('#13294B','#13294B','#13294B'))),
                 row_names_gp = gpar(fontsize=20)))
od <- column_order(ht)
cs <- comb_size(mat.new.egenes)
rs <- set_size(mat.new.egenes)
ro <- row_order(ht)
decorate_annotation("intersection_size", {grid.text(cs[od], x = seq_along(cs), 
                                                    y = unit(cs[od], "native") + unit(2, "pt"), default.units = "native", 
                                                    just = "bottom", gp = gpar(fontsize = 20, col = c('#7BAFD4', '#13294B', '#13294B', '#13294B', '#13294B', '#13294B', '#13294B')))})
decorate_annotation('set_size', {grid.text(rs[ro], x = unit(rs[ro], 'native') + unit(16, 'pt'),
                                           y = rev(seq_len(length(rs))), default.units = 'native', just = 'bottom', 
                                           gp = gpar(fontsize = 20, col = c('#13294B','#13294B','#13294B')))})
dev.off()
