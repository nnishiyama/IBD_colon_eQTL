#!/usr/bin/R

# Select eQTL model using kneedle and plot number of significant eQTL based on number of RUV covariates

args <- commandArgs(TRUE)

library(data.table)
library(reticulate)
library(qvalue)
use_python('/nas/longleaf/apps/python/3.6.6/bin/python')
kneed <- import('kneed')

# set QTL results directory
path <- args[1]
# create a list of files
file <- list.files(path)
# loop through files an quantify significant QTL hits
total <- c()
tested <- c()
sig.p <- c()
sig.q <- c()
for (i in 1:length(file)){
  tmp <- fread(paste(path, file[i], sep = '/'))
  tmp$q <- qvalue(tmp[[ncol(tmp)]])$qvalues
  total <- c(total, nrow(tmp))
  tested <- c(tested, length(na.omit(tmp$V20)))
  sig.p <- c(sig.p, nrow(tmp[tmp$V20<0.05]))
  sig.q <- c(sig.q, nrow(tmp[tmp$q<0.05]))
}
# create matrix of results
res <- as.data.frame(cbind(file,total,tested,sig.p,sig.q))
################################################################################
# write out results
write.table(res, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE,
            file = 'kneedle_eQTL.txt')
################################################################################
# reformat
res <- cbind(data.frame(do.call('rbind', strsplit(as.character(res$file), '.',
              fixed = TRUE))), res$total, res$tested, res$sig.p, res$sig.q)
# position of RUV
n <- ncol(res)-5
# reorder in ascending
res <- res[order(as.numeric(res[[n]])),]
# use kneedle to find the knee
ka <- kneed$KneeLocator(as.numeric(res[[n]]), res$`res$sig.q`, curve = 'concave', direction = 'increasing')
# print results
# number of RUV factors
ka$knee
# plot
png(filename = 'kneedle_eQTL.png', res = 300, units = 'in', height = 4, width = 4)
plot(res[[n]], res$`res$sig.q`, xlab='number of RUV factors',
     ylab='significant primary eQTL (FDR<0.05)',
     main='cis-eQTL results: RUV Factor Testing', col='black', pch=19)
abline(v=ka$knee, col='red')
legend(40,2200, paste('k =', ka$knee), col='red', lty = c(1), bty='n')
dev.off()
