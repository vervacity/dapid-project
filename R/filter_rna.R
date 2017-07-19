#!/usr/bin/env Rscript

# Description:
# removes media timepoints and
# filters out low expression genes
library(reshape)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
mat_gzfile <- args[1] # input RNA file
prefix <- args[2] # prefix for output files
data_mean_cutoff <- as.numeric(args[3]) # mean cutoff (empirically derived below)

# -------------------------------------------
# read in data and preprocess
# -------------------------------------------
data <- read.table(gzfile(mat_gzfile), header=TRUE, row.names=1)
sample_labels <- colnames(data)

# remove media timepoints and obvious zeros
data_nomedia <- data[,-c(15, 16, 17, 18, 23, 24)]
data_nonzero <- data_nomedia[rowSums(data_nomedia)>0, ] # remove zeros

# -------------------------------------------
# empirical threshold for low expr genes
# -------------------------------------------
data_dist_file <- paste(prefix, '.expression_distr.png', sep='')
if (!file.exists(data_dist_file)) {
    data_reshaped <- melt(data_nonzero)
    ggplot(data_reshaped, aes(x=value)) +
        geom_density() + geom_vline(xintercept=as.numeric(data_mean_cutoff))
    ggsave(data_dist_file)
}

data_row_max <- apply(data_nonzero, 1, max)
data_keep_indices <- which(data_row_max > data_mean_cutoff)
data_expressed <- data_nonzero[data_keep_indices,]

# change ENSG names
rownames(data_expressed) <- gsub("\\..*", "", rownames(data_expressed))

# write out to file
out_file = paste(prefix, ".expressed.mat.gz", sep='')
write.table(data_expressed, file=gzfile(out_file), quote=FALSE, sep='\t')


