library(DESeq2)
setwd("C:/Users/Nick/Desktop/GBM/raw_data")
count_matrix <- read.table("example_exp_mat.txt", sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
annotations <- read.table("example_annotations.txt", sep = "\t", stringsAsFactors = F, header = T, row.names = 1, comment.char = "")
#check annotations in proper order
all(colnames(count_matrix) == annotations$Short.ID)

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = annotations, design = ~ RNA.Batch.. + label)
seq <- DESeq(dds)
resultsNames(seq)
res <- results(seq, contrast = c("label", "vitro.fail", "vitro.form"))
res_lfc <- lfcShrink(seq, res = res, contrast = c("label", "vitro.fail", "vitro.form"), type = "ashr")
