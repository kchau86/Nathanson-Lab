library(dplyr)

########################
#### Nathanson Data ####
########################

setwd("C:/Users/Nick/Desktop/GBM/raw_data")

# read in expression matrix
nat_expression_data <-read.table("Nathanson_RSEM_all+bbseal_xenos_raw_counts_20210208.txt", sep = "\t", stringsAsFactors = F, row.names = 1, header = T)
# remove zero variance genes and round the values to whole numbers (necessary for DESeq2)
nat_mat <- data.matrix(round(nat_expression_data[apply(nat_expression_data,1,var) > 0,]))
# read in annotation information
annotation <- read.table(file ="Sequencing Metadata batch 10 fixed.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
rna_anno <- annotation[!is.na(annotation$"RNA.Batch..") & annotation$FLAG.RNA != "terminate",]
rna_anno$Short.ID <- gsub("[-,_]", "\\.", rna_anno$Short.ID)
# useful for accessing rows by sample name
rownames(rna_anno) <- rna_anno$Short.ID
# remove some bad samples
rna_anno <- rna_anno %>% 
	arrange(Line..,Sample.Type, Short.ID) %>%
	filter(!FLAG.RNA %in% c("terminate", "relabel", "remove"))
# reorder annotations to match expression matrix in case they don't match
pair_ind <- match(colnames(nat_mat), rna_anno$Short.ID)
nat_anno <- rna_anno[pair_ind,]
# remove some more bad samples/non-GBM samples
filt_nat_mat <- nat_mat[,(!nat_anno$Line.. %in% c(104, 154, 133, 199)) & !is.na(nat_anno$Line..)]
filt_nat_anno <- nat_anno[(!nat_anno$Line.. %in% c(104, 154, 133, 199)) & !is.na(nat_anno$Line..),]
# filter genes with very low expression, these thresholds are arbitrary
count_matrix <- filt_nat_mat[rowSums(filt_nat_mat > 250) >= 15 | rowMeans(filt_nat_mat) > 100,]
