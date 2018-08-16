library(clusterProfiler)
setwd("/data5/galaxy/project/tf_analysis/motif_analysis/fimo_result/")
data <- read.table("zscore_filter-by-RPKM_result.txt", sep="\t", header = TRUE)
# symbol <- data$motif_id
entrez <- data$ensembl_id
# trans_result = bitr(symbol, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")
trans_result = bitr(data$ensembl_id, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
head(trans_result)
write.table(trans_result, "id_trans.txt", sep = "\t", quote = FALSE, row.names = FALSE)



data <- read.table("genes_promoters.bed", sep = "\t", header = F)
trans_result = bitr(data$V6, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = "org.Mm.eg.db")
head(trans_result)
write.table(trans_result, file = "bioMart_transform.txt", sep = "\t", row.names = F, quote = F)
