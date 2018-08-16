library(clusterProfiler)
setwd("/data5/galaxy/project/lncRNA_analysis/m6a_expression/mRNA")

file_list <- list.files(".", "*.txt")
for (file in file_list){
  geneEnsemble_to_symbol(file)
  # sys.on.exit(0)
}
  


geneEnsemble_to_symbol <- function(in_file){
  data <- read.table(in_file, sep = "\t", header = TRUE)
  tissue = colnames(data)[2]
  print(tissue)
  # print(head(data))
  # names(data) <- c("chr", "start", "end", "length", "strand", "entrezid")
  eg = bitr(data$Gene_ID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
  # print(head(eg))
  df = merge(data, eg, by.x = "Gene_ID", by.y = "ENSEMBL")[c("SYMBOL", tissue)]
  colnames(df) <- c("Gene_ID", tissue)
  df = df[ordered(df$Gene_ID), ]
  write.table(unique(df), paste("symbol", in_file, sep = "_"), sep = "\t", row.names = FALSE, quote = FALSE)
  # print(head(df))
}
