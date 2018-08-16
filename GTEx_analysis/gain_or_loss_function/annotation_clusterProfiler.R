#annotation
library(clusterProfiler)
library(org.Hs.eg.db)
library(plyr)
library(pathview)
setwd("/data5/galaxy/shell_dir/2018_3_17/GTEx_analysis/gain_or_loss_function")
file_data <- read.table("kidney.txt", header = FALSE, sep = "\t")
gene_name <- file_data$V1
pos_gene_id <- bitr(gene_name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID

# GO enrichment
go_enri <- function(i_genes){
  for(an_type in c("BP", "MF", "CC")){
    ego <- enrichGO(gene = i_genes,
                    OrgDb = org.Hs.eg.db,
                    ont = an_type,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
    file_prefix <- paste("GO_ENRICH_", an_type, sep = "")
    pdf(file = paste(file_prefix, ".pdf", sep = ""), width = 8, height = 8)
    print(barplot(ego, showCategory=12))
    dev.off()
  }
}


go_enri(pos_gene_id)


# KEGG annotation
# KEGG Enrichment
# transform SYMBOL to ENTREZID

#
kegg_rich <- function(i_genes){
  kk <- enrichKEGG(gene = i_genes,
                   organism = "hsa",
                   pvalueCutoff = 0.05)
  print(head(kk))
  file_prefix <- paste("KEGG_ENRICH_", sep = "")
  pdf(file = paste(file_prefix, ".pdf", sep = ""), width = 8, height = 8)
  print(barplot(kk, showCategory = 12))
  dev.off()
  # plot each kegg id pathway
  kegg_id <- summary(kk)$ID
  for(i_id in kegg_id){
    pdf(file = paste("KEGG_pathway_", i_id, ".pdf", sep = ""))
    print(browseKEGG(kk, i_id))
    dev.off()
  }
}

kegg_rich(pos_gene_id)








