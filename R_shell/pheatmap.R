library(pheatmap)
setwd("/data5/galaxy/project/tf_analysis/motif_analysis/fimo_result")
data <- read.table("zscore_filter-by-RPKM_result.txt", header = T, sep = "\t", row.names = 1)
data = na.omit(data)
pheatmap(data, show_rownames = FALSE, fontsize = 15, cluster_cols = FALSE)

makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}
#########
# cutree_rows = 8,
# color = c("#F5F5F5", "#FA8072"), breaks = c(0, 1.5, 5)
cols <- makeColorRampPalette(c("white", "grey", "#ffdab9", "#F4A460"), 1.5/max(data), 100)
res <- pheatmap(data, cluster_cols = FALSE,  number_format = "%.4f", display_numbers = FALSE, show_rownames = FALSE, color = cols, cluster_rows = TRUE)
data.cluster <- cbind(data, cluster = cutree(res$tree_row, k = 7))
head(data.cluster)
write.table(data.cluster, file = "zscore_full_info.txt", sep = "\t", quote = FALSE)
###
res2 <- data[c(res$tree_row[["order"]]), ]
write.table(res2, "ordered_output.txt", sep = "\t", quote = FALSE)

# 
setwd("/data5/galaxy/project/GWAS_analysis/cluster_by_map/Fisher_exact_test")
data = read.table("Fisher_m6a_vs_input_matrix.txt", sep = "\t", header = T, row.names = 1)
pheatmap(data)

