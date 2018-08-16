library("ChIPseeker")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
setwd("/data5/galaxy/project/mettl3_enrich/result_parse")
setwd("/data5/galaxy/project/mettl3_enrich/macs2_peak")
query_bed = "/data5/galaxy/project/mettl3_enrich/macs2_peak/mettl3.bed"

tfbs_bed = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_hg19_to_GRCh38/RAD51.bed"
#
tfbs_dir = "/data5/galaxy/project/what/TF_narrowPeak/total_TF_narrowPeak/03_2_GRCh38_clean"
tfbs_list = list.files(path = tfbs_dir, pattern = "*.bed", full.names = TRUE)
#
name_list <- list()
i = 1
for(x in tfbs_list){
  print(x)
  x_name = strsplit(strsplit(x, split = "/")[[1]][9], split = ".bed")[[1]][1]
  name_list[[i]] <- x_name
  i = i+1
}
names(tfbs_list) <- name_list
# test
# test_dir = "/data5/galaxy/project/mettl3_enrich/macs2_peak/test_bed"
# test_list = list.files(path = test_dir, pattern = "*.bed", full.names = TRUE)
# names(test_list) = c("HNF4G", "MZF1", "ZNF511")
enrichPeakOverlap(queryPeak = tfbs_bed , targetPeak = query_bed, TxDb = txdb, pAdjustMethod = "BH", nShuffle = 5000)

result = enrichPeakOverlap(queryPeak = query_bed, 
                  targetPeak = unlist(tfbs_list),
                  TxDb = txdb,
                  pAdjustMethod = "BH",
                  nShuffle = 5000,
                  chainFile     = NULL,
                  verbose       = FALSE)

### download GEO data
downloadGEObedFiles(genome = "hg18", destDir = "/data5/galaxy/project/mettl3_enrich/GEO_data")
####  TSS plot ###########
setwd("/data5/galaxy/project/mettl3_enrich/mettl3_enrich_high_CpG")
promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000)
peak <- readPeakFile("low_CpG_peaks.txt")
tagMatrix <- getTagMatrix(peak, windows=promoter)
plotAvgProf(tagMatrix, xlim=c(-5000, 5000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

write.table(promoter, "/data5/galaxy/project/mettl3_enrich/clean_bam/plot_profile/current_version/promoter_3_3.bed", sep = "\t", row.names = FALSE, quote = F, col.names = F)


###################################################################
setwd("/data5/galaxy/project/data/promoter/mouse/2k_100")
library("ChIPseeker")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_txdb <- genes(txdb)
promoters_txdb <- promoters(genes_txdb, upstream = 2000, downstream = 100)
promoters_txdb
write.table(promoters_txdb, "genes_promoters.bed", sep = "\t", col.names = F, row.names = F, quote = F)

###################################################################
library("ChIPseeker")
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_txdb <- genes(txdb)
promoters_txdb <- promoters(genes_txdb, upstream = 2000, downstream = 100)
write.table(promoters_txdb, "test.bed", sep = "\t", col.names = F, row.names = F, quote = F)

promoters_txdb

##########
# files <- list(peak1 = "/myFolder/peak1.bed", peak2 = "/myFolder/peak2.bed")
###########################################################################
peak <- readPeakFile('C:/Users/74469/Desktop/mettl3.bed')
df = read.table('C:/Users/74469/Desktop/clean_promoter2.bed')
# 将原始promoter文件转化为Granges object
promoter = makeGRangesFromDataFrame(df,keep.extra.columns=FALSE,ignore.strand=TRUE,seqinfo=NULL,seqnames.field=c("seqnames", "seqname","chromosome", "chrom","chr", "chromosome_name","seqid"),start.field="start",end.field=c("end", "stop"),strand.field="strand",starts.in.df.are.0based=FALSE)
print(promoter)
tagMatrix <- getTagMatrix(peak, windows=promoter)
p = plotAvgProf(tagMatrix,cex.axis=15,cex.lab=15, xlim=c(-2000, 2000),cex.lab=50,xlab="Genomic Region (5'->3')", ylab = "Frequency of METTL3 \n occupation (%)")
# 更改纵坐标刻度值
p1 = p + scale_y_continuous(breaks = c(0e+00,3e-04,6e-04,9e-04),labels = c(0,0.03,0.06,0.09))
# 去除网格、边框及更改坐标轴字体
p2 = p1 + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16),axis.text = element_text(size = 15,face = "italic"))
p2
###




