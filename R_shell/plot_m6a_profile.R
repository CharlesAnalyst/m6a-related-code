library(Guitar)
library(GenomicFeatures)
setwd("/data5/galaxy/project/DNMT1_KO/macs2_peak")
m6a_wt <- narrowPeaktoGRanges("CGR8_peaks.narrowPeak")
m6a_ko <- narrowPeaktoGRanges("DTK_peaks.narrowPeak")
feature_hg38 <- list(m6a_wt, m6a_ko)
names(feature_hg38) <- c("WT", "DNMT1-KO")
# ~ 0.5 hours
txdb <- makeTxDbFromUCSC(genome = "hg38")  # Download
gc_txdb <- makeGuitarCoordsFromTxDb(txdb, noBins = 100)
GuitarPlot(gfeatures = feature_hg38, GuitarCoordsFromTxDb = gc_txdb, saveToPDFprefix = "WT_vs_KO")

