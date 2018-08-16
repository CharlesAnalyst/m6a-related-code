library("TxDb.Hsapiens.UCSC.hg38.knownGene")

setwd("/data5/galaxy/project/data/promoter")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_txdb <- genes(txdb)
promoters_txdb <- promoters(genes_txdb)
#By default, the promoters function will fetch the 2000 nucleotides before the transcription start site (TSS) and the 200 nucleotides after the TSS.
promoters_txdb
