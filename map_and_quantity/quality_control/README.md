#### A tool designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance which developed by [OpenGene](https://github.com/OpenGene).

##  fastp features

1. comprehensive quality profiling for both before and after filtering data (quality curves, base contents, KMER, Q20/Q30, GC Ratio, duplication, adapter contents...)
2. filter out bad reads (too low quality, too short, or too many N...)
3. cut low quality bases for per read in its 5' and 3' by evaluating the mean quality from a sliding window (like Trimmomatic but faster).
4. trim all reads in front and tail
5. cut adapters. Adapter sequences can be automatically detected,which means you don't have to input the adapter sequences to trim them.
6. correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality
7. trim polyG in 3' ends, which is commonly seen in NovaSeq/NextSeq data. Trim polyX in 3' ends to remove unwanted polyX tailing (i.e. polyA tailing for mRNA-Seq data)
8. preprocess unique molecular identifer (UMI) enabled data, shift UMI to sequence name.
report JSON format result for further interpreting.
9. visualize quality control and filtering results on a single HTML page (like FASTQC but faster and more informative).
10. split the output to multiple files (0001.R1.gz, 0002.R1.gz...) to support parallel processing. Two modes can be used, limiting the total split file number, or limitting the lines of each split file.
11. support long reads (data from PacBio / Nanopore devices).
12. support reading from STDIN and writing to STDOUT
13. support interleaved input
...
This tool is being intensively developed, and new features can be implemented soon if they are considered useful. If you have any additional requirement for fastp, please file an issue:https://github.com/OpenGene/fastp/issues/new
