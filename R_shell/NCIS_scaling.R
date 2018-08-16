library(NCIS)
library(ShortRead)

options(echo=TRUE) # set to FALSE if you not  want see commands in output 
args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandsArgs(trailingOnly=FALSE))

ChIP_bam<- args[1]
input_control_bam<- args[2]

# NCIS usese the Aligned Reads object from the shortRead package, however, it is recommended
# to use GenomicAignments package to read in the bam files
ga_ChIP<- readGAlignments(ChIP_bam)
ga_input<-readGAlignments(input_control_bam)
# However, the resulting GenomicAlignment object is not recognized by NCIS.
# I have to use the legacy readAligned function from the ShortRead package.
# it takes around 15mins to finish

# ga_ChIP<- readAligned(ChIP_bam, type="BAM")
# ga_input<-readAligned(input_control_bam, type="BAM")

res<- NCIS(ga_ChIP, ga_input, data.type="AlignedRead")
res
res$est
res$r.seq.depth