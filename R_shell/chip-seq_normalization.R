library(DiffBind)
# For example, if you have a DBA object called myDBA:
  
myDBA <- dba.count(myDBA)
counts <- dba.peakset(myDBA, bRetrieve=TRUE)
# This will return a GRanges object with the count scores. The default count score is TMM normalised using the edgeR package. If you want raw read counts, you can change the score first:
    
myDBA <- dba.count(myDBA, peaks=NULL, score=DBA_SCORE_READS)
counts <- dba.peakset(myDBA, bRetrieve=TRUE)
# If you would prefer a data frame to a GRanges object:
    
counts <- dba.peakset(myDBA, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
# And if you want to write out the counts to a text file:
    
counts <- dba.peakset(myDBA, bRetrieve=TRUE, writeFile="readscores.txt")