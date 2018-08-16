library(edgeR)
setwd("/data5/galaxy/project/methyl_m6a/data/featureCount")
data = read.table("format_counts.txt", sep = "\t", header = TRUE, row.names = "Geneid")
# df = data[, c("SRR1596086", "SRR1596088", "SRR1596098", "SRR1596100")]
# df = data[, c("heart_3_L6.bam", "heart_1_L6.bam", "heart_2_L7.bam", "brain_1_L6.bam", "brain_2_L7.bam", "brain_3_LX.bam")]
df = data
head(df)
group = c("brain", "brain", "brain", "heart", "heart", "heart")
# group = c("lung", "heart", "kid", "mus", "brain", "mus", "heart", "liver", "pla", "kid", "sto", "pla", "liver", "brain", "brain", "liver", "heart", "pla", "lung", "sto", "kid")
cds <- DGEList(df, group = group)
y <- calcNormFactors(cds)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)
# qlf.2vs1 <- glmQLFTest(fit, coef=2)
# topTags(qlf.2vs1)
qlf <- glmQLFTest(fit, coef = 1)
topTags(qlf)
write.table(topTags(qlf, n = length(qlf)), sep = "\t", file="DEgenes.txt", quote = FALSE)
