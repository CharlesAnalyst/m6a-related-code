library(ggplot2)
setwd("/data5/galaxy/project/trend_test")

#scale_fill_manual(values=c("skyblue", "seagreen3", "darkorange2"))

###########group_total_data.txt
total_data <- read.table("group_total_data.txt", sep = "\t", header = FALSE)
names(total_data) <- c("score", "type", "group")
total_data$group <- factor(total_data$group, levels = c("low", "medium", "high"), ordered = TRUE)
ggplot(total_data, aes(x=group, y=score)) + geom_violin(width=0.9, position = position_dodge(0.75), trim = TRUE, aes(fill=group)) + scale_fill_manual(values=c("skyblue", "skyblue", "skyblue")) + geom_boxplot(width=0.1, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("total.pdf", width = 8, height = 4)
###########group_mRNA_linc_data.txt
m_l_data <- read.table("group_mRNA_linc_data.txt", sep = "\t", header = FALSE)
names(m_l_data) <- c("score", "a", "degree", "group")
m_l_data$degree <- factor(m_l_data$degree, levels = c("low", "medium", "high"), ordered = TRUE)
# m_l_data$group <- factor(m_l_data$group, levels = c("mRNA", "linc"), ordered = TRUE)
m_l_data$group <- factor(m_l_data$group, levels = c("linc", "mRNA"), ordered = TRUE)
ggplot(m_l_data, aes(x=group, y=score, group=interaction(group, degree))) + geom_violin(width=0.9, position = position_dodge(0.75), trim = TRUE, aes(fill=degree)) + scale_fill_manual(values=c("skyblue", "skyblue", "skyblue")) + geom_boxplot(width=0.1, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("mRNA_linc.pdf", width = 16, height = 4)
##########group_5utr-cds-stop-3utr.txt
feature_data <- read.table("group_5utr-cds-stop-3utr.txt", sep = "\t", header = FALSE)
names(feature_data) <- c("score", "a", "degree", "group")
feature_data$degree <- factor(feature_data$degree, levels = c("low", "medium", "high"), ordered = TRUE)
feature_data$group <- factor(feature_data$group, levels = c("5utr", "cds", "stop", "3utr"), ordered = TRUE)
ggplot(feature_data, aes(x=group, y=score, group=interaction(group, degree))) + geom_violin(width=0.7, position = position_dodge(0.75), trim = TRUE, aes(fill=degree)) + scale_fill_manual(values=c("skyblue", "skyblue", "skyblue")) + geom_boxplot(width=0.08, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("5utr-cds-stop-3utr.pdf", width = 18, height = 4)
####################################################################################################
###########methy_cv.txt
setwd("/data5/galaxy/project/plot_figure")
total_data <- read.table("methy_cv.txt", sep = "\t", header = FALSE)
names(total_data) <- c("gene", "score", "degree")
total_data$degree <- factor(total_data$degree, levels = c("methy_low_cv", "methy_high_cv"), ordered = TRUE)
ggplot(total_data, aes(x=degree, y=score)) + geom_violin(width=0.9, position = position_dodge(0.75), trim = TRUE, aes(fill=degree)) + scale_fill_manual(values=c("skyblue", "skyblue")) + geom_boxplot(width=0.1, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("total.pdf", width = 8, height = 4)
##########methy_cv_mRNA_lincRNA
m_l_data <- read.table("methy_cv_mRNA_lincRNA", sep = "\t", header = FALSE)
names(m_l_data) <- c("gene", "score", "a", "group", "degree")
m_l_data$degree <- factor(m_l_data$degree, levels = c("low", "high"), ordered = TRUE)
m_l_data$group <- factor(m_l_data$group, levels = c("lincRNA", "mRNA"), ordered = TRUE)
ggplot(m_l_data, aes(x=group, y=score, group=interaction(group, degree))) + geom_violin(width=0.7, position = position_dodge(0.75), trim = TRUE, aes(fill=degree)) + scale_fill_manual(values=c("skyblue", "skyblue")) + geom_boxplot(width=0.08, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("mRNA_lincRNA.pdf", width = 16, height = 4)
###################################################################################################
setwd("/data5/galaxy/project/promoter_TF_enrich/m6a_vs_free/compare_CpG")
capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}
#
data <- read.table("m6a_vs_free.txt", sep = "\t", header = FALSE)
names(data) <- c("tissue", "degree", "score")
head(data)
#
data$tissue <- capFirst(data$tissue)
data$degree <- factor(data$degree, levels = c("free", "m6a"), ordered = TRUE)
data$tissue <- factor(data$tissue, ordered = TRUE)
p <- ggplot(data, aes(x=tissue, y=score, group=interaction(tissue, degree))) + geom_violin(width=0.7, position = position_dodge(0.75), trim = TRUE, aes(fill=degree)) + scale_fill_manual(values=c("skyblue", "darkorange2")) + geom_boxplot(width=0.08, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1 <- p + xlab("") + ylab("CpG O/E ratio") + theme(axis.title.y = element_text(family="Arial",face="plain", color = "black", size = 20))
p2 <- p1 + theme(axis.text.x=element_text(family="Arial",face="plain",color="black",size=rel(2.0), margin = margin(t=25))) + theme(axis.text.y=element_text(family="Times",face="plain",color="black",size=rel(2.0)))
p2 + guides(fill=F)
ggsave("CpG_OE_ratio.pdf", width = 20, height = 4)
####################################################################################################
setwd("/data5/galaxy/project/plot_figure/Rank_expCV")
data <- read.table("low_high_cv.txt", sep = "\t", header = FALSE)
names(data) <- c("gene", "score", "degree")
data$degree <- factor(data$degree, levels = c("low_cv", "high_cv"), ordered = TRUE)
ggplot(data, aes(x=degree, y=score)) + geom_violin(width=1, position = position_dodge(0.75), trim = TRUE, aes(fill=degree)) + scale_fill_manual(values=c("skyblue", "skyblue")) + geom_boxplot(width=0.03, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("low_high.pdf", width = 4, height = 4)

######
setwd("/data5/galaxy/project/promoter_TF_enrich/m6a_vs_free/3-5_UTR/compare_CpG")
data <- read.table("5UTR_vs_3UTR.txt", sep = "\t", header = FALSE)
names(data) <- c("tissue", "degree", "score")
data$degree <- factor(data$degree, levels = c("5UTR", "3UTR"), ordered = TRUE)
data$tissue <- factor(data$tissue, ordered = TRUE)
ggplot(data, aes(x=tissue, y=score, group=interaction(tissue, degree))) + geom_violin(width=0.7, position = position_dodge(0.75), trim = TRUE, aes(fill=degree)) + scale_fill_manual(values=c("skyblue", "darkorange2")) + geom_boxplot(width=0.08, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("CpG_OE_ratio.pdf", width = 20, height = 4)


#
setwd("/data5/xlj/new_tissue_m6A/7_10_2018_Figure1/DMR_Input0/Expression/Data_for_plot")
m_l_data <- read.table("data_forplot", sep = "\t", header = TRUE)
head(m_l_data)
names(m_l_data) <- c("group", "degree", "score")
# head(m_l_data)
m_l_data$degree <- factor(m_l_data$degree, levels = c("common", "DMR"), ordered = TRUE)
# m_l_data$group <- factor(m_l_data$group, levels = c("mRNA", "linc"), ordered = TRUE)
ggplot(m_l_data, aes(x=group, y=score, group=interaction(group, degree))) + geom_violin(width=0.9, position = position_dodge(0.75), trim = TRUE, aes(fill=degree)) + scale_fill_manual(values=c("skyblue", "darkorange2")) + geom_boxplot(width=0.1, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("/data5/galaxy/shell_dir/2018_3_17/jupyter_shell/plot_vi.pdf", width = 12, height = 4)

###
setwd("/data5/galaxy/project/lncRNA_analysis/m6a_expression")
data <- read.table("total-tissues_phastCons_vertebrate.txt", sep = "\t", header = TRUE)
head(data)
# names(data) <- c("group", "degree", "score")
# head(m_l_data)
# data = sub_data
for (x in tissue_list){
  print(x)
  sub_data = data[data$tissue == x, ]
  #
  m6a = sub_data[sub_data$type == "m6a", ]
  unm6a = sub_data[sub_data$type == "unm6a", ]
  m6a_score = m6a[m6a$mean > 0.3, ]
  unm6a_score = unm6a[unm6a$mean > 0.3, ]
  m6a_propor = nrow(m6a_score) / nrow(m6a)
  unm6a_propor = nrow(unm6a_score) / nrow(unm6a)
  print(paste(median(m6a_score$mean), median(unm6a_score$mean)))
  print(paste(median(m6a_propor), median(unm6a_propor)))
  # print(wilcox.test(m6a_score, unm6a_score))
}
data$type <- factor(data$type, levels = c("unm6a", "m6a"), ordered = TRUE)
# m_l_data$group <- factor(m_l_data$group, levels = c("mRNA", "linc"), ordered = TRUE)
ggplot(data, aes(x=tissue, y=mean0, group=interaction(tissue, type))) + geom_violin(width=0.9, position = position_dodge(0.75), trim = TRUE, aes(fill=type)) + scale_fill_manual(values=c("skyblue", "darkorange2")) + geom_boxplot(width=0.1, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("/data5/galaxy/project/lncRNA_analysis/m6a_expression/primate.pdf", width = 12, height = 4)


##########################################################
### peak phastCons
########################################
library(Hmisc)
setwd("/data5/galaxy/project/lncRNA_analysis/lncRNA_cons/formatted")
f_list  = list.files(".", pattern = "\\_discrete.tab")
df <- data.frame()
for (file in f_list){
  tissue = capitalize(strsplit(file, "_")[[1]][2])
  print(tissue)
  m6a_tab = paste("m6a_", tissue, ".tab", sep = "")
  df_input <- read_tab("input", tissue, file)
  df_m6a <- read_tab("m6a", tissue, m6a_tab)
  # filter
  # df_input <- df_input[df_input$mean >= 0.0, ]
  # df_m6a <- df_m6a[df_m6a$mean >= 0.0, ]
  #
  print(paste(nrow(df_input), nrow(df_m6a)))
  print(paste(median(df_m6a$mean), median(df_input$mean)))
  print(wilcox.test(df_m6a$mean, df_input$mean))
  df <- rbind(df_input, df)
  df <- rbind(df_m6a, df)
}
# head(df)
tail(df)
ggplot(df, aes(x=tissue, y=mean, group=interaction(tissue, type))) + geom_violin(width=0.9, position = position_dodge(0.75), trim = TRUE, aes(fill=type)) + scale_fill_manual(values=c("skyblue", "darkorange2")) + geom_boxplot(width=0.1, outlier.shape = NA, position = position_dodge(0.75), fill="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


read_tab <- function(file_type, tissue, x){
  data = read.table(x, sep = "\t", header = FALSE)
  # print(head(data))
  names(data) <- c("name", "size", "covered", "sum", "mean0", "mean")
  data["tissue"] = tissue
  data["type"] = file_type
  sub_data <- data[c("tissue", "type", "mean")]
  # print(head(sub_data))
  return(sub_data)
}


#######################################################
setwd("/data5/galaxy/project/tf_analysis/tf_and_PromoterCpG")
data <- read.table("result.bed", sep = "\t", header = FALSE)
head(data)
names(data) <- c("TFactor", "Type", "OERatio")
head(data)
data$Type <- factor(data$Type, levels = c("negative", "positive"), ordered = TRUE)
# m_l_data$group <- factor(m_l_data$group, levels = c("mRNA", "linc"), ordered = TRUE)
p <- ggplot(data, aes(x=TFactor, y=OERatio, group=interaction(TFactor, Type))) +
geom_violin(width=0.9, position = position_dodge(0.75), trim = TRUE, aes(fill=Type)) +
  scale_fill_manual(values=c("skyblue", "darkorange2")) +
  geom_boxplot(width=0.1, outlier.shape = NA, position = position_dodge(0.75), fill="white") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(axis.text.x=element_text(family="Arial",
                                 face="plain",
                                 color="black",
                                 size=rel(2.0), 
                                 margin = margin(t=5))) +
  theme(axis.text.y=element_text(family="Arial",
                                 face="plain",
                                 color="black",
                                 size=rel(2.0))) +
  theme(axis.title.y = element_text(family = "Arial", 
                                    face = "plain", 
                                    color = "black", 
                                    size=rel(1.8))) +
  theme(axis.title.x = element_blank())
p
ggsave("plot_vi.pdf", device = cairo_pdf, width = 8, height = 4)




