# library(coin)
# http://yatani.jp/teaching/doku.php?id=hcistats:mannwhitney
# How to report
# You can report the results of Mann-Whitney's U test as follows:
# 
# The medians of Group A and Group B were 2.5 and 3.5, respectively. 
# We ran a Mann-Whitney's U test to evaluate the difference in the responses of our 5-Likert scale question. 
# We found a significant effect of Group (The mean ranks of Group A and Group B were 7.8 and 13.2, respectively; 
#                                         U = 23, Z = -2.11, p < 0.05, r = 0.47).
# 
library(effsize)
setwd("/data5/galaxy/project/p-value/fig2")
total_df = read.table("data_forplot", sep = "\t", header = TRUE)
head(total_df)
#
tissue_list = list()
col_list = colnames(total_df)
i = 0
for (x in col_list){
  i = i + 1
  tissue_list[[i]] = strsplit(x, "_")[[1]][1]
}
uniq_tissue_list = unique(tissue_list)
#
for (tissue in uniq_tissue_list){
  print(tissue)
  calculate_effectiveSize(tissue)
}

calculate_effectiveSize <- function(tissue){
  df = total_df[total_df$Organ == "Stomach",]
  # print(paste(tissue, "_m6a", sep = ""))
  # GroupA <- total_df[, paste(tissue, "_m6a", sep = "")]
  GroupA <- df[df$Type == "DMR", ][, 3]
  GroupB <- df[df$Type == "Common", ][, 3]
  # class(GroupA)
  # GroupB <- total_df[, paste(tissue, "_no", sep = "")]
  # print(head(GroupB))
  g = factor(c(rep("GroupA", length(GroupA)), rep("GroupB", length(GroupB))))
  v = c(GroupA, GroupB)
  print(wilcox_test(v ~ g)) # exact test
  r = rank(v)
  data = data.frame(g, r)
  d <- (split(data, data$g))
  print(length(GroupA) + length(GroupB))
  }

#
#       small size	medium size	large size
# abs(r)	0.1	           0.3	      0.5
#####################################################################################
setwd("/data5/galaxy/project/p-value/fig6")
df = read.table("CpG_OE_ratio.txt", header = TRUE, sep = "\t")
head(df)
GroupA <- df[df$V1 == "high_cv", ][, 2]
class(GroupA)
GroupB <- df[df$V1 == "low_cv", ][, 2]
# length(GroupA)
# length(GroupB)
g = factor(c(rep("GroupA", length(GroupA)), rep("GroupB", length(GroupB))))
v = c(GroupA, GroupB)
print(wilcox_test(v ~ g)) # exact test
#
r = rank(v)
data = data.frame(g, r)
d <- (split(data, data$g))
# mean(d$GroupA$r)
# mean(d$GroupB$r)
print(length(GroupA) + length(GroupB))

#####################################################################################
################################## cohen's d ########################################
# The magnitude is assessed using the thresholds provided in (Cohen 1992), 
# i.e. |d|<0.2 "negligible", |d|<0.5 "small", |d|<0.8 "medium", otherwise "large"
#####################################################################################
library(effsize)
setwd("/data5/galaxy/project/p-value/fig2")
######################################
total_df = read.table("data_forplot", sep = "\t", header = TRUE)
head(total_df)

for (tissue in levels(total_df$Organ)){
  print(tissue)
  calculate_effectiveSize(tissue)
}

calculate_effectiveSize <- function(tissue){
  df = total_df[total_df$Organ == tissue,]
  GroupA <- df[df$Type == "DMR", ][, 3]
  GroupB <- df[df$Type == "Common", ][, 3]
  result = cohen.d(GroupA, GroupB, pooled = TRUE, conf.level = 0.95)
  print(result)
  # print(result$estimate)
  # print(result$conf.int)
}


#################################
library(effsize)
setwd("/data5/galaxy/project/p-value/fig4/b")
df = read.table("high_low_cv.txt", sep = "\t", header = FALSE)
head(df)
GroupA <- df[df$V1 == "high_cv", ][, 2]
GroupB <- df[df$V1 == "low_cv", ][, 2]
result = cohen.d(GroupA, GroupB, pooled = TRUE, conf.level = 0.95)
print(result)
################################



##################################
library(effsize)
setwd("/data5/galaxy/project/p-value/fig4/d e")
############
# d 
df = read.table("high_low_cv.txt", sep = "\t", header = FALSE)
head(df)
GroupA <- df[df$V1 == "high_cv", ][, 2]
GroupB <- df[df$V1 == "low_cv", ][, 2]
result = cohen.d(GroupA, GroupB, pooled = TRUE, conf.level = 0.95)
print(result)
############
# e  lincRNA
############
df = read.table("high_low_cv_lincRNA.txt", sep = "\t", header = FALSE)
head(df)
GroupA <- df[df$V1 == "high_cv", ][, 2]
GroupB <- df[df$V1 == "low_cv", ][, 2]
result = cohen.d(GroupA, GroupB, pooled = TRUE, conf.level = 0.95)
print(result)
############
# e mRNA
############
df = read.table("high_low_cv_mRNA.txt", sep = "\t", header = FALSE)
head(df)
GroupA <- df[df$V1 == "high_cv", ][, 2]
GroupB <- df[df$V1 == "low_cv", ][, 2]
result = cohen.d(GroupA, GroupB, pooled = TRUE, conf.level = 0.95)
print(result)
################################



############################################################
library(effsize)
setwd("/data5/galaxy/project/p-value/fig6")
######################################
total_df = read.table("CpG_OE_ratio.txt", sep = "\t", header = TRUE)
head(total_df)
tissue_list = list()
col_list = colnames(total_df)
i = 0
for (x in col_list){
  i = i + 1
  tissue_list[[i]] = strsplit(x, "_")[[1]][1]
}
uniq_tissue_list = unique(tissue_list)
for (tissue in uniq_tissue_list){
  print(tissue)
  calculate_effectiveSize_2(tissue)
}


calculate_effectiveSize_2 <- function(tissue){
  GroupA <- na.omit(total_df[, paste(tissue, "_m6a", sep = "")])
  # print(paste(tissue, "_m6a", sep = ""))
  GroupB <- na.omit(total_df[, paste(tissue, "_no", sep = "")])
  # print(head(GroupB))
  result = cohen.d(GroupA, GroupB, pooled = TRUE, conf.level = 0.95)
  print(result)
  # print(result$estimate)
  # print(result$conf.int)
}


########################
#####plot overlap ######
# http://rpsychologist.com/short-r-script-to-plot-effect-sizes-cohens-d-and-shade-overlapping-area
plot_overlap <- function(){
  require("ggplot2")
  
  # Standardized Mean Difference (Cohen's d)
  ES <- 0.8
  # get mean2 depending on value of ES from d = (u1 - u2)/sd
  mean1 <- ES*1 + 1
  # create x sequence
  x <- seq(1 - 3*1, mean1 + 3*1, .01)
  # generate normal dist #1
  y1 <- dnorm(x, 1, 1)
  # put in data frame
  df1 <- data.frame("x" = x, "y" = y1)
  # generate normal dist #2
  y2 <- dnorm(x, mean1, 1)
  # put in data frame
  df2 <- data.frame("x" = x, "y" = y2)
  # get y values under overlap
  y.poly <- pmin(y1,y2)
  # put in data frame
  poly <- data.frame("x" = x, "y" = y.poly)
  
  # Cohen's U3, proportion of control > 50th perc. treatment
  u3 <- 1 - pnorm(1, mean1,1)
  u3 <- round(u3,3)
  
  # plot with ggplot2
  ggplot(df1, aes(x,y, color="treatment")) +
    # add line for treatment group
    geom_line(size=1) + 
    # add line for control group
    geom_line(data=df2, aes(color="control"),size=1) +
    # shade overlap
    geom_polygon(aes(color=NULL), data=poly, fill="red", alpha=I(4/10),
                 show_guide=F) +
    # add vlines for group means
    geom_vline(xintercept = 1, linetype="dotted") + 
    geom_vline(xintercept = mean1, linetype="dotted") + 
    # add plot title
    opts(title=paste("Visualizing Effect Sizes 
      (Cohen's d = ",ES,"; U3 = ",u3,")", sep="")) +
    # change colors and legend annotation
    scale_color_manual("Group", 
                       values= c("treatment" = "black","control" = "red")) +
    # remove axis labels
    ylab(NULL) + xlab(NULL)
}

