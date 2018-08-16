library(car)
data = read.delim("modi_123_456_78.txt", header = TRUE)
# 正太性检验
qqPlot(lm(measurement ~ condition, data = data_long), simulate = T, labels=F)
# 数据转换 
data$squr = sqrt(data_long$measurement)
# plot density
plot(density(na.omit(data$tissue_78)))
# 方差齐性分析
bartlett.test(pingf ~ condition, data=data_long)

aov(measurement ~ condition, data_long)
# 数据转换
p1 <- powerTransform(data_123_sub)
bcPower(data_123_sub, p1$roundlam)

