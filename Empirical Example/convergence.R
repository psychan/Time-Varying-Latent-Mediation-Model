library(coda)

remain_ed<- read.csv(file='F:/course/statistic/tvm/S6/s6a/remain_ed_2.txt',header=TRUE,sep='\t')
remain_ed1<- read.csv(file='F:/course/statistic/tvm/S6/s6a/remain_ed_1.txt',header=TRUE,sep='\t')
df1 <- as.mcmc(remain_ed)
df2 <- as.mcmc(remain_ed1)
chains <- mcmc.list(df1,df2)
# 计算 Gelman-Rubin 统计量
gelman.diag(chains, multivariate = FALSE)


remain_ed2 <- read.csv(file='F:/course/statistic/tvm/S6/s6b/remain_ed.txt',header=TRUE,sep='\t')
geweke.diag(as.mcmc(remain_ed2))
