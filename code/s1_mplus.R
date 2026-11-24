library(MplusAutomation)
library(textreg)

n_sample = c(20,50,100,200)
m_measure = c(10,20)
for (n in 1:4) {
  for (m in 1:2) {
    wd = paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n_sample[n]),
                'm',as.character(m_measure[m]))
    setwd(wd)
    for (i in 1:200) {
      data0 <- read.table(paste0('data_valley',as.character(i),'.txt'),header=TRUE)
      nid <- unique(data0$id)
      nobs <- as.numeric(as.vector(table(data0$id)))
      t <- matrix(rep(1:nobs,length(nid)))
      data0$t <- t
      prepareMplusData(data0,paste0('multilevel/data',as.character(i),'.dat'))
    }
    createModels('multilevel/template_MLM.txt')
    runModels('multilevel',recursive=TRUE)
    mlm <- readModels('multilevel',recursive=TRUE)
    for (i in 1:200) {
      write.table(mlm[[i]]$parameters[[1]],paste0('multilevel/para.est',as.character(i),'.txt'))
    }
  }
}
