n_sample = c(20,50,100,200)
m_measure = c(10,20)
type = c("flat","decay","hump","valley")
source('F:/course/statistic/tvm/S1.1/DynMed.R')
for(it in 1:4){
  for(n in 1:4){
    for(m in 1:2){
      library(parallel)
      clnum <- detectCores()
      
      library(foreach)
      library(iterators)
      library(doParallel)
      cl <- makeCluster(clnum-2)
      registerDoParallel(cl)
      foreach(i=1:200,.errorhandling="remove") %dopar% {
        filename <- paste0('F:/course/statistic/tvm/S1/',type[it],'/n',
                           as.character(n_sample[n]),'m',as.character(m_measure[m]),
                           '/data_',type[it],as.character(i),'.txt')
        data0 <- read.table(file=filename,header=TRUE,sep='\t')
        wd <- paste0('F:/course/statistic/tvm/S1.1/',type[it],'/n',
                     as.character(n_sample[n]),'m',as.character(m_measure[m]),'/')
        results <- DynMed(data=data0,x=rowMeans(data0[,2:4]),
                          m=rowMeans(data0[,5:7]),y=rowMeans(data0[,8:10]),
                          cluster = data0[,1],t=data0[,11])
        f1 <- paste0(wd,'raw_nocfa',as.character(i),'.txt')
        f2 <- paste0(wd,'remain_nocfa',as.character(i),'.txt')
        f3 <- paste0(wd,'est_nocfa',as.character(i),'.txt')
        write.table(results$raw,file=f1,sep='\t',quote=FALSE,col.names = TRUE,row.names = FALSE)
        write.table(results$remain,file=f2,sep='\t',quote=FALSE,col.names = TRUE,row.names = FALSE)
        write.table(results$est,file=f3,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE)
      }
      stopCluster(cl)
    }
  }
}
