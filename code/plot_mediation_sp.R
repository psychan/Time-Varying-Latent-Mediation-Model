library(foreach)
library(iterators)
n <- c(20,50,100,200)
m <- c(10,20)
#valley
foreach (k=2:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n[k]),
              'm10/sp')
  setwd(wd)
  data1 <- array(dim=c(10,200))
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:10) {
        data1[j,i] <- data0[(420+j),3]
      }
    }
  }
  data <- data1[,!is.na(data1[1,])]
  #画图
  time <- c(1:10)*0.1
  m_data <- rowMeans(data)
  lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
  uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
  beta <- 1/(1+exp((-7)*(time-0.5)))-0.5
  true_effect <- beta*beta
  f_name <- paste0('F:/course/statistic/tvm/S1/valley/sp_valley_n',
                   as.character(n[k]),'m10.png')
  png(f_name)
  plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#00bfff',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  abline(h=0,col='#000000')
  dev.off()
  f_name <- paste0('F:/course/statistic/tvm/S1/valley/sp_valley_every_n',
                   as.character(n[k]),'m10.png')
  png(f_name)
  plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#cc0000',lwd=2)
  for (i in 1:length(data[1,])) {
    lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
  }
  lines(time,m_data,col='#00bfff',lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  abline(h=0,col='#000000')
  dev.off()
}
foreach (k=3:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n[k]),
              'm20/sp')
  setwd(wd)
  data1 <- array(dim=c(20,200))
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:20) {
        data1[j,i] <- data0[(1040+j),3]
      }
    }
  }
  data <- data1[,!is.na(data1[1,])]
  #画图
  time <- c(1:20)*0.05
  m_data <- rowMeans(data)
  lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
  uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
  beta <- 1/(1+exp((-7)*(time-0.5)))-0.5
  true_effect <- beta*beta
  f_name <- paste0('F:/course/statistic/tvm/S1/valley/sp_valley_n',
                   as.character(n[k]),'m20.png')
  png(f_name)
  plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#00bfff',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  abline(h=0,col='#000000')
  dev.off()
  f_name <- paste0('F:/course/statistic/tvm/S1/valley/sp_valley_every_n',
                   as.character(n[k]),'m20.png')
  png(f_name)
  plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#cc0000',lwd=2)
  for (i in 1:length(data[1,])) {
    lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
  }
  lines(time,m_data,col='#00bfff',lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  abline(h=0,col='#000000')
  dev.off()
}
#flat
foreach (k=2:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n[k]),
              'm10/sp')
  setwd(wd)
  data1 <- array(dim=c(10,200))
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:10) {
        data1[j,i] <- data0[(420+j),3]
      }
    }
  }
  data <- data1[,!is.na(data1[1,])]
  #画图
  time <- c(1:10)*0.1
  m_data <- rowMeans(data)
  lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
  uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
  beta <- rep(0.4,10)
  true_effect <- beta*beta
  f_name <- paste0('F:/course/statistic/tvm/S1/flat/sp_flat_n',
                   as.character(n[k]),'m10.png')
  png(f_name)
  plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#00bfff',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  abline(h=0,col='#000000')
  dev.off()
  f_name <- paste0('F:/course/statistic/tvm/S1/flat/sp_flat_every_n',
                   as.character(n[k]),'m10.png')
  png(f_name)
  plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#cc0000',lwd=2)
  for (i in 1:length(data[1,])) {
    lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
  }
  lines(time,m_data,col='#00bfff',lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  abline(h=0,col='#000000')
  dev.off()
}
foreach (k=4:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n[k]),
              'm20/sp')
  setwd(wd)
  data1 <- array(dim=c(20,200))
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:20) {
        data1[j,i] <- data0[(1040+j),3]
      }
    }
  }
  data <- data1[,!is.na(data1[1,])]
  #画图
  time <- c(1:20)*0.05
  m_data <- rowMeans(data)
  lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
  uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
  beta <- rep(0.4,20)
  true_effect <- beta*beta
  f_name <- paste0('F:/course/statistic/tvm/S1/flat/sp_flat_n',
                   as.character(n[k]),'m20.png')
  png(f_name)
  plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#00bfff',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  abline(h=0,col='#000000')
  dev.off()
  f_name <- paste0('F:/course/statistic/tvm/S1/flat/sp_flat_every_n',
                   as.character(n[k]),'m20.png')
  png(f_name)
  plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#cc0000',lwd=2)
  for (i in 1:length(data[1,])) {
    lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
  }
  lines(time,m_data,col='#00bfff',lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  abline(h=0,col='#000000')
  dev.off()
}
#hump
foreach (k=2:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n[k]),
              'm10/sp')
  setwd(wd)
  data1 <- array(dim=c(10,200))
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:10) {
        data1[j,i] <- data0[(420+j),3]
      }
    }
  }
  data <- data1[,!is.na(data1[1,])]
  #画图
  time <- c(1:10)*0.1
  m_data <- rowMeans(data)
  lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
  uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
  beta <- 0.3*time+0.3*exp(-(3*(time-0.5))^2)
  true_effect <- beta*beta
  f_name <- paste0('F:/course/statistic/tvm/S1/hump/sp_hump_n',
                   as.character(n[k]),'m10.png')
  png(f_name)
  plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#00bfff',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  abline(h=0,col='#000000')
  dev.off()
  f_name <- paste0('F:/course/statistic/tvm/S1/hump/sp_hump_every_n',
                   as.character(n[k]),'m10.png')
  png(f_name)
  plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#cc0000',lwd=2)
  for (i in 1:length(data[1,])) {
    lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
  }
  lines(time,m_data,col='#00bfff',lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  abline(h=0,col='#000000')
  dev.off()
}
foreach (k=4:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n[k]),
              'm20/sp')
  setwd(wd)
  data1 <- array(dim=c(20,200))
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:20) {
        data1[j,i] <- data0[(1040+j),3]
      }
    }
  }
  data <- data1[,!is.na(data1[1,])]
  #画图
  time <- c(1:20)*0.05
  m_data <- rowMeans(data)
  lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
  uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
  beta <- 0.3*time+0.3*exp(-(3*(time-0.5))^2)
  true_effect <- beta*beta
  f_name <- paste0('F:/course/statistic/tvm/S1/hump/sp_hump_n',
                   as.character(n[k]),'m20.png')
  png(f_name)
  plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#00bfff',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  abline(h=0,col='#000000')
  dev.off()
  f_name <- paste0('F:/course/statistic/tvm/S1/hump/sp_hump_every_n',
                   as.character(n[k]),'m20.png')
  png(f_name)
  plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
       ylab = "indirect effect",col='#cc0000',lwd=2)
  for (i in 1:length(data[1,])) {
    lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
  }
  lines(time,m_data,col='#00bfff',lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  abline(h=0,col='#000000')
  dev.off()
}
#decay
foreach (k=2:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n[k]),
              'm10/sp')
  setwd(wd)
  data1 <- array(dim=c(10,200))
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:10) {
        data1[j,i] <- data0[(420+j),3]
      }
    }
  }
  data <- data1[,!is.na(data1[1,])]
  #画图
  time <- c(1:10)*0.1
  m_data <- rowMeans(data)
  lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
  uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
  beta <- -0.1+0.5*exp(-time^0.5)
  true_effect <- beta*beta
  f_name <- paste0('F:/course/statistic/tvm/S1/decay/sp_decay_n',
                   as.character(n[k]),'m10.png')
  png(f_name)
  plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.2),"l", xlab="time",
       ylab = "indirect effect",col='#00bfff',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  abline(h=0,col='#000000')
  dev.off()
  f_name <- paste0('F:/course/statistic/tvm/S1/decay/sp_decay_every_n',
                   as.character(n[k]),'m10.png')
  png(f_name)
  plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.2),"l", xlab="time",
       ylab = "indirect effect",col='#cc0000',lwd=2)
  for (i in 1:length(data[1,])) {
    lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
  }
  lines(time,m_data,col='#00bfff',lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  abline(h=0,col='#000000')
  dev.off()
}
foreach (k=3:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n[k]),
              'm20/sp')
  setwd(wd)
  data1 <- array(dim=c(20,200))
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:20) {
        data1[j,i] <- data0[(1040+j),3]
      }
    }
  }
  data <- data1[,!is.na(data1[1,])]
  #画图
  time <- c(1:20)*0.05
  m_data <- rowMeans(data)
  lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
  uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
  beta <- -0.1+0.5*exp(-time^0.5)
  true_effect <- beta*beta
  f_name <- paste0('F:/course/statistic/tvm/S1/decay/sp_decay_n',
                   as.character(n[k]),'m20.png')
  png(f_name)
  plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.2),"l", xlab="time",
       ylab = "indirect effect",col='#00bfff',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  abline(h=0,col='#000000')
  dev.off()
  f_name <- paste0('F:/course/statistic/tvm/S1/decay/sp_decay_every_n',
                   as.character(n[k]),'m20.png')
  png(f_name)
  plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.2),"l", xlab="time",
       ylab = "indirect effect",col='#cc0000',lwd=2)
  for (i in 1:length(data[1,])) {
    lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
  }
  lines(time,m_data,col='#00bfff',lwd=2)
  lines(time,true_effect,col='#cc0000',lwd=2)
  lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
  lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
  abline(h=0,col='#000000')
  dev.off()
}
