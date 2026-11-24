library(foreach)
library(iterators)
n <- c(20,50,100,200)
m <- c(10,20)
#valley
foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd = paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n[j]),
                'm',as.character(m[k]),'/multilevel')
    setwd(wd)
    data <- array(dim=c(9,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      data[,i] <- data0[62,3]
    }
    #画图
    time <- c(1:9)*0.1
    m_data <- rowMeans(data)
    lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
    uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
    beta <- 1/(1+exp((-7)*(time-0.5)))-0.5
    true_effect <- beta*beta
    f_name <- paste0('F:/course/statistic/tvm/S1/valley/mlm_valley_n',
                     as.character(n[j]),'m',as.character(m[k]),'.png')
    png(f_name)
    plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
         ylab = "indirect effect",col='#00bfff',lwd=2)
    lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,true_effect,col='#cc0000',lwd=2)
    abline(h=0,col='#000000')
    dev.off()
    f_name <- paste0('F:/course/statistic/tvm/S1/valley/mlm_valley_every_n',
                     as.character(n[j]),'m',as.character(m[k]),'.png')
    png(f_name)
    plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
         ylab = "indirect effect",col='#cc0000',lwd=2)
    for (i in 1:200) {
      lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
    }
    lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,m_data,col='#00bfff',lwd=2)
    lines(time,true_effect,col='#cc0000',lwd=2)
    abline(h=0,col='#000000')
    dev.off()
  }
}

foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd = paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n[j]),
                'm',as.character(m[k]),'/multilevel')
    setwd(wd)
    data <- array(dim=c(9,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      data[,i] <- data0[62,3]
    }
    #画图
    time <- c(1:9)*0.1
    m_data <- rowMeans(data)
    lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
    uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
    beta <- rep(0.4,9)
    true_effect <- beta*beta
    f_name <- paste0('F:/course/statistic/tvm/S1/flat/mlm_flat_n',
                     as.character(n[j]),'m',as.character(m[k]),'.png')
    png(f_name)
    plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
         ylab = "indirect effect",col='#00bfff',lwd=2)
    lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,true_effect,col='#cc0000',lwd=2)
    abline(h=0,col='#000000')
    dev.off()
    f_name <- paste0('F:/course/statistic/tvm/S1/flat/mlm_flat_every_n',
                     as.character(n[j]),'m',as.character(m[k]),'.png')
    png(f_name)
    plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
         ylab = "indirect effect",col='#cc0000',lwd=2)
    for (i in 1:200) {
      lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
    }
    lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,m_data,col='#00bfff',lwd=2)
    lines(time,true_effect,col='#cc0000',lwd=2)
    abline(h=0,col='#000000')
    dev.off()
  }
}

foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd = paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n[j]),
                'm',as.character(m[k]),'/multilevel')
    setwd(wd)
    data <- array(dim=c(9,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      data[,i] <- data0[62,3]
    }
    #画图
    time <- c(1:9)*0.1
    m_data <- rowMeans(data)
    lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
    uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
    beta <- 0.3*time+0.3*exp(-(3*(time-0.5))^2)
    true_effect <- beta*beta
    f_name <- paste0('F:/course/statistic/tvm/S1/hump/mlm_hump_n',
                     as.character(n[j]),'m',as.character(m[k]),'.png')
    png(f_name)
    plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
         ylab = "indirect effect",col='#00bfff',lwd=2)
    lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,true_effect,col='#cc0000',lwd=2)
    abline(h=0,col='#000000')
    dev.off()
    f_name <- paste0('F:/course/statistic/tvm/S1/hump/mlm_hump_every_n',
                     as.character(n[j]),'m',as.character(m[k]),'.png')
    png(f_name)
    plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.4),"l", xlab="time",
         ylab = "indirect effect",col='#cc0000',lwd=2)
    for (i in 1:200) {
      lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
    }
    lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,m_data,col='#00bfff',lwd=2)
    lines(time,true_effect,col='#cc0000',lwd=2)
    abline(h=0,col='#000000')
    dev.off()
  }
}

foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd = paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n[j]),
                'm',as.character(m[k]),'/multilevel')
    setwd(wd)
    data <- array(dim=c(9,200))
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      data[,i] <- data0[62,3]
    }
    #画图
    time <- c(1:9)*0.1
    m_data <- rowMeans(data)
    lci_data <- apply(data,1,quantile,probs = c(0.025),na.rm = TRUE)
    uci_data <- apply(data,1,quantile,probs = c(0.975),na.rm = TRUE)
    beta <- -0.1+0.5*exp(-time^0.5)
    true_effect <- beta*beta
    f_name <- paste0('F:/course/statistic/tvm/S1/decay/mlm_decay_n',
                     as.character(n[j]),'m',as.character(m[k]),'.png')
    png(f_name)
    plot(time,m_data,xlim=c(0,1),ylim=c(-0.1,0.2),"l", xlab="time",
         ylab = "indirect effect",col='#00bfff',lwd=2)
    lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,true_effect,col='#cc0000',lwd=2)
    abline(h=0,col='#000000')
    dev.off()
    f_name <- paste0('F:/course/statistic/tvm/S1/decay/mlm_decay_every_n',
                     as.character(n[j]),'m',as.character(m[k]),'.png')
    png(f_name)
    plot(time,true_effect,xlim=c(0,1),ylim=c(-0.1,0.2),"l", xlab="time",
         ylab = "indirect effect",col='#cc0000',lwd=2)
    for (i in 1:200) {
      lines(time,data[,i],col='#bbbb00',lty=2,lwd=0.5)
    }
    lines(time,lci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,uci_data,col='#00bfff',lty=2,lwd=2)
    lines(time,m_data,col='#00bfff',lwd=2)
    lines(time,true_effect,col='#cc0000',lwd=2)
    abline(h=0,col='#000000')
    dev.off()
  }
}
