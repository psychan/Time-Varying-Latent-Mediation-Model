library(foreach)
library(iterators)
library(ggplot2)
library(patchwork)
library(gridExtra)

n <- c(20,50,100,200)
m <- c(10,20)

tvm_cfa_valley <- tvm_ncfa_valley <- mlm_valley <- array(dim=c(9,200,4,2))
foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd <- paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data2 <- data1[68:76,1]
      tvm_cfa_valley[,i,j,k] <- data2
    }
    
    wd <- paste0('F:/course/statistic/tvm/S1.1/valley/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1.1/valley/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est_nocfa',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data2 <- data1[50:58,1]
      tvm_ncfa_valley[,i,j,k] <- data2
    }
    
    wd = paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n[j]),
                'm',as.character(m[k]),'/multilevel')
    setwd(wd)
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      mlm_valley[,i,j,k] <- data0[62,3]
    }
  }
}
sp_valley_m10 <- array(dim=c(10,200,4))
foreach (k=2:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n[k]),
              'm10/sp')
  setwd(wd)
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:10) {
        sp_valley_m10[j,i,k] <- data0[(420+j),3]
      }
    }
  }
}

sp_valley_m20 <- array(dim=c(20,200,4))
foreach (k=3:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/valley/n',as.character(n[k]),
              'm20/sp')
  setwd(wd)
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:20) {
        sp_valley_m20[j,i,k] <- data0[(1040+j),3]
      }
    }
  }
}
m_tvm_cfa_valley <- apply(tvm_cfa_valley,1,mean,na.rm = TRUE)
m_tvm_ncfa_valley <- apply(tvm_ncfa_valley,1,mean,na.rm = TRUE)
m_mlm_valley <- apply(mlm_valley,1,mean,na.rm = TRUE)
m_sp_valley_m10 <- apply(sp_valley_m10,1,mean,na.rm = TRUE)
m_sp_valley_m20 <- apply(sp_valley_m20,1,mean,na.rm = TRUE)
time_m10 <- c(1:10)*0.1
time_m20 <- c(1:20)*0.05
time_9 <- c(1:9)*0.1
beta <- 1/(1+exp((-7)*(time_m20-0.5)))-0.5
true_effect <- beta*beta
df0 <- data.frame(time=time_m20,value=true_effect,model='True Effect')
df1 <- data.frame(time=time_9,value=m_tvm_cfa_valley,model='TVM_CFA')
df2 <- data.frame(time=time_9,value=m_tvm_ncfa_valley,model='TVM_NCFA')
df3 <- data.frame(time=time_9,value=m_mlm_valley,model='MLM')
df4 <- data.frame(time=time_m10,value=m_sp_valley_m10,model='SP_M10')
df5 <- data.frame(time=time_m20,value=m_sp_valley_m20,model='SP_M20')
combined_df <- rbind(df0,df1,df2,df3,df4,df5)
combined_df$model <- factor(combined_df$model,levels = c("True Effect","TVM_CFA","TVM_NCFA","MLM","SP_M10","SP_M20"))

p1 <- ggplot(data = combined_df,mapping = aes(x=time,y=value,group=model,
                                              colour=model,linetype=model)) +
  xlab("Time")+ylab("Indirect Effect")+
  scale_color_manual(values=c("grey23","#0070C0", "#FF5733", "#90d743", "#e9b012", 
                              "#d17c0a"),name="model")+
  scale_linetype_manual(values=c("solid","dashed","dashed","dashed","dashed",
                                 "dashed"))+
  geom_line(lineend="butt",linewidth=1.25)+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))


tvm_cfa_hump <- tvm_ncfa_hump <- mlm_hump <- array(dim=c(9,200,4,2))
foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd <- paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data2 <- data1[68:76,1]
      tvm_cfa_hump[,i,j,k] <- data2
    }
    
    wd <- paste0('F:/course/statistic/tvm/S1.1/hump/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1.1/hump/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est_nocfa',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data2 <- data1[50:58,1]
      tvm_ncfa_hump[,i,j,k] <- data2
    }
    
    wd = paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n[j]),
                'm',as.character(m[k]),'/multilevel')
    setwd(wd)
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      mlm_hump[,i,j,k] <- data0[62,3]
    }
  }
}
sp_hump_m10 <- array(dim=c(10,200,4))
foreach (k=2:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n[k]),
              'm10/sp')
  setwd(wd)
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:10) {
        sp_hump_m10[j,i,k] <- data0[(420+j),3]
      }
    }
  }
}

sp_hump_m20 <- array(dim=c(20,200,4))
foreach (k=3:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/hump/n',as.character(n[k]),
              'm20/sp')
  setwd(wd)
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:20) {
        sp_hump_m20[j,i,k] <- data0[(1040+j),3]
      }
    }
  }
}
m_tvm_cfa_hump <- apply(tvm_cfa_hump,1,mean,na.rm = TRUE)
m_tvm_ncfa_hump <- apply(tvm_ncfa_hump,1,mean,na.rm = TRUE)
m_mlm_hump <- apply(mlm_hump,1,mean,na.rm = TRUE)
m_sp_hump_m10 <- apply(sp_hump_m10,1,mean,na.rm = TRUE)
m_sp_hump_m20 <- apply(sp_hump_m20,1,mean,na.rm = TRUE)
time_m10 <- c(1:10)*0.1
time_m20 <- c(1:20)*0.05
time_9 <- c(1:9)*0.1
beta <- 0.3*time_m20+0.3*exp(-(3*(time_m20-0.5))^2)
true_effect <- beta*beta
df0 <- data.frame(time=time_m20,value=true_effect,model='True Effect')
df1 <- data.frame(time=time_9,value=m_tvm_cfa_hump,model='TVM_CFA')
df2 <- data.frame(time=time_9,value=m_tvm_ncfa_hump,model='TVM_NCFA')
df3 <- data.frame(time=time_9,value=m_mlm_hump,model='MLM')
df4 <- data.frame(time=time_m10,value=m_sp_hump_m10,model='SP_M10')
df5 <- data.frame(time=time_m20,value=m_sp_hump_m20,model='SP_M20')
combined_df <- rbind(df0,df1,df2,df3,df4,df5)
combined_df$model <- factor(combined_df$model,levels = c("True Effect","TVM_CFA","TVM_NCFA","MLM","SP_M10","SP_M20"))

p2 <- ggplot(data = combined_df,mapping = aes(x=time,y=value,group=model,
                                              colour=model,linetype=model)) +
  xlab("Time")+ylab("Indirect Effect")+
  scale_color_manual(values=c("grey23","#0070C0", "#FF5733", "#90d743", "#e9b012", 
                              "#d17c0a"),name="model")+
  scale_linetype_manual(values=c("solid","dashed","dashed","dashed","dashed",
                                 "dashed"))+
  geom_line(lineend="butt",linewidth=1.25)+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))

tvm_cfa_decay <- tvm_ncfa_decay <- mlm_decay <- array(dim=c(9,200,4,2))
foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd <- paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data2 <- data1[68:76,1]
      tvm_cfa_decay[,i,j,k] <- data2
    }
    
    wd <- paste0('F:/course/statistic/tvm/S1.1/decay/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1.1/decay/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est_nocfa',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data2 <- data1[50:58,1]
      tvm_ncfa_decay[,i,j,k] <- data2
    }
    
    wd = paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n[j]),
                'm',as.character(m[k]),'/multilevel')
    setwd(wd)
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      mlm_decay[,i,j,k] <- data0[62,3]
    }
  }
}
sp_decay_m10 <- array(dim=c(10,200,4))
foreach (k=2:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n[k]),
              'm10/sp')
  setwd(wd)
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:10) {
        sp_decay_m10[j,i,k] <- data0[(420+j),3]
      }
    }
  }
}

sp_decay_m20 <- array(dim=c(20,200,4))
foreach (k=3:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/decay/n',as.character(n[k]),
              'm20/sp')
  setwd(wd)
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:20) {
        sp_decay_m20[j,i,k] <- data0[(1040+j),3]
      }
    }
  }
}
m_tvm_cfa_decay <- apply(tvm_cfa_decay,1,mean,na.rm = TRUE)
m_tvm_ncfa_decay <- apply(tvm_ncfa_decay,1,mean,na.rm = TRUE)
m_mlm_decay <- apply(mlm_decay,1,mean,na.rm = TRUE)
m_sp_decay_m10 <- apply(sp_decay_m10,1,mean,na.rm = TRUE)
m_sp_decay_m20 <- apply(sp_decay_m20,1,mean,na.rm = TRUE)
time_m10 <- c(1:10)*0.1
time_m20 <- c(1:20)*0.05
time_9 <- c(1:9)*0.1
beta <- -0.1+0.5*exp(-time_m20^0.5)
true_effect <- beta*beta
df0 <- data.frame(time=time_m20,value=true_effect,model='True Effect')
df1 <- data.frame(time=time_9,value=m_tvm_cfa_decay,model='TVM_CFA')
df2 <- data.frame(time=time_9,value=m_tvm_ncfa_decay,model='TVM_NCFA')
df3 <- data.frame(time=time_9,value=m_mlm_decay,model='MLM')
df4 <- data.frame(time=time_m10,value=m_sp_decay_m10,model='SP_M10')
df5 <- data.frame(time=time_m20,value=m_sp_decay_m20,model='SP_M20')
combined_df <- rbind(df0,df1,df2,df3,df4,df5)
combined_df$model <- factor(combined_df$model,levels = c("True Effect","TVM_CFA","TVM_NCFA","MLM","SP_M10","SP_M20"))

p3 <- ggplot(data = combined_df,mapping = aes(x=time,y=value,group=model,
                                              colour=model,linetype=model)) +
  xlab("Time")+ylab("Indirect Effect")+
  scale_color_manual(values=c("grey23","#0070C0", "#FF5733", "#90d743", "#e9b012", 
                              "#d17c0a"),name="model")+
  scale_linetype_manual(values=c("solid","dashed","dashed","dashed","dashed",
                                 "dashed"))+
  geom_line(lineend="butt",linewidth=1.25)+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))

tvm_cfa_flat <- tvm_ncfa_flat <- mlm_flat <- array(dim=c(9,200,4,2))
foreach (j=1:4,.errorhandling="remove") %do% {
  foreach (k=1:2,.errorhandling="remove") %do% {
    wd <- paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data2 <- data1[68:76,1]
      tvm_cfa_flat[,i,j,k] <- data2
    }
    
    wd <- paste0('F:/course/statistic/tvm/S1.1/flat/n',as.character(n[j]),'m',as.character(m[k]))
    setwd(wd)
    for (i in 1:200) {
      filename <- paste0('F:/course/statistic/tvm/S1.1/flat/n',as.character(n[j]),'m',as.character(m[k]),
                         '/est_nocfa',as.character(i),'.txt')
      data0 <- read.csv(file=filename,header=TRUE,sep='\t')
      data1 <- as.matrix(data0)
      data2 <- data1[50:58,1]
      tvm_ncfa_flat[,i,j,k] <- data2
    }
    
    wd = paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n[j]),
                'm',as.character(m[k]),'/multilevel')
    setwd(wd)
    for (i in 1:200) {
      data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
      mlm_flat[,i,j,k] <- data0[62,3]
    }
  }
}
sp_flat_m10 <- array(dim=c(10,200,4))
foreach (k=2:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n[k]),
              'm10/sp')
  setwd(wd)
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:10) {
        sp_flat_m10[j,i,k] <- data0[(420+j),3]
      }
    }
  }
}

sp_flat_m20 <- array(dim=c(20,200,4))
foreach (k=3:4,.errorhandling="remove") %do% {
  wd = paste0('F:/course/statistic/tvm/S1/flat/n',as.character(n[k]),
              'm20/sp')
  setwd(wd)
  for (i in 1:200) {
    data0 <- read.table(paste0('para.est',as.character(i),'.txt'),header=TRUE)
    if(ncol(data0)==8){
      for (j in 1:20) {
        sp_flat_m20[j,i,k] <- data0[(1040+j),3]
      }
    }
  }
}
m_tvm_cfa_flat <- apply(tvm_cfa_flat,1,mean,na.rm = TRUE)
m_tvm_ncfa_flat <- apply(tvm_ncfa_flat,1,mean,na.rm = TRUE)
m_mlm_flat <- apply(mlm_flat,1,mean,na.rm = TRUE)
m_sp_flat_m10 <- apply(sp_flat_m10,1,mean,na.rm = TRUE)
m_sp_flat_m20 <- apply(sp_flat_m20,1,mean,na.rm = TRUE)
time_m10 <- c(1:10)*0.1
time_m20 <- c(1:20)*0.05
time_9 <- c(1:9)*0.1
beta <- rep(0.4,20)
true_effect <- beta*beta
df0 <- data.frame(time=time_m20,value=true_effect,model='True Effect')
df1 <- data.frame(time=time_9,value=m_tvm_cfa_flat,model='TVM_CFA')
df2 <- data.frame(time=time_9,value=m_tvm_ncfa_flat,model='TVM_NCFA')
df3 <- data.frame(time=time_9,value=m_mlm_flat,model='MLM')
df4 <- data.frame(time=time_m10,value=m_sp_flat_m10,model='SP_M10')
df5 <- data.frame(time=time_m20,value=m_sp_flat_m20,model='SP_M20')
combined_df <- rbind(df0,df1,df2,df3,df4,df5)
combined_df$model <- factor(combined_df$model,levels = c("True Effect","TVM_CFA","TVM_NCFA","MLM","SP_M10","SP_M20"))

p4 <- ggplot(data = combined_df,mapping = aes(x=time,y=value,group=model,
                                              colour=model,linetype=model)) +
  xlab("Time")+ylab("Indirect Effect")+
  scale_color_manual(values=c("grey23","#0070C0", "#FF5733", "#90d743", "#e9b012", 
                              "#d17c0a"),name="model")+
  scale_linetype_manual(values=c("solid","dashed","dashed","dashed","dashed",
                                 "dashed"))+
  geom_line(lineend="butt",linewidth=1.25)+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))

p5 <- (p1 + p2)/(p3 + p4)
ggsave(filename ="F:/course/statistic/tvm/p_mediation.png",
       dpi=300,width = 12,height = 12 )

d_f <- "F:/course/statistic/tvm/number/compare.txt"
data0 <- read.csv(file = d_f,header = TRUE,sep='\t')
data0 <- data0[data0$coefficient=='indirect effect',]
data0$sample.size <- factor(data0$sample.size,levels = c(10,15,20,25,50,100,200))
data0$measurement.occasion <- factor(data0$measurement.occasion,
                                     levels = c(10,15,20,25,50,100,200))
data0$type <- factor(data0$type,levels = c("decay","hump","valley"))
p_cr <- ggplot(data = data0,
               mapping = aes(x=sample.size,y=cover.rate,group=measurement.occasion,
                             colour=measurement.occasion,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(0,1,2,5,6,7,9),name="Number of Measurement")+
  scale_color_manual(values=c("#90d743", "#5ec962", "#35b779", "#21918c", 
                              "#3b528b", "#482878", "#440154"),
                     name="Number of Measurement")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type))+ylim(0,1)+
  geom_hline(yintercept = 0.9,colour="grey23",linetype="dashed")+
  labs(x="Sample Size",y="95% Cover Rate")+
  theme(legend.position =  "none",
        legend.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
p_mse <- ggplot(data = data0,
                mapping = aes(x=sample.size,y=MSE,group=measurement.occasion,
                              colour=measurement.occasion,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(0,1,2,5,6,7,9),name="Number of Measurement")+
  scale_color_manual(values=c("#90d743", "#5ec962", "#35b779", "#21918c", 
                              "#3b528b", "#482878", "#440154"),
                     name="Number of Measurement")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type))+ylim(0,0.04)+
  geom_hline(yintercept = 0,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("MSE")+
  theme(legend.position =  "none",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
p_rb <- ggplot(data=data0,
               mapping=aes(x=sample.size,y=relative.bias,group=measurement.occasion,
                           colour=measurement.occasion,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(0,1,2,5,6,7,9),name="Number of Measurement")+
  scale_color_manual(values=c("#90d743", "#5ec962", "#35b779", "#21918c", 
                              "#3b528b", "#482878", "#440154"),
                     name="Number of Measurement")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type))+ylim(-0.25,0.25)+
  geom_hline(yintercept = 0,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("Relative Bias")+
  theme(legend.position =  "none",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
p_tpr <- ggplot(data = data0,
                mapping = aes(x=sample.size,y=true.positive.rate,group=measurement.occasion,
                              colour=measurement.occasion,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(0,1,2,5,6,7,9),name="Number of Measurement")+
  scale_color_manual(values=c("#90d743", "#5ec962", "#35b779", "#21918c", 
                              "#3b528b", "#482878", "#440154"),
                     name="Number of Measurement")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type))+ylim(0,1)+
  geom_hline(yintercept = 0.8,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("True Positive Rate")+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
p_sum <- p_cr/p_mse/p_rb/p_tpr
ggsave(filename ="F:/course/statistic/tvm/p_number1.png",
       dpi=300,width = 9,height = 12 )

data1 <- data0[!is.na(data0$true.negative.rate),]
data1$fnr <- 1-data1$true.negative.rate
p_fpr <- ggplot(data = data1,
                mapping = aes(x=sample.size,y=fnr,group=measurement.occasion,
                              colour=measurement.occasion,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(0,1,2,5,6,7,9),name="Number of Measurement")+
  scale_color_manual(values=c("#90d743", "#5ec962", "#35b779", "#21918c", 
                              "#3b528b", "#482878", "#440154"),
                     name="Number of Measurement")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type))+ylim(0,0.2)+
  geom_hline(yintercept = 0.05,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("False Positive Rate")+
  theme(legend.position =  "right",legend.box = "vertical",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
ggsave(filename ="F:/course/statistic/tvm/p_fpr.png",
       dpi=300,width = 5,height = 3 )

library(ggplot2)
library(gridExtra)
d_f <- "F:/course/statistic/tvm/S1.1/results.txt"
data0 <- read.csv(file = d_f,header = TRUE,sep='\t')
data0 <- data0[data0$coefficient=='indirect effect',]
data0$methods <- factor(data0$methods,levels = c('without measurement error',
                                                 'considering measurement error'))
data0$sample.size <- factor(data0$sample.size,levels = c(20,50,100,200))
data0$measurement.occasion <- factor(data0$measurement.occasion,levels = c(10,20))
data0$type <- factor(data0$type,levels = c("flat","decay","hump","valley"))
p_cr <- ggplot(data = data0,
               mapping = aes(x=sample.size,y=cover.rate,
                             group=interaction(methods,measurement.occasion),
                             colour=methods,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(2,7),name="Number of Measurement")+
  scale_color_manual(values=c("#FF5733","#0070C0"),name="Measurement Error")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type))+ylim(0,1)+
  geom_hline(yintercept = 0.9,colour="grey23",linetype="dashed")+
  labs(x="Sample Size",y="95% Cover Rate")+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
p_mse <- ggplot(data = data0,
                mapping = aes(x=sample.size,y=MSE,
                              group=interaction(methods,measurement.occasion),
                              colour=methods,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(2,7),name="Number of Measurement")+
  scale_color_manual(values=c("#FF5733","#0070C0"),name="Measurement Error")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type))+ylim(0,0.015)+
  geom_hline(yintercept = 0,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("MSE")+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
p_tpr <- ggplot(data = data0,
                mapping = aes(x=sample.size,y=true.positive.rate,
                              group=interaction(methods,measurement.occasion),
                              colour=methods,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(2,7),name="Number of Measurement")+
  scale_color_manual(values=c("#FF5733","#0070C0"),name="Measurement Error")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type))+ylim(0,1)+
  geom_hline(yintercept = 0.8,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("True Positive Rate")+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
p_rb <- ggplot(data=data0,
               mapping=aes(x=sample.size,y=relative.bias,
                           group=interaction(methods,measurement.occasion),
                           colour=methods,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(2,7),name="Number of Measurement")+
  scale_color_manual(values=c("#FF5733","#0070C0"),name="Measurement Error")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type))+ylim(-0.3,0.1)+
  geom_hline(yintercept = 0,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("Relative Bias")+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
p_sum <- p_cr/p_mse/p_rb/p_tpr
ggsave(filename ="F:/course/statistic/tvm/p_tvm1.png",
       dpi=300,width = 9,height = 12 )

data1 <- data0[!is.na(data0$true.negative.rate),]
data1$fnr <- 1-data1$true.negative.rate
p_fpr <- ggplot(data = data1,
                mapping = aes(x=sample.size,y=fnr,group=measurement.occasion,
                              colour=methods,shape=measurement.occasion))+
  geom_point(size=1.8)+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(2,7),name="Number of Measurement")+
  scale_color_manual(values=c("#FF5733","#0070C0"),name="Measurement Error")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type))+ylim(0,0.2)+
  geom_hline(yintercept = 0.05,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("False Positive Rate")+
  theme(legend.position =  "right",legend.box = "vertical",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
ggsave(filename ="F:/course/statistic/tvm/p_tvm2.png",
       dpi=300,width = 5,height = 3 )
