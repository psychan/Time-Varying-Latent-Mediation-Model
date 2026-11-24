#plot
library(ggplot2)
library(gridExtra)
d_f <- "F:/course/statistic/tvm/number/compare.txt"
data0 <- read.csv(file = d_f,header = TRUE,sep='\t')
data0$sample.size <- factor(data0$sample.size,levels = c(10,15,20,25,50,100,200))
data0$measurement.occasion <- factor(data0$measurement.occasion,
                                     levels = c(10,15,20,25,50,100,200))
data0$type <- factor(data0$type,levels = c("decay","hump","valley"))
data0$coefficient <- factor(data0$coefficient,
                            levels = c("loadings","measurement error","intercept",
                                       "path coefficients","regression residual",
                                       "indirect effect"))
p_cr <- ggplot(data = data0,
               mapping = aes(x=sample.size,y=cover.rate,group=measurement.occasion,
                             colour=measurement.occasion,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(0,1,2,5,6,7,9),name="Number of Measurement")+
  scale_color_manual(values=c("#BFD9F0","#8FBCE6","#5F9FDC","#3F82D2","#2665C8",
                              "#1648BE","#0A2BA4"),name="Number of Measurement")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type),rows = vars(coefficient))+ylim(0,1)+
  geom_hline(yintercept = 0.9,colour="grey23",linetype="dashed")+
  labs(x="Sample Size",y="95% Cover Rate")+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
ggsave(filename ="F:/course/statistic/tvm/number/p_cr.png",
       dpi=300,width = 12,height = 9 )


p_mse <- ggplot(data = data0,
                mapping = aes(x=sample.size,y=MSE,group=measurement.occasion,
                              colour=measurement.occasion,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(0,1,2,5,6,7,9),name="Number of Measurement")+
  scale_color_manual(values=c("#BFD9F0","#8FBCE6","#5F9FDC","#3F82D2","#2665C8",
                              "#1648BE","#0A2BA4"),name="Number of Measurement")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type),rows = vars(coefficient))+ylim(0,0.04)+
  geom_hline(yintercept = 0,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("MSE")+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
ggsave(filename ="F:/course/statistic/tvm/number/p_mse.png",
       dpi=300,width = 12,height = 9 )

p_tpr <- ggplot(data = data0,
                mapping = aes(x=sample.size,y=true.positive.rate,group=measurement.occasion,
                              colour=measurement.occasion,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(0,1,2,5,6,7,9),name="Number of Measurement")+
  scale_color_manual(values=c("#BFD9F0","#8FBCE6","#5F9FDC","#3F82D2","#2665C8",
                              "#1648BE","#0A2BA4"),name="Number of Measurement")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type),rows = vars(coefficient))+ylim(0,1)+
  geom_hline(yintercept = 0.8,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("True Positive Rate")+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
ggsave(filename ="F:/course/statistic/tvm/number/p_tpr.png",
       dpi=300,width = 12,height = 9 )

p_rb <- ggplot(data=data0,
               mapping=aes(x=sample.size,y=relative.bias,group=measurement.occasion,
                           colour=measurement.occasion,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(0,1,2,5,6,7,9),name="Number of Measurement")+
  scale_color_manual(values=c("#BFD9F0","#8FBCE6","#5F9FDC","#3F82D2","#2665C8",
                              "#1648BE","#0A2BA4"),name="Number of Measurement")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type),rows = vars(coefficient))+ylim(-0.25,0.25)+
  geom_hline(yintercept = 0,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("Relative Bias")+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
ggsave(filename ="F:/course/statistic/tvm/number/p_rb.png",
       dpi=300,width = 12,height = 9 )

data1 <- data0[!is.na(data0$true.negative.rate),]
p_tnr <- ggplot(data = data1,
                mapping = aes(x=sample.size,y=true.negative.rate,group=measurement.occasion,
                              colour=measurement.occasion,shape=measurement.occasion))+
  geom_point(size=1.8)+
  scale_shape_manual(values = c(0,1,2,5,6,7,9),name="Number of Measurement")+
  scale_color_manual(values=c("#BFD9F0","#8FBCE6","#5F9FDC","#3F82D2","#2665C8",
                              "#1648BE","#0A2BA4"),name="Number of Measurement")+
  geom_line(lineend="butt",linewidth=0.5)+
  facet_grid(cols = vars(type),rows = vars(coefficient))+ylim(0.8,1)+
  geom_hline(yintercept = 0.95,colour="grey23",linetype="dashed")+
  xlab("Sample Size")+ylab("True Negative Rate")+
  theme(legend.position =  "bottom",legend.box = "horizontal",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray84",linetype = 1,linewidth = 0.08),
        panel.border = element_rect(fill = NA,color = "grey40",linetype = 1),
        strip.background = element_rect(fill = "gray80",color = "grey40",linetype = 1))
ggsave(filename ="F:/course/statistic/tvm/number/p_tnr.png",
       dpi=300,width = 4,height = 6 )
