
# TL - 200623
library(ggplot2)
library("ggrepel")
library("grid")
library("gridExtra")
setwd("/home/thibaultleroy/Rose/Papier/Fig4")


#### PCA with info TN
pca_TN=read.table("res_PCA_071021_withinfo_withTN.txt",h=T)

plotA<-ggplot(pca_TN,aes(x=Axis1,y=Axis2,colour=as.factor(genetic_group)))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(ID), fill = as.factor(genetic_group)), color = 'black', size = 2)+
  labs(colour="Genetic Groups (SSR Liorzou)",fill="Genetic Groups (SSR Liorzou)")

plotB<-ggplot(pca_TN,aes(x=Axis1,y=Axis2,colour=TN2015))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = TN2014, fill = TN2014), color = 'black', size = 3.5)+
  scale_fill_gradient2(name = "Score A1 (BS2014)",low = "navyblue", mid="lightgoldenrod1",high = "darkred")+
  scale_colour_gradient2(name = "Score A1 (BS2014)",low = "navyblue", mid="lightgoldenrod1",high = "darkred")


plotC<-ggplot(pca_TN,aes(x=Axis1,y=Axis2,colour=TN2016))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = TN2015, fill = TN2015), color = 'black', size = 3.5)+
  scale_fill_gradient2(name = "Score A2 (BS2015)",low = "navyblue", mid="lightgoldenrod1",high = "darkred")+
  scale_colour_gradient2(name = "Score A2 (BS2015)",low = "navyblue", mid="lightgoldenrod1",high = "darkred")

plotD<-ggplot(pca_TN,aes(x=Axis1,y=Axis2,colour=TN2014))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = TN2016, fill = TN2016), color = 'black', size = 3.5)+
  scale_fill_gradient2(name = "Score A3 (BS2016)",low = "navyblue", mid="lightgoldenrod1",high = "darkred")+
  scale_colour_gradient2(name = "Score A3 (BS2016)",low = "navyblue", mid="lightgoldenrod1",high = "darkred")

plotE<-ggplot(pca_TN,aes(x=Axis1,y=TN2014))+
  stat_smooth(method='lm', formula = y~poly(x,2))+
  geom_point(size=2.5)+
  xlab("Axis 1 (PCA)")+
  ylab("BS suscept. score\nA1 (BS2014)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  annotate(geom="text", x=20, y=3.5, label="model:score~poly(Axis1,2),p-value=1.02e-06,R2=0.134", color="navyblue")

plotF<-ggplot(pca_TN,aes(x=Axis1,y=TN2015))+
  stat_smooth(method='lm', formula = y~poly(x,2))+
  geom_point(size=2.5)+
  xlab("Axis 1 (PCA)")+
  ylab("BS suscept. score\nA2 (BS2015)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  annotate(geom="text", x=20, y=3.5, label="model:score~poly(Axis1,2),p-value=1.03e-09,R2=0.195", color="navyblue")


plotG<-ggplot(pca_TN,aes(x=Axis1,y=TN2016))+
  stat_smooth(method='lm', formula = y~poly(x,2))+
  geom_point(size=2.5)+
  xlab("Axis 1 (PCA)")+
  ylab("BS suscept. score\nA3 (BS2016)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  annotate(geom="text", x=20, y=3.5, label="model:score~poly(Axis1,2),p-value=3.19e-06,R2=0.124", color="navyblue")




pdf(height=20,width=16,file="PCA_TN_204rosesCC_200623.pdf")
#grid.arrange(plotA, plotB, plotC, widths=c(2,1,1))

pushViewport(viewport(layout=grid.layout(9,10)))
vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col = y)
#print(plotA, vp=vplayout(1,1:3))
print(plotE, vp=vplayout(1,1:9))
print(plotB, vp=vplayout(2:3,1:10))
print(plotF, vp=vplayout(4,1:9))
print(plotC, vp=vplayout(5:6,1:10))
print(plotG, vp=vplayout(7,1:9))
print(plotD, vp=vplayout(8:9,1:10))

dev.off()

summary(lm(formula=pca_TN$TN2014~poly(pca_TN$Axis1,degree=2)))
summary(lm(formula=pca_TN$TN2015~poly(pca_TN$Axis1,degree=2)))
summary(lm(formula=pca_TN$TN2016~poly(pca_TN$Axis1,degree=2)))


