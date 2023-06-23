# TL - 171022 # MAJ: 060523
library(ggplot2)
library(grid)
library(gridExtra)
library(agricolae)
setwd("~/Rose/Papier/Fig1/")

### Evolution taches noires Diplo

data=read.table("dataFigure1_171022_Petals_RecFlowering_Flagr_BlackSp_periodsbtwinf1800and1910.txt",h=T)

data$periode_corCOP_2decades_reordered <- with(data, relevel(periode_corCOP_2decades,"Botanical"))

# compute mean for the 3 years of Black Spot scoring
data$meanBlackSpot <- apply(data[,33:35],1,mean,na.rm=TRUE)

colfunc <- colorRampPalette(c("deeppink3", "white","yellow"))
colfunc(7)
#[1] "#CD1076" "#DD5FA3" "#EEAFD1" "#FFFFFF" "#FFFFAA" "#FFFF55" "#FFFF00"

# before : khaki2, then #e6ab02
plotA <- ggplot(data,aes(x=periode_corCOP_2decades_reordered,y=mean_petals))+geom_violin(scale="width",adjust=0.75,fill="#CD1076",colour="white") + geom_boxplot(colour="navyblue",outlier.shape=NA,width=0.2) + coord_flip()+
  stat_summary(fun.y=mean, geom="point", shape="|", size=15, color="cornflowerblue", fill="cornflowerblue") +
  geom_jitter(shape=20,width=0.04,colour="olivedrab3",size=1.1)+#,colour=c(rep(c("violetred1","forestgreen","darkred"),49),"forestgreen","violetred1"))+
  xlab("")+ylab("Petal counts")+
  theme_bw()+
  theme(panel.background = element_rect(fill="transparent",colour = NA), plot.background = element_rect(fill="transparent",colour=NA),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'white', size = 1.25), axis.ticks = element_line(colour = 'white', size = 1.25), 
        axis.text.x = element_text(colour="white",size=8,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="white",size=0,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="white",size=12,angle=0,hjust=.5,vjust=.2,face="bold"),
        axis.title.y = element_text(colour="white",size=14,angle=90,hjust=.5,vjust=.5,face="bold"))+
  annotate("text",x =1.2,y =165,label="c",color="white",size=6,fontface=2)+
  annotate("text",x =2.2,y =165,label="b",color="white",size=6,fontface=2)+
  annotate("text",x =3.2,y =165,label="ab",color="white",size=6,fontface=2)+
  annotate("text",x =4.2,y =165,label="a",color="white",size=6,fontface=2)+
  annotate("text",x =5.2,y =165,label="ab",color="white",size=6,fontface=2)+
  annotate("text",x =6.2,y =165,label="ab",color="white",size=6,fontface=2)

# before darkorchid1
plotB <- ggplot(data,aes(x=periode_corCOP_2decades_reordered,y=codenumremonteeAvg3years))+geom_violin(scale="width",adjust=0.75,fill="#EEAFD1",colour="white") + geom_boxplot(colour="navyblue",outlier.shape=NA,width=0.2) + coord_flip()+
  stat_summary(fun.y=mean, geom="point", shape="|", size=15, color="cornflowerblue", fill="cornflowerblue") +
  geom_jitter(shape=20,width=0.04,colour="olivedrab3",size=1.1)+#,colour=c(rep(c("violetred1","forestgreen","darkred"),49),"forestgreen","violetred1"))+
  xlab("")+ylab("Recurrent flowering")+
  theme_bw()+
  theme(panel.background = element_rect(fill="transparent",colour = NA), plot.background = element_rect(fill="transparent",colour=NA),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'white', size = 1.25), axis.ticks = element_line(colour = 'white', size = 1.25), 
        axis.text.x = element_text(colour="white",size=8,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="white",size=0,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="white",size=12,angle=0,hjust=.5,vjust=.2,face="bold"),
        axis.title.y = element_text(colour="white",size=14,angle=90,hjust=.5,vjust=.5,face="bold"))+  annotate("text",x =1.15,y =2.7,label="bc",color="white",size=6,fontface=2)+
  annotate("text",x =2.2,y =2.7,label="c",color="white",size=6,fontface=2)+
  annotate("text",x =3.2,y =2.7,label="bc",color="white",size=6,fontface=2)+
  annotate("text",x =4.2,y =2.7,label="abc",color="white",size=6,fontface=2)+
  annotate("text",x =5.2,y =2.7,label="ab",color="white",size=6,fontface=2)+
  annotate("text",x =6.2,y =2.7,label="a",color="white",size=6,fontface=2)

# before lightpink2
plotC <- ggplot(data,aes(x=periode_corCOP_2decades_reordered,y=X2phenylethanol))+geom_violin(scale="width",adjust=0.75,fill="antiquewhite",colour="white") + geom_boxplot(colour="navyblue",outlier.shape=NA,width=0.2) + coord_flip()+
  stat_summary(fun.y=mean, geom="point", shape="|", size=15, color="cornflowerblue", fill="cornflowerblue") +
  geom_jitter(shape=20,width=0.04,colour="olivedrab3",size=1.1)+#,colour=c(rep(c("violetred1","forestgreen","darkred"),49),"forestgreen","violetred1"))+
  xlab("")+ylab("2-PE content")+
  theme_bw()+
  theme(panel.background = element_rect(fill="transparent",colour = NA), plot.background = element_rect(fill="transparent",colour=NA),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'white', size = 1.25), axis.ticks = element_line(colour = 'white', size = 1.25), 
        axis.text.x = element_text(colour="white",size=8,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="white",size=0,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="white",size=12,angle=0,hjust=.5,vjust=.2,face="bold"),
        axis.title.y = element_text(colour="white",size=14,angle=90,hjust=.5,vjust=.5,face="bold"))+
  annotate("text",x =1.15,y =1100,label="ab",color="white",size=6,fontface=2)+
  annotate("text",x =2.15,y =1100,label="a",color="white",size=6,fontface=2)+
  annotate("text",x =3.15,y =1100,label="a",color="white",size=6,fontface=2)+
  annotate("text",x =4.15,y =1100,label="ab",color="white",size=6,fontface=2)+
  annotate("text",x =5.15,y =1100,label="ab",color="white",size=6,fontface=2)+
  annotate("text",x =6.15,y =1100,label="b",color="white",size=6,fontface=2)

# before: olivedrab3
plotD <- ggplot(data,aes(x=periode_corCOP_2decades_reordered,y=geraniol))+geom_violin(scale="width",adjust=0.75,fill="#FFFFAA",colour="white") + geom_boxplot(colour="navyblue",outlier.shape=NA,width=0.2) + coord_flip()+
  stat_summary(fun.y=mean, geom="point", shape="|", size=15, color="cornflowerblue", fill="cornflowerblue") +
  geom_jitter(shape=20,width=0.04,colour="olivedrab3",size=1.1)+#,colour=c(rep(c("violetred1","forestgreen","darkred"),49),"forestgreen","violetred1"))+
  xlab("")+ylab("Geraniol content")+
  theme_bw()+
  theme(panel.background = element_rect(fill="transparent",colour = NA), plot.background = element_rect(fill="transparent",colour=NA),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'white', size = 1.25), axis.ticks = element_line(colour = 'white', size = 1.25), 
        axis.text.x = element_text(colour="white",size=8,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="white",size=0,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="white",size=12,angle=0,hjust=.5,vjust=.2,face="bold"),
        axis.title.y = element_text(colour="white",size=14,angle=90,hjust=.5,vjust=.5,face="bold"))+
  annotate("text",x =1.15,y =265,label="a",color="white",size=6,fontface=2)+
  annotate("text",x =2.15,y =265,label="a",color="white",size=6,fontface=2)+
  annotate("text",x =3.15,y =265,label="a",color="white",size=6,fontface=2)+
  annotate("text",x =4.15,y =265,label="a",color="white",size=6,fontface=2)+
  annotate("text",x =5.15,y =265,label="a",color="white",size=6,fontface=2)+
  annotate("text",x =6.15,y =265,label="a",color="white",size=6,fontface=2)

# before: burlywood4
plotE <- ggplot(data,aes(x=periode_corCOP_2decades_reordered,y=meanBlackSpot))+geom_violin(scale="width",adjust=0.75,fill="gold1",colour="white") + geom_boxplot(colour="navyblue",outlier.shape=NA,width=0.2) + coord_flip()+
  stat_summary(fun.y=mean, geom="point", shape="|", size=15, color="cornflowerblue", fill="cornflowerblue") +
  geom_jitter(shape=20,width=0.04,colour="olivedrab3",size=1.1)+#,colour=c(rep(c("violetred1","forestgreen","darkred"),49),"forestgreen","violetred1"))+
  xlab("")+ylab("Suscept. to blackspot")+
  theme_bw()+
  theme(panel.background = element_rect(fill="transparent",colour = NA), plot.background = element_rect(fill="transparent",colour=NA),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'white', size = 1.25), axis.ticks = element_line(colour = 'white', size = 1.25), 
        axis.text.x = element_text(colour="white",size=8,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="white",size=0,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="white",size=12,angle=0,hjust=.5,vjust=.2,face="bold"),
        axis.title.y = element_text(colour="white",size=14,angle=90,hjust=.5,vjust=.5,face="bold"))+
  annotate("text",x =1.15,y =5.2,label="c",color="white",size=6,fontface=2)+
  annotate("text",x =2.15,y =5.2,label="c",color="white",size=6,fontface=2)+
  annotate("text",x =3.15,y =4.9,label="bc",color="white",size=6,fontface=2)+
  annotate("text",x =4.15,y =4.9,label="ab",color="white",size=6,fontface=2)+
  annotate("text",x =5.15,y =5.2,label="a",color="white",size=6,fontface=2)+
  annotate("text",x =6.15,y =4.62,label="abc",color="white",size=6,fontface=2)

#pdf("Fig2_131019",14,7)
# Move to a new page
grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 11)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plotA, vp = define_region(row = 1, col = 1:3))   # Span over two columns
print(plotB, vp = define_region(row = 1, col = 4:5))
print(plotC, vp = define_region(row = 1, col = 6:7))
print(plotD, vp = define_region(row = 1, col = 8:9))
print(plotE, vp = define_region(row = 1, col = 10:11))

ggsave("test_figure1_bplotA_v180623.png",plotA,bg="transparent",width = 2,height = 8)
ggsave("test_figure1_bplotB_v180623.png",plotB,bg="transparent",width = 2,height = 8)
ggsave("test_figure1_bplotC_v180623.png",plotC,bg="transparent",width = 2,height = 8)
ggsave("test_figure1_bplotD_v180623.png",plotD,bg="transparent",width = 2,height = 8)
ggsave("test_figure1_bplotE_v180623.png",plotE,bg="transparent",width = 2,height = 8)


datapheno=read.table("dataFigure1_171022_Petals_RecFlowering_Flagr_BlackSp_periodsbtwinf1800and1910.txt",h=T)

datapheno$periode_corCOP_2decades_reordered <- with(datapheno, relevel(periode_corCOP_2decades,"Botanical"))

# compute mean for the 3 years of Black Spot scoring
datapheno$meanBlackSpot <- apply(data[,33:35],1,mean,na.rm=TRUE)

modelp1<-aov(data=datapheno,mean_petals~periode_corCOP_2decades_reordered)
outp1 <- HSD.test(modelp1,"periode_corCOP_2decades_reordered", group=TRUE,console=TRUE)


modelp2<-aov(data=datapheno,codenumremonteeAvg3years~periode_corCOP_2decades_reordered)
outp2 <- HSD.test(modelp2,"periode_corCOP_2decades_reordered", group=TRUE,console=TRUE)

modelp3<-aov(data=datapheno,X2phenylethanol~periode_corCOP_2decades_reordered)
outp3 <- HSD.test(modelp3,"periode_corCOP_2decades_reordered", group=TRUE,console=TRUE)

modelp4<-aov(data=datapheno,geraniol~periode_corCOP_2decades_reordered)
outp4 <- HSD.test(modelp4,"periode_corCOP_2decades_reordered", group=TRUE,console=TRUE)

modelp5<-aov(data=datapheno,meanBlackSpot~periode_corCOP_2decades_reordered)
outp5 <- HSD.test(modelp5,"periode_corCOP_2decades_reordered", group=TRUE,console=TRUE)



### petal colors vs. blackspot for rose flowers with a unique color only

petalcolBS=read.table("dataFigure1_204cc_colorPetals_BlackSp_070523_unicoloronly.txt",h=T,sep="\t")
petalcolBS$meanBlackSpot <- apply(petalcolBS[,33:35],1,mean,na.rm=TRUE)
aggregate(data=petalcolBS, meanBlackSpot~main_petal_color, summary)

petalcolBS$periode_corCOP_2decades_reordered <- with(petalcolBS, relevel(periode_corCOP_2decades,"Botanical"))

plotA<- ggplot(petalcolBS,aes(x=main_petal_color, y=meanBlackSpot))+geom_violin(scale="width",adjust=0.75,fill="#d95f02",colour="white") + geom_boxplot(colour="navyblue",outlier.shape=NA,width=0.2) + 
  stat_summary(fun.y=mean, geom="point", shape="-", size=15, color="red", fill="red") +
  geom_jitter(shape=20,width=0.05,aes(colour=petalcolBS$periode_corCOP_2decades_reordered),size=3)+#,colour=c(rep(c("violetred1","forestgreen","darkred"),49),"forestgreen","violetred1"))+
  xlab("Main petal color")+ylab("Susceptibility to blackspot disease")+
  scale_colour_manual(values = c("deeppink","#56B4E9","forestgreen","goldenrod1", "red1","darkred"))+
  theme_bw()+
  theme(legend.text=element_text(color="black",size=18),legend.title=element_text(color="white",size=0),legend.spacing.x = unit(0.15 ,'cm'),legend.spacing.y = unit(0.05 ,'cm'),legend.box.spacing = unit(0.1 ,'cm'))+
  theme(panel.background = element_rect(fill="transparent",colour = NA), plot.background = element_rect(fill="transparent",colour=NA),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="bold"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="bold"))+
  geom_segment(x=0.7,xend=5.3,y=5.1,yend=5.1)+
annotate("text",x =3,y =5.2,label="No significant differences",color="black",size=3,fontface=2)


modelp1<-aov(data=petalcolBS,meanBlackSpot~main_petal_color)
outp1 <- HSD.test(modelp1,"main_petal_color", group=TRUE,console=TRUE)


petalcolBS2=read.table("dataFigure1_204cc_colorPetals_BlackSp_070523_unicoloronly_botanicalandbefore1850.txt",h=T,sep="\t")
petalcolBS2$meanBlackSpot <- apply(petalcolBS2[,33:35],1,mean,na.rm=TRUE)
aggregate(data=petalcolBS2, meanBlackSpot~main_petal_color, summary)

petalcolBS2$periode_corCOP_2decades_reordered <- with(petalcolBS2, relevel(periode_corCOP_2decades,"Botanical"))

plotB<- ggplot(petalcolBS2,aes(x=main_petal_color, y=meanBlackSpot))+geom_violin(scale="width",adjust=0.75,fill="#d95f02",colour="white") + geom_boxplot(colour="navyblue",outlier.shape=NA,width=0.2) + 
  stat_summary(fun.y=mean, geom="point", shape="-", size=15, color="red", fill="red") +
  geom_jitter(shape=20,width=0.05,aes(colour=petalcolBS2$periode_corCOP_2decades_reordered),size=3)+#,colour=c(rep(c("violetred1","forestgreen","darkred"),49),"forestgreen","violetred1"))+
  xlab("Main petal color (botanical + cultivars <1850)")+ylab("Susceptibility to blackspot disease")+
  scale_colour_manual(values = c("deeppink","#56B4E9","forestgreen"))+
  theme_bw()+
  theme(legend.text=element_text(color="black",size=18),legend.title=element_text(color="white",size=0),legend.spacing.x = unit(0.15 ,'cm'),legend.spacing.y = unit(0.05 ,'cm'),legend.box.spacing = unit(0.1 ,'cm'))+
  theme(panel.background = element_rect(fill="transparent",colour = NA), plot.background = element_rect(fill="transparent",colour=NA),panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="bold"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="bold"))+
  geom_segment(x=0.7,xend=3.3,y=5.1,yend=5.1)+
  annotate("text",x =2,y =5.2,label="No significant differences",color="black",size=3,fontface=2)



modelp2<-aov(data=petalcolBS,meanBlackSpot~main_petal_color)
outp2 <- HSD.test(modelp2,"main_petal_color", group=TRUE,console=TRUE)

pdf("Avgblackspotdisease_petalcolors.pdf",14,9)

grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plotA, vp = define_region(row = 1, col = 1))   # Span over two columns
print(plotB, vp = define_region(row = 1, col = 2))

dev.off()


### correlation between number of of metals, remontée et sensibilité
datalatesel=read.table("dataFigure1_204cc_colorPetals_BlackSp_070523_unicoloronly_only1890-1909.txt",h=TRUE,sep="\t")
datalatesel$meanBlackSpot <- apply(datalatesel[,30:32],1,mean,na.rm=TRUE)

plot_ly(x=datalatesel$mean_petals, y=datalatesel$codenumremonteeAvg3years, z=datalatesel$meanBlackSpot, type="scatter3d", mode="markers", color=datalatesel$mean_petals) %>%
  layout(scene = list(xaxis = list(title = "Petal counts"), yaxis = list(title = "Rec. flowering"),zaxis = list(title = "Susceptibility BS")))

summary(lm(formula=datalatesel$codenumremonteeAvg3years~datalatesel$mean_petals))
summary(lm(formula=datalatesel$meanBlackSpot~datalatesel$mean_petals))
summary(lm(formula=datalatesel$meanBlackSpot~datalatesel$codenumremonteeAvg3years))
summary(lm(formula=datalatesel$codenumremonteeAvg3years~datalatesel$PCo1))
