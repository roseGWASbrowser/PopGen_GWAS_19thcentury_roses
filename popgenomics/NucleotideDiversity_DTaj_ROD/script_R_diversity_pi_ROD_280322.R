# TL - 280322
# TL - 220322
library(ggplot2)
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/second_round_calling/diversity/final_280322")

sumstatspi=read.table("rose_round2.jointvcf_nucleotidediversity.AllChr.stats_220623.prR",h=T,sep="\t")

### NUCLEOTIDE DIVERSITY:  pi
plotA <- ggplot(sumstatspi, aes(x=group, y=pi_site)) +
  geom_jitter(aes(color=sumstatspi$chr),shape=20, position=position_jitter(0.15),cex=0.2)+#geom_jitter(shape=20, position=position_jitter(0.3))+
  geom_boxplot(width=0.4,outlier.size=0,alpha=0.5)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab("Rose groups")+ ylab(expression(pi))+
  #ylim(0,0.03)+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


ggplot(sumstatspi, aes(x=pi_site)) +
  geom_density(aes(color=group),lwd=1.4)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab(expression(pi))+
  #ylim(0,0.03)+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  #theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

# chr per chr
plotB <- ggplot(sumstatspi, aes(x=pi_site)) +
  geom_density(aes(color=group),lwd=1.4)+ facet_wrap(~chr)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab(expression(pi))+
  #ylim(0,0.03)+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  #theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

# TAJIMA'S D
plotC <- ggplot(sumstatspi, aes(x=group, y=D_Pop)) +
  geom_jitter(aes(color=sumstatspi$chr),shape=20, position=position_jitter(0.15),cex=0.2)+#geom_jitter(shape=20, position=position_jitter(0.3))+
  geom_boxplot(width=0.4,outlier.size=0,alpha=0.5)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab("Rose groups")+ ylab("Tajima's D")+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


ggplot(sumstatspi, aes(x=D_Pop)) +
  geom_density(aes(color=group),lwd=1.4)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab("Tajima's D")+
  #ylim(0,0.03)+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  #theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

# chr per chr
plotD <- ggplot(sumstatspi, aes(x=D_Pop)) +
  geom_density(aes(color=group),lwd=1.4)+ facet_wrap(~chr)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab("Tajima's D")+
  #ylim(0,0.03)+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  #theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

### Reduction of diversity, computed as 1 - (pi new / pi ancien), the value ranges between 1 (100% of the diversity loss, and - infiinite). Positive values indicate a net loss of diversity, negative a gain!
dataROD=read.table("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.ROD.prR",h=TRUE)
ggplot(dataROD, aes(x=comparison, y=RoD_pi)) + geom_hline(yintercept = 0,col="red",lwd=1.1)+
  geom_jitter(aes(color=dataROD$chr),shape=20, position=position_jitter(0.15),cex=0.2)+#geom_jitter(shape=20, position=position_jitter(0.3))+
  geom_boxplot(width=0.4,outlier.size=0,alpha=0.5)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab("Rose groups")+ ylab("RoD")+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

ggplot(dataROD, aes(x=comparison, y=RoD_pi)) +  geom_hline(yintercept = 0,col="red",lwd=1.1)+
  geom_jitter(aes(color=dataROD$chr),shape=20, position=position_jitter(0.15),cex=0.2)+#geom_jitter(shape=20, position=position_jitter(0.3))+
  geom_boxplot(width=0.4,outlier.size=0,alpha=0.5)+
  ylim(-1,1)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab("Rose group comparisons")+ ylab("RoD")+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


ggplot(dataROD, aes(x=comparison, y=RoD_pi)) +   facet_wrap(~chr)+ geom_hline(yintercept = 0,col="red",lwd=1.1)+
  geom_jitter(aes(color=sumstatspi$chr),shape=20, position=position_jitter(0.15),cex=0.2)+#geom_jitter(shape=20, position=position_jitter(0.3))+
  geom_boxplot(width=0.4,outlier.size=0,alpha=0.5)+
  coord_flip()+
  ylim(-1,1)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab("Rose groups")+ ylab("RoD")+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


### reduction of diversity in Europe between European and Tea Hybrids
dataRODhybTeavsEurope=subset(dataROD,dataROD$comparison=="HybTea_vs_Europe")
mean(dataRODhybTeavsEurope$RoD_pi,na.rm=TRUE)
# 0.2689011
# If all the european old cultivars are lost, the net loss of genetic diversity in Europe would be around 27% -> conservation

# step by step:
dataRODHyb1stvsEurope=subset(dataROD,dataROD$comparison=="HybFirst_vs_Europe")
mean(dataRODHyb1stvsEurope$RoD_pi,na.rm=TRUE)
# 0.08786848
dataRODHybTeavsHyb1st=subset(dataROD,dataROD$comparison=="HybTea_vs_HybFirst")
mean(dataRODHybTeavsHyb1st$RoD_pi,na.rm=TRUE)
#  0.1869149
# This means that the diversity reduced more in the tea hybrids than in the 1st hybrids of our dataset (8.8% between European and 1st hybrids and then the remaining diversity was reduced by 18.7% between the 1st hybrids and the tea hybrids)


dataRODHybTeavsAsia=subset(dataROD,dataROD$comparison=="HybTea_vs_Asia")
mean(dataRODHybTeavsAsia$RoD_pi,na.rm=TRUE)
# -0.07115337 
#> median(dataRODHybTeavsAsia$RoD_pi,na.rm=TRUE)
#[1] -0.01590499
# Slight improvment of the diversity if compared to Asian, but far less than the loss for the European compartment

dataRODHyb1stvsAsia=subset(dataROD,dataROD$comparison=="HybFirst_vs_Asia")
mean(dataRODHyb1stvsAsia$RoD_pi,na.rm=TRUE)

### Absolute tajima's D deviation between groups

ggplot(dataROD, aes(x=comparison, y=deltaD)) + 
  geom_jitter(aes(color=sumstatspi$chr),shape=20, position=position_jitter(0.15),cex=0.2)+#geom_jitter(shape=20, position=position_jitter(0.3))+
  geom_boxplot(width=0.4,outlier.size=0,alpha=0.5)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab("Rose groups")+ ylab("Absolute Tajima's D difference between groups")+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


ggplot(dataROD, aes(x=comparison, y=deltaD)) + facet_wrap(~chr)+  
  geom_jitter(aes(color=sumstatspi$chr),shape=20, position=position_jitter(0.15),cex=0.2)+#geom_jitter(shape=20, position=position_jitter(0.3))+
  geom_boxplot(width=0.4,outlier.size=0,alpha=0.5)+
  coord_flip()+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab("Rose groups")+ ylab("Absolute Tajima's D difference between groups")+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


