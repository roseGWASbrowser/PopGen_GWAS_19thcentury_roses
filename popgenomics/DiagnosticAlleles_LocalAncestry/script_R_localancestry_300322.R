# TL - 300322$
library(ggplot2)
library(grid)
library(gridExtra)
# generate estimates of the local ancestry using diagnostic markers
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/second_round_calling/localancestry")
diagnostic=read.table("rose_round2.jointvcf.clean.AllFreqRefPerGroup.diagnostics.txt",h=TRUE)

# plot boxplot
dataset=cbind("Ancient Asian",diagnostic$chr,diagnostic$pos,diagnostic$freqasia)
dataset=rbind(dataset,cbind("Ancient European",diagnostic$chr,diagnostic$pos,diagnostic$freqeuro))
dataset=rbind(dataset,cbind("Early European x Asian",diagnostic$chr,diagnostic$pos,diagnostic$freqfirsthyb))
dataset=rbind(dataset,cbind("Hybrid Tea roses",diagnostic$chr,diagnostic$pos,diagnostic$freqhybtea))
dataset=as.data.frame(dataset)
colnames(dataset)<- c("Group","chr","pos","freq")
dataset$freq=as.numeric(as.character(dataset$freq))
dataset$Group=as.factor(as.character(dataset$Group))


cols <- c("Ancient Asian"="goldenrod2","Ancient European"="navyblue",
          "Early European x Asian" ="grey Tea roses" ="grey")

ggplot(dataset)+
  geom_hline(yintercept = 0.5,lty=2,color="red")+
  geom_hline(yintercept = 0.25,lty=2,color="orange")+
  geom_hline(yintercept = 0.75,lty=2,color="orange")+
  geom_violin(aes(x=Group,y=freq),draw_quantiles =0.5)+
  ylim(0,1)+
  xlab("Groups")+  ylab("Allele frequency at diagnostic markers")+
  scale_colour_manual(values = cols)+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

chr3=subset(diagnostic,diagnostic$chr=="Chr03")
ggplot(chr3,aes(x=pos,y=freqhybtea))+
  geom_point(size=0.1)+
  geom_smooth(adjust=0.2)+
  xlab("Position on Chr3")+ylab("Allele frequency (Hyb Tea)")+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

chr3fixed1sthybasia=subset(diagnostic,diagnostic$chr=="Chr03" & diagnostic$freqfirsthyb==1)
chr3fixed1sthybeuro=subset(diagnostic,diagnostic$chr=="Chr03" & diagnostic$freqfirsthyb==0)
chr3fixedhybteaasia=subset(diagnostic,diagnostic$chr=="Chr03" & diagnostic$freqhybtea==1)
chr3fixedhybteaeuro=subset(diagnostic,diagnostic$chr=="Chr03" & diagnostic$freqhybtea==0)

ggplot()+
  geom_density(aes(chr3fixed1sthybasia$pos),adjust=0.2,colour="goldenrod2")+
  geom_density(aes(chr3fixed1sthybeuro$pos),adjust=0.2,colour="navyblue")+  
  xlab("Position on Chr3")+ylab("Density of fixed diagnostic allele")+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

ggplot()+
  geom_density(aes(chr3fixedhybteaasia$pos),adjust=0.2,colour="goldenrod2")+
  #geom_density(aes(chr3fixedhybteaeuro$pos),adjust=0.2,colour="navyblue")+  
  xlab("Position on Chr3")+ylab("Density of fixed diagnostic allele")+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))



windowsize=100000
sumstats=cbind("chr","startwin","endwin","NbDiagnoSNPs","medianfreqfirsthyb","meanfreqfirsthyb","medianfreqhybtea","meanfreqhybtea")
# for each chromosome
for (chr in c(0:7)){
  chromo=paste("Chr0",chr,sep="")
  subsetdiagnochr=subset(diagnostic,diagnostic$chr==chromo)
  print(paste(chromo,nrow(subsetdiagnochr)))
  # for each 100 kb window
  for (win in c(1:ceiling(max(subsetdiagnochr$pos)/windowsize))){
    startwin=win*windowsize - (windowsize-1)
    endwin=win*windowsize
    subsetregion=subset(subsetdiagnochr,(subsetdiagnochr$pos >= startwin & subsetdiagnochr$pos <= endwin))
    if (nrow(subsetregion) >= 5) {
      medianfreqfirsthybsubsetregion=median(subsetregion$freqfirsthyb,na.rm=TRUE)
      meanfreqfirsthybsubsetregion=mean(subsetregion$freqfirsthyb,na.rm=TRUE)
      medianfreqhybteasubsetregion=median(subsetregion$freqhybtea,na.rm=TRUE)
      meanfreqhybteasubsetregion=mean(subsetregion$freqhybtea,na.rm=TRUE)
    }else {
      medianfreqfirsthybsubsetregion=NA
      meanfreqfirsthybsubsetregion=NA
      medianfreqhybteasubsetregion=NA
      meanfreqhybteasubsetregion=NA
    }
    sumstats=rbind(sumstats,cbind(chromo,startwin,endwin,nrow(subsetregion),medianfreqfirsthybsubsetregion,meanfreqfirsthybsubsetregion,medianfreqhybteasubsetregion,meanfreqhybteasubsetregion))
  }
}
write.table(sumstats,file="localancestryestimates_diagnostic_300322.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


### generate plots

# median
pdf(file="localancestry_median_slidwin100kb_300322.pdf",height=9,width=10)
quantilesmedianfirsthyb100kb=quantile(read.table("localancestryestimates_diagnostic_300322.txt",header =TRUE)$medianfreqfirsthyb,probs=c(0.01,0.05,0.5,0.95,0.99),na.rm = TRUE)
quantilesmedianhybtea100kb=quantile(read.table("localancestryestimates_diagnostic_300322.txt",header =TRUE)$medianfreqhybtea,probs=c(0.01,0.05,0.5,0.95,0.99),na.rm = TRUE)
for (chr in c(0:7)){
  dataancestry=read.table(paste("localancestryestimates_diagnostic_300322.txt.Chr0",chr,sep=""),h=TRUE)
  quantilesmedianfirsthyb100kbchr=quantile(dataancestry$meanfreqfirsthyb,probs=c(0.01,0.05,0.5,0.95,0.99),na.rm=TRUE)
  quantilesmedianhybtea100kbchr=quantile(dataancestry$meanfreqhybtea,probs=c(0.01,0.05,0.5,0.95,0.99),na.rm=TRUE)
  plotA<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=medianfreqfirsthyb),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=medianfreqfirsthyb),color="navyblue")+
    geom_hline(yintercept = quantilesmedianfirsthyb100kb[1],color="black",lty=2)+
    geom_hline(yintercept = quantilesmedianfirsthyb100kb[2],color="grey40",lty=2)+
    geom_hline(yintercept = quantilesmedianfirsthyb100kbchr[3],color="grey",lty=1)+ # local median (chr)
    geom_hline(yintercept = quantilesmedianfirsthyb100kb[3],color="grey",lty=2)+
    geom_hline(yintercept = quantilesmedianfirsthyb100kb[4],color="grey40",lty=2)+
    geom_hline(yintercept = quantilesmedianfirsthyb100kb[5],color="black",lty=2)+
    geom_segment(x=-2000000,y=quantilesmedianfirsthyb100kb[3],xend=-2000000,yend=quantilesmedianfirsthyb100kbchr[3],arrow=arrow(length=unit(0.10,"cm")),size=1.1,color=ifelse(quantilesmedianfirsthyb100kbchr[3]>=quantilesmedianfirsthyb100kb[3],"goldenrod2","navyblue"))+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  plotB<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=medianfreqhybtea),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=medianfreqhybtea),color="navyblue")+
    geom_hline(yintercept = quantilesmedianhybtea100kb[1],color="black",lty=2)+
    geom_hline(yintercept = quantilesmedianhybtea100kb[2],color="grey40",lty=2)+
    geom_hline(yintercept = quantilesmedianhybtea100kbchr[3],color="grey",lty=1)+ # local median (chr)
    geom_hline(yintercept = quantilesmedianhybtea100kb[3],color="grey",lty=2)+
    geom_hline(yintercept = quantilesmedianhybtea100kb[4],color="grey40",lty=2)+
    geom_hline(yintercept = quantilesmedianhybtea100kb[5],color="black",lty=2)+
    geom_segment(x=-2000000,y=quantilesmedianhybtea100kb[3],xend=-2000000,yend=quantilesmedianhybtea100kbchr[3],arrow=arrow(length=unit(0.10,"cm")),size=1.1,color=ifelse(quantilesmedianhybtea100kbchr[3]>=quantilesmedianhybtea100kb[3],"goldenrod2","navyblue"))+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

  grid.arrange(plotA,plotB,nrow=2)
}
dev.off()


# mean

pdf(file="localancestry_mean_slidwin100kb_300322.pdf",height=9,width=10)
quantilesmeanfirsthyb100kb=quantile(read.table("localancestryestimates_diagnostic_300322.txt",header =TRUE)$meanfreqfirsthyb,probs=c(0.01,0.05,0.5,0.95,0.99),na.rm = TRUE)
quantilesmeanhybtea100kb=quantile(read.table("localancestryestimates_diagnostic_300322.txt",header =TRUE)$meanfreqhybtea,probs=c(0.01,0.05,0.5,0.95,0.99),na.rm = TRUE)
for (chr in c(0:7)){
  dataancestry=read.table(paste("localancestryestimates_diagnostic_300322.txt.Chr0",chr,sep=""),h=TRUE)
  quantilesmeanfirsthyb100kbchr=quantile(dataancestry$meanfreqfirsthyb,probs=c(0.01,0.05,0.5,0.95,0.99),na.rm=TRUE)
  quantilesmeanhybtea100kbchr=quantile(dataancestry$meanfreqhybtea,probs=c(0.01,0.05,0.5,0.95,0.99),na.rm=TRUE)
  plotA<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=meanfreqfirsthyb),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=meanfreqfirsthyb),color="navyblue")+
    geom_hline(yintercept = quantilesmeanfirsthyb100kb[1],color="black",lty=2)+ # bottom 1% values (whole genome)
    geom_hline(yintercept = quantilesmeanfirsthyb100kb[2],color="grey40",lty=2)+ # bottom 5%
    geom_hline(yintercept = quantilesmeanfirsthyb100kbchr[3],color="grey",lty=1)+ # local median (chr)
    geom_hline(yintercept = quantilesmeanfirsthyb100kb[3],color="grey",lty=2)+ # median whole genome
    geom_hline(yintercept = quantilesmeanfirsthyb100kb[4],color="grey40",lty=2)+ # top 5%
    geom_hline(yintercept = quantilesmeanfirsthyb100kb[5],color="black",lty=2)+ # top 1% values (whole genome)
    geom_segment(x=-2000000,y=quantilesmeanfirsthyb100kb[3],xend=-2000000,yend=quantilesmeanfirsthyb100kbchr[3],arrow=arrow(length=unit(0.10,"cm")),size=1.1,color=ifelse(quantilesmeanfirsthyb100kbchr[3]>=quantilesmeanfirsthyb100kb[3],"goldenrod2","navyblue"))+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  plotB<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=meanfreqhybtea),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=meanfreqhybtea),color="navyblue")+
    geom_hline(yintercept = quantilesmeanhybtea100kb[1],color="black",lty=2)+
    geom_hline(yintercept = quantilesmeanhybtea100kb[2],color="grey40",lty=2)+
    geom_hline(yintercept = quantilesmeanhybtea100kbchr[3],color="grey",lty=1)+ # local median (chr)
    geom_hline(yintercept = quantilesmeanhybtea100kb[3],color="grey",lty=2)+
    geom_hline(yintercept = quantilesmeanhybtea100kb[4],color="grey40",lty=2)+
    geom_hline(yintercept = quantilesmeanhybtea100kb[5],color="black",lty=2)+
    geom_segment(x=-2000000,y=quantilesmeanhybtea100kb[3],xend=-2000000,yend=quantilesmeanhybtea100kbchr[3],arrow=arrow(length=unit(0.10,"cm")),size=1.1,color=ifelse(quantilesmeanhybtea100kbchr[3]>=quantilesmeanhybtea100kb[3],"goldenrod2","navyblue"))+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  grid.arrange(plotA,plotB,nrow=2)
}
dev.off()






### computation with slidingwindows


windowsize=300000
shift=50000
sumstats=cbind("chr","startwin","medianwin","endwin","NbDiagnoSNPs","medianfreqfirsthyb","meanfreqfirsthyb","medianfreqhybtea","meanfreqhybtea")
# for each chromosome
for (chr in c(0:7)){
  chromo=paste("Chr0",chr,sep="")
  subsetdiagnochr=subset(diagnostic,diagnostic$chr==chromo)
  print(paste(chromo,nrow(subsetdiagnochr)))
  # for each 100 kb window
  for (win in c(1:ceiling(max(subsetdiagnochr$pos)/shift))){
    startwin=win*shift - ((windowsize/2)-1)
    midwin=win*shift
    endwin=win*shift + (windowsize/2)
    subsetregion=subset(subsetdiagnochr,(subsetdiagnochr$pos >= startwin & subsetdiagnochr$pos <= endwin))
    if (nrow(subsetregion) >= 30) {
      medianfreqfirsthybsubsetregion=median(subsetregion$freqfirsthyb,na.rm=TRUE)
      meanfreqfirsthybsubsetregion=mean(subsetregion$freqfirsthyb,na.rm=TRUE)
      medianfreqhybteasubsetregion=median(subsetregion$freqhybtea,na.rm=TRUE)
      meanfreqhybteasubsetregion=mean(subsetregion$freqhybtea,na.rm=TRUE)
    }else {
      medianfreqfirsthybsubsetregion=NA
      meanfreqfirsthybsubsetregion=NA
      medianfreqhybteasubsetregion=NA
      meanfreqhybteasubsetregion=NA
    }
    sumstats=rbind(sumstats,cbind(chromo,startwin,midwin,endwin,nrow(subsetregion),medianfreqfirsthybsubsetregion,meanfreqfirsthybsubsetregion,medianfreqhybteasubsetregion,meanfreqhybteasubsetregion))
  }
}
write.table(sumstats,file="localancestryestimates_diagnostic_win300kb_step50k_310322.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

# median
pdf(file="localancestry_median_slidwin300kb_step50kb_310322.pdf",height=9,width=10)
for (chr in c(0:7)){
  dataancestry=read.table(paste("localancestryestimates_diagnostic_win300kb_step50k_310322.txt.Chr0",chr,sep=""),h=TRUE)
  plotA<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=medianfreqfirsthyb),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=medianfreqfirsthyb),color="navyblue")+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  plotB<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=medianfreqhybtea),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=medianfreqhybtea),color="navyblue")+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  grid.arrange(plotA,plotB,nrow=2)
}
dev.off()


# mean

pdf(file="localancestry_mean_slidwin300kb_step50kb_310322.pdf",height=9,width=10)
for (chr in c(0:7)){
  dataancestry=read.table(paste("localancestryestimates_diagnostic_win300kb_step50k_310322.txt.Chr0",chr,sep=""),h=TRUE)
  plotA<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=meanfreqfirsthyb),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=meanfreqfirsthyb),color="navyblue")+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  plotB<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=meanfreqhybtea),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=meanfreqhybtea),color="navyblue")+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  grid.arrange(plotA,plotB,nrow=2)
}
dev.off()


## 1 Mb


### computation with slidingwindows


windowsize=1000000
shift=100000
sumstats=cbind("chr","startwin","medianwin","endwin","NbDiagnoSNPs","medianfreqfirsthyb","meanfreqfirsthyb","medianfreqhybtea","meanfreqhybtea")
# for each chromosome
for (chr in c(0:7)){
  chromo=paste("Chr0",chr,sep="")
  subsetdiagnochr=subset(diagnostic,diagnostic$chr==chromo)
  print(paste(chromo,nrow(subsetdiagnochr)))
  # for each 100 kb window
  for (win in c(1:ceiling(max(subsetdiagnochr$pos)/shift))){
    startwin=win*shift - ((windowsize/2)-1)
    midwin=win*shift
    endwin=win*shift + (windowsize/2)
    subsetregion=subset(subsetdiagnochr,(subsetdiagnochr$pos >= startwin & subsetdiagnochr$pos <= endwin))
    if (nrow(subsetregion) >= 30) {
      medianfreqfirsthybsubsetregion=median(subsetregion$freqfirsthyb,na.rm=TRUE)
      meanfreqfirsthybsubsetregion=mean(subsetregion$freqfirsthyb,na.rm=TRUE)
      medianfreqhybteasubsetregion=median(subsetregion$freqhybtea,na.rm=TRUE)
      meanfreqhybteasubsetregion=mean(subsetregion$freqhybtea,na.rm=TRUE)
    }else {
      medianfreqfirsthybsubsetregion=NA
      meanfreqfirsthybsubsetregion=NA
      medianfreqhybteasubsetregion=NA
      meanfreqhybteasubsetregion=NA
    }
    sumstats=rbind(sumstats,cbind(chromo,startwin,midwin,endwin,nrow(subsetregion),medianfreqfirsthybsubsetregion,meanfreqfirsthybsubsetregion,medianfreqhybteasubsetregion,meanfreqhybteasubsetregion))
  }
}
write.table(sumstats,file="localancestryestimates_diagnostic_win1Mb_step100k_310322.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

# median
pdf(file="localancestry_median_slidwin1Mb_step100kb_310322.pdf",height=9,width=10)
for (chr in c(0:7)){
  dataancestry=read.table(paste("localancestryestimates_diagnostic_win300kb_step50k_310322.txt.Chr0",chr,sep=""),h=TRUE)
  plotA<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=medianfreqfirsthyb),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=medianfreqfirsthyb),color="navyblue")+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  plotB<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=medianfreqhybtea),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=medianfreqhybtea),color="navyblue")+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  grid.arrange(plotA,plotB,nrow=2)
}
dev.off()


# mean

pdf(file="localancestry_mean_slidwin1Mb_step100kb_310322.pdf",height=9,width=10)
for (chr in c(0:7)){
  dataancestry=read.table(paste("localancestryestimates_diagnostic_win300kb_step50k_310322.txt.Chr0",chr,sep=""),h=TRUE)
  plotA<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=meanfreqfirsthyb),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=meanfreqfirsthyb),color="navyblue")+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  plotB<- ggplot(dataancestry)+
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=meanfreqhybtea),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=meanfreqhybtea),color="navyblue")+
    #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    ylim(0,1)+
    xlab(paste("Chromosome",chr))+  ylab("Local ancestry")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  grid.arrange(plotA,plotB,nrow=2)
}
dev.off()

