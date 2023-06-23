# TL - 300322
library(grid)
library(gridExtra)
library(ggplot2)
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/second_round_calling/diversity/final_280322")

ROD=read.table("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.ROD",h=TRUE)
quantileRODHyb1stEuro=quantile(ROD$RoD_pi_Hyb1stvsEurope,probs=c(0.05,0.5,0.95),na.rm=TRUE)
quantileRODHybTeaEuro=quantile(ROD$RoD_pi_HybTeavsEurope,probs=c(0.05,0.5,0.95),na.rm=TRUE)
quantileRODHyb1stAsia=quantile(ROD$RoD_pi_Hyb1stvsAsia,probs=c(0.05,0.5,0.95),na.rm=TRUE)
quantileRODHybTeaAsia=quantile(ROD$RoD_pi_HybTeavsAsia,probs=c(0.05,0.5,0.95),na.rm=TRUE)

pdf(file="diversity_ROD_KperK_050422_withcandidategenes.pdf",height=20,width=10)
for (chr in c(0:7)){
  datapi=read.table(paste("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.prR.Chr0",chr,sep=""),h=TRUE)
  plotA<- ggplot(datapi,aes(x=startwin,y=pi_site,group=group))+
    geom_point(aes(color=datapi$group),size=1)+
    geom_smooth(aes(color=datapi$group),span=0.1,fill="gray84",alpha=0.5)+
    scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    #ylim(0,max(datapi$pi_persite))+
  xlab(paste("Chromosome",chr))+  ylab(expression(pi))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",28500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",30750000,0),y=ifelse(chr=="3",max(datapi$pi_site),0),size=ifelse(chr=="3",2.5,0),label="roKSN\n(28.5-33Mb)",color="forestgreen")+
  # Number of petals
  geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="1",41514796,0),y=ifelse(chr=="1",max(datapi$pi_site),0),size=ifelse(chr=="1",2.5,0),label="AGL6",color="forestgreen")+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="1",52096899,0),y=ifelse(chr=="1",max(datapi$pi_site),0),size=ifelse(chr=="1",2.5,0),label="TOE2\neuAP1",color="forestgreen")+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",11833234,0),y=ifelse(chr=="2",max(datapi$pi_site),0),size=ifelse(chr=="2",2.5,0),label="TM6",color="forestgreen")+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",17461203,0),y=ifelse(chr=="2",max(datapi$pi_site),0),size=ifelse(chr=="2",2.5,0),label="AP2",color="forestgreen")+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",57253197,0),y=ifelse(chr=="2",max(datapi$pi_site),0),size=ifelse(chr=="2",2.5,0),label="PLE\nAgD",color="forestgreen")+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",67016526,0),y=ifelse(chr=="2",max(datapi$pi_site),0),size=ifelse(chr=="2",2.5,0),label="euFUL\nSEP1/2",color="forestgreen")+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",33238136,0),y=ifelse(chr=="3",max(datapi$pi_site),0),size=ifelse(chr=="3",2.5,0),label="AP2\nTOE",color="forestgreen")+    
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="4",46474063,0),y=ifelse(chr=="4",max(datapi$pi_site),0),size=ifelse(chr=="4",2.5,0),label="SEP3",color="forestgreen")+    
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="5",8425314,0),y=ifelse(chr=="5",max(datapi$pi_site),0),size=ifelse(chr=="5",2.5,0),label="AgC",color="forestgreen")+    
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="5",61446161,0),y=ifelse(chr=="5",max(datapi$pi_site),0),size=ifelse(chr=="5",2.5,0),label="TOE1",color="forestgreen")+    
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="6",58611899,0),y=ifelse(chr=="6",max(datapi$pi_site),0),size=ifelse(chr=="6",2.5,0),label="PI/\nGLO",color="forestgreen")+    
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="6",60224075,0),y=ifelse(chr=="6",max(datapi$pi_site),0),size=ifelse(chr=="6",2.5,0),label="AP3",color="forestgreen")+    
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="7",4414048,0),y=ifelse(chr=="7",max(datapi$pi_site),0),size=ifelse(chr=="7",2.5,0),label="SEP4\neuFUL",color="forestgreen")+    
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="7",50555186,0),y=ifelse(chr=="7",max(datapi$pi_site),0),size=ifelse(chr=="7",2.5,0),label="D",color="forestgreen")+    
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="7",52218507,0),y=ifelse(chr=="7",max(datapi$pi_site),0),size=ifelse(chr=="7",2.5,0),label="AGL12",color="forestgreen")+    
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",41500000,0),y=ifelse(chr=="3",max(datapi$pi_site),0),size=ifelse(chr=="3",2.5,0),label="S\nlocus",color="forestgreen") + 
    # scent
    #2PE
   geom_vline(xintercept = ifelse(chr=="6",5220887,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="6",5220887,0),y=ifelse(chr=="6",max(datapi$pi_site),0),size=ifelse(chr=="6",2.5,0),label="RhPAAS",color="forestgreen")+  
   geom_vline(xintercept = ifelse(chr=="0",23612237,-500000),lwd=ifelse(chr=="0",1,0),color=ifelse(chr=="0","forestgreen","white"))+
     annotate("text",x=ifelse(chr=="0",23612237,0),y=ifelse(chr=="0",max(datapi$pi_site),0),size=ifelse(chr=="0",2.5,0),label="RhPAAS\n-like",color="forestgreen")+
   geom_vline(xintercept = ifelse(chr=="5",10508387,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
     annotate("text",x=ifelse(chr=="5",10508387,0),y=ifelse(chr=="5",max(datapi$pi_site),0),size=ifelse(chr=="5",2.5,0),label="RhPAAS\n-like",color="forestgreen")+
    #NUDIX
    geom_vline(xintercept = ifelse(chr=="4",51340000,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="4",51340000,0),y=ifelse(chr=="4",max(datapi$pi_site),0),size=ifelse(chr=="4",2.5,0),label="NUDIX1\n(anc)",color="forestgreen")+  
    geom_vline(xintercept = ifelse(chr=="0",3210000,-500000),lwd=ifelse(chr=="0",1,0),color=ifelse(chr=="0","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="0",3210000,0),y=ifelse(chr=="0",max(datapi$pi_site),0),size=ifelse(chr=="0",2.5,0),label="NUDIX1-1a",color="forestgreen")+    
    geom_vline(xintercept = ifelse(chr=="6",2816735,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="6",2816735,0),y=ifelse(chr=="6",max(datapi$pi_site),0),size=ifelse(chr=="6",2.5,0),label="NUDIX1-2b",color="forestgreen")+    
    geom_vline(xintercept = ifelse(chr=="7",37345000,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="7",37345000,0),y=ifelse(chr=="7",max(datapi$pi_site),0),size=ifelse(chr=="7",2.5,0),label="NUDIX1-2c",color="forestgreen")  
    
    
  plotB<- ggplot(datapi,aes(x=startwin,y=D_Pop,group=group))+
    geom_point(aes(color=datapi$group),size=1)+
    geom_smooth(aes(color=datapi$group),span=0.1,fill="gray84",alpha=0.5)+
    scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
    #ylim(0,max(datapi$D_Pop))+
    xlab(paste("Chromosome",chr))+  ylab("Tajima's D")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",28500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",30750000,0),y=ifelse(chr=="3",max(datapi$D_Pop),0),size=ifelse(chr=="3",2.5,0),label="roKSN\n(28.5-33Mb)",color="forestgreen")
  
  

  datadiv=read.table(paste("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.ROD.Chr",chr,sep=""),h=TRUE)
  
  plotC<- ggplot(datadiv, aes(x=startwin, y=RoD_pi_Hyb1stvsEurope)) +  
    geom_segment(aes(x=startwin,y=0,xend=startwin,yend=RoD_pi_Hyb1stvsEurope,col=ifelse(RoD_pi_Hyb1stvsEurope<quantileRODHyb1stEuro[1],"a",ifelse(RoD_pi_Hyb1stvsEurope<0,"b",ifelse(RoD_pi_Hyb1stvsEurope>quantileRODHyb1stEuro[3],"d","c")))))+
    ylim(-1,1)+
    scale_color_manual(values=c("forestgreen", "chartreuse3","orange","red3"))+
    xlab(paste("Chromosome",chr))+ ylab("RoD First Hybrid vs. Europe")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",28500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",30750000,0),y=ifelse(chr=="3",1,0),size=ifelse(chr=="3",2.5,0),label="roKSN\n(28.5-33Mb)",color="forestgreen")
  
  
  plotD<- ggplot(datadiv, aes(x=startwin, y=RoD_pi_HybTeavsEurope)) +  
    geom_segment(aes(x=startwin,y=0,xend=startwin,yend=RoD_pi_HybTeavsEurope,col=ifelse(RoD_pi_HybTeavsEurope<0,"a",ifelse(RoD_pi_HybTeavsEurope>quantileRODHybTeaEuro[3],"c","b"))))+
    ylim(-1,1)+
    scale_color_manual(values=c("forestgreen","orange","red3"))+ # specific case where very few negative RoD values = 5% quantile is already positive
    xlab(paste("Chromosome",chr))+ ylab("RoD Hybrid Tea vs. Europe")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",28500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",30750000,0),y=ifelse(chr=="3",1,0),size=ifelse(chr=="3",2.5,0),label="roKSN\n(28.5-33Mb)",color="forestgreen")
  
  
  plotE<- ggplot(datadiv, aes(x=startwin, y=RoD_pi_Hyb1stvsAsia)) +  
    geom_segment(aes(x=startwin,y=0,xend=startwin,yend=RoD_pi_Hyb1stvsAsia,col=ifelse(RoD_pi_Hyb1stvsAsia<quantileRODHyb1stAsia[1],"a",ifelse(RoD_pi_Hyb1stvsAsia<0,"b",ifelse(RoD_pi_Hyb1stvsAsia>quantileRODHyb1stAsia[3],"d","c")))))+
    ylim(-1,1)+
    ylim(-1,1)+
    scale_color_manual(values=c("forestgreen", "chartreuse3","orange","red3"))+
    xlab(paste("Chromosome",chr))+ ylab("RoD First Hybrid vs. Asia")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",28500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",30750000,0),y=ifelse(chr=="3",1,0),size=ifelse(chr=="3",2.5,0),label="roKSN\n(28.5-33Mb)",color="forestgreen")
  
  plotF<- ggplot(datadiv, aes(x=startwin, y=RoD_pi_HybTeavsAsia)) +  
    geom_segment(aes(x=startwin,y=0,xend=startwin,yend=RoD_pi_HybTeavsAsia,col=ifelse(RoD_pi_HybTeavsAsia<quantileRODHybTeaAsia[1],"a",ifelse(RoD_pi_HybTeavsAsia<0,"b",ifelse(RoD_pi_HybTeavsAsia>quantileRODHybTeaAsia[3],"d","c")))))+
    ylim(-1,1)+
    scale_color_manual(values=c("forestgreen", "chartreuse3","orange","red3"))+
    xlab(paste("Chromosome",chr))+ ylab("RoD Hybrid Tea vs. Asia")+
    theme_bw()+
    theme(legend.position="none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",28500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",30750000,0),y=ifelse(chr=="3",1,0),size=ifelse(chr=="3",2.5,0),label="roKSN\n(28.5-33Mb)",color="forestgreen")
  
  
  grid.arrange(plotA,plotB,plotC,plotD,plotE,plotF,nrow=6)
 }
dev.off()




#### chr 7 
pdf(file="diversity_ROD_KperK_080422_withcandidategenes_chr7.pdf",height=14,width=10)

datapi=read.table("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.prR.Chr07",h=TRUE)

plotA<- ggplot(datapi,aes(x=startwin,y=pi_site,group=group))+
geom_point(aes(color=datapi$group),size=1)+
geom_smooth(aes(color=datapi$group),span=0.1,fill="gray84",alpha=0.5)+
scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
#ylim(0,max(datapi$pi_persite))+
xlab("Chromosome 7")+  ylab(expression(pi))+
theme_bw()+
theme(legend.position="none")+
theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
      axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
      axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
      axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
      axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))#+
#geom_vline(xintercept =4414048,y=max(datapi$pi_site),lwd=1,color="forestgreen")+
#annotate("text",x=4414048,y=max(datapi$pi_site),size=2.5,label="SEP4\neuFUL",color="forestgreen")+    
#geom_vline(xintercept =50555186,y=max(datapi$pi_site),lwd=1,color="forestgreen")+
#annotate("text",x=50555186,y=max(datapi$pi_site),size=2.5,label="D",color="forestgreen")+    
#geom_vline(xintercept =52218507,y=max(datapi$pi_site),lwd=1,color="forestgreen")+
#annotate("text",x=52218507,y=max(datapi$pi_site),size=2.5,label="AGL12",color="forestgreen")+    
# scent
#geom_vline(xintercept = 37345000,lwd=1,color="forestgreen")+
#annotate("text",x=37345000,y=max(datapi$pi_site),size=2.5,label="NUDIX1-2c",color="forestgreen")  


plotB<- ggplot(datapi,aes(x=startwin,y=D_Pop,group=group))+
  geom_point(aes(color=datapi$group),size=1)+
  geom_smooth(aes(color=datapi$group),span=0.1,fill="gray84",alpha=0.5)+
  scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
  #ylim(0,max(datapi$D_Pop))+
  xlab("Chromosome 7")+  ylab("Tajima's D")+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

dataancestry=read.table("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/second_round_calling/localancestry/localancestryestimates_diagnostic_300322.txt.Chr07",h=TRUE)

plotC<- ggplot(dataancestry)+
  geom_segment(aes(x=startwin,xend=startwin,y=0,yend=meanfreqhybtea),color="goldenrod2")+
  geom_segment(aes(x=startwin,xend=startwin,y=1,yend=meanfreqhybtea),color="navyblue")+
  #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
  ylim(0,1)+
  xlab(paste("Chromosome 7"))+  ylab("Local ancestry")+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

GWASgeraniolquali50=read.table("~/Rose/documents_transfert_autreordi/GWAS_204cc/scentAvg2years/geraniol_binary50_gwas_results/geraniol_binary50_scores.txt",h=TRUE)
GWASgeraniolquali50chr7=subset(GWASgeraniolquali50,GWASgeraniolquali50$Chrom=="7")

plotD<- ggplot(GWASgeraniolquali50chr7,aes(x=Position,y=-log10(geraniol_binary50.1.dom.alt_P_value)))+
  geom_point(size=1,color=ifelse(GWASgeraniolquali50chr7$geraniol_binary50.1.dom.alt_P_value<0.001,"darkred",ifelse(GWASgeraniolquali50chr7$geraniol_binary50.1.dom.alt_P_value<0.01,"orange","navyblue")))+
  #scale_color_manual(values=c("darkgrey","goldenrod2", "navyblue","hotpink1","mediumvioletred"))+
  #ylim(0,max(datapi$D_Pop))+
  xlab(paste("Chromosome 7"))+  ylab("-log10(p-values)")+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

grid.arrange(plotA,plotB,plotC,plotD,nrow=4)

dev.off()