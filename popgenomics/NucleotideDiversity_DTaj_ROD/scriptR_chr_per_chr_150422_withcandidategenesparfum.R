# TL - 150422
library(grid)
library(gridExtra)
library(ggplot2)
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/second_round_calling/diversity/final_280322")

ROD=read.table("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.ROD",h=TRUE)
candidategeneparfum=read.table("candidate_genes_parfum_avril2022_fromBenoit.txt.withgenomicpositions",h=FALSE)
limitchrzero=read.table("chr0_limit_scaffold_syntheticchr0.txt",h=F)
quantileRODHyb1stEuro=quantile(ROD$RoD_pi_Hyb1stvsEurope,probs=c(0.05,0.5,0.95),na.rm=TRUE)
quantileRODHybTeaEuro=quantile(ROD$RoD_pi_HybTeavsEurope,probs=c(0.05,0.5,0.95),na.rm=TRUE)
quantileRODHyb1stAsia=quantile(ROD$RoD_pi_Hyb1stvsAsia,probs=c(0.05,0.5,0.95),na.rm=TRUE)
quantileRODHybTeaAsia=quantile(ROD$RoD_pi_HybTeavsAsia,probs=c(0.05,0.5,0.95),na.rm=TRUE)

pdf(file="diversity_ROD_KperK_150422_withcandidategenesScent.pdf",height=20,width=10)
for (chr in c(0:7)){
  datapi=read.table(paste("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.prR.Chr0",chr,sep=""),h=TRUE)
  chrfullname=paste("Chr0",chr,sep="")
  subsetgeneannotationparfum=subset(candidategeneparfum,candidategeneparfum$V4==chrfullname)
  subsetchrzeroannot=subset(limitchrzero,limitchrzero$V1==chrfullname)
  plotA<- ggplot(datapi,aes(x=startwin,y=pi_site,group=group))+
    # limitchr0
    geom_vline(xintercept = subsetchrzeroannot$V3,lwd=1,color="lightgrey",lwd=0.5)+
    # gene info parfum
    geom_vline(xintercept = subsetgeneannotationparfum$V7,lwd=1,color="cornflowerblue")+
    annotate("text",x= subsetgeneannotationparfum$V7,y=max(datapi$pi_site)*0.9,label=paste(subsetgeneannotationparfum$V1," (",subsetgeneannotationparfum$V3,")",sep=""),color="black",size=1.4,angle=90)+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",29000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",31000000,0),y=ifelse(chr=="3",max(datapi$pi_site),0),label="roKSN\n(29-33Mb)",color=ifelse(chr=="3","forestgreen","white"),size=ifelse(chr=="3",2.5,0))+
    # Number of petals
    geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="1",41514796,0),y=ifelse(chr=="1",max(datapi$pi_site),0),label="AGL6",color=ifelse(chr=="1","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="1",52096899,0),y=ifelse(chr=="1",max(datapi$pi_site),0),label="TOE2\neuAP1",color=ifelse(chr=="1","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",11833234,0),y=ifelse(chr=="2",max(datapi$pi_site),0),label="TM6",color=ifelse(chr=="2","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",17461203,0),y=ifelse(chr=="2",max(datapi$pi_site),0),label="AP2",color=ifelse(chr=="2","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",57253197,0),y=ifelse(chr=="2",max(datapi$pi_site),0),label="PLE AgD",color=ifelse(chr=="2","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",67016526,0),y=ifelse(chr=="2",max(datapi$pi_site),0),label="euFUL SEP1/2",color=ifelse(chr=="2","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",33238136,0),y=ifelse(chr=="3",max(datapi$pi_site),0),label="AP2 TOE",color=ifelse(chr=="3","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="4",46474063,0),y=ifelse(chr=="4",max(datapi$pi_site),0),label="SEP3",color=ifelse(chr=="4","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="5",8425314,0),y=ifelse(chr=="5",max(datapi$pi_site),0),label="AgC",color=ifelse(chr=="5","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="5",61446161,0),y=ifelse(chr=="5",max(datapi$pi_site),0),label="TOE1",color=ifelse(chr=="5","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="6",58611899,0),y=ifelse(chr=="6",max(datapi$pi_site),0),label="PI/GLO",color=ifelse(chr=="6","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="6",60224075,0),y=ifelse(chr=="6",max(datapi$pi_site),0),label="AP3",color=ifelse(chr=="6","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="7",4414048,0),y=ifelse(chr=="7",max(datapi$pi_site),0),label="SEP4\neuFUL",color=ifelse(chr=="7","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="7",50555186,0),y=ifelse(chr=="7",max(datapi$pi_site),0),label="D",color=ifelse(chr=="7","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="7",52218507,0),y=ifelse(chr=="7",max(datapi$pi_site),0),label="AGL12",color=ifelse(chr=="7","black","white"),size=1.4,angle=90)+    
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",41500000,0),y=ifelse(chr=="3",max(datapi$pi_site),0),label="S locus",color=ifelse(chr=="3","black","white"),size=1.4,angle=90) + 
    # generate plot
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
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
    
  plotB<- ggplot(datapi,aes(x=startwin,y=D_Pop,group=group))+
    # limitchr0
    geom_vline(xintercept = subsetchrzeroannot$V3,lwd=1,color="lightgrey",lwd=0.5)+
    # gene info parfum
    geom_vline(xintercept = subsetgeneannotationparfum$V7,lwd=1,color="cornflowerblue")+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",29000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",31000000,0),y=ifelse(chr=="3",max(datapi$D_Pop),0),label="roKSN\n(29-33Mb)",color=ifelse(chr=="3","forestgreen","white"),size=ifelse(chr=="3",2.5,0))+
    # Number of petals
    geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    # generate plot
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
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

  datadiv=read.table(paste("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.ROD.Chr",chr,sep=""),h=TRUE)
  
  plotC<- ggplot(datadiv, aes(x=startwin, y=RoD_pi_Hyb1stvsEurope)) +   
    # limitchr0
    geom_vline(xintercept = subsetchrzeroannot$V3,lwd=1,color="lightgrey",lwd=0.5)+
    # gene info parfum
    geom_vline(xintercept = subsetgeneannotationparfum$V7,lwd=1,color="cornflowerblue")+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",29000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",31000000,0),y=ifelse(chr=="3",max(datadiv$RoD_pi_Hyb1stvsEurope),0),label="roKSN\n(29-33Mb)",color="forestgreen",size=ifelse(chr=="3",2.5,0))+
    # Number of petals
    geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    # generate plot
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
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  plotD<- ggplot(datadiv, aes(x=startwin, y=RoD_pi_HybTeavsEurope)) +  
    # limitchr0
    geom_vline(xintercept = subsetchrzeroannot$V3,lwd=1,color="lightgrey",lwd=0.5)+
    # gene info parfum
    geom_vline(xintercept = subsetgeneannotationparfum$V7,lwd=1,color="cornflowerblue")+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",29000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",31000000,0),y=ifelse(chr=="3",max(datadiv$RoD_pi_HybTeavsEurope),0),label="roKSN\n(29-33Mb)",color="forestgreen",size=ifelse(chr=="3",2.5,0))+
    # Number of petals
    geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    # generate plot
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
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  plotE<- ggplot(datadiv, aes(x=startwin, y=RoD_pi_Hyb1stvsAsia)) +  
    # limitchr0
    geom_vline(xintercept = subsetchrzeroannot$V3,lwd=1,color="lightgrey",lwd=0.5)+
    # gene info parfum
    geom_vline(xintercept = subsetgeneannotationparfum$V7,lwd=1,color="cornflowerblue")+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",29000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",31000000,0),y=ifelse(chr=="3",max(datadiv$RoD_pi_Hyb1stvsAsia),0),label="roKSN\n(29-33Mb)",color="forestgreen",size=ifelse(chr=="3",2.5,0))+
    # Number of petals
    geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    # generate plot
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
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  plotF<- ggplot(datadiv, aes(x=startwin, y=RoD_pi_HybTeavsAsia)) +  
    # limitchr0
    geom_vline(xintercept = subsetchrzeroannot$V3,lwd=1,color="lightgrey",lwd=0.5)+
    # gene info parfum
    geom_vline(xintercept = subsetgeneannotationparfum$V7,lwd=1,color="cornflowerblue")+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",29000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",31000000,0),y=ifelse(chr=="3",max(datadiv$RoD_pi_HybTeavsAsia),0),label="roKSN\n(29-33Mb)",color="forestgreen",size=ifelse(chr=="3",2.5,0))+
    # Number of petals
    geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    # generate plot
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
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  grid.arrange(plotA,plotB,plotC,plotD,plotE,plotF,nrow=6)
 }
dev.off()

#### with localancestry 
pdf(file="diversity_ROD_ancestry_KperK_150422_withcandidategenesScent.pdf",height=20,width=10)
for (chr in c(0:7)){
  datapi=read.table(paste("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.prR.Chr0",chr,sep=""),h=TRUE)
  chrfullname=paste("Chr0",chr,sep="")
  subsetgeneannotationparfum=subset(candidategeneparfum,candidategeneparfum$V4==chrfullname)
  subsetchrzeroannot=subset(limitchrzero,limitchrzero$V1==chrfullname)
  plotA<- ggplot(datapi,aes(x=startwin,y=pi_site,group=group))+
    # limitchr0
    geom_vline(xintercept = subsetchrzeroannot$V3,lwd=1,color="lightgrey",lwd=0.5)+
    # gene info parfum
    geom_vline(xintercept = subsetgeneannotationparfum$V7,lwd=1,color="cornflowerblue")+
    annotate("text",x= subsetgeneannotationparfum$V7,y=max(datapi$pi_site)*0.9,label=paste(subsetgeneannotationparfum$V1," (",subsetgeneannotationparfum$V3,")",sep=""),color="black",size=1.4,angle=90)+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",29000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",31000000,0),y=ifelse(chr=="3",max(datapi$pi_site),0),label="roKSN\n(29-33Mb)",color=ifelse(chr=="3","forestgreen","white"),size=ifelse(chr=="3",2.5,0))+
    # Number of petals
    geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="1",41514796,0),y=ifelse(chr=="1",max(datapi$pi_site),0),label="AGL6",color=ifelse(chr=="1","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="1",52096899,0),y=ifelse(chr=="1",max(datapi$pi_site),0),label="TOE2\neuAP1",color=ifelse(chr=="1","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",11833234,0),y=ifelse(chr=="2",max(datapi$pi_site),0),label="TM6",color=ifelse(chr=="2","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",17461203,0),y=ifelse(chr=="2",max(datapi$pi_site),0),label="AP2",color=ifelse(chr=="2","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",57253197,0),y=ifelse(chr=="2",max(datapi$pi_site),0),label="PLE AgD",color=ifelse(chr=="2","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="2",67016526,0),y=ifelse(chr=="2",max(datapi$pi_site),0),label="euFUL SEP1/2",color=ifelse(chr=="2","black","white"),size=1.4,angle=90)+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",33238136,0),y=ifelse(chr=="3",max(datapi$pi_site),0),label="AP2 TOE",color=ifelse(chr=="3","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="4",46474063,0),y=ifelse(chr=="4",max(datapi$pi_site),0),label="SEP3",color=ifelse(chr=="4","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="5",8425314,0),y=ifelse(chr=="5",max(datapi$pi_site),0),label="AgC",color=ifelse(chr=="5","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="5",61446161,0),y=ifelse(chr=="5",max(datapi$pi_site),0),label="TOE1",color=ifelse(chr=="5","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="6",58611899,0),y=ifelse(chr=="6",max(datapi$pi_site),0),label="PI/GLO",color=ifelse(chr=="6","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="6",60224075,0),y=ifelse(chr=="6",max(datapi$pi_site),0),label="AP3",color=ifelse(chr=="6","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="7",4414048,0),y=ifelse(chr=="7",max(datapi$pi_site),0),label="SEP4\neuFUL",color=ifelse(chr=="7","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="7",50555186,0),y=ifelse(chr=="7",max(datapi$pi_site),0),label="D",color=ifelse(chr=="7","black","white"),size=1.4,angle=90)+    
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="7",52218507,0),y=ifelse(chr=="7",max(datapi$pi_site),0),label="AGL12",color=ifelse(chr=="7","black","white"),size=1.4,angle=90)+    
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",41500000,0),y=ifelse(chr=="3",max(datapi$pi_site),0),label="S locus",color=ifelse(chr=="3","black","white"),size=1.4,angle=90) + 
    # generate plot
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
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  plotB<- ggplot(datapi,aes(x=startwin,y=D_Pop,group=group))+
    # limitchr0
    geom_vline(xintercept = subsetchrzeroannot$V3,lwd=1,color="lightgrey",lwd=0.5)+
    # gene info parfum
    geom_vline(xintercept = subsetgeneannotationparfum$V7,lwd=1,color="cornflowerblue")+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",29000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",31000000,0),y=ifelse(chr=="3",max(datapi$D_Pop),0),label="roKSN\n(29-33Mb)",color=ifelse(chr=="3","forestgreen","white"),size=ifelse(chr=="3",2.5,0))+
    # Number of petals
    geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    # generate plot
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
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  
  dataancestry=read.table(paste("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/second_round_calling/localancestry/localancestryestimates_diagnostic_300322.txt.Chr0",chr,sep=""),h=TRUE)
  
  plotC<- ggplot(dataancestry)+
    # limitchr0
    geom_vline(xintercept = subsetchrzeroannot$V3,lwd=1,color="lightgrey",lwd=0.5)+
    # gene info parfum
    geom_vline(xintercept = subsetgeneannotationparfum$V7,lwd=1,color="cornflowerblue")+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",29000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",31000000,0),y=ifelse(chr=="3",1,0),label="roKSN\n(29-33Mb)",color="forestgreen",size=ifelse(chr=="3",2.5,0))+
    # Number of petals
    geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    # generate plot
    geom_segment(aes(x=startwin,xend=startwin,y=0,yend=meanfreqhybtea),color="goldenrod2")+
    geom_segment(aes(x=startwin,xend=startwin,y=1,yend=meanfreqhybtea),color="navyblue")+
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
  
  datadiv=read.table(paste("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.ROD.Chr",chr,sep=""),h=TRUE)
  plotD<- ggplot(datadiv, aes(x=startwin, y=RoD_pi_HybTeavsEurope)) +  
    # limitchr0
    geom_vline(xintercept = subsetchrzeroannot$V3,lwd=1,color="lightgrey",lwd=0.5)+
    # gene info parfum
    geom_vline(xintercept = subsetgeneannotationparfum$V7,lwd=1,color="cornflowerblue")+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",29000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",31000000,0),y=ifelse(chr=="3",max(datadiv$RoD_pi_HybTeavsEurope),0),label="roKSN\n(29-33Mb)",color="forestgreen",size=ifelse(chr=="3",2.5,0))+
    # Number of petals
    geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    # generate plot
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
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  

  plotE<- ggplot(datadiv, aes(x=startwin, y=RoD_pi_HybTeavsAsia)) +  
    # limitchr0
    geom_vline(xintercept = subsetchrzeroannot$V3,lwd=1,color="lightgrey",lwd=0.5)+
    # gene info parfum
    geom_vline(xintercept = subsetgeneannotationparfum$V7,lwd=1,color="cornflowerblue")+
    # roKSN locus
    geom_vline(xintercept = ifelse(chr=="3",29000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33000000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    annotate("text",x=ifelse(chr=="3",31000000,0),y=ifelse(chr=="3",max(datadiv$RoD_pi_HybTeavsAsia),0),label="roKSN\n(29-33Mb)",color="forestgreen",size=ifelse(chr=="3",2.5,0))+
    # Number of petals
    geom_vline(xintercept = ifelse(chr=="1",41514796,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="1",52096899,-500000),lwd=ifelse(chr=="1",1,0),color=ifelse(chr=="1","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",11833234,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",17461203,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",57253197,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="2",67016526,-500000),lwd=ifelse(chr=="2",1,0),color=ifelse(chr=="2","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="3",33238136,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="4",46474063,-500000),lwd=ifelse(chr=="4",1,0),color=ifelse(chr=="4","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",8425314,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="5",61446161,-500000),lwd=ifelse(chr=="5",1,0),color=ifelse(chr=="5","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",58611899,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="6",60224075,-500000),lwd=ifelse(chr=="6",1,0),color=ifelse(chr=="6","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",4414048,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",50555186,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    geom_vline(xintercept = ifelse(chr=="7",52218507,-500000),lwd=ifelse(chr=="7",1,0),color=ifelse(chr=="7","forestgreen","white"))+
    # S locus
    geom_vline(xintercept = ifelse(chr=="3",41500000,-500000),lwd=ifelse(chr=="3",1,0),color=ifelse(chr=="3","forestgreen","white"))+
    # generate plot
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
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  
  grid.arrange(plotA,plotB,plotC,plotD,plotE,nrow=5)
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