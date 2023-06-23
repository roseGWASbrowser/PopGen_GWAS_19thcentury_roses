# TL - 130121
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/first_round_calling")
coverage=read.table(file="rose32.jointvcf.clean.PASSonly.SNPonly.vcf.2pcsites.vcf.depth.hzposonly")

# filter per indivdual all sites with cov higher than 20
coverage_to_use=subset(coverage,(coverage$V1=="175580")&(coverage$V5>=20))
coverage_to_use$V9<- ifelse((coverage_to_use$V7/coverage_to_use$V5)<=0.5,(coverage_to_use$V7/coverage_to_use$V5),(coverage_to_use$V8/coverage_to_use$V5))
length(coverage_to_use$V9)
ggplot(coverage_to_use, aes(x=V9))+ 
  geom_histogram(aes(y=..density..), colour="black", fill="white") + 
  geom_density(alpha=.2, fill="#FF6666",adjust=1.5) +
  xlim(0,0.51)+ 
  ylab("Density")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


ggplot(coverage_to_use, aes(x=V8/V5)) + geom_density()  + 
  ylab("Density")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


# filter per indivdual all sites with cov higher than 20
coverage_to_use=subset(coverage,(coverage$V5>=20))
coverage_to_use$V9<- ifelse((coverage_to_use$V7/coverage_to_use$V5)<=0.5,(coverage_to_use$V7/coverage_to_use$V5),(coverage_to_use$V8/coverage_to_use$V5))

# mode fonction
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}
# loop compute mode
mode_val=NULL
for (ind in sort(unique(coverage_to_use$V1))){
  covind=subset(coverage_to_use,coverage_to_use$V1==ind)
  mode_val=rbind(mode_val,cbind(ind,Mode(covind$V9,na.rm=TRUE)))
}
mode_val<-as.data.frame.matrix(mode_val)

pdf(file = "allelicbalance_perind.pdf",width = 10,height = 8)
for (ind in sort(unique(coverage_to_use$V1))){
 print(ind) 
 p1<- ggplot(coverage_to_use, aes(x=V9,color=V1)) + geom_density(adjust=1.5)  + 
   scale_colour_manual(values = ifelse(sort(unique(coverage_to_use$V1))==ind,"red2","grey70"))+  # "#ffe458","#cb745b"))+
   xlab("Allelic balance at heterozygous calls (minor allele)")+
    ylab("Density")+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
 print(p1)
 #Sys.sleep(2)
}
dev.off()



pdf(file = "allelicbalance_perind_interchromosomicvar.pdf",width = 10,height = 8)
for (ind in sort(unique(coverage_to_use$V1))){
  coverage_to_use_ind=subset(coverage,(coverage$V1==ind)&(coverage$V5>=20))
  coverage_to_use_ind$V9<- ifelse((coverage_to_use_ind$V7/coverage_to_use_ind$V5)<=0.5,(coverage_to_use_ind$V7/coverage_to_use_ind$V5),(coverage_to_use_ind$V8/coverage_to_use_ind$V5))
  print(ind) 
  p1<- ggplot(coverage_to_use_ind, aes(x=V9,color=V2)) + geom_density(adjust=1.5)  + 
    xlab("Allelic balance at heterozygous calls (minor allele)")+
    ylab("Density")+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  print(p1+ggtitle(ind) )
  #Sys.sleep(2)
}
dev.off()




ggplot(coverage_to_use, aes(x=V9,color=V1)) + geom_density(adjust=1.5)  + 
  ylab("Density")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
