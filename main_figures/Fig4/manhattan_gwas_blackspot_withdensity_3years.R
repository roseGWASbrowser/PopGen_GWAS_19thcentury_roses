## TL - 160324
library(qqman)

setwd("/home/thibaultleroy/Rose/Papier/Fig4/")
densitygene=read.table("densitygene_slidwin500kb.txt",h=T)

candidateRgene=read.table("list_candidateRgenes_Laurence.txt",h=F)


#gwasres=read.table("TN2014_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="TN2014"

gwasres1=read.table("TN2014_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait1="TN2014"
gwasres2=read.table("TN2015_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait2="TN2015"
gwasres3=read.table("TN2016_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait3="TN2016"

# remove chr00 -> reviewer wants chr0
#densitygene=subset(densitygene,densitygene$chr!=0)
#candidateRgene=subset(candidateRgene,candidateRgene$V2!=0)

#"navyblue","dodgerblue3"))
#"darkolivegreen","darkolivegreen4"))
#"darkmagenta","darkorchid1"))

par(mfrow=c(3,1))
gwasres1_summary=as.data.frame(cbind(gwasres1$Marker,gwasres1$Chrom,gwasres1$Position,gwasres1$TN2014.general_P_value))
colnames(gwasres1_summary)<-c("SNP","CHR","BP","P")
is.numeric(gwasres1_summary$P)
gwasres1_summary<-na.omit(gwasres1_summary)
manhattan(gwasres1_summary,ylim=c(0,8),main="Susceptibility to blackspot disease - 2014",col = c("navyblue","dodgerblue3"),suggestiveline = F,genomewideline = F)

gwasres2_summary=as.data.frame(cbind(gwasres2$Marker,gwasres2$Chrom,gwasres2$Position,gwasres2$TN2015.general_P_value))
colnames(gwasres2_summary)<-c("SNP","CHR","BP","P")
is.numeric(gwasres2_summary$P)
gwasres2_summary<-na.omit(gwasres2_summary)
manhattan(gwasres2_summary,ylim=c(0,8),main="Susceptibility to blackspot disease - 2015",col = c("darkolivegreen","darkolivegreen4"),suggestiveline = F,genomewideline = F)

gwasres3_summary=as.data.frame(cbind(gwasres3$Marker,gwasres3$Chrom,gwasres3$Position,gwasres3$TN2016.general_P_value))
colnames(gwasres3_summary)<-c("SNP","CHR","BP","P")
is.numeric(gwasres3_summary$P)
gwasres3_summary<-na.omit(gwasres3_summary)
manhattan(gwasres3_summary,ylim=c(0,8),main="Susceptibility to blackspot disease - 2016",col = c("darkmagenta","darkorchid1"),suggestiveline = F,genomewideline = F)

