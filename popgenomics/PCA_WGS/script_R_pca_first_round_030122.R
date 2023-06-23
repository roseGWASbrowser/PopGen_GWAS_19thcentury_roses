# TL - 030121
library(ggplot2)
library("ggrepel")
library("grid")
library("gridExtra")
library(ggtree)
library(ape)

# perform PCA
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/first_round_calling/pca_first_round/")
library("SNPRelate")
vcf.fn<-"rose32.jointvcf.clean.PASSonly.varonly.biallelic50kwithheader.vcf" # finalVCF 0.0151 of all sites among 69,676,016 PASS sites (33,332,595 SNPs!)  subsampling: awk '$4 == "A" || $4 == "C" || $4 == "G" || $4 == "T" {print $0}' rose32.jointvcf.clean.PASSonly.varonly.vcf | awk '5 == "A" || $5 == "C" || $5 == "G" || $5 == "T" {print $0}' | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .00151) print $0}' > rose32.jointvcf.clean.PASSonly.varonly.biallelic50k.vcf
#vcf.fn<-"rose32.jointvcf.clean.PASSonly.SNPonly.vcf.2pcsites.vcf" # finalVCF 2% of all sites among 33,332,595 SNPs!
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- openfn.gds("ccm.gds")


# tree
#dendogram
snpgdsVCF2GDS(vcf.fn, "my.gds")
snpgdsSummary("my.gds")
genofile <- openfn.gds("my.gds")

dissMatrix  =  snpgdsIBS(genofile , sample.id=NULL, snp.id=NULL, autosome.only=TRUE, 
                         remove.monosnp=TRUE,  maf=NaN, missing.rate=NaN, num.thread=2, verbose=TRUE)
snpgdsClose(genofile)

snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)

cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL, 
                        col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, label.Z=TRUE, 
                        verbose=TRUE)

snpgdsDrawTree(cutTree, main = "32 roses",edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),y.label.kinship=F,leaflab="perpendicular")
snpgdsDrawTree(cutTree, type = "dendrogram")

cutTreedendro = cutTree$dendrogram

tree2 = plot(as.phylo(as.hclust(cutTreedendro)))  # , layout="circular",color='darkgreen', branch.length="branch.length")+ geom_tiplab(size=2.5, aes(angle=angle))+ ggtitle("32 roses")
tree2

# pca

ccm_pca<-snpgdsPCA(genofile)

# Generate plots
plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2] ,col=as.numeric(substr(ccm_pca$sample, 1,3) == 'CCM')+3, pch=20)

# with ggplot2, first transform file
pca_summary = noquote(cbind(ccm_pca$sample,ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2],ccm_pca$eigenvect[,3],ccm_pca$eigenvect[,4],ccm_pca$eigenvect[,5]))
pca_summary = as.data.frame.matrix(pca_summary)
colnames(pca_summary)<-c("sampleID","axis1","axis2","axis3","axis4","axis5")
#write.table(pca_summary, file = "PCA_1stround_filtered50kSNPs_061221.txt",sep = "\t",quote = FALSE)

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/first_round_calling/pca_first_round/")
# just add index for the first column + add groups ID
pca_summary=read.table("PCA_1stround_filtered50kSNPs_061221.txt",header = TRUE)
# then plot

pdf(file = "PCA_firstround_filtered50kSNPs.pdf",height = 10, width = 12)
ggplot(pca_summary,aes(x=axis1,axis2))+geom_point(size=2.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 2 (8.9%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(sampleID2),fill=groupID), color = 'black', size = 3)


ggplot(pca_summary,aes(x=axis1,axis3))+geom_point(size=2.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 3 (7.5%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(sampleID2),fill=groupID), color = 'black', size = 3)


ggplot(pca_summary,aes(x=axis1,axis4))+geom_point(size=2.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 4 (7.2%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(sampleID2),fill=groupID), color = 'black', size = 3)

ggplot(pca_summary,aes(x=axis1,axis5))+geom_point(size=2.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 5 (5.0%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(sampleID2),fill=groupID), color = 'black', size = 3)


ggplot(pca_summary,aes(x=axis2,axis3))+geom_point(size=2.5)+
  xlab("Axis 2 (8.9%)")+
  ylab("Axis 3 (7.5%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(sampleID2),fill=groupID), color = 'black', size = 3)

dev.off()


plotminorA <- ggplot(pca_summary,aes(x=axis1,axis2))+geom_point(size=2.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 2 (8.9%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=11,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=11,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.position = "none")+ 
  geom_label_repel(aes(label = factor(substring(sampleID2,0,8)),fill=groupID), color = 'black', size = 1.8)

plotminorB <- ggplot(pca_summary,aes(x=axis1,axis4))+geom_point(size=2.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 4 (7.2%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=11,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=11,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.position = "none")+ 
  geom_label_repel(aes(label = factor(substring(sampleID2,0,8)),fill=groupID), color = 'black', size = 1.8)

plotminorC <- ggplot(pca_summary,aes(x=axis2,axis3))+geom_point(size=2.5)+
  xlab("Axis 2 (8.9%)")+
  ylab("Axis 3 (7.5%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=11,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=11,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.position = "none")+ 
  geom_label_repel(aes(label = factor(substring(sampleID2,0,8)),fill=groupID), color = 'black', size = 1.8)


plotmain <- ggplot(pca_summary,aes(x=axis1,axis3))+geom_point(size=2.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 3 (7.5%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.position='bottom')  + 
  geom_label_repel(aes(label = factor(sampleID2),fill=groupID), color = 'black', size = 3)



pushViewport(viewport(layout = grid.layout(3, 3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(plotmain, vp = vplayout(1:3, 1:2))
print(plotminorA, vp = vplayout(1, 3))
print(plotminorB, vp = vplayout(2, 3))
print(plotminorC, vp = vplayout(3, 3))



pca_summary=read.table("PCA_1stround_filtered50kSNPs_061221_withmappingres.txt",header = TRUE)
ggplot(pca_summary,aes(x=axis1,axis3))+geom_point(size=2.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 3 (7.5%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.position='bottom')  + 
  geom_label_repel(aes(label = factor(sampleID2),fill=groupID), color = 'black', size = 3)


#% mapped
ggplot(pca_summary,aes(x=axis1,axis3,colour=pc_mapped))+geom_point(size=2.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 3 (7.5%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(IDshort),fill=pc_mapped), color = 'black', size = 3)+
scale_fill_gradient2(midpoint=80,low  = "navyblue",mid = "lightblue", high = "red",na.value = "grey50")+
scale_colour_gradient2(midpoint=80,low  = "navyblue", mid = "lightblue", high = "red",na.value = "grey50")


ggplot(pca_summary,aes(x=axis1,axis3,colour=pc_properly_paired))+geom_point(size=2.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 3 (7.5%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(IDshort),fill=pc_properly_paired), color = 'black', size = 3)+
  scale_fill_gradient2(midpoint=80,low  = "navyblue",mid = "lightblue", high = "red",na.value = "grey50")+
  scale_colour_gradient2(midpoint=80,low  = "navyblue", mid = "lightblue", high = "red",na.value = "grey50")







ggplot(pca_summary,aes(x=axis1,axis3,color=groupID))+geom_point(size=3.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 3 (7.5%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))#+
  #geom_label_repel(aes(label = factor(sampleID2),fill=groupID), color = 'black', size = 3)

