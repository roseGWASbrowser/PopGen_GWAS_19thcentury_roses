# TL - 051124
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ggpubr)
#plot res from PCA with poolfstat (Mathieu Gautier, PCI)
# PCA performed on allele counts at 65 million SNPs (no ploidy assumption)

setwd("/home/tleroy/Papiers/Rose_19eme/PCIEvolBiol/First_round_PCI")
pca_result=read.table("rose32.jointvcf.clean.PASSonly.varonly.vcf.allele_counts_table.tsv.withoutSNPinfo.65M.PCAres.randomallele.txt.PCArespoolfstat_formatted.txt",header=TRUE)
pca_summary=read.table("PCA_1stround_filtered50kSNPs_061221_withmappingres_120324.txt",h=T,sep = "\t") # PCA main text

# Color code
cols <- c("1" = "#000080", "2" = "#000080", "3" = "#000080", "4" = "#000080",
          "5" = "#000080", "6" = "#000080", "7" = "#EEB422", "8" = "#EEB422",
          "9" = "#EEB422", "10" = "#EEB422", "11" = "#EEB422", "12" = "#EEB422",
          "13" = "#EEB422", "14" = "#EEB422", "15" = "#228B22", "16" = "#228B22", 
          "17"= "#228B22", "18" = "#228B22", "19"="#889F22", "20" = "#889F22",
          "21" = "#889F22", "22" = "#889F22", "23" = "#889F22", "24" = "#889F22", 
          "25" = "#889F22", "26" = "#FF0000", "27" = "#FF0000", "28" = "#FF0000",
          "29" = "#FF0000", "30" = "#FF0000", "31" = "#FF0000", "32" ="#FF0000")


cols_restricted <- c("Ancient European" = "#000080", "Early European x Asian" = "#228B22", "Hybrid tea roses" = "#889F22", 
                     "Ancient Asian" = "#EEB422", "Botanical" = "#FF0000")

# contributions PCs (i.e. pca_result)
#PC1 8.16021544280086
#PC2 7.07883052940982
#PC3 6.32576493561204
#PC4 6.13326332834651
#PC5 4.50775398677962

pdf(width=16,height=8,file="SupFigureX_PCA_allelecounts_vs_PCAcallsFig2B_131224.pdf")

plot1 <- ggplot(pca_summary,aes(x=axis1,axis3,colour=as.factor(sampleID4)))+geom_point(size=2.5)+
  xlab("SNPRelate - PC1 (9.7%)")+
  ylab("SNPRelate - PC3 (6.3%)")+
  scale_colour_manual(values = cols )+
  scale_fill_manual(values = cols )+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=18,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=18,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.position="none") +
  geom_label_repel(aes(label = as.factor(sampleID4),fill=as.factor(sampleID4)), color = 'white', size = 3,segment.color = 'black')+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the "a" in the legend



plot2 <- ggplot(pca_result,aes(x=Axis1,y=Axis3,colour=as.factor(ID_number)))+geom_point(size=4)+
  xlab("Poolfstat - PC1 (8.2%)")+
  ylab("Poolfstat - PC3 (6.3%)")+
  theme_bw()+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=18,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=18,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(ID_number),fill=as.factor(ID_number)), color = 'white', size = 3,segment.color = 'black')+
  theme(legend.position="none") +
  labs(col="",fill="")+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the "a" in the legend


ggarrange(plot1, plot2,ncol=2)

dev.off()






ggplot(pca_result,aes(x=Axis1,y=Axis2,colour=as.factor(ID_number)))+geom_point(size=4)+
  xlab("PC1 (8.2%)")+
  ylab("PC2 (7.1%)")+
  theme_bw()+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=18,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=18,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(ID_number),fill=as.factor(ID_number)), color = 'white', size = 3,segment.color = 'black')+
  theme(legend.position="none") +
  labs(col="",fill="")+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the "a" in the legend

















# white bg version
pdf("Figure2A_PCA_211122_SNPpruned130324.pdf",width = 12, height = 9)
ggplot(pca_faststr2,aes(x=Axis1,y=Axis2,colour=as.factor(genetic_group)))+geom_point(size=4)+
  xlab("PC1 (21.4%)")+
  ylab("PC2 (11.6%)")+
  theme_bw()+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=18,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=18,angle=90,hjust=.5,vjust=.5,face="italic"))+
  labs(colour="Genetic Groups\n(Liorzou et al. 2016)",fill="Genetic Groups\n(Liorzou et al. 2016)")
dev.off()










plot(pca_result$Axis1, pca_result$Axis2,
xlab = "PC1", ylab = "PC2",
#main = "PCA of Samples with Variable Ploidy",
col = "blue", pch = 20)
text(pca_result$Axis1, pca_result$Axis2, labels = pca_result$ID_number, pos = 3, cex = 0.8)

plot(pca_result$Axis1, pca_result$Axis3,
     xlab = "PC1", ylab = "PC3",
     #main = "PCA of Samples with Variable Ploidy",
     col = "blue", pch = 20)
text(pca_result$Axis1, pca_result$Axis3, labels = pca_result$ID_number, pos = 3, cex = 0.8)


plot(pca_result$Axis1, pca_result$Axis2,
     xlab = "PC1", ylab = "PC2",
     #main = "PCA of Samples with Variable Ploidy",
     col = "blue", pch = 20)
text(pca_result$Axis1, pca_result$Axis2, labels = pca_result$ID_long, pos = 3, cex = 0.8)

plot(pca_result$Axis1, pca_result$Axis3,
     xlab = "PC1", ylab = "PC3",
     #main = "PCA of Samples with Variable Ploidy",
     col = "blue", pch = 20)
text(pca_result$Axis1, pca_result$Axis3, labels = pca_result$ID_long, pos = 3, cex = 0.8)
