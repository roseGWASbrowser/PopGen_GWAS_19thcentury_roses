#Qplots pour pop NATIVES
library(ggplot2)
library(cowplot)
library("gridExtra")
library(grid)
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/first_round_calling/faststructure_first_round/")

png("barplots_faststructure_220623.png",900,1200)

cols_restricted <- c("Ancient European" = "#000080", "Early European x Asian" = "#228B22", "Hybrid tea roses" = "#889F22", 
                     "Ancient Asian" = "#EEB422", "Botanical" = "#FF0000")

admix1=read.table("rose32.jointvcf.clean.PASSonly.varonly.biallelic50kwithheader.2.meanQ.final", header=TRUE)
admix2=read.table("rose32.jointvcf.clean.PASSonly.varonly.biallelic50kwithheader.3.meanQ.final", header=TRUE)
admix3=read.table("rose32.jointvcf.clean.PASSonly.varonly.biallelic50kwithheader.4.meanQ.final", header=TRUE)
admix4=read.table("rose32.jointvcf.clean.PASSonly.varonly.biallelic50kwithheader.5.meanQ.final", header=TRUE)
admix5=read.table("rose32.jointvcf.clean.PASSonly.varonly.biallelic50kwithheader.6.meanQ.final", header=TRUE)

a1 <- ggplot(admix1, aes(Name3, qval, fill=groupK)) +  geom_bar(width=0.9, stat="identity") + 
  scale_fill_manual(values = c("#000080","#EEB422"))+
  xlab("")+ ylab("Ancestry K2")+theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=0,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=0,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

a2 <- ggplot(admix2, aes(Name3, qval, fill=groupK)) +  geom_bar(width=0.9, stat="identity") + 
  scale_fill_manual(values = c("#EEB422","grey","#000080"))+
  xlab("")+ ylab("Ancestry K3")+theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=0,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=0,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

a3 <- ggplot(admix3, aes(Name3, qval, fill=groupK)) +  geom_bar(width=0.9, stat="identity") + 
  scale_fill_manual(values = c("#EEB422","grey","#000080","#FF0000"))+
  xlab("")+ ylab("Ancestry K4")+theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=0,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=0,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

a4 <- ggplot(admix4, aes(Name3, qval, fill=groupK)) +  geom_bar(width=0.9, stat="identity") + 
  scale_fill_manual(values = c("grey","#EEB422","#FF0000","purple","#000080"))+
  xlab("")+ ylab("Ancestry K5")+theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=0,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=0,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

a5 <- ggplot(admix5, aes(Name3, qval, fill=groupK)) +  geom_bar(width=0.9, stat="identity") + 
  scale_fill_manual(values = c("#000080","#EEB422","#FF0000","purple","grey","orange"))+
  xlab("")+ ylab("Ancestry K6")+theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

#grid.arrange(a1,a2,a3,a4,a5, ncol=1)

pushViewport(viewport(layout = grid.layout(11, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(a1, vp = vplayout(1:2, 1))
print(a2, vp = vplayout(3:4, 1))
print(a3, vp = vplayout(5:6, 1))
print(a4, vp = vplayout(7:8, 1))
print(a5, vp = vplayout(9:11, 1))


dev.off()