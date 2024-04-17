# TL - 110921
library(ggplot2)
library(RColorBrewer)
library("ggrepel")
library("grid")
library("gridExtra")
setwd("/home/thibaultleroy/Rose/Papier/Fig2/")

#write.table(pca.cows$li,file="res_PCA_270genotypes.txt",sep="\t",quote=FALSE)

pca_faststr=read.table("dataFigure2_211122_Petals_RecFlowering_Flagr_BlackSp.txt",h=T,sep = "\t")
pca_faststr2=read.table("dataFigure2_211122_Petals_RecFlowering_Flagr_BlackSp_withPCASNPpruning.txt",h=T,sep = "\t")


# Color code
cols <- c("1" = "#000080", "2" = "#082268", "3" = "#114551", "4" = "#196839",
          "5" = "#228B22", "6" = "#559522", "7" = "#889F22", "8" = "#BBA922",
          "9" = "#EEB422", "10" = "#FF0000", "11" = "#F93539", "12" = "#F36B73",
          "13" = "#EEA2AD", "14" = "#CD6F8E", "15" = "#AC3C6F", "16" = "#8B0A50", "NA" = "grey")

# color palette between 1 and 9 <-> Europe <-> Asia
colfunc <- colorRampPalette(c("navyblue", "forestgreen","goldenrod2"))
colfunc(9)

colfunc <- colorRampPalette(c("red","lightpink2","deeppink4"))
colfunc(7)


# PCA group with names -> Supplementary BEFORE SNP PRUNING

pdf("Figure2A_PCA_180623_withnames_blackbg.pdf",width = 12, height = 9)
ggplot(pca_faststr,aes(x=Axis1,y=Axis2,colour=as.factor(genetic_group)))+geom_point(size=2.5)+
  xlab("PC1 (29.3%)")+
  ylab("PC2 (10.6%)")+
  theme_bw()+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(plot.background = element_rect(fill = "black"),panel.background = element_rect(fill = "black"),legend.background = element_rect(fill = "gray30",color = "white"))+
  theme(legend.text=element_text(color="white",size=10),legend.title=element_text(color="white",size=11))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'white', size = 1.25), axis.ticks = element_line(colour = 'white', size = 1.25), 
        axis.text.x = element_text(colour="white",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="white",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="white",size=16,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="white",size=16,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(Name_simpl), fill = as.factor(genetic_group)), color = 'snow2', size = 1.75)+
  labs(colour="Genetic Groups\n(Liorzou et al. 2016)",fill="Genetic Groups\n(Liorzou et al. 2016)")+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the A
dev.off()

# PCA group without names -> Panel 1A BEFORE SNP PRUNING

align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
}

library(cowplot)
pdf("Figure2A_PCA_180623_blackbg.pdf",width = 12, height = 9)
p<- ggplot(pca_faststr,aes(x=Axis1,y=Axis2,fill=as.factor(genetic_group)))+geom_point(size=6,shape=21,colour="white")+
  xlab("PC1 (29.3%)")+
  ylab("PC2 (10.6%)")+
  theme_bw()+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(plot.background = element_rect(fill = "black"),panel.background = element_rect(fill = "black"),legend.background = element_rect(fill = "gray30",color = "white"),legend.key = element_rect(fill = "gray30"))+
  theme(legend.text=element_text(color="white",size=20),legend.title=element_text(color="white",size=20),legend.spacing.x = unit(0.15 ,'cm'),legend.spacing.y = unit(0.05 ,'cm'),legend.box.spacing = unit(0.1 ,'cm'))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'white', size = 1.25), axis.ticks = element_line(colour = 'white', size = 1.25), 
        axis.text.x = element_text(colour="white",size=18,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="white",size=18,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="white",size=20,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="white",size=20,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.title.align = 0.5,legend.key.width = unit(15, 'mm'),
      legend.direction = "vertical",
      legend.box.just = "center")+
   labs(colour="Genetic Groups\n(Liorzou et\nal. 2016)",fill="Genetic Groups\n(Liorzou et\nal. 2016)")+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the A
ggdraw(align_legend(p, hjust = 0.5))
dev.off()

# white bg version
pdf("Figure2A_PCA_211122.pdf",width = 12, height = 9)
ggplot(pca_faststr,aes(x=Axis1,y=Axis2,colour=as.factor(genetic_group)))+geom_point(size=4)+
  xlab("PC1 (29.3%)")+
  ylab("PC2 (10.6%)")+
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





# PCA group with names -> Supplementary AFTER SNP PRUNING

pdf("Figure2A_PCA_180623_withnames_blackbg_SNPpruned130324_whitebg.pdf",width = 12, height = 9)
ggplot(pca_faststr2,aes(x=Axis1_pruned,y=Axis2_pruned,colour=as.factor(genetic_group)))+geom_point(size=2.5)+
  xlab("PC1 (21.4%)")+
  ylab("PC2 (11.6%)")+
  theme_bw()+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(plot.background = element_rect(fill = "white"),panel.background = element_rect(fill = "white"),legend.background = element_rect(fill = "white",color = "black"))+
  theme(legend.text=element_text(color="black",size=10),legend.title=element_text(color="black",size=11))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=16,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(Name_simpl), fill = as.factor(genetic_group)), color = 'snow2', size = 1.75)+
  labs(colour="Genetic Groups\n(Liorzou et al. 2016)",fill="Genetic Groups\n(Liorzou et al. 2016)")+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the A
dev.off()

# PCA group without names -> Panel 1A AFTER SNP PRUNING

align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
}

library(cowplot)
pdf("Figure2A_PCA_180623_blackbg_SNPpruned130324.pdf",width = 12, height = 9)
p<- ggplot(pca_faststr2,aes(x=Axis1_pruned,y=Axis2_pruned,fill=as.factor(genetic_group)))+geom_point(size=6,shape=21,colour="white")+
  xlab("PC1 (21.4%)")+
  ylab("PC2 (11.6%)")+
  theme_bw()+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(plot.background = element_rect(fill = "black"),panel.background = element_rect(fill = "black"),legend.background = element_rect(fill = "gray30",color = "white"),legend.key = element_rect(fill = "gray30"))+
  theme(legend.text=element_text(color="white",size=20),legend.title=element_text(color="white",size=20),legend.spacing.x = unit(0.15 ,'cm'),legend.spacing.y = unit(0.05 ,'cm'),legend.box.spacing = unit(0.1 ,'cm'))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'white', size = 1.25), axis.ticks = element_line(colour = 'white', size = 1.25), 
        axis.text.x = element_text(colour="white",size=18,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="white",size=18,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="white",size=20,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="white",size=20,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.title.align = 0.5,legend.key.width = unit(15, 'mm'),
        legend.direction = "vertical",
        legend.box.just = "center")+
  labs(colour="Genetic Groups\n(Liorzou et\nal. 2016)",fill="Genetic Groups\n(Liorzou et\nal. 2016)")+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the A
ggdraw(align_legend(p, hjust = 0.5))
dev.off()

# white bg version
pdf("Figure2A_PCA_211122_SNPpruned130324.pdf",width = 12, height = 9)
ggplot(pca_faststr2,aes(x=Axis1_pruned,y=Axis2_pruned,colour=as.factor(genetic_group)))+geom_point(size=4)+
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

# comparison

pdf("Comparaison_SNPwithout_SNPpruned130324.pdf",width = 18, height = 7)
plot1 <- ggplot(pca_faststr2,aes(x=Axis1,y=Axis1_pruned,colour=as.factor(genetic_group)))+geom_point(size=4)+
  xlab("PC1 (all SNPs, 29.3%)")+
  ylab("PC1 (pruned SNPs, 21.4%)")+
  geom_abline(slope = 3.546e-01, intercept = -8.219e-11,lty=2,lwd=1.5)+
  theme_bw()+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=18,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=18,angle=90,hjust=.5,vjust=.5,face="italic"))

plot2<- ggplot(pca_faststr2,aes(x=Axis2,y=Axis2_pruned,colour=as.factor(genetic_group)))+geom_point(size=4)+
  xlab("PC2 (all SNPs, 10.6%)")+
  ylab("PC2 (pruned SNPs, 11.6%)")+
  geom_abline(slope = 3.862e-01, intercept = 6.248e-11,lty=2,lwd=1.5)+
  theme_bw()+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=18,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=18,angle=90,hjust=.5,vjust=.5,face="italic"))

pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, 2))
dev.off()



## 2B

cols_restricted <- c("Ancient European" = "#000080", "Early European x Asian" = "#228B22", "Hybrid tea roses" = "#889F22", 
          "Ancient Asian" = "#EEB422", "Botanical" = "#FF0000")

#PCA_1stround_filtered50kSNPs_061221_withmappingres_mod221122.txt # initiallement puis changé first par early European x Asian
pca_summary=read.table("PCA_1stround_filtered50kSNPs_061221_withmappingres_120324.txt",header = TRUE,sep = "\t")

pdf("Figure2B_PCA_120324_blackbg_fullnames.pdf",width = 9, height = 9)
ggplot(pca_summary,aes(x=-(axis1),axis3,colour=groupID))+geom_point(size=2.5)+
  xlab("Axis 1 (9.7%)")+
  ylab("Axis 3 (6.3%)")+
  #xlab("PC1 (10.9%)")+ # before SNP pruning
  #ylab("PC3 (7.5%)")+# before SNP pruning
  scale_colour_manual(values = cols_restricted )+
  scale_fill_manual(values = cols_restricted )+
  theme_bw()+
  theme(plot.background = element_rect(fill = "black"),panel.background = element_rect(fill = "black"),legend.background = element_rect(fill = "gray30",color = "white"),legend.key = element_rect(fill = "gray30"))+
  theme(legend.text=element_text(color="white",size=9),legend.title=element_text(color="white",size=10))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'white', size = 1.25), axis.ticks = element_line(colour = 'white', size = 1.25), 
        axis.text.x = element_text(colour="white",size=17,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="white",size=17,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="white",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="white",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  labs(colour="Groups",fill="Groups")+
  theme(legend.position='top') +
  geom_label_repel(aes(label = factor(sampleID),fill=groupID), color = 'white', size = 3,segment.color = 'white')+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the "a" in the legend
dev.off()

pdf("Figure2B_PCA_120324_blackbg_justnumbers.pdf",width = 9, height = 9)
ggplot(pca_summary,aes(x=-(axis1),axis3,colour=groupID))+geom_point(size=2.5)+
  xlab("Axis 1 (9.7%)")+
  ylab("Axis 3 (6.3%)")+
  #xlab("PC1 (10.9%)")+ # before SNP pruning
  #ylab("PC3 (7.5%)")+# before SNP pruning
  scale_colour_manual(values = cols_restricted )+
  scale_fill_manual(values = cols_restricted )+
  theme_bw()+
  theme(plot.background = element_rect(fill = "black"),panel.background = element_rect(fill = "black"),legend.background = element_rect(fill = "gray30",color = "white"),legend.key = element_rect(fill = "gray30"))+
  theme(legend.text=element_text(color="white",size=9),legend.title=element_text(color="white",size=10))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'white', size = 1.25), axis.ticks = element_line(colour = 'white', size = 1.25), 
        axis.text.x = element_text(colour="white",size=17,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="white",size=17,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="white",size=18,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="white",size=18,angle=90,hjust=.5,vjust=.5,face="italic"))+
  labs(colour="Groups",fill="Groups")+
  theme(legend.position=c(0.35,0.88),legend.title=element_text(size=0), legend.text=element_text(size=20),legend.spacing.x = unit(0.15 ,'cm'),legend.spacing.y = unit(0.05 ,'cm'),legend.box.spacing = unit(0.1 ,'cm')) +
  geom_label_repel(aes(label = factor(sampleID4),fill=groupID), color = 'white', size = 7,segment.color = 'white')+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the "a" in the legend
dev.off()


# white bg version for supplementary


plotminorA <- ggplot(pca_summary,aes(x=axis1,axis2,colour=groupID))+geom_point(size=2.5)+
  xlab("Axis 1 (9.7%)")+
  ylab("Axis 2 (8.2%)")+
  scale_colour_manual(values = cols_restricted )+
  scale_fill_manual(values = cols_restricted )+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.position="none") +
  geom_label_repel(aes(label = factor(sampleID4),fill=groupID), color = 'white', size = 3,segment.color = 'black')+
  labs(col="",fill="")+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the "a" in the legend

plotminorB <- ggplot(pca_summary,aes(x=axis1,axis4,colour=groupID))+geom_point(size=2.5)+
  xlab("Axis 1 (9.7%)")+
  ylab("Axis 4 (5.9%)")+
  
  scale_colour_manual(values = cols_restricted )+
  scale_fill_manual(values = cols_restricted )+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.position="none") +
  geom_label_repel(aes(label = factor(sampleID4),fill=groupID), color = 'white', size = 3,segment.color = 'black')+
  labs(col="",fill="")+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the "a" in the legend

plotminorC <- ggplot(pca_summary,aes(x=axis2,axis3,colour=groupID))+geom_point(size=2.5)+
  xlab("Axis 2 (7.2%)")+
  ylab("Axis 3 (6.3%)")+
  scale_colour_manual(values = cols_restricted )+
  scale_fill_manual(values = cols_restricted )+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.position="none") +
  geom_label_repel(aes(label = factor(sampleID4),fill=groupID), color = 'white', size = 3,segment.color = 'black')+
  labs(col="",fill="")+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the "a" in the legend


plotmain <- ggplot(pca_summary,aes(x=axis1,axis3,colour=groupID))+geom_point(size=2.5)+
  xlab("Axis 1 (9.7%)")+
  ylab("Axis 3 (6.3%)")+
  scale_colour_manual(values = cols_restricted )+
  scale_fill_manual(values = cols_restricted )+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.position='top') +
  geom_label_repel(aes(label = factor(sampleID4),fill=groupID), color = 'white', size = 3,segment.color = 'black')+
  labs(col="",fill="")+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the "a" in the legend

pdf("SupFigureX_PCA_WGS_afterSNPpruning_130324.pdf",width = 16, height = 12)
pushViewport(viewport(layout = grid.layout(3, 3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(plotmain, vp = vplayout(1:3, 1:2))
print(plotminorA, vp = vplayout(1, 3))
print(plotminorB, vp = vplayout(2, 3))
print(plotminorC, vp = vplayout(3, 3))
dev.off()

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
  geom_label_repel(aes(label = factor(sampleID3),fill=pc_mapped), color = 'black', size = 3)+
  scale_fill_gradient2(midpoint=90,low  = "navyblue",mid = "lightblue", high = "red",na.value = "grey50")+
  scale_colour_gradient2(midpoint=90,low  = "navyblue", mid = "lightblue", high = "red",na.value = "grey50")+
  labs(col="Mapping rate (%)",fill="Mapping rate (%)")

cols_restricted <- c("Ancient_European" = "#000080", "First_European_x_Asian" = "#228B22", "Tea_Hybrids" = "#889F22", 
                     "Ancient_Asian" = "#EEB422", "Botanical" = "#FF0000")


pdf("Figure2B_PCA_211122.pdf",width = 9, height = 9)
ggplot(pca_summary,aes(x=-(axis1),axis3,colour=groupID))+geom_point(size=2.5)+
  xlab("Axis 1 (10.9%)")+
  ylab("Axis 3 (7.5%)")+
  scale_colour_manual(values = cols_restricted )+
  scale_fill_manual(values = cols_restricted )+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  theme(legend.position='top') +
  geom_label_repel(aes(label = factor(sampleID3),fill=groupID), color = 'white', size = 3,segment.color = 'black')+
  labs(col="",fill="")+
  guides(fill = guide_legend(override.aes = aes(label = ""))) # to remove the "a" in the legend
dev.off()



### Panel 2C network
IDs=read.table("IDs4kinshipnetwork",h=T)
relationships=read.table("relationships4network.txt",h=T)
require(igraph)

par(mar=c(0,0,0,0))
g=graph.data.frame(relationships, directed=FALSE)
graph_layout=layout.fruchterman.reingold(g)
plot(g,layout=graph_layout)

plot(g,layout=graph_layout,edge.width=4-relationships$Weight)  # ,vertex.color=IDs$Type


net = graph_from_data_frame(relationships,vertices = IDs,directed = F)
colrs.v = c(AncientAsia = "#EEB422", AncientEU = "#000080", Botanical = "#FF0000", FirstHyb = "#228B22", HybTea = "#889F22",Ambiguous="white") #node colours
V(net)$color = colrs.v[V(net)$Type]

colrs.e = c(level1 = "firebrick1", level2 = "darkorange",level3="gold1") #edge colours
E(net)$color = colrs.e[E(net)$Weight] 

#https://r-graph-gallery.com/248-igraph-plotting-parameters.html
par(bg="black")
#plot(net, edge.curved=0.2,vertex.size=8,edge.width=4,vertex.label=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32"),vertex.label.cex=1.1,vertex.shape = shapes()[8],vertex.label.color="white",vertex.frame.color = "white",vertex.label.font=2,vertex.label.family="Times")
plot(net, edge.curved=0.2,vertex.size=8,edge.width=4,vertex.label=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32"),vertex.label.cex=1.1,vertex.shape = shapes()[3],vertex.label.color="white",vertex.frame.color = "white",vertex.label.font=2,vertex.label.family="Times")


# panel 2D - diversity

sumstatspi=read.table("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.modif180623.prR",h=T,sep="\t")

# version 180623 - without chr0 # with early 
sumstatspiwithoutchr0=subset(sumstatspi,sumstatspi$chr!="Chr00")
ggplot(sumstatspiwithoutchr0, aes(x=group, y=pi_site)) +
  geom_jitter(aes(color=sumstatspiwithoutchr0$chr),shape=20, position=position_jitter(0.15),cex=0.2)+#geom_jitter(shape=20, position=position_jitter(0.3))+
  geom_boxplot(width=0.4,outlier.size=0,alpha=0.5)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  #scale_colour_manual(values=group.colors)+
  xlab("Rose groups")+ ylab(expression(pi))+
  #ylim(0,0.03)+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  

# sort the file
library(dplyr)
sumstatspiwithoutchr0 <- sumstatspiwithoutchr0 %>%
  mutate( group=factor(group,levels=c("All accessions","Hybrid tea roses","Early European x Asian","Ancient Asian","Ancient European")) )

pdf("Figure2D_diversity_120324.pdf",width = 9, height = 9)
ggplot(sumstatspiwithoutchr0, aes(x=group, y=pi_site)) +
  # 12/03/24: no dots #geom_jitter(aes(color=sumstatspiwithoutchr0$chr),shape=20, position=position_jitter(0.40),cex=0.5)+#geom_jitter(shape=20, position=position_jitter(0.3))+
  geom_hline(yintercept = 0.0167325040,color="#000080",linetype="longdash",lwd=1.1)+ # median values
  geom_hline(yintercept = 0.0120770664,color="#EEB422",linetype="longdash",lwd=1.1)+
  geom_hline(yintercept = 0.0154966073,color="#228B22",linetype="longdash",lwd=1.1)+
  geom_hline(yintercept = 0.0123191425+0.00005,color="#889F22",linetype="longdash",lwd=1.1,lty=2)+ #very slight shift to increase the readability
  geom_violin(alpha=0.5,color="white")+
  geom_boxplot(width=0.1,outlier.size=0,outlier.shape = NA,alpha=0.5,lwd=1.2,color=c("black","#889F22","#228B22","#EEB422","#000080"))+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7,aes(fill=pinpis$Taxa))+
  scale_color_brewer(palette = rev("Set1"))+
  coord_flip()+ 
  xlab("")+ ylab(expression(paste("Nucleotide diversity (",pi,")")))+
  #ylim(0,0.03)+
  #theme(legend.key.size = unit(1.5,"line"))+
  theme_bw()+
  labs(colour="Chromosomes",fill="Chromosomes")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(legend.position='bottom') +
  theme(plot.background = element_rect(fill = "black"),panel.background = element_rect(fill = "black"),legend.background = element_rect(fill = "gray30",color = "white"),legend.key = element_rect(fill = "gray30"))+
  theme(legend.text=element_text(color="white",size=13),legend.title=element_text(color="white",size=12))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
  theme(axis.line = element_line(colour = 'white', size = 1.25), axis.ticks = element_line(colour = 'white', size = 1.25), 
        axis.text.x = element_text(colour="white",size=18,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="white",size=18,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="white",size=18,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="white",size=18,angle=90,hjust=.5,vjust=.5,face="italic"))+
  annotate("text", y=0.033,x=5,label = "a",color="white",size=7,fontface=2)+
  annotate("text", y=0.03,x=4,label = "d",color="white",size=7,fontface=2)+
  annotate("text", y=0.03,x=3,label = "b",color="white",size=7,fontface=2)+
  annotate("text", y=0.026,x=2,label = "d",color="white",size=7,fontface=2)+
  annotate("text", y=0.027,x=1,label = "c",color="white",size=7,fontface=2)
dev.off()

# median values per group
aggregate(data=sumstatspi, pi_site~group, summary)

library(agricolae)
modelp1<-aov(data=sumstatspi,pi_site~group)
HSD.test(modelp1,"group", group=TRUE,console=TRUE)


#Study: modelp1 ~ "group"
#HSD Test for pi_site 
#Mean Square Error:  1.139683e-05 
#group,  means

#pi_site         std    r          Min        Max
#All samples            0.01465169 0.003689961 4629 0.0005271455 0.02504209
#Ancient Asian          0.01201257 0.003592339 4629 0.0000000000 0.02787033
#Ancient European       0.01673758 0.002803474 4629 0.0000000000 0.03108274
#First European x Asian 0.01520233 0.003433318 4629 0.0000000000 0.02793215
#Tea hybrids            0.01213776 0.003288816 4629 0.0000000000 0.02393039

#Alpha: 0.05 ; DF Error: 23140 
#Critical Value of Studentized Range: 3.857957 

#Minimun Significant Difference: 0.0001914282 

#Treatments with the same letter are not significantly different.

#pi_site groups
#Ancient European       0.01673758      a
#First European x Asian 0.01520233      b
#All samples            0.01465169      c
#Tea hybrids            0.01213776      d
#Ancient Asian          0.01201257      d


### Supplementary ROD per chromosome

### Reduction of diversity, computed as 1 - (pi new / pi ancien), the value ranges between 1 (100% of the diversity loss, and - infiinite). Positive values indicate a net loss of diversity, negative a gain!
dataROD=read.table("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.ROD2_220623.prR",h=TRUE,sep="\t")

ggplot(dataROD, aes(x=comparison, y=RoD_pi)) +   facet_wrap(~chr, ncol = 4)+ 
  geom_hline(yintercept = -0.25,col="grey",lwd=0.2,lty=2)+
  geom_hline(yintercept = -0.5,col="grey",lwd=0.5)+
  geom_hline(yintercept = -0.75,col="grey",lwd=0.2)+
  geom_hline(yintercept = -1,col="grey",lwd=0.5)+
  geom_hline(yintercept = 0,col="red",lwd=1.1)+
  geom_hline(yintercept = 0.25,col="grey",lwd=0.2)+
  geom_hline(yintercept = 0.5,col="grey",lwd=0.5)+
  geom_hline(yintercept = 0.75,col="grey",lwd=0.2)+
  geom_hline(yintercept = 1,col="grey",lwd=0.5)+
  geom_jitter(aes(color=dataROD$chr),shape=20, position=position_jitter(0.25),cex=0.2)+#geom_jitter(shape=20, position=position_jitter(0.3))+
  geom_violin(width=0.6,alpha=0.2)+
  geom_boxplot(width=0.1,alpha=0.2,outlier.size=0)+
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
        axis.text.x = element_text(colour="black",size=9,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=7,angle=0,hjust=1,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
