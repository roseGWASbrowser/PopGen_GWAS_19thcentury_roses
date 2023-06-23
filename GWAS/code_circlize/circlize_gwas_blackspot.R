## TL - 100621

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/")

gwasres1=read.table("./TN2016_gwas_results/TN2016_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
gwasres2=read.table("./TN2015_gwas_results/TN2015_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
gwasres3=read.table("./TN2014_gwas_results/TN2014_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE

# 1st circlize: general model
library("circlize")
circos.clear()
par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
circos.par("track.height"=0.18,start.degree=90,gap.degree=4,points.overflow.warning=FALSE,cell.padding=c(0, 0, 0, 0),track.margin=c(0.005,0.005))
circos.initialize(factors = gwasres1$Chrom, x = gwasres1$Position)

# 1st track (external): 2016
circos.trackPlotRegion(factors=  gwasres1$Chrom, y=-log10(gwasres1$TN2016.general_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2016.general_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres1$Chrom, gwasres1$Position, -log10(gwasres1$TN2016.general_P_value),pch=20,col=ifelse(-log10(gwasres1$TN2016.general_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres1$TN2016.general_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres1$TN2016.general_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres1$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres1$TN2016.general_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres1$TN2016.general_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres1$TN2016.general_P_value_FDR)>=-log10(0.20),1.5,0.75))))

# 2nd track : 2015
circos.trackPlotRegion(factors=  gwasres2$Chrom, y=-log10(gwasres2$TN2015.general_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2015.general_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres2$Chrom, gwasres2$Position, -log10(gwasres2$TN2015.general_P_value),pch=20,col=ifelse(-log10(gwasres2$TN2015.general_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres2$TN2015.general_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres2$TN2015.general_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres2$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres2$TN2015.general_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres2$TN2015.general_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres2$TN2015.general_P_value_FDR)>=-log10(0.20),1.5,0.75))))


# 3rd track: 2014
circos.trackPlotRegion(factors=  gwasres3$Chrom, y=-log10(gwasres3$TN2014.general_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2014.general_P_value){
  circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.5)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), sector.index, facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
#circos.trackPoints(gwasres$Chrom, gwasres$Position, -log10(gwasres$mean_petals.additive_P_value),pch=20,col=ifelse(-log10(gwasres$mean_petals.additive_P_value_FDR)>=-log10(0.05),"red",ifelse(-log10(gwasres$mean_petals.additive_P_value_FDR)>=-log10(0.10),"orange",ifelse(-log10(gwasres$mean_petals.additive_P_value_FDR)>=-log10(0.20),"goldenrod",ifelse(gwasres$Chrom%%2==0,"navyblue","dodgerblue3")))))
circos.trackPoints(gwasres3$Chrom, gwasres3$Position, -log10(gwasres3$TN2014.general_P_value),pch=20,col=ifelse(-log10(gwasres3$TN2014.general_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres3$TN2014.general_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres3$TN2014.general_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres3$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres3$TN2014.general_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres3$TN2014.general_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres3$TN2014.general_P_value_FDR)>=-log10(0.20),1.5,0.75))))


text(x=0.4,y=-0.05,font=2,labels="2014",cex=2,col="purple")
text(x=0.68,y=-0.05,font=2,labels="2015",cex=2,col="purple")
text(x=0.95,y=-0.05,font=2,labels="2016",cex=2,col="purple")

text(x=0,y=0,font=2,labels="General\nModel",cex=2.5)



### 2nd circlize : additive
library("circlize")
circos.clear()
par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
circos.par("track.height"=0.18,start.degree=90,gap.degree=4,points.overflow.warning=FALSE,cell.padding=c(0, 0, 0, 0),track.margin=c(0.005,0.005))
circos.initialize(factors = gwasres1$Chrom, x = gwasres1$Position)

# 1st track (external): 2016
circos.trackPlotRegion(factors=  gwasres1$Chrom, y=-log10(gwasres1$TN2016.additive_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2016.additive_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres1$Chrom, gwasres1$Position, -log10(gwasres1$TN2016.additive_P_value),pch=20,col=ifelse(-log10(gwasres1$TN2016.additive_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres1$TN2016.additive_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres1$TN2016.additive_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres1$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres1$TN2016.additive_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres1$TN2016.additive_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres1$TN2016.additive_P_value_FDR)>=-log10(0.20),1.5,0.75))))

# 2nd track : 2015
circos.trackPlotRegion(factors=  gwasres2$Chrom, y=-log10(gwasres2$TN2015.additive_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2015.additive_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres2$Chrom, gwasres2$Position, -log10(gwasres2$TN2015.additive_P_value),pch=20,col=ifelse(-log10(gwasres2$TN2015.additive_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres2$TN2015.additive_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres2$TN2015.additive_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres2$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres2$TN2015.additive_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres2$TN2015.additive_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres2$TN2015.additive_P_value_FDR)>=-log10(0.20),1.5,0.75))))


# 3rd track: 2014
circos.trackPlotRegion(factors=  gwasres3$Chrom, y=-log10(gwasres3$TN2014.additive_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2014.additive_P_value){
  circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.5)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), sector.index, facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
#circos.trackPoints(gwasres$Chrom, gwasres$Position, -log10(gwasres$mean_petals.additive_P_value),pch=20,col=ifelse(-log10(gwasres$mean_petals.additive_P_value_FDR)>=-log10(0.05),"red",ifelse(-log10(gwasres$mean_petals.additive_P_value_FDR)>=-log10(0.10),"orange",ifelse(-log10(gwasres$mean_petals.additive_P_value_FDR)>=-log10(0.20),"goldenrod",ifelse(gwasres$Chrom%%2==0,"navyblue","dodgerblue3")))))
circos.trackPoints(gwasres3$Chrom, gwasres3$Position, -log10(gwasres3$TN2014.additive_P_value),pch=20,col=ifelse(-log10(gwasres3$TN2014.additive_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres3$TN2014.additive_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres3$TN2014.additive_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres3$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres3$TN2014.additive_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres3$TN2014.additive_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres3$TN2014.additive_P_value_FDR)>=-log10(0.20),1.5,0.75))))


text(x=0.4,y=-0.05,font=2,labels="2014",cex=2,col="purple")
text(x=0.68,y=-0.05,font=2,labels="2015",cex=2,col="purple")
text(x=0.95,y=-0.05,font=2,labels="2016",cex=2,col="purple")

text(x=0,y=0,font=2,labels="Additive\nModel",cex=2.5)




### 3rd circlize : diplo-additive
library("circlize")
circos.clear()
par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
circos.par("track.height"=0.18,start.degree=90,gap.degree=4,points.overflow.warning=FALSE,cell.padding=c(0, 0, 0, 0),track.margin=c(0.005,0.005))
circos.initialize(factors = gwasres1$Chrom, x = gwasres1$Position)

# 1st track (external): 2016
circos.trackPlotRegion(factors=  gwasres1$Chrom, y=-log10(gwasres1$TN2016.diplo.additive_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2016.diplo.additive_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres1$Chrom, gwasres1$Position, -log10(gwasres1$TN2016.diplo.additive_P_value),pch=20,col=ifelse(-log10(gwasres1$TN2016.diplo.additive_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres1$TN2016.diplo.additive_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres1$TN2016.diplo.additive_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres1$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres1$TN2016.diplo.additive_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres1$TN2016.diplo.additive_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres1$TN2016.diplo.additive_P_value_FDR)>=-log10(0.20),1.5,0.75))))

# 2nd track : 2015
circos.trackPlotRegion(factors=  gwasres2$Chrom, y=-log10(gwasres2$TN2015.diplo.additive_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2015.diplo.additive_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres2$Chrom, gwasres2$Position, -log10(gwasres2$TN2015.diplo.additive_P_value),pch=20,col=ifelse(-log10(gwasres2$TN2015.diplo.additive_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres2$TN2015.diplo.additive_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres2$TN2015.diplo.additive_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres2$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres2$TN2015.diplo.additive_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres2$TN2015.diplo.additive_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres2$TN2015.diplo.additive_P_value_FDR)>=-log10(0.20),1.5,0.75))))


# 3rd track: 2014
circos.trackPlotRegion(factors=  gwasres3$Chrom, y=-log10(gwasres3$TN2014.diplo.additive_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2014.diplo.additive_P_value){
  circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.5)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), sector.index, facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
#circos.trackPoints(gwasres$Chrom, gwasres$Position, -log10(gwasres$mean_petals.diplo.additive_P_value),pch=20,col=ifelse(-log10(gwasres$mean_petals.diplo.additive_P_value_FDR)>=-log10(0.05),"red",ifelse(-log10(gwasres$mean_petals.diplo.additive_P_value_FDR)>=-log10(0.10),"orange",ifelse(-log10(gwasres$mean_petals.diplo.additive_P_value_FDR)>=-log10(0.20),"goldenrod",ifelse(gwasres$Chrom%%2==0,"navyblue","dodgerblue3")))))
circos.trackPoints(gwasres3$Chrom, gwasres3$Position, -log10(gwasres3$TN2014.diplo.additive_P_value),pch=20,col=ifelse(-log10(gwasres3$TN2014.diplo.additive_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres3$TN2014.diplo.additive_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres3$TN2014.diplo.additive_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres3$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres3$TN2014.diplo.additive_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres3$TN2014.diplo.additive_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres3$TN2014.diplo.additive_P_value_FDR)>=-log10(0.20),1.5,0.75))))

text(x=0.4,y=-0.05,font=2,labels="2014",cex=2,col="purple")
text(x=0.68,y=-0.05,font=2,labels="2015",cex=2,col="purple")
text(x=0.95,y=-0.05,font=2,labels="2016",cex=2,col="purple")

text(x=0,y=0,font=2,labels="Diploidized\nAdditive\nModel",cex=2.5)


### 5th circlize : 1.dom.ref
library("circlize")
circos.clear()
par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
circos.par("track.height"=0.18,start.degree=90,gap.degree=4,points.overflow.warning=FALSE,cell.padding=c(0, 0, 0, 0),track.margin=c(0.005,0.005))
circos.initialize(factors = gwasres1$Chrom, x = gwasres1$Position)

# 1st track (external): 2016
circos.trackPlotRegion(factors=  gwasres1$Chrom, y=-log10(gwasres1$TN2016.1.dom.ref_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2016.1.dom.ref_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres1$Chrom, gwasres1$Position, -log10(gwasres1$TN2016.1.dom.ref_P_value),pch=20,col=ifelse(-log10(gwasres1$TN2016.1.dom.ref_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres1$TN2016.1.dom.ref_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres1$TN2016.1.dom.ref_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres1$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres1$TN2016.1.dom.ref_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres1$TN2016.1.dom.ref_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres1$TN2016.1.dom.ref_P_value_FDR)>=-log10(0.20),1.5,0.75))))

# 2nd track : 2015
circos.trackPlotRegion(factors=  gwasres2$Chrom, y=-log10(gwasres2$TN2015.1.dom.ref_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2015.1.dom.ref_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres2$Chrom, gwasres2$Position, -log10(gwasres2$TN2015.1.dom.ref_P_value),pch=20,col=ifelse(-log10(gwasres2$TN2015.1.dom.ref_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres2$TN2015.1.dom.ref_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres2$TN2015.1.dom.ref_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres2$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres2$TN2015.1.dom.ref_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres2$TN2015.1.dom.ref_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres2$TN2015.1.dom.ref_P_value_FDR)>=-log10(0.20),1.5,0.75))))


# 3rd track: 2014
circos.trackPlotRegion(factors=  gwasres3$Chrom, y=-log10(gwasres3$TN2014.1.dom.ref_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2014.1.dom.ref_P_value){
  circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.5)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), sector.index, facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
#circos.trackPoints(gwasres$Chrom, gwasres$Position, -log10(gwasres$mean_petals.1.dom.ref_P_value),pch=20,col=ifelse(-log10(gwasres$mean_petals.1.dom.ref_P_value_FDR)>=-log10(0.05),"red",ifelse(-log10(gwasres$mean_petals.1.dom.ref_P_value_FDR)>=-log10(0.10),"orange",ifelse(-log10(gwasres$mean_petals.1.dom.ref_P_value_FDR)>=-log10(0.20),"goldenrod",ifelse(gwasres$Chrom%%2==0,"navyblue","dodgerblue3")))))
circos.trackPoints(gwasres3$Chrom, gwasres3$Position, -log10(gwasres3$TN2014.1.dom.ref_P_value),pch=20,col=ifelse(-log10(gwasres3$TN2014.1.dom.ref_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres3$TN2014.1.dom.ref_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres3$TN2014.1.dom.ref_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres3$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres3$TN2014.1.dom.ref_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres3$TN2014.1.dom.ref_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres3$TN2014.1.dom.ref_P_value_FDR)>=-log10(0.20),1.5,0.75))))

text(x=0.4,y=-0.05,font=2,labels="2014",cex=2,col="purple")
text(x=0.68,y=-0.05,font=2,labels="2015",cex=2,col="purple")
text(x=0.95,y=-0.05,font=2,labels="2016",cex=2,col="purple")

text(x=0,y=0,font=2,labels="Simplex\nDominant\nRef\nModel",cex=2.5)


### 6th circlize : 1.dom.alt
library("circlize")
circos.clear()
par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
circos.par("track.height"=0.18,start.degree=90,gap.degree=4,points.overflow.warning=FALSE,cell.padding=c(0, 0, 0, 0),track.margin=c(0.005,0.005))
circos.initialize(factors = gwasres1$Chrom, x = gwasres1$Position)

# 1st track (external): 2016
circos.trackPlotRegion(factors=  gwasres1$Chrom, y=-log10(gwasres1$TN2016.1.dom.alt_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2016.1.dom.alt_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres1$Chrom, gwasres1$Position, -log10(gwasres1$TN2016.1.dom.alt_P_value),pch=20,col=ifelse(-log10(gwasres1$TN2016.1.dom.alt_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres1$TN2016.1.dom.alt_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres1$TN2016.1.dom.alt_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres1$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres1$TN2016.1.dom.alt_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres1$TN2016.1.dom.alt_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres1$TN2016.1.dom.alt_P_value_FDR)>=-log10(0.20),1.5,0.75))))

# 2nd track : 2015
circos.trackPlotRegion(factors=  gwasres2$Chrom, y=-log10(gwasres2$TN2015.1.dom.alt_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2015.1.dom.alt_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres2$Chrom, gwasres2$Position, -log10(gwasres2$TN2015.1.dom.alt_P_value),pch=20,col=ifelse(-log10(gwasres2$TN2015.1.dom.alt_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres2$TN2015.1.dom.alt_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres2$TN2015.1.dom.alt_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres2$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres2$TN2015.1.dom.alt_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres2$TN2015.1.dom.alt_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres2$TN2015.1.dom.alt_P_value_FDR)>=-log10(0.20),1.5,0.75))))


# 3rd track: 2014
circos.trackPlotRegion(factors=  gwasres3$Chrom, y=-log10(gwasres3$TN2014.1.dom.alt_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2014.1.dom.alt_P_value){
  circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.5)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), sector.index, facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
#circos.trackPoints(gwasres$Chrom, gwasres$Position, -log10(gwasres$mean_petals.1.dom.alt_P_value),pch=20,col=ifelse(-log10(gwasres$mean_petals.1.dom.alt_P_value_FDR)>=-log10(0.05),"red",ifelse(-log10(gwasres$mean_petals.1.dom.alt_P_value_FDR)>=-log10(0.10),"orange",ifelse(-log10(gwasres$mean_petals.1.dom.alt_P_value_FDR)>=-log10(0.20),"goldenrod",ifelse(gwasres$Chrom%%2==0,"navyblue","dodgerblue3")))))
circos.trackPoints(gwasres3$Chrom, gwasres3$Position, -log10(gwasres3$TN2014.1.dom.alt_P_value),pch=20,col=ifelse(-log10(gwasres3$TN2014.1.dom.alt_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres3$TN2014.1.dom.alt_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres3$TN2014.1.dom.alt_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres3$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres3$TN2014.1.dom.alt_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres3$TN2014.1.dom.alt_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres3$TN2014.1.dom.alt_P_value_FDR)>=-log10(0.20),1.5,0.75))))

text(x=0.4,y=-0.05,font=2,labels="2014",cex=2,col="purple")
text(x=0.68,y=-0.05,font=2,labels="2015",cex=2,col="purple")
text(x=0.95,y=-0.05,font=2,labels="2016",cex=2,col="purple")

text(x=0,y=0,font=2,labels="Simplex\nDominant\nAlt\nModel",cex=2.5)



### 7th circlize : 2.dom.ref
library("circlize")
circos.clear()
par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
circos.par("track.height"=0.18,start.degree=90,gap.degree=4,points.overflow.warning=FALSE,cell.padding=c(0, 0, 0, 0),track.margin=c(0.005,0.005))
circos.initialize(factors = gwasres1$Chrom, x = gwasres1$Position)

# 1st track (external): 2016
circos.trackPlotRegion(factors=  gwasres1$Chrom, y=-log10(gwasres1$TN2016.2.dom.ref_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2016.2.dom.ref_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres1$Chrom, gwasres1$Position, -log10(gwasres1$TN2016.2.dom.ref_P_value),pch=20,col=ifelse(-log10(gwasres1$TN2016.2.dom.ref_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres1$TN2016.2.dom.ref_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres1$TN2016.2.dom.ref_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres1$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres1$TN2016.2.dom.ref_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres1$TN2016.2.dom.ref_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres1$TN2016.2.dom.ref_P_value_FDR)>=-log10(0.20),1.5,0.75))))

# 2nd track : 2015
circos.trackPlotRegion(factors=  gwasres2$Chrom, y=-log10(gwasres2$TN2015.2.dom.ref_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2015.2.dom.ref_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres2$Chrom, gwasres2$Position, -log10(gwasres2$TN2015.2.dom.ref_P_value),pch=20,col=ifelse(-log10(gwasres2$TN2015.2.dom.ref_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres2$TN2015.2.dom.ref_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres2$TN2015.2.dom.ref_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres2$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres2$TN2015.2.dom.ref_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres2$TN2015.2.dom.ref_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres2$TN2015.2.dom.ref_P_value_FDR)>=-log10(0.20),1.5,0.75))))


# 3rd track: 2014
circos.trackPlotRegion(factors=  gwasres3$Chrom, y=-log10(gwasres3$TN2014.2.dom.ref_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2014.2.dom.ref_P_value){
  circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.5)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), sector.index, facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
#circos.trackPoints(gwasres$Chrom, gwasres$Position, -log10(gwasres$mean_petals.2.dom.ref_P_value),pch=20,col=ifelse(-log10(gwasres$mean_petals.2.dom.ref_P_value_FDR)>=-log10(0.05),"red",ifelse(-log10(gwasres$mean_petals.2.dom.ref_P_value_FDR)>=-log10(0.10),"orange",ifelse(-log10(gwasres$mean_petals.2.dom.ref_P_value_FDR)>=-log10(0.20),"goldenrod",ifelse(gwasres$Chrom%%2==0,"navyblue","dodgerblue3")))))
circos.trackPoints(gwasres3$Chrom, gwasres3$Position, -log10(gwasres3$TN2014.2.dom.ref_P_value),pch=20,col=ifelse(-log10(gwasres3$TN2014.2.dom.ref_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres3$TN2014.2.dom.ref_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres3$TN2014.2.dom.ref_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres3$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres3$TN2014.2.dom.ref_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres3$TN2014.2.dom.ref_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres3$TN2014.2.dom.ref_P_value_FDR)>=-log10(0.20),1.5,0.75))))

text(x=0.4,y=-0.05,font=2,labels="2014",cex=2,col="purple")
text(x=0.68,y=-0.05,font=2,labels="2015",cex=2,col="purple")
text(x=0.95,y=-0.05,font=2,labels="2016",cex=2,col="purple")

text(x=0,y=0,font=2,labels="Duplex\nDominant\nRef\nModel",cex=2.5)


### 8th circlize : 2.dom.alt
library("circlize")
circos.clear()
par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
circos.par("track.height"=0.18,start.degree=90,gap.degree=4,points.overflow.warning=FALSE,cell.padding=c(0, 0, 0, 0),track.margin=c(0.005,0.005))
circos.initialize(factors = gwasres1$Chrom, x = gwasres1$Position)

# 1st track (external): 2016
circos.trackPlotRegion(factors=  gwasres1$Chrom, y=-log10(gwasres1$TN2016.2.dom.alt_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2016.2.dom.alt_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres1$Chrom, gwasres1$Position, -log10(gwasres1$TN2016.2.dom.alt_P_value),pch=20,col=ifelse(-log10(gwasres1$TN2016.2.dom.alt_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres1$TN2016.2.dom.alt_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres1$TN2016.2.dom.alt_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres1$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres1$TN2016.2.dom.alt_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres1$TN2016.2.dom.alt_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres1$TN2016.2.dom.alt_P_value_FDR)>=-log10(0.20),1.5,0.75))))

# 2nd track : 2015
circos.trackPlotRegion(factors=  gwasres2$Chrom, y=-log10(gwasres2$TN2015.2.dom.alt_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2015.2.dom.alt_P_value){
  #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
circos.trackPoints(gwasres2$Chrom, gwasres2$Position, -log10(gwasres2$TN2015.2.dom.alt_P_value),pch=20,col=ifelse(-log10(gwasres2$TN2015.2.dom.alt_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres2$TN2015.2.dom.alt_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres2$TN2015.2.dom.alt_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres2$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres2$TN2015.2.dom.alt_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres2$TN2015.2.dom.alt_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres2$TN2015.2.dom.alt_P_value_FDR)>=-log10(0.20),1.5,0.75))))


# 3rd track: 2014
circos.trackPlotRegion(factors=  gwasres3$Chrom, y=-log10(gwasres3$TN2014.2.dom.alt_P_value), track.height=0.2, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,6.1), panel.fun=function(Position,TN2014.2.dom.alt_P_value){
  circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.5)
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), sector.index, facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.8)
})
#circos.trackPoints(gwasres$Chrom, gwasres$Position, -log10(gwasres$mean_petals.2.dom.alt_P_value),pch=20,col=ifelse(-log10(gwasres$mean_petals.2.dom.alt_P_value_FDR)>=-log10(0.05),"red",ifelse(-log10(gwasres$mean_petals.2.dom.alt_P_value_FDR)>=-log10(0.10),"orange",ifelse(-log10(gwasres$mean_petals.2.dom.alt_P_value_FDR)>=-log10(0.20),"goldenrod",ifelse(gwasres$Chrom%%2==0,"navyblue","dodgerblue3")))))
circos.trackPoints(gwasres3$Chrom, gwasres3$Position, -log10(gwasres3$TN2014.2.dom.alt_P_value),pch=20,col=ifelse(-log10(gwasres3$TN2014.2.dom.alt_P_value_FDR)>=-log10(0.05),"firebrick",ifelse(-log10(gwasres3$TN2014.2.dom.alt_P_value_FDR)>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres3$TN2014.2.dom.alt_P_value_FDR)>=-log10(0.20),"goldenrod2",ifelse(gwasres3$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres3$TN2014.2.dom.alt_P_value_FDR)>=-log10(0.05),2.5,ifelse(-log10(gwasres3$TN2014.2.dom.alt_P_value_FDR)>=-log10(0.10),2,ifelse(-log10(gwasres3$TN2014.2.dom.alt_P_value_FDR)>=-log10(0.20),1.5,0.75))))

text(x=0.4,y=-0.05,font=2,labels="2014",cex=2,col="purple")
text(x=0.68,y=-0.05,font=2,labels="2015",cex=2,col="purple")
text(x=0.95,y=-0.05,font=2,labels="2016",cex=2,col="purple")

text(x=0,y=0,font=2,labels="Duplex\nDominant\nAlt\nModel",cex=2.5)


