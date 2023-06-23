
# TL - 070523
library(ggplot2)
library(grid)
library(gridExtra)
# generate estimates of the local ancestry using diagnostic markers
setwd("/home/thibaultleroy/Rose/Papier/Fig3/")
diagnostic=read.table("rose_round2.jointvcf.clean.AllFreqRefPerGroup.diagnostics.txt",h=TRUE)


### for supplementary figure diagnostic see below !

### circlize fig3
div=read.table("rose_round2.jointvcf_nucleotidediversity.AllChr.stats.ROD.withoutchr0",h=TRUE)
localancestry=read.table("localancestryestimates_diagnostic_300322.withoutchr0.withoutNA_190623.txt",h=TRUE)
# first detect windows of interest
# 5% lower pi excepted on chr0
# 5% top negative Tajima's D
# 5% highest RoD compared for Hybrid Tea as compared to European

# compute quantiles
quantile(div$pi_persite_HybTea,probs = c(0.01,0.05,0.10,0.5),na.rm = TRUE)
quantile(div$D_Pop_HybTea,probs = c(0.01,0.05,0.10,0.5),na.rm = TRUE)
quantile(div$RoD_pi_HybTeavsEurope,probs = c(0.5,0.9,0.95,0.99),na.rm = TRUE)

# identify windows
divfootprints5pc <- subset(div,(div$pi_persite_HybTea < quantile(div$pi_persite_HybTea,probs =0.05,na.rm = TRUE))&(div$D_Pop_HybTea < quantile(div$D_Pop_HybTea,probs =0.05,na.rm = TRUE))&((div$RoD_pi_HybTeavsEurope > quantile(div$RoD_pi_HybTeavsEurope,probs =0.95,na.rm = TRUE))|(div$RoD_pi_HybTeavsAsia > quantile(div$RoD_pi_HybTeavsAsia,probs =0.99,na.rm = TRUE))))
nrow(divfootprints5pc) # 99 windows  # 98 only considering in the direction of Europeans
divfootprints1pc <- subset(div,(div$pi_persite_HybTea < quantile(div$pi_persite_HybTea,probs =0.01,na.rm = TRUE))&(div$D_Pop_HybTea < quantile(div$D_Pop_HybTea,probs =0.01,na.rm = TRUE))&((div$RoD_pi_HybTeavsEurope > quantile(div$RoD_pi_HybTeavsEurope,probs =0.99,na.rm = TRUE))|(div$RoD_pi_HybTeavsAsia > quantile(div$RoD_pi_HybTeavsAsia,probs =0.99,na.rm = TRUE))))
nrow(divfootprints1pc) # 18 windows

write.table(divfootprints1pc,file = "windows_detected_footprintselection_top1pc_pi_D_ROD_070523.txt",row.names = FALSE,quote = FALSE,col.names = TRUE,sep = "\t")
write.table(divfootprints5pc,file = "windows_detected_footprintselection_top5pc_pi_D_ROD_070523.txt",row.names = FALSE,quote = FALSE,col.names = TRUE,sep = "\t")
# a few edits/information parsed on excel
detectedwindows=read.table("windows_detected_footprintselection_top1or5pc_pi_D_ROD_070523.txt",h=TRUE)
# Fist of all, we should compute a loess function per chromosome, for readability
detectedwindowschr1=subset(detectedwindows,detectedwindows$chr=="Chr01")
detectedwindowschr2=subset(detectedwindows,detectedwindows$chr=="Chr02")
detectedwindowschr3=subset(detectedwindows,detectedwindows$chr=="Chr03")
detectedwindowschr4=subset(detectedwindows,detectedwindows$chr=="Chr04")
detectedwindowschr5=subset(detectedwindows,detectedwindows$chr=="Chr05")
detectedwindowschr6=subset(detectedwindows,detectedwindows$chr=="Chr06")
detectedwindowschr7=subset(detectedwindows,detectedwindows$chr=="Chr07")

# Loess Pi
PredictionPiEurope=NULL
PredictionPiEuropeTotal=NULL
PredictionPiAsia=NULL
PredictionPiAsiaTotal=NULL
PredictionPiHyb1st=NULL
PredictionPiHyb1stTotal=NULL
PredictionPiHybTea=NULL
PredictionPiHybTeaTotal=NULL
for (chromo in c(1:7)){
  chromo2=paste0("Chr0",chromo)
  divchr=subset(div,div$chr==chromo2)
  # Europe
  loess10Europe <- loess(pi_persite_Europe ~ startwin, data=divchr, span=.1)
  smooth10Europe <- predict(loess10Europe) 
  print(paste(chromo,length(smooth10Europe)))
  PredictionPiEurope=cbind(chromo2,divchr$startwin,smooth10Europe)
  print(head(PredictionPiEurope))
  print(tail(PredictionPiEurope))
  print(nrow(PredictionPiEurope))
  PredictionPiEuropeTotal=rbind(PredictionPiEuropeTotal,PredictionPiEurope)
  # Asia
  loess10Asia <- loess(pi_persite_Asia ~ startwin, data=divchr, span=.1)
  smooth10Asia <- predict(loess10Asia) 
  print(paste(chromo,length(smooth10Asia)))
  PredictionPiAsia=cbind(chromo2,divchr$startwin,smooth10Asia)
  print(head(PredictionPiAsia))
  print(tail(PredictionPiAsia))
  print(nrow(PredictionPiAsia))
  PredictionPiAsiaTotal=rbind(PredictionPiAsiaTotal,PredictionPiAsia)
  # Hyb1st
  loess10Hyb1st <- loess(pi_persite_HybFirst ~ startwin, data=divchr, span=.1)
  smooth10Hyb1st <- predict(loess10Hyb1st) 
  print(paste(chromo,length(smooth10Hyb1st)))
  PredictionPiHyb1st=cbind(chromo2,divchr$startwin,smooth10Hyb1st)
  print(head(PredictionPiHyb1st))
  print(tail(PredictionPiHyb1st))
  print(nrow(PredictionPiHyb1st))
  PredictionPiHyb1stTotal=rbind(PredictionPiHyb1stTotal,PredictionPiHyb1st)
  # HybTea
  loess10HybTea <- loess(pi_persite_HybTea ~ startwin, data=divchr, span=.1)
  smooth10HybTea <- predict(loess10HybTea) 
  print(paste(chromo,length(smooth10HybTea)))
  PredictionPiHybTea=cbind(chromo2,divchr$startwin,smooth10HybTea)
  print(head(PredictionPiHybTea))
  print(tail(PredictionPiHybTea))
  print(nrow(PredictionPiHybTea))
  PredictionPiHybTeaTotal=rbind(PredictionPiHybTeaTotal,PredictionPiHybTea)
} 
# Dataset Europe
PredictionPiEuropeTotal=as.data.frame(PredictionPiEuropeTotal)
colnames(PredictionPiEuropeTotal)<- c("chr","pos","pipersiteEurope")
PredictionPiEuropeTotal$chr <- as.factor(as.character(PredictionPiEuropeTotal$chr))
PredictionPiEuropeTotal$pos <- as.numeric(as.character(PredictionPiEuropeTotal$pos))
PredictionPiEuropeTotal$pipersiteEurope <- as.numeric(as.character(PredictionPiEuropeTotal$pipersiteEurope))
# Dataset Asia
PredictionPiAsiaTotal=as.data.frame(PredictionPiAsiaTotal)
colnames(PredictionPiAsiaTotal)<- c("chr","pos","pipersiteAsia")
PredictionPiAsiaTotal$chr <- as.factor(as.character(PredictionPiAsiaTotal$chr))
PredictionPiAsiaTotal$pos <- as.numeric(as.character(PredictionPiEuropeTotal$pos))
PredictionPiAsiaTotal$pipersiteAsia <- as.numeric(as.character(PredictionPiAsiaTotal$pipersiteAsia))
# Dataset Hyb1st
PredictionPiHyb1stTotal=as.data.frame(PredictionPiHyb1stTotal)
colnames(PredictionPiHyb1stTotal)<- c("chr","pos","pipersiteHyb1st")
PredictionPiHyb1stTotal$chr <- as.factor(as.character(PredictionPiHyb1stTotal$chr))
PredictionPiHyb1stTotal$pos <- as.numeric(as.character(PredictionPiEuropeTotal$pos))
PredictionPiHyb1stTotal$pipersiteHyb1st <- as.numeric(as.character(PredictionPiHyb1stTotal$pipersiteHyb1st))
# Dataset HybTea
PredictionPiHybTeaTotal=as.data.frame(PredictionPiHybTeaTotal)
colnames(PredictionPiHybTeaTotal)<- c("chr","pos","pipersiteHybTea")
PredictionPiHybTeaTotal$chr <- as.factor(as.character(PredictionPiHybTeaTotal$chr))
PredictionPiHybTeaTotal$pos <- as.numeric(as.character(PredictionPiEuropeTotal$pos))
PredictionPiHybTeaTotal$pipersiteHybTea <- as.numeric(as.character(PredictionPiHybTeaTotal$pipersiteHybTea))



### Loess Tajima's D 
# no windows with NA
PredictionDEurope=NULL
PredictionDEuropeTotal=NULL
PredictionDAsia=NULL
PredictionDAsiaTotal=NULL
PredictionDHyb1st=NULL
PredictionDHyb1stTotal=NULL
PredictionDHybTea=NULL
PredictionDHybTeaTotal=NULL
for (chromo in c(1:7)){
  chromo2=paste0("Chr0",chromo)
  divchr=subset(div,div$chr==chromo2)
  # Europe
  loess10Europe <- loess(D_Pop_Europe ~ startwin, data=divchr, span=.1)
  smooth10Europe <- predict(loess10Europe) 
  print(paste(chromo,length(smooth10Europe)))
  PredictionDEurope=cbind(chromo2,divchr$startwin,smooth10Europe)
  print(head(PredictionDEurope))
  print(tail(PredictionDEurope))
  print(nrow(PredictionDEurope))
  PredictionDEuropeTotal=rbind(PredictionDEuropeTotal,PredictionDEurope)
  # Asia
  loess10Asia <- loess(D_Pop_Asia ~ startwin, data=divchr, span=.1)
  smooth10Asia <- predict(loess10Asia) 
  print(paste(chromo,length(smooth10Asia)))
  PredictionDAsia=cbind(chromo2,divchr$startwin,smooth10Asia)
  print(head(PredictionDAsia))
  print(tail(PredictionDAsia))
  print(nrow(PredictionDAsia))
  PredictionDAsiaTotal=rbind(PredictionDAsiaTotal,PredictionDAsia)
  # Hyb1st
  loess10Hyb1st <- loess(D_Pop_HybFirst ~ startwin, data=divchr, span=.1)
  smooth10Hyb1st <- predict(loess10Hyb1st) 
  print(paste(chromo,length(smooth10Hyb1st)))
  PredictionDHyb1st=cbind(chromo2,divchr$startwin,smooth10Hyb1st)
  print(head(PredictionDHyb1st))
  print(tail(PredictionDHyb1st))
  print(nrow(PredictionDHyb1st))
  PredictionDHyb1stTotal=rbind(PredictionDHyb1stTotal,PredictionDHyb1st)
  # HybTea
  loess10HybTea <- loess(D_Pop_HybTea ~ startwin, data=divchr, span=.1)
  smooth10HybTea <- predict(loess10HybTea) 
  print(paste(chromo,length(smooth10HybTea)))
  PredictionDHybTea=cbind(chromo2,divchr$startwin,smooth10HybTea)
  print(head(PredictionDHybTea))
  print(tail(PredictionDHybTea))
  print(nrow(PredictionDHybTea))
  PredictionDHybTeaTotal=rbind(PredictionDHybTeaTotal,PredictionDHybTea)
} 
# Dataset Europe
PredictionDEuropeTotal=as.data.frame(PredictionDEuropeTotal)
colnames(PredictionDEuropeTotal)<- c("chr","pos","DpopEurope")
PredictionDEuropeTotal$chr <- as.factor(as.character(PredictionDEuropeTotal$chr))
PredictionDEuropeTotal$pos <- as.numeric(as.character(PredictionDEuropeTotal$pos))
PredictionDEuropeTotal$DpopEurope <- as.numeric(as.character(PredictionDEuropeTotal$DpopEurope))
# Dataset Asia
PredictionDAsiaTotal=as.data.frame(PredictionDAsiaTotal)
colnames(PredictionDAsiaTotal)<- c("chr","pos","DpopAsia")
PredictionDAsiaTotal$chr <- as.factor(as.character(PredictionDAsiaTotal$chr))
PredictionDAsiaTotal$pos <- as.numeric(as.character(PredictionDEuropeTotal$pos))
PredictionDAsiaTotal$DpopAsia <- as.numeric(as.character(PredictionDAsiaTotal$DpopAsia))
# Dataset Hyb1st
PredictionDHyb1stTotal=as.data.frame(PredictionDHyb1stTotal)
colnames(PredictionDHyb1stTotal)<- c("chr","pos","DpopHyb1st")
PredictionDHyb1stTotal$chr <- as.factor(as.character(PredictionDHyb1stTotal$chr))
PredictionDHyb1stTotal$pos <- as.numeric(as.character(PredictionDEuropeTotal$pos))
PredictionDHyb1stTotal$DpopHyb1st <- as.numeric(as.character(PredictionDHyb1stTotal$DpopHyb1st))
# Dataset HybTea
PredictionDHybTeaTotal=as.data.frame(PredictionDHybTeaTotal)
colnames(PredictionDHybTeaTotal)<- c("chr","pos","DpopHybTea")
PredictionDHybTeaTotal$chr <- as.factor(as.character(PredictionDHybTeaTotal$chr))
PredictionDHybTeaTotal$pos <- as.numeric(as.character(PredictionDEuropeTotal$pos))
PredictionDHybTeaTotal$DpopHybTea <- as.numeric(as.character(PredictionDHybTeaTotal$DpopHybTea))

library(circlize)
#pdf(file = "Fig3_circlize_footprints_detectedwindows.pdf",height = 16,width = 16)
#png(file = "Fig3_circlize_footprints_detectedwindows.png",height = 800,width = 800)
circos.clear()
par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
circos.par("track.height"=0.12,start.degree=90,gap.degree=4,points.overflow.warning=FALSE,cell.padding=c(0, 0, 0, 0),track.margin=c(0.005,0.005))
circos.initialize(factors = div$chr, x = div$startwin)

# 1st track density pi per site
circos.trackPlotRegion(factors= div$chr, y=div$pi_persite_Europe, track.height=0.21, track.margin=c(0.01,0.01), bg.border = NA, panel.fun=function(startwin,pi_persite_Europe){
  #circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), max(ylim), paste(sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.6)
})

# add detected windows
for (chromo in c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07")){
  detectedwindowschr=subset(detectedwindows,detectedwindows$chr==chromo)
  circos.segments(sector.index = detectedwindowschr$chr,x0=detectedwindowschr$startwin,x1=detectedwindowschr$startwin,y0=rep(0,nrow(detectedwindowschr)),y1=rep(1,nrow(detectedwindowschr)),col="darkred")
}
circos.segments(sector.index = "Chr01",x0=detectedwindowschr1$startwin,x1=detectedwindowschr1$startwin,y0=rep(-0.085,nrow(detectedwindowschr1)),y1=rep(max(div$pi_persite_Europe),nrow(detectedwindowschr1))*0.8,col=ifelse(detectedwindowschr1$type=="top1pc",add_transparency("darkviolet",0.5),add_transparency("salmon",0.5)),lwd = ifelse(detectedwindowschr1$type=="top1pc",0.6,0.4))
circos.segments(sector.index = "Chr02",x0=detectedwindowschr2$startwin,x1=detectedwindowschr2$startwin,y0=rep(-0.085,nrow(detectedwindowschr2)),y1=rep(max(div$pi_persite_Europe),nrow(detectedwindowschr2))*0.8,col=ifelse(detectedwindowschr2$type=="top1pc",add_transparency("darkviolet",0.5),add_transparency("salmon",0.5)),lwd = ifelse(detectedwindowschr2$type=="top1pc",0.6,0.4))
circos.segments(sector.index = "Chr03",x0=detectedwindowschr3$startwin,x1=detectedwindowschr3$startwin,y0=rep(-0.085,nrow(detectedwindowschr3)),y1=rep(max(div$pi_persite_Europe),nrow(detectedwindowschr3))*0.8,col=ifelse(detectedwindowschr3$type=="top1pc",add_transparency("darkviolet",0.5),add_transparency("salmon",0.5)),lwd = ifelse(detectedwindowschr3$type=="top1pc",0.6,0.4))
circos.segments(sector.index = "Chr04",x0=detectedwindowschr4$startwin,x1=detectedwindowschr4$startwin,y0=rep(-0.085,nrow(detectedwindowschr4)),y1=rep(max(div$pi_persite_Europe),nrow(detectedwindowschr4))*0.8,col=ifelse(detectedwindowschr4$type=="top1pc",add_transparency("darkviolet",0.5),add_transparency("salmon",0.5)),lwd = ifelse(detectedwindowschr4$type=="top1pc",0.6,0.4))
circos.segments(sector.index = "Chr05",x0=detectedwindowschr5$startwin,x1=detectedwindowschr5$startwin,y0=rep(-0.085,nrow(detectedwindowschr5)),y1=rep(max(div$pi_persite_Europe),nrow(detectedwindowschr5))*0.8,col=ifelse(detectedwindowschr5$type=="top1pc",add_transparency("darkviolet",0.5),add_transparency("salmon",0.5)),lwd = ifelse(detectedwindowschr5$type=="top1pc",0.6,0.4))
circos.segments(sector.index = "Chr06",x0=detectedwindowschr6$startwin,x1=detectedwindowschr6$startwin,y0=rep(-0.085,nrow(detectedwindowschr6)),y1=rep(max(div$pi_persite_Europe),nrow(detectedwindowschr6))*0.8,col=ifelse(detectedwindowschr6$type=="top1pc",add_transparency("darkviolet",0.5),add_transparency("salmon",0.5)),lwd = ifelse(detectedwindowschr6$type=="top1pc",0.6,0.4))
circos.segments(sector.index = "Chr07",x0=detectedwindowschr7$startwin,x1=detectedwindowschr7$startwin,y0=rep(-0.085,nrow(detectedwindowschr7)),y1=rep(max(div$pi_persite_Europe),nrow(detectedwindowschr7))*0.8,col=ifelse(detectedwindowschr7$type=="top1pc",add_transparency("darkviolet",0.5),add_transparency("salmon",0.5)),lwd = ifelse(detectedwindowschr7$type=="top1pc",0.6,0.4))


circos.trackPoints(div$chr, div$startwin, div$pi_persite_Europe, col=add_transparency("#000080",0.3),pch=20,cex=0.3)
circos.trackPoints(div$chr, div$startwin, div$pi_persite_Asia, col=add_transparency("#EEB422",0.3),pch=20,cex=0.3)
circos.trackPoints(div$chr, div$startwin, div$pi_persite_HybFirst, col=add_transparency("#228B22",0.3),pch=20,cex=0.3)
circos.trackPoints(div$chr, div$startwin, div$pi_persite_HybTea, col=add_transparency("#889F22",0.3),pch=20,cex=0.3)


circos.trackLines(PredictionPiAsiaTotal$chr,PredictionPiAsiaTotal$pos,PredictionPiAsiaTotal$pipersiteAsia,col="#EEB422",lwd=2,type="l")
circos.trackLines(PredictionPiEuropeTotal$chr,PredictionPiEuropeTotal$pos,PredictionPiEuropeTotal$pipersiteEurope,col="#000080",lwd=2,type="l")
circos.trackLines(PredictionPiHyb1stTotal$chr,PredictionPiHyb1stTotal$pos,PredictionPiHyb1stTotal$pipersiteHyb1st,col="#228B22",lwd=2,type="l")
circos.trackLines(PredictionPiHybTeaTotal$chr,PredictionPiHybTeaTotal$pos,PredictionPiHybTeaTotal$pipersiteHybTea,col="#889F22",lwd=2.5,type="l")


# 2nd track density Tajima's D
circos.trackPlotRegion(factors= div$chr, y=div$D_Pop_Asia, track.height=0.21, track.margin=c(0.01,0.01), bg.border = NA, panel.fun=function(startwin,D_Pop_Asia){
  #circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), paste(sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.6)
})

circos.trackPoints(div$chr, div$startwin, div$D_Pop_Europe, col=add_transparency("#000080",0.3),pch=20,cex=0.3)
circos.trackPoints(div$chr, div$startwin, div$D_Pop_Asia, col=add_transparency("#EEB422",0.3),pch=20,cex=0.3)
circos.trackPoints(div$chr, div$startwin, div$D_Pop_HybFirst, col=add_transparency("#228B22",0.3),pch=20,cex=0.3)
circos.trackPoints(div$chr, div$startwin, div$D_Pop_HybTea, col=add_transparency("#889F22",0.3),pch=20,cex=0.3)


circos.trackLines(PredictionDAsiaTotal$chr,PredictionDAsiaTotal$pos,PredictionDAsiaTotal$DpopAsia,col="#EEB422",lwd=2,type="l")
circos.trackLines(PredictionDEuropeTotal$chr,PredictionDEuropeTotal$pos,PredictionDEuropeTotal$DpopEurope,col="#000080",lwd=2,type="l")
circos.trackLines(PredictionDHyb1stTotal$chr,PredictionDHyb1stTotal$pos,PredictionDHyb1stTotal$DpopHyb1st,col="#228B22",lwd=2,type="l")
circos.trackLines(PredictionDHybTeaTotal$chr,PredictionDHybTeaTotal$pos,PredictionDHybTeaTotal$DpopHybTea,col="#889F22",lwd=2.5,type="l")

# 3rd track density RoD
circos.trackPlotRegion(factors= div$chr, y=div$RoD_pi_HybTeavsEurope, track.height=0.15, track.margin=c(0.01,0.01), bg.border = NA, panel.fun=function(startwin,RoD_pi_HybTeavsEurope){
  #circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), paste(sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  circos.yaxis("left", lwd=1,labels.cex=0.6)
})
circos.trackLines(div$chr,div$startwin,div$RoD_pi_HybTeavsEurope,baseline=0,col=ifelse(div$RoD_pi_HybTeavsEurope<0,"cornflowerblue",ifelse(div$RoD_pi_HybTeavsEurope>quantile(div$RoD_pi_HybTeavsEurope,probs =0.99,na.rm = TRUE),"darkred",ifelse(div$RoD_pi_HybTeavsEurope>quantile(div$RoD_pi_HybTeavsEurope,probs =0.95,na.rm = TRUE),"orange","#889F22"))),lwd=ifelse(div$RoD_pi_HybTeavsEurope<0,0.7,ifelse(div$RoD_pi_HybTeavsEurope>quantile(div$RoD_pi_HybTeavsEurope,probs =0.99,na.rm = TRUE),1.5,ifelse(div$RoD_pi_HybTeavsEurope>quantile(div$RoD_pi_HybTeavsEurope,probs =0.95,na.rm = TRUE),0.7,0.4))),type = "h")

# 4th local ancestry
circos.trackPlotRegion(factors= localancestry$chr, y=localancestry$medianfreqhybtea, track.height=0.15, track.margin=c(0.01,0.01), bg.border = NA, panel.fun=function(startwin,medianfreqhybtea){
  circos.axis(h=-0.2, major.at=c(10000,25000000,50000000,75000000), labels=c(0,25,50,75), major.tick=TRUE, minor.ticks = 4, lwd=2,labels.cex=0.6,col="black",direction= "inside")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  #circos.text(mean(xlim), max(ylim), paste(sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
  #circos.yaxis("left", lwd=1,labels.cex=0.6)
})
#circos.trackLines(localancestry$chr,localancestry$startwin,(1-localancestry$medianfreqhybtea),baseline=1,col="#000080",type = "h")
#circos.trackLines(localancestry$chr,localancestry$startwin,localancestry$medianfreqhybtea,col="#EEB422",type = "h")



for (chromolocal in c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07")){
  localancestrychr=subset(localancestry,localancestry$chr==chromolocal) #  & is.na(localancestry)==FALSE
  localancestrychr=localancestrychr[order(localancestrychr$startwin),]
  #circos.rect(sector.index = chromolocal ,xleft=localancestrychr$startwin[1],xright=localancestrychr$endwin[nrow(localancestrychr)],ybottom=0,ytop=1,color="black")
  circos.segments(sector.index = chromolocal ,x0=localancestrychr$startwin,x1=localancestrychr$endwin,y0=rep(0,nrow(localancestrychr)),y1=(1-localancestrychr$medianfreqhybtea),col="#000080")
  circos.segments(sector.index = chromolocal ,x0=localancestrychr$startwin,x1=localancestrychr$endwin,y0=rep(1,nrow(localancestrychr)),y1=(1-localancestrychr$medianfreqhybtea),col="#EEB422")
  #for (line in nrow(localancestrychr)){
  #  circos.rect(xleft=localancestrychr$startwin[line],ybottom=0,xright=localancestrychr$endwin[line],ytop=(1-localancestrychr$medianfreqhybtea[line]),sector.index=localancestrychr$chr[line],track.index = 4, col="#000080")
  #  #circos.rect(localancestry$chr,localancestry$startwin,localancestry$medianfreqhybtea,col="#EEB422",type = "h")
  #}
}

#dev.off()

### supplementary figure diagnostic 

# plot boxplot
dataset=cbind("Asia",diagnostic$chr,diagnostic$pos,diagnostic$freqasia)
dataset=rbind(dataset,cbind("Europe",diagnostic$chr,diagnostic$pos,diagnostic$freqeuro))
dataset=rbind(dataset,cbind("Hyb1st",diagnostic$chr,diagnostic$pos,diagnostic$freqfirsthyb))
dataset=rbind(dataset,cbind("HybTea",diagnostic$chr,diagnostic$pos,diagnostic$freqhybtea))
dataset=as.data.frame(dataset)
colnames(dataset)<- c("Group","chr","pos","freq")
dataset$freq=as.numeric(as.character(dataset$freq))


ggplot(dataset)+
  geom_hline(yintercept = 0.5,lty=2,color="red")+
  geom_hline(yintercept = 0.25,lty=2,color="orange")+
  geom_hline(yintercept = 0.75,lty=2,color="orange")+
  geom_violin(aes(x=Group,y=freq,colour=Group),draw_quantiles =0.5)+
  scale_color_manual(values=c("#EEB422","#000080","#228B22","#889F22"),labels=c("Ancient Asian","Ancient European","First European x Asian","Hybrid Tea"))+
  ylim(0,1)+
  xlab("Groups")+  ylab("Allele frequency at diagnostic markers\n(n=170,637 SNPs)\n")+
  scale_x_discrete(labels=c("Ancient Asian","Ancient European","First European x Asian","Hybrid Tea")) +
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


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



## Supplementary variation Tajima's D
dataset=cbind("Acient Asia",div$chr,div$startwin,div$pi_persite_Asia)
dataset=rbind(dataset,cbind("Ancient Europe",div$chr,div$startwin,div$pi_persite_Europe))
dataset=rbind(dataset,cbind("First Europe x Asia",div$chr,div$startwin,div$pi_persite_HybFirst))
dataset=rbind(dataset,cbind("Hybrid Tea",div$chr,div$startwin,div$pi_persite_HybTea))
dataset=as.data.frame(dataset,sep="\t")
colnames(dataset)<- c("Group","chr","pos","pi")
dataset$pi=as.numeric(as.character(dataset$pi))

ggplot(dataset)+
  geom_violin(aes(x=Group,y=pi,colour=Group),draw_quantiles =0.5)+
  scale_color_manual(values=c("#EEB422","#000080","#228B22","#889F22"),labels=c("Ancient Asian","Ancient European","First European x Asian","Hybrid Tea"))+
  #ylim(0,1)+
  xlab("Groups")+  ylab("Nucleotide diversity")+
  scale_x_discrete(labels=c("Ancient Asian","Ancient European","First European x Asian","Hybrid Tea")) +
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

dataset=cbind("Acient Asia",div$chr,div$startwin,div$D_Pop_Asia)
dataset=rbind(dataset,cbind("Ancient Europe",div$chr,div$startwin,div$D_Pop_Europe))
dataset=rbind(dataset,cbind("First Europe x Asia",div$chr,div$startwin,div$D_Pop_HybFirst))
dataset=rbind(dataset,cbind("Hybrid Tea",div$chr,div$startwin,div$D_Pop_HybTea))
dataset=as.data.frame(dataset,sep="\t")
colnames(dataset)<- c("Group","chr","pos","D")
dataset$D=as.numeric(as.character(dataset$D))

ggplot(dataset, aes(x=Group, y=D)) +
  geom_violin(aes(colour=Group),fill="transparent",draw_quantiles =0.5)+
  geom_boxplot(width=0.4,outlier.size=0,alpha=0.5)+
  scale_color_manual(values=c("#EEB422","#000080","#228B22","#889F22"),labels=c("Ancient Asian","Ancient European","First European x Asian","Hybrid Tea"))+
  #ylim(0,1)+
  xlab("Groups")+  ylab("Tajima's D")+
  scale_x_discrete(labels=c("Ancient Asian","Ancient European","First European x Asian","Hybrid Tea")) +
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

# Correlation matrix (Pearson)
data4cor <- div[, c(11,13,15,17,12,14,16,18)]
rescor <- cor(data4cor,method = "Pearson")
round(rescor, 3)
colnames(rescor)<- c("pi_Asia","pi_Europe","pi_First_Europe_x_Asia","pi_Hybrid_Tea","D_Asia","D_Europe","D_First_Europe_x_Asia","D_Hybrid_Tea")
rownames(rescor)<- c("pi_Asia","pi_Europe","pi_First_Europe_x_Asia","pi_Hybrid_Tea","D_Asia","D_Europe","D_First_Europe_x_Asia","D_Hybrid_Tea")
library(corrplot)

corrplot(rescor, method= "number",type = "upper", sig.level = 0.05, insig = "label_sig",
         col=colorRampPalette(c("navyblue", "lightgrey", "darkred"))(100), 
         tl.col = "black", tl.srt = 90)

#
corrplot(rescor, type="upper", method= "number",
         col=colorRampPalette(c("navyblue", "lightgrey", "darkred"))(100), 
         tl.col = "black", tl.srt = 90,number.cex=0.85)
pval <- psych::corr.test(data4cor, adjust="none")$p
