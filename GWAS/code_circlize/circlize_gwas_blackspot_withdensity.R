## TL - 021121

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/")
densitygene=read.table("densitygene_slidwin500kb.txt",h=T)

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/")
candidateRgene=read.table("list_candidateRgenes_Laurence.txt",h=F)

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/TN2014_gwas_results/")
#gwasres=read.table("TN2014_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="TN2014"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/TN2015_gwas_results/")
#gwasres=read.table("TN2015_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="TN2015"

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/TN2016_gwas_results/")
gwasres=read.table("TN2016_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait="TN2016"
#gwasres1=read.table("./TN2016_gwas_results/TN2016_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#gwasres2=read.table("./TN2015_gwas_results/TN2015_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#gwasres3=read.table("./TN2014_gwas_results/TN2014_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE


#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/mean_petals_binary_less10_gwas_results")
#gwasres=read.table("mean_petals_binary_less10_scores.txt.withoutchrNA.sed",h=T)
#trait="mean_petals_binary_less10"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/mean_petals_binary_more25_gwas_results")
#gwasres=read.table("mean_petals_binary_more25_scores.txt.withoutchrNA.sed",h=T)
#trait="mean_petals_binary_more25"

# for all models
for (model in c("additive","diplo.additive","general","diplo.general","1.dom.alt","1.dom.ref","2.dom.alt","2.dom.ref")){
  print(paste("Currently working on model: ",model))
  ### slidingwindows R
  ## fist loop by chromosome
  dataset <- data.frame(matrix(ncol = 12, nrow = 0))
  colnames(dataset)=c("chr","startwin","endwin","Nb_SNPs_Win","Prop_Win_Total","NbSNPWinPvalueInf05","NbSNPWinPvalueInf01","NbSNPWinPvalueInf001","SNPswithoutpvalues","PropSNPWinPvalueInf05","PropSNPWinPvalueInf01","PropSNPWinPvalueInf001")
  window=1000000
  for (chr in c(0:7)){
    print(paste("currently working on chr", chr))
    number_of_snp=nrow(gwasres)
    currentchr=subset(gwasres,gwasres$Chrom==chr)
    SNPcurrentchr=nrow(currentchr)
    print(paste("proportion of SNPs on chr",chr,"=",SNPcurrentchr/number_of_snp))
    ### 
    for (slidwin in c(0:floor(max(currentchr$Position)/window))){  # or 1:ceiling = but here 0 to have something easier with position
      start=slidwin*window+1
      end=slidwin*window+window
      currentwindow=subset(currentchr,currentchr$Position>=start & currentchr$Position<=end)
      SNPwindows=nrow(currentwindow)
      columntosubsample=paste0(trait,".",model,"_P_value")
      SNPwindowswithouttests=nrow(subset(currentwindow,is.na(currentwindow[[columntosubsample]])))
      SNPwindowswithtests=nrow(subset(currentwindow,!is.na(currentwindow[[columntosubsample]])))
      #SNPwindowswithouttests=nrow(subset(currentwindow,is.na(currentwindow$mean_petals.general_P_value)))
      #SNPwindowswithtests=nrow(subset(currentwindow,!is.na(currentwindow$mean_petals.general_P_value)))
      ## outliers
      #currentwindowoutlier05=subset(currentwindow,currentwindow$mean_petals.general_P_value <= 0.05)
      #currentwindowoutlier01=subset(currentwindow,currentwindow$mean_petals.general_P_value <= 0.01)
      #currentwindowoutlier001=subset(currentwindow,currentwindow$mean_petals.general_P_value <= 0.001)
      currentwindowoutlier05=subset(currentwindow,currentwindow[[columntosubsample]] <= 0.05)
      currentwindowoutlier01=subset(currentwindow,currentwindow[[columntosubsample]] <= 0.01)
      currentwindowoutlier001=subset(currentwindow,currentwindow[[columntosubsample]] <= 0.001)
      if (SNPwindowswithtests>0){
        propcurrentwindowoutlier05=nrow(currentwindowoutlier05)/SNPwindowswithtests
        propcurrentwindowoutlier01=nrow(currentwindowoutlier01)/SNPwindowswithtests
        propcurrentwindowoutlier001=nrow(currentwindowoutlier001)/SNPwindowswithtests
      }else {
        propcurrentwindowoutlier05=0
        propcurrentwindowoutlier01=0
        propcurrentwindowoutlier001=0
      }
      # print
      currentline=data.frame(cbind(chr,start,end,SNPwindows,SNPwindows/number_of_snp,nrow(currentwindowoutlier05),nrow(currentwindowoutlier01),nrow(currentwindowoutlier001),SNPwindowswithouttests,propcurrentwindowoutlier05,propcurrentwindowoutlier01,propcurrentwindowoutlier001))
      colnames(currentline)=c("chr","startwin","endwin","Nb_SNPs_Win","Prop_Win_Total","NbSNPWinPvalueInf05","NbSNPWinPvalueInf01","NbSNPWinPvalueInf001","SNPswithoutpvalues","PropSNPWinPvalueInf05","PropSNPWinPvalueInf01","PropSNPWinPvalueInf001")
      dataset=rbind(dataset,currentline)
    }
  }
  write.csv(dataset, file=paste0("res_densityassociated_",trait,"_",model,"_slidwin_",format(window,scientific=FALSE),".csv"),row.names=FALSE)

  #### Circlize with density curves # general model
  print(paste("Currently plotting results for model: ",model))
  columntoplot=paste0(trait,".",model,"_P_value")
  columntofilter=paste0(trait,".",model,"_P_value_FDR")
  
  pdf(file=paste0("circlize_",trait,"_",model,"slidwin",format(window,scientific=FALSE),".pdf"),width=12,height=12)
  library("circlize")
  circos.clear()
  par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
  circos.par("track.height"=0.18,start.degree=90,gap.degree=4,points.overflow.warning=FALSE,cell.padding=c(0, 0, 0, 0),track.margin=c(0.005,0.005))
  circos.initialize(factors = gwasres$Chrom, x = gwasres$Position)

  # 1st track density SNP
  circos.trackPlotRegion(factors= dataset$chr, y=dataset$Prop_Win_Total, track.height=0.08, track.margin=c(0.03,0.03), bg.border = NA, panel.fun=function(startwin,Prop_Win_Total){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.6)
  })

  # metaQTL Diana
  circos.rect(sector.index=3,xleft=21605249,ybottom=max(dataset$Prop_Win_Total)*-0.1,xright=24567362,ytop=max(dataset$Prop_Win_Total)*0.9,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=3,xleft=34220246,ybottom=max(dataset$Prop_Win_Total)*-0.1,xright=37772912,ytop=max(dataset$Prop_Win_Total)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=2414969,ybottom=max(dataset$Prop_Win_Total)*-0.1,xright=4219224,ytop=max(dataset$Prop_Win_Total)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=18827666,ybottom=max(dataset$Prop_Win_Total)*-0.1,xright=24889549,ytop=max(dataset$Prop_Win_Total)*1.1,lwd=0.01,lty=1,col="pink")
  
  for (gene in 1:nrow(candidateRgene)){
    #print(candidateRgene[gene,1])
    circos.segments(sector.index=as.factor(candidateRgene[gene,2]),x0=candidateRgene[gene,3],y0=0,x1=candidateRgene[gene,3],y1=max(dataset$Prop_Win_Total)*0.85,lwd=0.4,lty=1,col="darkgrey")
  }

  # density gene

  circos.trackLines(densitygene$chr, densitygene$startwin+(densitygene$endwin-densitygene$startwin), densitygene$nb_genes/20000, col="forestgreen", lwd=1.5, type="l")

 # density SNP
  circos.trackLines(dataset$chr, as.numeric(dataset$startwin)+window*0.5, dataset$Prop_Win_Total,lwd=2)



  # 2nd track = GWAS general Petals
  circos.trackPlotRegion(factors=  gwasres$Chrom, y=-log10(gwasres[[columntoplot]]), track.height=0.15, track.margin=c(0.03,0.03), bg.border = NA, ylim = c(0,max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)+0.3), panel.fun=function(Position,columntoplot){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.8)
  })

  # metaQTL Diana
  circos.rect(sector.index=3,xleft=21605249,ybottom=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*-0.1,xright=24567362,ytop=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=3,xleft=34220246,ybottom=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*-0.1,xright=37772912,ytop=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=2414969,ybottom=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*-0.1,xright=4219224,ytop=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=18827666,ybottom=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*-0.1,xright=24889549,ytop=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  for (gene in 1:nrow(candidateRgene)){
    #print(candidateRgene[gene,1])
    circos.segments(sector.index=as.factor(candidateRgene[gene,2]),x0=candidateRgene[gene,3],y0=0,x1=candidateRgene[gene,3],y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE),lwd=0.4,lty=1,col="darkgrey")
  }
  
  
  # values
  circos.trackPoints(gwasres$Chrom, gwasres$Position, -log10(gwasres[[columntoplot]]),pch=20,col=ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.05),"firebrick",ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.20),"goldenrod2",ifelse(gwasres$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.05),2.5,ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.10),2,ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.20),1.5,0.75))))

  # 3rd track = density significant p < 0.05

  circos.trackPlotRegion(factors= dataset$chr, y=dataset$PropSNPWinPvalueInf05, track.height=0.08, track.margin=c(0.03,0.03), bg.border = NA, panel.fun=function(startwin,PropSNPWinPvalueInf05){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.6)
  })
 
  # metaQTL Diana
  circos.rect(sector.index=3,xleft=21605249,ybottom=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE)*-0.1,xright=24567362,ytop=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=3,xleft=34220246,ybottom=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE)*-0.1,xright=37772912,ytop=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=2414969,ybottom=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE)*-0.1,xright=4219224,ytop=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=18827666,ybottom=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE)*-0.1,xright=24889549,ytop=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  for (gene in 1:nrow(candidateRgene)){
    #print(candidateRgene[gene,1])
    circos.segments(sector.index=as.factor(candidateRgene[gene,2]),x0=candidateRgene[gene,3],y0=0,x1=candidateRgene[gene,3],y1=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE),lwd=0.4,lty=1,col="darkgrey")
  }
  
  
  circos.trackLines(dataset$chr, dataset$startwin+window*0.5, dataset$PropSNPWinPvalueInf05,lwd=2,col="goldenrod2")


  # p < 0.01
  circos.trackPlotRegion(factors= dataset$chr, y=dataset$PropSNPWinPvalueInf01, track.height=0.08, track.margin=c(0.03,0.03), bg.border = NA, panel.fun=function(startwin,PropSNPWinPvalueInf01){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.6)
  })
  #
  # metaQTL Diana
  circos.rect(sector.index=3,xleft=21605249,ybottom=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE)*-0.1,xright=24567362,ytop=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=3,xleft=34220246,ybottom=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE)*-0.1,xright=37772912,ytop=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=2414969,ybottom=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE)*-0.1,xright=4219224,ytop=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=18827666,ybottom=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE)*-0.1,xright=24889549,ytop=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  for (gene in 1:nrow(candidateRgene)){
    #print(candidateRgene[gene,1])
    circos.segments(sector.index=as.factor(candidateRgene[gene,2]),x0=candidateRgene[gene,3],y0=0,x1=candidateRgene[gene,3],y1=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=0.4,lty=1,col="darkgrey")
  }
  
  circos.trackLines(dataset$chr, dataset$startwin+window*0.5, dataset$PropSNPWinPvalueInf01,lwd=2,col="darkorange2")

  ### p < 0.001
  circos.trackPlotRegion(factors= dataset$chr, y=dataset$PropSNPWinPvalueInf001, track.height=0.08, track.margin=c(0.03,0.03), bg.border = NA, panel.fun=function(startwin,PropSNPWinPvalueInf001){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    circos.axis(h=-max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE)*0.15, direction="inside",labels.facing="inside", lwd=0.5,labels.cex=0.4)
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.6)
  })
  
  # metaQTL Diana
  circos.rect(sector.index=3,xleft=21605249,ybottom=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE)*-0.1,xright=24567362,ytop=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=3,xleft=34220246,ybottom=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE)*-0.1,xright=37772912,ytop=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=2414969,ybottom=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE)*-0.1,xright=4219224,ytop=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=18827666,ybottom=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE)*-0.1,xright=24889549,ytop=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  for (gene in 1:nrow(candidateRgene)){
    #print(candidateRgene[gene,1])
    circos.segments(sector.index=as.factor(candidateRgene[gene,2]),x0=candidateRgene[gene,3],y0=0,x1=candidateRgene[gene,3],y1=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE),lwd=0.4,lty=1,col="darkgrey")
  }
  
  circos.trackLines(dataset$chr, dataset$startwin+window*0.5, dataset$PropSNPWinPvalueInf001,lwd=2,col="firebrick")

 dev.off() 
}




