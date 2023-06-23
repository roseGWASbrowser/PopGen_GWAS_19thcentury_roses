## TL - 021121

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/")
densitygene=read.table("densitygene_slidwin500kb.txt",h=T)

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/")
chr0=read.table("chr0_limit_scaffold_syntheticchr0.txt",h=F)

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/Nb_prickles_branch_gwas_results/")
gwasres=read.table("Nb_prickles_branch_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait="Nb_prickles_branch"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/Nb_prickles_branch_binaryTatiana_gwas_results/")
#gwasres=read.table("Nb_prickles_branch_binaryTatiana_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="Nb_prickles_branch_binaryTatiana"


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

  
  circos.segments(sector.index="0",x0=0,y0=max(dataset$Prop_Win_Total)*1.305,x1=52084813,y1=max(dataset$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=1,x0=0,y0=max(dataset$Prop_Win_Total)*1.305,x1=64743805,y1=max(dataset$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=2,x0=0,y0=max(dataset$Prop_Win_Total)*1.305,x1=75091153,y1=max(dataset$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=3,x0=0,y0=max(dataset$Prop_Win_Total)*1.305,x1=46750189,y1=max(dataset$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=4,x0=0,y0=max(dataset$Prop_Win_Total)*1.305,x1=58981928,y1=max(dataset$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=5,x0=0,y0=max(dataset$Prop_Win_Total)*1.305,x1=85844054,y1=max(dataset$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=6,x0=0,y0=max(dataset$Prop_Win_Total)*1.305,x1=67387004,y1=max(dataset$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=7,x0=0,y0=max(dataset$Prop_Win_Total)*1.305,x1=67077172,y1=max(dataset$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  
  # QTLs
  circos.rect(sector.index=1,xleft=6412287,ybottom=max(dataset$Prop_Win_Total)*1.21,xright=7764439,ytop=max(dataset$Prop_Win_Total)*1.4,lwd=0.01,lty=1,col="#FFA500") # lightcoral = 2017
  circos.rect(sector.index=3,xleft=27934327,ybottom=max(dataset$Prop_Win_Total)*1.21,xright=46440369,ytop=max(dataset$Prop_Win_Total)*1.3,lwd=0.01,lty=1,col="#FFA500") #  = 2019
  circos.rect(sector.index=3,xleft=41648024,ybottom=max(dataset$Prop_Win_Total)*1.31,xright=42317122,ytop=max(dataset$Prop_Win_Total)*1.4,lwd=0.01,lty=1,col="#FFA500") #  = 2019
  circos.rect(sector.index=4,xleft=34803638,ybottom=max(dataset$Prop_Win_Total)*1.21,xright=56107784,ytop=max(dataset$Prop_Win_Total)*1.4,lwd=0.01,lty=1,col="#FFA500") #  = 2019
  circos.rect(sector.index=6,xleft=1518964,ybottom=max(dataset$Prop_Win_Total)*1.21,xright=44264630,ytop=max(dataset$Prop_Win_Total)*1.4,lwd=0.01,lty=1,col="#FFA500") #  = 2019
  
  # annotations
  
  circos.segments(sector.index=1,x0=44468298,y0=0,x1=44468298,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=47708266,y0=0,x1=47708266,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=62070383,y0=0,x1=62070383,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=63982095,y0=0,x1=63982095,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=2,x0=2470719,y0=0,x1=2470719,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=47908413,y0=0,x1=47908413,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=54366345,y0=0,x1=54366345,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=3,x0=23331984,y0=0,x1=23331984,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33399696,y0=0,x1=33399696,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=39896892,y0=0,x1=39896892,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=3,x0=52096899,y0=0,x1=52096899,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=57125905,y0=0,x1=57125905,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=50315805,y0=0,x1=50315805,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=6,x0=52002793,y0=0,x1=52002793,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=55856328,y0=0,x1=55856328,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=7,x0=11958961,y0=0,x1=11958961,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=15536877,y0=0,x1=15536877,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=1,x0=63982095,y0=0,x1=63982095,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33399696,y0=0,x1=33399696,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=3,x0=52096899,y0=0,x1=52096899,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=57125905,y0=0,x1=57125905,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=50315805,y0=0,x1=50315805,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=57253197,y0=0,x1=57253197,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=57254457,y0=0,x1=57254457,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=67016526,y0=0,x1=67016526,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=67036489,y0=0,x1=67036489,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=3,x0=33238136,y0=0,x1=33238136,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=4,x0=46474063,y0=0,x1=46474063,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=5,x0=8425314,y0=0,x1=8425314,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=5,x0=8430126,y0=0,x1=8430126,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=5,x0=61446161,y0=0,x1=61446161,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=6,x0=58611899,y0=0,x1=58611899,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=6,x0=60224075,y0=0,x1=60224075,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=4414048,y0=0,x1=4414048,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=4468385,y0=0,x1=4468385,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=50555186,y0=0,x1=50555186,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=52218507,y0=0,x1=52218507,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.rect(sector.index=3,xleft=41280000,ybottom=0,xright=41990000,ytop=max(dataset$Prop_Win_Total,na.rm=TRUE),lwd=0.01,lty=1,col="purple",border=0)
  # WAG
  #circos.segments(sector.index=2,x0=59116541,y0=0,x1=59116541,y1=max(dataset$Prop_Win_Total),lwd=2,lty=2,col="forestgreen")
  #circos.segments(sector.index=3,x0=28893385,y0=0,x1=28893385,y1=max(dataset$Prop_Win_Total),lwd=2,lty=2,col="forestgreen")
  #circos.segments(sector.index=3,x0=33215504,y0=0,x1=33215504,y1=max(dataset$Prop_Win_Total),lwd=2,lty=2,col="forestgreen")
  #circos.segments(sector.index=4,x0=58205451,y0=0,x1=58205451,y1=max(dataset$Prop_Win_Total),lwd=2,lty=2,col="forestgreen")
  #circos.segments(sector.index=5,x0=8125816,y0=0,x1=8125816,y1=max(dataset$Prop_Win_Total),lwd=2,lty=2,col="forestgreen")
  
  
  # circos.segments(sector.index=2,x0=65290839,y0=0,x1=65290839,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="pink")
  # circos.segments(sector.index=6,x0=53520115,y0=0,x1=53520115,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="pink")
  # circos.segments(sector.index=5,x0=5255677,y0=0,x1=5255677,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="pink")
  
  # limits chr0
  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(rbind(max(dataset$Prop_Win_Total,na.rm=TRUE),max(densitygene$nb_genes/20000,na.rm=TRUE))*0.8),lwd=0.2,col = "lightgrey")
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

  # annotations
  circos.segments(sector.index=1,x0=44468298,y0=0,x1=44468298,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=47708266,y0=0,x1=47708266,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=62070383,y0=0,x1=62070383,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=63982095,y0=0,x1=63982095,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")

  circos.segments(sector.index=2,x0=2470719,y0=0,x1=2470719,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=47908413,y0=0,x1=47908413,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=54366345,y0=0,x1=54366345,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=3,x0=23331984,y0=0,x1=23331984,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33399696,y0=0,x1=33399696,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=39896892,y0=0,x1=39896892,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=3,x0=52096899,y0=0,x1=52096899,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=57125905,y0=0,x1=57125905,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=50315805,y0=0,x1=50315805,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")

  circos.segments(sector.index=6,x0=52002793,y0=0,x1=52002793,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=55856328,y0=0,x1=55856328,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=7,x0=11958961,y0=0,x1=11958961,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=15536877,y0=0,x1=15536877,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  
  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE),lwd=0.2,col = "lightgrey")
  }
  
  # values
  circos.trackPoints(gwasres$Chrom, gwasres$Position, -log10(gwasres[[columntoplot]]),pch=20,col=ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.05),"firebrick",ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.20),"goldenrod2",ifelse(gwasres$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.05),2.5,ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.10),2,ifelse(-log10(gwasres[[columntofilter]])>=-log10(0.20),1.5,0.75))))

  circos.text(sector.index=1, 44468298, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "MYC1",cex=0.75,col="purple")
  circos.text(sector.index=1, 47708266, max(-log10(gwasres[[columntoplot]])*0.9,na.rm=TRUE), "CPC",cex=0.75,col="purple")
  circos.text(sector.index=1, 62070383, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "TRY",cex=0.75,col="purple")
  circos.text(sector.index=1, 63982095, max(-log10(gwasres[[columntoplot]])*0.9,na.rm=TRUE), "TTG1",cex=0.75,col="purple")
 
  circos.text(sector.index=2, 2470719, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "MYB82",cex=0.75,col="purple") 
  circos.text(sector.index=2, 47908413, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "ZFP1\n-like1",cex=0.75,col="purple")
  circos.text(sector.index=2, 54366345, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "GL2",cex=0.75,col="purple")
  
  circos.text(sector.index=3, 23331984, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "GIS2",cex=0.75,col="purple")
  circos.text(sector.index=3, 33399696, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "TTG2",cex=0.75,col="purple")
  circos.text(sector.index=3, 39896892, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "MYB61",cex=0.75,col="purple")
  
  circos.text(sector.index=4, 50315805, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "GIS3",cex=0.75,col="purple")
  circos.text(sector.index=4, 57125905, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "ZFP5",cex=0.75,col="purple") 

  circos.text(sector.index=6, 52002793, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "TT8",cex=0.75,col="purple") 
  circos.text(sector.index=6, 55856328, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "ZFP1\n-like2",cex=0.75,col="purple") 
  
  circos.text(sector.index=7, 11958961, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "GL1",cex=0.75,col="purple") 
  circos.text(sector.index=7, 15536877, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "GL3",cex=0.75,col="purple") 
  
    # 3rd track = density significant p < 0.05

  circos.trackPlotRegion(factors= dataset$chr, y=dataset$PropSNPWinPvalueInf05, track.height=0.08, track.margin=c(0.03,0.03), bg.border = NA, panel.fun=function(startwin,PropSNPWinPvalueInf05){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.6)
  })
  
  circos.segments(sector.index=1,x0=44468298,y0=0,x1=44468298,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=47708266,y0=0,x1=47708266,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=62070383,y0=0,x1=62070383,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=63982095,y0=0,x1=63982095,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=2,x0=2470719,y0=0,x1=2470719,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=47908413,y0=0,x1=47908413,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=54366345,y0=0,x1=54366345,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=3,x0=23331984,y0=0,x1=23331984,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33399696,y0=0,x1=33399696,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=39896892,y0=0,x1=39896892,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=3,x0=52096899,y0=0,x1=52096899,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=57125905,y0=0,x1=57125905,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=50315805,y0=0,x1=50315805,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=6,x0=52002793,y0=0,x1=52002793,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=55856328,y0=0,x1=55856328,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=7,x0=11958961,y0=0,x1=11958961,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=15536877,y0=0,x1=15536877,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  
  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE),lwd=0.2,col = "lightgrey")
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

  circos.segments(sector.index=1,x0=44468298,y0=0,x1=44468298,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=47708266,y0=0,x1=47708266,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=62070383,y0=0,x1=62070383,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=63982095,y0=0,x1=63982095,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=2,x0=2470719,y0=0,x1=2470719,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=47908413,y0=0,x1=47908413,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=54366345,y0=0,x1=54366345,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=3,x0=23331984,y0=0,x1=23331984,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33399696,y0=0,x1=33399696,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=39896892,y0=0,x1=39896892,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=3,x0=52096899,y0=0,x1=52096899,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=57125905,y0=0,x1=57125905,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=50315805,y0=0,x1=50315805,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=6,x0=52002793,y0=0,x1=52002793,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=55856328,y0=0,x1=55856328,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=7,x0=11958961,y0=0,x1=11958961,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=15536877,y0=0,x1=15536877,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")

  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=0.2,col = "lightgrey")
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
  
  circos.segments(sector.index=1,x0=44468298,y0=0,x1=44468298,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=47708266,y0=0,x1=47708266,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=62070383,y0=0,x1=62070383,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=63982095,y0=0,x1=63982095,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=2,x0=2470719,y0=0,x1=2470719,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=47908413,y0=0,x1=47908413,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=54366345,y0=0,x1=54366345,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=3,x0=23331984,y0=0,x1=23331984,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33399696,y0=0,x1=33399696,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=39896892,y0=0,x1=39896892,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=3,x0=52096899,y0=0,x1=52096899,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=57125905,y0=0,x1=57125905,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=50315805,y0=0,x1=50315805,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=6,x0=52002793,y0=0,x1=52002793,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=55856328,y0=0,x1=55856328,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  
  circos.segments(sector.index=7,x0=11958961,y0=0,x1=11958961,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=15536877,y0=0,x1=15536877,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  
  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE),lwd=0.2,col = "lightgrey")
  }
  
  circos.trackLines(dataset$chr, dataset$startwin+window*0.5, dataset$PropSNPWinPvalueInf001,lwd=2,col="firebrick")

 dev.off() 
}



