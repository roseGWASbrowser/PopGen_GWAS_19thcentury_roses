## TL - 021121

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/")
densitygene=read.table("densitygene_slidwin500kb.txt",h=T)

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/")
chr0=read.table("chr0_limit_scaffold_syntheticchr0.txt",h=F)

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/Height_mean_gwas_results/")
#gwasres=read.table("Height_mean_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="Height_mean"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/Height_sd_gwas_results/")
#gwasres=read.table("Height_sd_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="Height_sd"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/Height_2012_gwas_results/")
#gwasres=read.table("Height_2012_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="Height_2012"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/Height_2013_gwas_results/")
#gwasres=read.table("Height_2013_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="Height_2013"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/Height_2014_gwas_results/")
#gwasres=read.table("Height_2014_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="Height_2014"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/Height_NA_out_mean_gwas_results/")
#gwasres=read.table("Height_NA_out_mean_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="Height_NA_out_mean"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/Height_NA_out_sd_gwas_results/")
#gwasres=read.table("Height_NA_out_sd_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="Height_NA_out_sd"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/BLUP_height_gwas_results/")
#gwasres=read.table("BLUP_height_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="BLUP_height"

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

  
  #circos.rect(sector.index=2,xleft=61857659,ybottom=max(dataset$Prop_Win_Total)*-0.1,xright=73734182,ytop=max(dataset$Prop_Win_Total)*1.1,lwd=0.01,lty=1,col="pink")
  #circos.rect(sector.index=6,xleft=46446433,ybottom=max(dataset$Prop_Win_Total)*-0.1,xright=56572104,ytop=max(dataset$Prop_Win_Total)*1.1,lwd=0.01,lty=1,col="pink")
  #circos.rect(sector.index=5,xleft=1946573,ybottom=max(dataset$Prop_Win_Total)*-0.1,xright=14484359,ytop=max(dataset$Prop_Win_Total)*1.1,lwd=0.01,lty=1,col="pink")
  
  
  # annotations
  #circos.segments(sector.index=1,x0=41514796,y0=0,x1=41514796,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=1,x0=52096899,y0=0,x1=52096899,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=1,x0=52554849,y0=0,x1=52554849,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=11833234,y0=0,x1=11833234,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=17461203,y0=0,x1=17461203,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
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

  # annotation QTL
  #circos.rect(sector.index=2,xleft=61857659,ybottom=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*-0.1,xright=73734182,ytop=max(-log10(gwasres[[columntoplot]])*1.1,na.rm=TRUE),lwd=0.01,lty=1,col="pink")
  #circos.rect(sector.index=6,xleft=46446433,ybottom=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*-0.1,xright=56572104,ytop=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  #circos.rect(sector.index=5,xleft=1946573,ybottom=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*-0.1,xright=14484359,ytop=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  
  # annotation genes
  #circos.segments(sector.index=1,x0=41514796,y0=0,x1=41514796,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=1,x0=52096899,y0=0,x1=52096899,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=1,x0=52554849,y0=0,x1=52554849,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=11833234,y0=0,x1=11833234,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=17461203,y0=0,x1=17461203,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=57253197,y0=0,x1=57253197,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=57254457,y0=0,x1=57254457,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=67016526,y0=0,x1=67016526,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=67036489,y0=0,x1=67036489,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=3,x0=33238136,y0=0,x1=33238136,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=4,x0=46474063,y0=0,x1=46474063,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=5,x0=8425314,y0=0,x1=8425314,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=5,x0=8430126,y0=0,x1=8430126,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=5,x0=61446161,y0=0,x1=61446161,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=6,x0=58611899,y0=0,x1=58611899,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=6,x0=60224075,y0=0,x1=60224075,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=4414048,y0=0,x1=4414048,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=4468385,y0=0,x1=4468385,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=50555186,y0=0,x1=50555186,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=52218507,y0=0,x1=52218507,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  
  #circos.rect(sector.index=3,xleft=41280000,ybottom=0,xright=41990000,ytop=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.01,lty=1,col="purple",border=0)
  
  # WAG
  #circos.segments(sector.index=2,x0=59116541,y0=0,x1=59116541,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.75,lwd=2,lty=2,col="forestgreen")
  #circos.segments(sector.index=3,x0=28893385,y0=0,x1=28893385,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.75,lwd=2,lty=2,col="forestgreen")
  #circos.segments(sector.index=3,x0=33215575,y0=0,x1=33215575,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.75,lwd=2,lty=2,col="forestgreen")
  #circos.segments(sector.index=4,x0=58205451,y0=0,x1=58205451,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.75,lwd=2,lty=2,col="forestgreen")
  #circos.segments(sector.index=5,x0=8125816,y0=0,x1=8125816,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.75,lwd=2,lty=2,col="forestgreen")
  
  
  #circos.text(sector.index=1, 41514796, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "AGL6",cex=0.75,col="purple")
  #circos.text(sector.index=1, 52096899, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "TOE2",cex=0.75,col="purple")
  #circos.text(sector.index=1, 52554849, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "euAP1",cex=0.75,col="purple")  
  
  #circos.text(sector.index=2, 11833234, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "TM6",cex=0.75,col="purple")  
  #circos.text(sector.index=2, 17461203, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "AP2",cex=0.75,col="purple")
  #circos.text(sector.index=2, 57253197, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "PLE",cex=0.75,col="purple") 
  #circos.text(sector.index=2, 57254457, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "Agamous D",cex=0.75,col="purple")
  #circos.text(sector.index=2, 67016526, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "euFUL",cex=0.75,col="purple") 
  #circos.text(sector.index=2, 67036489, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "SEP1/2",cex=0.75,col="purple") 

  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE),lwd=0.2,col = "lightgrey")
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
  
  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE),lwd=0.2,col = "lightgrey")
  }
  
  circos.trackLines(dataset$chr, dataset$startwin+window*0.5, dataset$PropSNPWinPvalueInf001,lwd=2,col="firebrick")

 dev.off() 
}



