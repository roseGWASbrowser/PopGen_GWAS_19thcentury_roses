## TL - 021121

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/")
densitygene=read.table("densitygene_slidwin500kb.txt",h=T)

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/")
chr0=read.table("chr0_limit_scaffold_syntheticchr0.txt",h=F)

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/")
gene=read.table("candidate_genes_parfum_avril2022_fromBenoit.txt.withgenomicpositions",h=F)

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/scent2017/X35dimethoxytoluene_gwas_results")
#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/scent2018/X35dimethoxytoluene_gwas_results")
#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/scentAvg2years/X35dimethoxytoluene_gwas_results")
#gwasres=read.table("X35dimethoxytoluene_scores.withoutchrNA.sed",h=T)
#trait="X35dimethoxytoluene"
#label_legend="GWAS DMT (2017)"
#label_legend="GWAS DMT (2018)"
#label_legend="GWAS DMT (2 years)"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/scent2017/X35dimethoxytoluene_binary0p1_gwas_results")
#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/scent2018/X35dimethoxytoluene_binary0p1_gwas_results")
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/scentAvg2years/X35dimethoxytoluene_binary0p1_gwas_results")
gwasres=read.table("X35dimethoxytoluene_binary0p1_scores.withoutchrNA.sed",h=T)
trait="X35dimethoxytoluene_binary0p1"
#label_legend="GWAS DMT (binary 2017 - cutoff: 10)"
#label_legend="GWAS DMT (binary 2018 - cutoff: 10)"
label_legend="GWAS DMT (binary 2 years - cutoff: 10)"

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
  
  pdf(file=paste0("circlize_v2_",trait,"_",model,"slidwin",format(window,scientific=FALSE),".pdf"),width=12,height=12)
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

  # limits chr0
  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(rbind(max(dataset$Prop_Win_Total,na.rm=TRUE),max(densitygene$nb_genes/20000,na.rm=TRUE))*0.8),lwd=0.2,col = "lightgrey")
  }
  
  for (line in 1:nrow(gene)){
    circos.segments(sector.index=as.character(gene[line,4]),x0=gene[line,8],y0=0,x1=gene[line,8],y1=max(dataset$Prop_Win_Total),col=ifelse(gene[line,1]=="OMT","purple","lightcyan2"),lwd=ifelse(gene[line,1]=="OMT",2,1))
  }
  # OOMT1: RC2GO374900 & RC2G0375300
  #circos.segments(sector.index=2,x0=42110000,y0=0,x1=42110000,y1=max(dataset$Prop_Win_Total)*0.8,lwd=2,lty=1,col="purple")
  # OOMT2: RC2GO348900
  #circos.segments(sector.index=2,x0=39094532,y0=0,x1=39094532,y1=max(dataset$Prop_Win_Total)*0.8,lwd=2,lty=1,col="purple")  
  # OOMT2: RC2GO34970
  #circos.segments(sector.index=2,x0=39326245,y0=0,x1=39326245,y1=max(dataset$Prop_Win_Total)*0.8,lwd=2,lty=1,col="purple")  

   # density gene
  circos.trackLines(densitygene$chr, densitygene$startwin+(densitygene$endwin-densitygene$startwin), densitygene$nb_genes/20000, col="forestgreen", lwd=1.5, type="l")


  # annotations
  #circos.segments(sector.index=2,x0=3297631,y0=0,x1=3297631,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=4652707,y0=0,x1=4652707,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=39778447,y0=0,x1=39778447,y1=max(dataset$Prop_Win_Total)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=46406190,y0=0,x1=46406190,y1=max(dataset$Prop_Win_Total)*0.8,lwd=2,lty=1,col="purple")
  #circos.segments(sector.index=2,x0=69535318,y0=0,x1=69535318,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="purple")
  
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
  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE),lwd=0.2,col = "lightgrey")
  }
  # annotations
  for (line in 1:nrow(gene)){
    circos.text(sector.index=as.character(gene$V4[line]),x=gene$V8[line],y=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*.85,labels=paste(gene$V2[line],"(",gene$V3[line],")",sep=""),col=ifelse(gene$V1[line]=="OMT","purple","lightcyan2"),cex=ifelse(gene$V1[line]=="OMT",0.6,0.5),facing = "reverse.clockwise")
  }
  # OOMT1: RC2GO374900 & RC2G0375300
  #circos.segments(sector.index=2,x0=42110000,y0=0,x1=42110000,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  # OOMT2: RC2GO348900
  #circos.segments(sector.index=2,x0=39094532,y0=0,x1=39094532,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")  
  # OOMT2: RC2GO34970
  #circos.segments(sector.index=2,x0=39326245,y0=0,x1=39326245,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")  
  
  #circos.text(sector.index=2, 41000000, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), labels="OOMT2 & OOMT1",cex=0.75,col="purple")
 
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
  
  for (line in 1:nrow(gene)){
    circos.segments(sector.index=as.character(gene$V4[line]),x0=gene$V8[line],y0=0,x1=gene$V8[line],y1=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE),col=ifelse(gene$V1[line]=="OMT","purple","lightcyan2"),lwd=ifelse(gene$V1[line]=="OMT",2,1))
  }
  
  #
  # OOMT1: RC2GO374900 & RC2G0375300
  #circos.segments(sector.index=2,x0=42110000,y0=0,x1=42110000,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")
  # OOMT2: RC2GO348900
  #circos.segments(sector.index=2,x0=39094532,y0=0,x1=39094532,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")  
  # OOMT2: RC2GO34970
  #circos.segments(sector.index=2,x0=39326245,y0=0,x1=39326245,y1=max(dataset$PropSNPWinPvalueInf05),lwd=2,lty=1,col="purple")  
  
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
  
  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=0.2,col = "lightgrey")
  }
  
  for (line in 1:nrow(gene)){
    circos.segments(sector.index=as.character(gene$V4[line]),x0=gene$V8[line],y0=0,x1=gene$V8[line],y1=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE),col=ifelse(gene$V1[line]=="OMT","purple","lightcyan2"),lwd=ifelse(gene$V1[line]=="OMT",2,1))
  }
  #
  # OOMT1: RC2GO374900 & RC2G0375300
  #circos.segments(sector.index=2,x0=42110000,y0=0,x1=42110000,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")
  # OOMT2: RC2GO348900
  #circos.segments(sector.index=2,x0=39094532,y0=0,x1=39094532,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")  
  # OOMT2: RC2GO34970
  #circos.segments(sector.index=2,x0=39326245,y0=0,x1=39326245,y1=max(dataset$PropSNPWinPvalueInf01),lwd=2,lty=1,col="purple")  
  

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
  
  for (line in 1:nrow(gene)){
    circos.segments(sector.index=as.character(gene$V4[line]),x0=gene$V8[line],y0=0,x1=gene$V8[line],y1=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE),col=ifelse(gene$V1[line]=="OMT","purple","lightcyan2"),lwd=ifelse(gene$V1[line]=="OMT",2,1))
  }
  #
  # OOMT1: RC2GO374900 & RC2G0375300
  #circos.segments(sector.index=2,x0=42110000,y0=0,x1=42110000,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")
  # OOMT2: RC2GO348900
  #circos.segments(sector.index=2,x0=39094532,y0=0,x1=39094532,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")  
  # OOMT2: RC2GO34970
  #circos.segments(sector.index=2,x0=39326245,y0=0,x1=39326245,y1=max(dataset$PropSNPWinPvalueInf001),lwd=2,lty=1,col="purple")  
  
  circos.trackLines(dataset$chr, dataset$startwin+window*0.5, dataset$PropSNPWinPvalueInf001,lwd=2,col="firebrick")

  # track numbers
  points(x=-0.015, y = 0.99 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.99, labels="1",cex=1.1,col="white")
  
  points(x=-0.015, y = 0.845 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.845, labels="2",cex=1.1,col="white")
  
  points(x=-0.015, y = 0.65 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.65, labels="3",cex=1.1,col="white") 
  
  points(x=-0.015, y = 0.505 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.505, labels="4",cex=1.1,col="white") 
  
  points(x=-0.015, y = 0.36 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.36, labels="5",cex=1.1,col="white") 
  
  # Key
  text(x=0, y = 0.175, labels="KEY",cex=1.2,font=2)
  
  segments(x0=-0.13,y0=0.13,x1=-0.13,y1=0.09,col="darkgrey",lwd=2)
  segments(x0=-0.115,y0=0.12,x1=-0.085,y1=0.12,col="black",lwd=2)
  segments(x0=-0.115,y0=0.10,x1=-0.085,y1=0.10,col="forestgreen",lwd=1.5)
  text(x=-0.08, y = 0.118, labels=expression("SNP density (MAF" >= "0.05)"),cex=0.7,col="black",pos=4)
  text(x=-0.08, y = 0.098, labels="Gene density",cex=0.7,col="forestgreen",pos=4)
  points(x=-0.15, y = 0.11 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.15, y = 0.11, labels="1",cex=1,col="white")
  
  segments(x0=-0.13,y0=0.07,x1=-0.13,y1=-0.03,col="darkgrey",lwd=2)
  points(x=-0.10, y = 0.06 , type = "p", cex=2,col="navyblue",pch=20)
  points(x=-0.10, y = 0.04 , type = "p", cex=2,col="dodgerblue3",pch=20)
  text(x=-0.08, y = 0.05, labels=label_legend,cex=0.7,col="blue",pos=4)
  points(x=-0.10, y = 0.02 , type = "p", cex=2,col="goldenrod2",pch=20)
  text(x=-0.08, y = 0.02, labels=expression("SNP with q-value"<= "0.20"),cex=0.7,col="goldenrod2",pos=4)
  points(x=-0.10, y = 0 , type = "p", cex=2,col="darkorange2",pch=20)
  text(x=-0.08, y = 0, labels=expression("SNP with q-value"<= "0.10"),cex=0.7,col="darkorange2",pos=4)
  points(x=-0.10, y = -0.02 , type = "p", cex=2,col="firebrick",pch=20)
  text(x=-0.08, y = -0.02, labels=expression("SNP with q-value"<= "0.05"),cex=0.7,col="firebrick",pos=4)
  points(x=-0.15, y = 0.02 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.15, y = 0.02, labels="2",cex=1,col="white")
  
  segments(x0=-0.13,y0=-0.06,x1=-0.13,y1=-0.08,col="darkgrey",lwd=2)
  segments(x0=-0.115,y0=-0.07,x1=-0.085,y1=-0.07,col="goldenrod2",lwd=1)
  text(x=-0.08, y = -0.075, labels=expression("Density of SNPs with p-value"<="0.05"),cex=0.7,col="goldenrod2",pos=4)
  points(x=-0.15, y = -0.07, type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.15, y = -0.07, labels="3",cex=1,col="white")
  
  segments(x0=-0.13,y0=-0.09,x1=-0.13,y1=-0.11,col="darkgrey",lwd=2)
  segments(x0=-0.115,y0=-0.10,x1=-0.085,y1=-0.10,col="darkorange2",lwd=1)
  text(x=-0.08, y = -0.105, labels=expression("Density of SNPs with p-value"<="0.01"),cex=0.7,col="darkorange2",pos=4)
  points(x=-0.15, y = -0.10, type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.15, y = -0.10, labels="4",cex=1,col="white")
  
  segments(x0=-0.13,y0=-0.12,x1=-0.13,y1=-0.14,col="darkgrey",lwd=2)
  segments(x0=-0.115,y0=-0.13,x1=-0.085,y1=-0.13,col="firebrick",lwd=1)
  text(x=-0.08, y = -0.135, labels=expression("Density of SNPs with p-value"<="0.001"),cex=0.7,col="firebrick",pos=4)
  points(x=-0.15, y = -0.13, type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.15, y = -0.13, labels="5",cex=1,col="white")
  
  #segments(x0=-0.115,y0=-0.16,x1=-0.085,y1=-0.16,col="purple",lwd=1)
  #text(x=-0.08, y = -0.165, labels=expression("Candidate genes"),cex=0.7,col="purple",pos=4)
  
 dev.off() 
}
