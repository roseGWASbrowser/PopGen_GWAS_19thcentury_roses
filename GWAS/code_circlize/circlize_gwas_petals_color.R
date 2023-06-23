## TL - 021121

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/")
densitygene=read.table("densitygene_slidwin500kb.txt",h=T)

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/")
chr0=read.table("chr0_limit_scaffold_syntheticchr0.txt",h=F)

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/main_petal_color_white_vs_other_gwas_results/")
#gwasres=read.table("main_petal_color_white_vs_other_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="main_petal_color_white_vs_other"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/main_petal_color_rose_vs_other_gwas_results/")
#gwasres=read.table("main_petal_color_rose_vs_other_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="main_petal_color_rose_vs_other"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/main_petal_color_purple_vs_other_gwas_results/")
#gwasres=read.table("main_petal_color_purple_vs_other_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="main_petal_color_purple_vs_other"

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/main_petal_color_redpinkpurple_vs_other_gwas_results/")
#gwasres=read.table("main_petal_color_redpinkpurple_vs_other_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="main_petal_color_redpinkpurple_vs_other"

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige/color_distribution_UnicolorOrNot_gwas_results/")
gwasres=read.table("color_distribution_UnicolorOrNot_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait="color_distribution_UnicolorOrNot"

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
  
  
  # fun_color_range <- colorRampPalette(c("orangered4","orange"))
  #un_color_range(3)
  #[1] "#8B2500" "#C46500" "#FFA500"
  circos.segments(sector.index="0",x0=0,y0=max(dataset$Prop_Win_Total)*1.405,x1=52084813,y1=max(dataset$Prop_Win_Total)*1.405,lwd=0.3,lty=1,col="#C46500")
  circos.segments(sector.index=1,x0=0,y0=max(dataset$Prop_Win_Total)*1.405,x1=64743805,y1=max(dataset$Prop_Win_Total)*1.405,lwd=0.3,lty=1,col="#C46500")
  circos.segments(sector.index=2,x0=0,y0=max(dataset$Prop_Win_Total)*1.405,x1=75091153,y1=max(dataset$Prop_Win_Total)*1.405,lwd=0.3,lty=1,col="#C46500")
  circos.segments(sector.index=3,x0=0,y0=max(dataset$Prop_Win_Total)*1.405,x1=46750189,y1=max(dataset$Prop_Win_Total)*1.405,lwd=0.3,lty=1,col="#C46500")
  circos.segments(sector.index=4,x0=0,y0=max(dataset$Prop_Win_Total)*1.405,x1=58981928,y1=max(dataset$Prop_Win_Total)*1.405,lwd=0.3,lty=1,col="#C46500")
  circos.segments(sector.index=5,x0=0,y0=max(dataset$Prop_Win_Total)*1.405,x1=85844054,y1=max(dataset$Prop_Win_Total)*1.405,lwd=0.3,lty=1,col="#C46500")
  circos.segments(sector.index=6,x0=0,y0=max(dataset$Prop_Win_Total)*1.405,x1=67387004,y1=max(dataset$Prop_Win_Total)*1.405,lwd=0.3,lty=1,col="#C46500")
  circos.segments(sector.index=7,x0=0,y0=max(dataset$Prop_Win_Total)*1.405,x1=67077172,y1=max(dataset$Prop_Win_Total)*1.405,lwd=0.3,lty=1,col="#C46500")
  
  # QTL OW
  
  # QTL 95%
  # 2021
  circos.rect(sector.index=2,xleft=9101761,ybottom=max(dataset$Prop_Win_Total)*1.31,xright=41938556,ytop=max(dataset$Prop_Win_Total)*1.5,lwd=0.01,lty=1,col="#C46500") # lightcoral = 2017
  circos.rect(sector.index=6,xleft=45439950,ybottom=max(dataset$Prop_Win_Total)*1.31,xright=51413466,ytop=max(dataset$Prop_Win_Total)*1.5,lwd=0.01,lty=1,col="#C46500") #  = 2019
  circos.rect(sector.index=7,xleft=344677,ybottom=max(dataset$Prop_Win_Total)*1.31,xright=13464107,ytop=max(dataset$Prop_Win_Total)*1.5,lwd=0.01,lty=1,col="#C46500") #  = 2019
  circos.rect(sector.index=7,xleft=54348019,ybottom=max(dataset$Prop_Win_Total)*1.31,xright=60580661,ytop=max(dataset$Prop_Win_Total)*1.5,lwd=0.01,lty=1,col="#C46500") #  = 2019
  
  
  # limits chr0
  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(rbind(max(dataset$Prop_Win_Total,na.rm=TRUE),max(densitygene$nb_genes/20000,na.rm=TRUE))*0.8),lwd=0.2,col = "lightgrey")
  }
  
  
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

  # annotations
  #cluster ancestral nudix (RC4G0402400,RC4G0402000,RC4G0402100,RC4G0401900)
  circos.segments(sector.index=4,x0=51340000,y0=0,x1=51340000,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-1a (RC0G0033800,RC0G0033900,RC0G0034100,RC0G0034200,RC0G0034300,RC0G0034400)
  circos.segments(sector.index="0",x0=3210000,y0=0,x1=3210000,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-2b (RC6G0027500)
  circos.segments(sector.index=6,x0=2816735,y0=0,x1=2816735,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-2c (not on OB but btw RC7G0379900 & RC7GO380000 assuming local synteny)
  circos.segments(sector.index=7,x0=37345000,y0=0,x1=37345000,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  
  # annotations Laurence
  circos.segments(sector.index=4,x0=46189442,y0=0,x1=46189442,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=6,x0=46729697,y0=0,x1=46729697,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=1028902,y0=0,x1=1028902,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=60240306,y0=0,x1=60240306,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  
  #SPL : squamosa promoter-binding protein
  # SP9 important gene according to https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0110-3/MediaObjects/41588_2018_110_MOESM1_ESM.pdf p30
  circos.segments(sector.index=1,x0=38104574,y0=0,x1=38104574,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=38066005,y0=0,x1=38066005,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=38428370,y0=0,x1=38428370,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=74455144,y0=0,x1=74455144,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=3,x0=22150938,y0=0,x1=22150938,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=38159209,y0=0,x1=38159209,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=44438170,y0=0,x1=44438170,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=44447115,y0=0,x1=44447115,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46455803,y0=0,x1=46455803,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46559516,y0=0,x1=46559516,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=52613548,y0=0,x1=52613548,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=5908590,y0=0,x1=5908590,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=25209212,y0=0,x1=25209212,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=4363642,y0=0,x1=4363642,y1=max(dataset$Prop_Win_Total),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4377780,y0=0,x1=4377780,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=9170060,y0=0,x1=9170060,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=34735777,y0=0,x1=34735777,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  
  # other putative Raymond 2018 https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0110-3/MediaObjects/41588_2018_110_MOESM1_ESM.pdf p30
  circos.segments(sector.index=1,x0=3298655,y0=0,x1=3298655,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=51114814,y0=0,x1=51114814,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52684502,y0=0,x1=52684502,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=61663883,y0=0,x1=61663883,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11399754,y0=0,x1=11399754,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=72308803,y0=0,x1=72308803,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=3,x0=45914297,y0=0,x1=45914297,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=31326504,y0=0,x1=31326504,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=6,x0=57312795,y0=0,x1=57312795,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4048486,y0=0,x1=4048486,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=17887444,y0=0,x1=17887444,y1=max(dataset$Prop_Win_Total)*0.8,lwd=0.3,lty=1,col="purple")
  
  # circos.segments(sector.index=2,x0=65290839,y0=0,x1=65290839,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="pink")
  # circos.segments(sector.index=6,x0=53520115,y0=0,x1=53520115,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="pink")
  # circos.segments(sector.index=5,x0=5255677,y0=0,x1=5255677,y1=max(dataset$Prop_Win_Total),lwd=2,lty=1,col="pink")

  
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

  for (line in 1:nrow(chr0)){
    circos.segments(sector.index = "0",x0=chr0[line,2],y0=0,x1=chr0[line,2],y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE),lwd=0.2,col = "lightgrey")
  }
  
  
  # annotation QTL
  #circos.rect(sector.index=2,xleft=61857659,ybottom=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*-0.1,xright=73734182,ytop=max(-log10(gwasres[[columntoplot]])*1.1,na.rm=TRUE),lwd=0.01,lty=1,col="pink")
  #circos.rect(sector.index=6,xleft=46446433,ybottom=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*-0.1,xright=56572104,ytop=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  #circos.rect(sector.index=5,xleft=1946573,ybottom=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*-0.1,xright=14484359,ytop=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  #
  # annotation roAP2
  # annotations
  #cluster ancestral nudix (RC4G0402400,RC4G0402000,RC4G0402100,RC4G0401900)
  circos.segments(sector.index=4,x0=51340000,y0=0,x1=51340000,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-1a (RC0G0033800,RC0G0033900,RC0G0034100,RC0G0034200,RC0G0034300,RC0G0034400)
  circos.segments(sector.index="0",x0=3210000,y0=0,x1=3210000,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-2b (RC6G0027500)
  circos.segments(sector.index=6,x0=2816735,y0=0,x1=2816735,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-2c (not on OB but btw RC7G0379900 & RC7GO380000 assuming local synteny)
  circos.segments(sector.index=7,x0=37345000,y0=0,x1=37345000,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")

  
  # annotations Laurence
  circos.segments(sector.index=4,x0=46189442,y0=0,x1=46189442,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=6,x0=46729697,y0=0,x1=46729697,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=1028902,y0=0,x1=1028902,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=60240306,y0=0,x1=60240306,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  
  #SPL : squamosa promoter-binding protein
  # SP9 important gene according to https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0110-3/MediaObjects/41588_2018_110_MOESM1_ESM.pdf p30
  circos.segments(sector.index=1,x0=38104574,y0=0,x1=38104574,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=38066005,y0=0,x1=38066005,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=38428370,y0=0,x1=38428370,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=74455144,y0=0,x1=74455144,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=3,x0=22150938,y0=0,x1=22150938,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=38159209,y0=0,x1=38159209,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=44438170,y0=0,x1=44438170,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=44447115,y0=0,x1=44447115,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46455803,y0=0,x1=46455803,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46559516,y0=0,x1=46559516,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=52613548,y0=0,x1=52613548,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=5908590,y0=0,x1=5908590,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=25209212,y0=0,x1=25209212,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=4363642,y0=0,x1=4363642,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4377780,y0=0,x1=4377780,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=9170060,y0=0,x1=9170060,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=34735777,y0=0,x1=34735777,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  
  # other putative Raymond 2018 https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0110-3/MediaObjects/41588_2018_110_MOESM1_ESM.pdf p30
  circos.segments(sector.index=1,x0=3298655,y0=0,x1=3298655,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=51114814,y0=0,x1=51114814,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52684502,y0=0,x1=52684502,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=61663883,y0=0,x1=61663883,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11399754,y0=0,x1=11399754,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=72308803,y0=0,x1=72308803,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=3,x0=45914297,y0=0,x1=45914297,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=31326504,y0=0,x1=31326504,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=6,x0=57312795,y0=0,x1=57312795,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4048486,y0=0,x1=4048486,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=17887444,y0=0,x1=17887444,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  
  circos.text(sector.index=4, 51340000, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*1.2, "ancestral\nNUDIX1",cex=0.75,col="purple")
  circos.text(sector.index="0", 3210000, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "NUDIX1-1a",cex=0.75,col="purple")
  circos.text(sector.index=6, 2816735, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "NUDIX1-2b",cex=0.75,col="purple")  
  circos.text(sector.index=7, 37345000, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "NUDIX1-2c",cex=0.75,col="purple")  
  
  circos.text(sector.index=4, 46189442, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*1.1, "Histidine\nKinase",cex=0.75,col="purple")
  circos.text(sector.index=6, 46729697, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "DENN",cex=0.75,col="purple")
  circos.text(sector.index=7, 60240306, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "PPR",cex=0.75,col="purple")  
  
  circos.text(sector.index=1, 38100000, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "SPL8&\nSPL10",cex=0.75,col="purple")  
  circos.text(sector.index=2, 38428370, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "SPL4",cex=0.75,col="purple")  
  circos.text(sector.index=2, 74455144, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "SPL",cex=0.75,col="purple")  
  circos.text(sector.index=3, 22150938, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "SPL9",cex=0.75,col="purple")  
  circos.text(sector.index=4, 38159209, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "SPL7",cex=0.75,col="purple")
  circos.text(sector.index=4, 44440000, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*1.3, "SPL12&\nSPL1",cex=0.75,col="purple")
  circos.text(sector.index=4, 46455803, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "+SPL",cex=0.75,col="purple")
  circos.text(sector.index=4, 52613548, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "SPL9",cex=0.75,col="purple") 
  circos.text(sector.index=5, 5908590, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "SPL3",cex=0.75,col="purple") 
  circos.text(sector.index=5, 25209212, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "SPL",cex=0.75,col="purple") 
  circos.text(sector.index=7, 4377780, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "SPL3",cex=0.75,col="purple")
  circos.text(sector.index=7, 9170060, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "SPL",cex=0.75,col="purple") 
  circos.text(sector.index=7, 34735777, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "SPL14",cex=0.75,col="purple") 
  
  circos.text(sector.index=1, 3298655, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "CHS",cex=0.75,col="purple") 
  circos.text(sector.index=1, 51114814, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "PAL4",cex=0.75,col="purple") 
  circos.text(sector.index=1, 52684502, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "CHI",cex=0.75,col="purple") 
  circos.text(sector.index=1, 61663883, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "GT1",cex=0.75,col="purple") 
  circos.text(sector.index=2, 11399754, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "F3H",cex=0.75,col="purple") 
  circos.text(sector.index=2, 72308803, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "MYB21",cex=0.75,col="purple") 
  circos.text(sector.index=3, 45914297, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "MYB10",cex=0.75,col="purple") 
  circos.text(sector.index=5, 31326504, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "GDS",cex=0.75,col="purple") 
  circos.text(sector.index=6, 57312795, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "DFR",cex=0.75,col="purple") 
  circos.text(sector.index=7, 17887444, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE), "ANS",cex=0.75,col="purple") 
  circos.text(sector.index=7, 4048486, max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.9, "F3'H",cex=0.75,col="purple") 
  
  # other putative Raymond 2018 https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0110-3/MediaObjects/41588_2018_110_MOESM1_ESM.pdf p30
  circos.segments(sector.index=1,x0=3298655,y0=0,x1=3298655,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=51114814,y0=0,x1=51114814,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52684502,y0=0,x1=52684502,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=61663883,y0=0,x1=61663883,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11399754,y0=0,x1=11399754,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=72308803,y0=0,x1=72308803,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=3,x0=45914297,y0=0,x1=45914297,y1=max(-log10(gwasres[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.3,lty=1,col="purple")
  
  
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
  
  # NUDIX
  #cluster ancestral nudix (RC4G0402400,RC4G0402000,RC4G0402100,RC4G0401900)
  circos.segments(sector.index=4,x0=51340000,y0=0,x1=51340000,y1=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-1a (RC0G0033800,RC0G0033900,RC0G0034100,RC0G0034200,RC0G0034300,RC0G0034400)
  circos.segments(sector.index="0",x0=3210000,y0=0,x1=3210000,y1=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-2b (RC6G0027500)
  circos.segments(sector.index=6,x0=2816735,y0=0,x1=2816735,y1=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-2c (not on OB but btw RC7G0379900 & RC7GO380000 assuming local synteny)
  circos.segments(sector.index=7,x0=37345000,y0=0,x1=37345000,y1=max(dataset$PropSNPWinPvalueInf05,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  
  # annotations Laurence
  circos.segments(sector.index=4,x0=46189442,y0=0,x1=46189442,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=6,x0=46729697,y0=0,x1=46729697,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=1028902,y0=0,x1=1028902,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=60240306,y0=0,x1=60240306,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  
  
  #SPL : squamosa promoter-binding protein
  # SP9 important gene according to https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0110-3/MediaObjects/41588_2018_110_MOESM1_ESM.pdf p30
  circos.segments(sector.index=1,x0=38104574,y0=0,x1=38104574,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=38066005,y0=0,x1=38066005,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=38428370,y0=0,x1=38428370,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=74455144,y0=0,x1=74455144,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=3,x0=22150938,y0=0,x1=22150938,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=38159209,y0=0,x1=38159209,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=44438170,y0=0,x1=44438170,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=44447115,y0=0,x1=44447115,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46455803,y0=0,x1=46455803,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46559516,y0=0,x1=46559516,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=52613548,y0=0,x1=52613548,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=5908590,y0=0,x1=5908590,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=25209212,y0=0,x1=25209212,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=4363642,y0=0,x1=4363642,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4377780,y0=0,x1=4377780,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=9170060,y0=0,x1=9170060,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=34735777,y0=0,x1=34735777,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  
  # other putative Raymond 2018 https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0110-3/MediaObjects/41588_2018_110_MOESM1_ESM.pdf p30
  circos.segments(sector.index=1,x0=3298655,y0=0,x1=3298655,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=51114814,y0=0,x1=51114814,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52684502,y0=0,x1=52684502,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=61663883,y0=0,x1=61663883,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11399754,y0=0,x1=11399754,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=72308803,y0=0,x1=72308803,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=3,x0=45914297,y0=0,x1=45914297,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=31326504,y0=0,x1=31326504,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=6,x0=57312795,y0=0,x1=57312795,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4048486,y0=0,x1=4048486,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=17887444,y0=0,x1=17887444,y1=max(dataset$PropSNPWinPvalueInf05),lwd=0.3,lty=1,col="purple")
  
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
  
  # NUDIX
  #cluster ancestral nudix (RC4G0402400,RC4G0402000,RC4G0402100,RC4G0401900)
  circos.segments(sector.index=4,x0=51340000,y0=0,x1=51340000,y1=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-1a (RC0G0033800,RC0G0033900,RC0G0034100,RC0G0034200,RC0G0034300,RC0G0034400)
  circos.segments(sector.index="0",x0=3210000,y0=0,x1=3210000,y1=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-2b (RC6G0027500)
  circos.segments(sector.index=6,x0=2816735,y0=0,x1=2816735,y1=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-2c (not on OB but btw RC7G0379900 & RC7GO380000 assuming local synteny)
  circos.segments(sector.index=7,x0=37345000,y0=0,x1=37345000,y1=max(dataset$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  
  # annotations Laurence
  circos.segments(sector.index=4,x0=46189442,y0=0,x1=46189442,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=6,x0=46729697,y0=0,x1=46729697,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=1028902,y0=0,x1=1028902,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=60240306,y0=0,x1=60240306,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  
  #SPL : squamosa promoter-binding protein
  # SP9 important gene according to https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0110-3/MediaObjects/41588_2018_110_MOESM1_ESM.pdf p30
  circos.segments(sector.index=1,x0=38104574,y0=0,x1=38104574,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=38066005,y0=0,x1=38066005,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=38428370,y0=0,x1=38428370,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=74455144,y0=0,x1=74455144,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=3,x0=22150938,y0=0,x1=22150938,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=38159209,y0=0,x1=38159209,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=44438170,y0=0,x1=44438170,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=44447115,y0=0,x1=44447115,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46455803,y0=0,x1=46455803,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46559516,y0=0,x1=46559516,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=52613548,y0=0,x1=52613548,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=5908590,y0=0,x1=5908590,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=25209212,y0=0,x1=25209212,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=4363642,y0=0,x1=4363642,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4377780,y0=0,x1=4377780,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=9170060,y0=0,x1=9170060,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=34735777,y0=0,x1=34735777,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  
  # other putative Raymond 2018 https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0110-3/MediaObjects/41588_2018_110_MOESM1_ESM.pdf p30
  circos.segments(sector.index=1,x0=3298655,y0=0,x1=3298655,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=51114814,y0=0,x1=51114814,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52684502,y0=0,x1=52684502,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=61663883,y0=0,x1=61663883,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11399754,y0=0,x1=11399754,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=72308803,y0=0,x1=72308803,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=3,x0=45914297,y0=0,x1=45914297,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=31326504,y0=0,x1=31326504,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=6,x0=57312795,y0=0,x1=57312795,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4048486,y0=0,x1=4048486,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=17887444,y0=0,x1=17887444,y1=max(dataset$PropSNPWinPvalueInf01),lwd=0.3,lty=1,col="purple")
  
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
  
  # NUDIX
  #cluster ancestral nudix (RC4G0402400,RC4G0402000,RC4G0402100,RC4G0401900)
  circos.segments(sector.index=4,x0=51340000,y0=0,x1=51340000,y1=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-1a (RC0G0033800,RC0G0033900,RC0G0034100,RC0G0034200,RC0G0034300,RC0G0034400)
  circos.segments(sector.index="0",x0=3210000,y0=0,x1=3210000,y1=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-2b (RC6G0027500)
  circos.segments(sector.index=6,x0=2816735,y0=0,x1=2816735,y1=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  # cluster nudix 1-2c (not on OB but btw RC7G0379900 & RC7GO380000 assuming local synteny)
  circos.segments(sector.index=7,x0=37345000,y0=0,x1=37345000,y1=max(dataset$PropSNPWinPvalueInf001,na.rm=TRUE),lwd=0.3,lty=1,col="purple")
  
  # annotations Laurence
  circos.segments(sector.index=4,x0=46189442,y0=0,x1=46189442,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=6,x0=46729697,y0=0,x1=46729697,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=1028902,y0=0,x1=1028902,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=60240306,y0=0,x1=60240306,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  
  #SPL : squamosa promoter-binding protein
  # SP9 important gene according to https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0110-3/MediaObjects/41588_2018_110_MOESM1_ESM.pdf p30
  circos.segments(sector.index=1,x0=38104574,y0=0,x1=38104574,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=38066005,y0=0,x1=38066005,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=38428370,y0=0,x1=38428370,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=74455144,y0=0,x1=74455144,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=3,x0=22150938,y0=0,x1=22150938,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=38159209,y0=0,x1=38159209,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=44438170,y0=0,x1=44438170,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=44447115,y0=0,x1=44447115,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46455803,y0=0,x1=46455803,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46559516,y0=0,x1=46559516,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=4,x0=52613548,y0=0,x1=52613548,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=5908590,y0=0,x1=5908590,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=25209212,y0=0,x1=25209212,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  #circos.segments(sector.index=7,x0=4363642,y0=0,x1=4363642,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4377780,y0=0,x1=4377780,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=9170060,y0=0,x1=9170060,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=34735777,y0=0,x1=34735777,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  
  # other putative Raymond 2018 https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0110-3/MediaObjects/41588_2018_110_MOESM1_ESM.pdf p30
  circos.segments(sector.index=1,x0=3298655,y0=0,x1=3298655,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=51114814,y0=0,x1=51114814,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52684502,y0=0,x1=52684502,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=1,x0=61663883,y0=0,x1=61663883,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11399754,y0=0,x1=11399754,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=2,x0=72308803,y0=0,x1=72308803,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=3,x0=45914297,y0=0,x1=45914297,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=5,x0=31326504,y0=0,x1=31326504,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=6,x0=57312795,y0=0,x1=57312795,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4048486,y0=0,x1=4048486,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  circos.segments(sector.index=7,x0=17887444,y0=0,x1=17887444,y1=max(dataset$PropSNPWinPvalueInf001),lwd=0.3,lty=1,col="purple")
  
  circos.trackLines(dataset$chr, dataset$startwin+window*0.5, dataset$PropSNPWinPvalueInf001,lwd=2,col="firebrick")

 dev.off() 
}



