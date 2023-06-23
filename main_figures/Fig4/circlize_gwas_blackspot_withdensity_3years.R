## TL - 021121 # MODIF 190623

setwd("/home/thibaultleroy/Rose/Papier/Fig4/")
densitygene=read.table("densitygene_slidwin500kb.txt",h=T)

candidateRgene=read.table("list_candidateRgenes_Laurence.txt",h=F)

#gwasres=read.table("TN2014_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="TN2014"

gwasres1=read.table("TN2014_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait1="TN2014"
gwasres2=read.table("TN2015_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait2="TN2015"
gwasres3=read.table("TN2016_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait3="TN2016"

# remove chr00
densitygene=subset(densitygene,densitygene$chr!=0)
candidateRgene=subset(candidateRgene,candidateRgene$V2!=0)
gwasres1=subset(gwasres1,gwasres1$Chrom!=0)
gwasres2=subset(gwasres2,gwasres2$Chrom!=0)
gwasres3=subset(gwasres3,gwasres3$Chrom!=0)

# for all models gwasres1
for (model in c("additive","diplo.additive","general","diplo.general","1.dom.alt","1.dom.ref","2.dom.alt","2.dom.ref")){
  print(paste("Currently working on model: ",model))
  # 2014
  ### slidingwindows R
  ## fist loop by chromosome
  dataset1 <- data.frame(matrix(ncol = 12, nrow = 0))
  colnames(dataset1)=c("chr","startwin","endwin","Nb_SNPs_Win","Prop_Win_Total","NbSNPWinPvalueInf05","NbSNPWinPvalueInf01","NbSNPWinPvalueInf001","SNPswithoutpvalues","PropSNPWinPvalueInf05","PropSNPWinPvalueInf01","PropSNPWinPvalueInf001")
  window=1000000
  for (chr in c(1:7)){
    print(paste("currently working on chr", chr))
    number_of_snp=nrow(gwasres1)
    currentchr=subset(gwasres1,gwasres1$Chrom==chr)
    SNPcurrentchr=nrow(currentchr)
    print(paste("proportion of SNPs on chr",chr,"=",SNPcurrentchr/number_of_snp))
    ### 
    for (slidwin in c(0:floor(max(currentchr$Position)/window))){  # or 1:ceiling = but here 0 to have something easier with position
      start=slidwin*window+1
      end=slidwin*window+window
      currentwindow=subset(currentchr,currentchr$Position>=start & currentchr$Position<=end)
      SNPwindows=nrow(currentwindow)
      columntosubsample=paste0(trait1,".",model,"_P_value")
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
      dataset1=rbind(dataset1,currentline)
    }
  }
  # 2015
    dataset2 <- data.frame(matrix(ncol = 12, nrow = 0))
    colnames(dataset2)=c("chr","startwin","endwin","Nb_SNPs_Win","Prop_Win_Total","NbSNPWinPvalueInf05","NbSNPWinPvalueInf01","NbSNPWinPvalueInf001","SNPswithoutpvalues","PropSNPWinPvalueInf05","PropSNPWinPvalueInf01","PropSNPWinPvalueInf001")
    window=1000000
    for (chr in c(1:7)){
      print(paste("currently working on chr", chr))
      number_of_snp=nrow(gwasres2)
      currentchr=subset(gwasres2,gwasres2$Chrom==chr)
      SNPcurrentchr=nrow(currentchr)
      print(paste("proportion of SNPs on chr",chr,"=",SNPcurrentchr/number_of_snp))
      ### 
      for (slidwin in c(0:floor(max(currentchr$Position)/window))){  # or 1:ceiling = but here 0 to have something easier with position
        start=slidwin*window+1
        end=slidwin*window+window
        currentwindow=subset(currentchr,currentchr$Position>=start & currentchr$Position<=end)
        SNPwindows=nrow(currentwindow)
        columntosubsample=paste0(trait2,".",model,"_P_value")
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
        dataset2=rbind(dataset2,currentline)
      }
    }
    # 2016
    dataset3 <- data.frame(matrix(ncol = 12, nrow = 0))
    colnames(dataset3)=c("chr","startwin","endwin","Nb_SNPs_Win","Prop_Win_Total","NbSNPWinPvalueInf05","NbSNPWinPvalueInf01","NbSNPWinPvalueInf001","SNPswithoutpvalues","PropSNPWinPvalueInf05","PropSNPWinPvalueInf01","PropSNPWinPvalueInf001")
    window=1000000
    for (chr in c(1:7)){
      print(paste("currently working on chr", chr))
      number_of_snp=nrow(gwasres3)
      currentchr=subset(gwasres3,gwasres3$Chrom==chr)
      SNPcurrentchr=nrow(currentchr)
      print(paste("proportion of SNPs on chr",chr,"=",SNPcurrentchr/number_of_snp))
      ### 
      for (slidwin in c(0:floor(max(currentchr$Position)/window))){  # or 1:ceiling = but here 0 to have something easier with position
        start=slidwin*window+1
        end=slidwin*window+window
        currentwindow=subset(currentchr,currentchr$Position>=start & currentchr$Position<=end)
        SNPwindows=nrow(currentwindow)
        columntosubsample=paste0(trait3,".",model,"_P_value")
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
        dataset3=rbind(dataset3,currentline)
      }
    }
    
  #write.csv(dataset, file=paste0("res_densityassociated_",trait,"_",model,"_slidwin_",format(window,scientific=FALSE),".csv"),row.names=FALSE)

  #### Circlize with density curves # general model
  print(paste("Currently plotting results for model: ",model))
  
  pdf(file=paste0("circlize_TN3yrs_",model,"slidwin",format(window,scientific=FALSE),".pdf"),width=12,height=12)
  library("circlize")
  circos.clear()
  par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
  circos.par("track.height"=0.18,start.degree=90,gap.degree=4,points.overflow.warning=FALSE,cell.padding=c(0, 0, 0, 0),track.margin=c(0.005,0.005))
  circos.initialize(factors = gwasres1$Chrom, x = gwasres1$Position)

  # 1st track density SNP
  circos.trackPlotRegion(factors= dataset1$chr, y=dataset1$Prop_Win_Total, track.height=0.06, track.margin=c(0.045,0.025), bg.border = NA, panel.fun=function(startwin,Prop_Win_Total){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), max(ylim)*1, paste0("Chr0",sector.index), facing = "inside", cex = 1.5,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.6)
  })

  # metaQTL Diana
  circos.rect(sector.index=3,xleft=21605249,ybottom=max(dataset1$Prop_Win_Total)*-0.1,xright=24567362,ytop=max(dataset1$Prop_Win_Total)*0.9,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=3,xleft=34220246,ybottom=max(dataset1$Prop_Win_Total)*-0.1,xright=37772912,ytop=max(dataset1$Prop_Win_Total)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=2414969,ybottom=max(dataset1$Prop_Win_Total)*-0.1,xright=4219224,ytop=max(dataset1$Prop_Win_Total)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=18827666,ybottom=max(dataset1$Prop_Win_Total)*-0.1,xright=24889549,ytop=max(dataset1$Prop_Win_Total)*1.1,lwd=0.01,lty=1,col="pink")

  
  # line per chr
  #circos.segments(sector.index="0",x0=0,y0=max(dataset1$Prop_Win_Total)*2.705,x1=52084813,y1=max(dataset1$Prop_Win_Total)*2.705,lwd=0.3,lty=1,col="#8B2500")
  circos.segments(sector.index=1,x0=0,y0=max(dataset1$Prop_Win_Total)*2.705,x1=64743805,y1=max(dataset1$Prop_Win_Total)*2.705,lwd=0.3,lty=1,col="#8B2500")
  circos.segments(sector.index=2,x0=0,y0=max(dataset1$Prop_Win_Total)*2.705,x1=75091153,y1=max(dataset1$Prop_Win_Total)*2.705,lwd=0.3,lty=1,col="#8B2500")
  circos.segments(sector.index=3,x0=0,y0=max(dataset1$Prop_Win_Total)*2.705,x1=46750189,y1=max(dataset1$Prop_Win_Total)*2.705,lwd=0.3,lty=1,col="#8B2500")
  circos.segments(sector.index=4,x0=0,y0=max(dataset1$Prop_Win_Total)*2.705,x1=58981928,y1=max(dataset1$Prop_Win_Total)*2.705,lwd=0.3,lty=1,col="#8B2500")
  circos.segments(sector.index=5,x0=0,y0=max(dataset1$Prop_Win_Total)*2.705,x1=85844054,y1=max(dataset1$Prop_Win_Total)*2.705,lwd=0.3,lty=1,col="#8B2500")
  circos.segments(sector.index=6,x0=0,y0=max(dataset1$Prop_Win_Total)*2.705,x1=67387004,y1=max(dataset1$Prop_Win_Total)*2.705,lwd=0.3,lty=1,col="#8B2500")
  circos.segments(sector.index=7,x0=0,y0=max(dataset1$Prop_Win_Total)*2.705,x1=67077172,y1=max(dataset1$Prop_Win_Total)*2.705,lwd=0.3,lty=1,col="#8B2500")
  
  #circos.segments(sector.index="0",x0=0,y0=max(dataset1$Prop_Win_Total)*2.505,x1=52084813,y1=max(dataset1$Prop_Win_Total)*2.505,lwd=0.3,lty=1,col="#9B3700")
  circos.segments(sector.index=1,x0=0,y0=max(dataset1$Prop_Win_Total)*2.505,x1=64743805,y1=max(dataset1$Prop_Win_Total)*2.505,lwd=0.3,lty=1,col="#9B3700")
  circos.segments(sector.index=2,x0=0,y0=max(dataset1$Prop_Win_Total)*2.505,x1=75091153,y1=max(dataset1$Prop_Win_Total)*2.505,lwd=0.3,lty=1,col="#9B3700")
  circos.segments(sector.index=3,x0=0,y0=max(dataset1$Prop_Win_Total)*2.505,x1=46750189,y1=max(dataset1$Prop_Win_Total)*2.505,lwd=0.3,lty=1,col="#9B3700")
  circos.segments(sector.index=4,x0=0,y0=max(dataset1$Prop_Win_Total)*2.505,x1=58981928,y1=max(dataset1$Prop_Win_Total)*2.505,lwd=0.3,lty=1,col="#9B3700")
  circos.segments(sector.index=5,x0=0,y0=max(dataset1$Prop_Win_Total)*2.505,x1=85844054,y1=max(dataset1$Prop_Win_Total)*2.505,lwd=0.3,lty=1,col="#9B3700")
  circos.segments(sector.index=6,x0=0,y0=max(dataset1$Prop_Win_Total)*2.505,x1=67387004,y1=max(dataset1$Prop_Win_Total)*2.505,lwd=0.3,lty=1,col="#9B3700")
  circos.segments(sector.index=7,x0=0,y0=max(dataset1$Prop_Win_Total)*2.505,x1=67077172,y1=max(dataset1$Prop_Win_Total)*2.505,lwd=0.3,lty=1,col="#9B3700")
  
  #circos.segments(sector.index="0",x0=0,y0=max(dataset1$Prop_Win_Total)*2.305,x1=52084813,y1=max(dataset1$Prop_Win_Total)*2.305,lwd=0.3,lty=1,col="#AC4900")
  circos.segments(sector.index=1,x0=0,y0=max(dataset1$Prop_Win_Total)*2.305,x1=64743805,y1=max(dataset1$Prop_Win_Total)*2.305,lwd=0.3,lty=1,col="#AC4900")
  circos.segments(sector.index=2,x0=0,y0=max(dataset1$Prop_Win_Total)*2.305,x1=75091153,y1=max(dataset1$Prop_Win_Total)*2.305,lwd=0.3,lty=1,col="#AC4900")
  circos.segments(sector.index=3,x0=0,y0=max(dataset1$Prop_Win_Total)*2.305,x1=46750189,y1=max(dataset1$Prop_Win_Total)*2.305,lwd=0.3,lty=1,col="#AC4900")
  circos.segments(sector.index=4,x0=0,y0=max(dataset1$Prop_Win_Total)*2.305,x1=58981928,y1=max(dataset1$Prop_Win_Total)*2.305,lwd=0.3,lty=1,col="#AC4900")
  circos.segments(sector.index=5,x0=0,y0=max(dataset1$Prop_Win_Total)*2.305,x1=85844054,y1=max(dataset1$Prop_Win_Total)*2.305,lwd=0.3,lty=1,col="#AC4900")
  circos.segments(sector.index=6,x0=0,y0=max(dataset1$Prop_Win_Total)*2.305,x1=67387004,y1=max(dataset1$Prop_Win_Total)*2.305,lwd=0.3,lty=1,col="#AC4900")
  circos.segments(sector.index=7,x0=0,y0=max(dataset1$Prop_Win_Total)*2.305,x1=67077172,y1=max(dataset1$Prop_Win_Total)*2.305,lwd=0.3,lty=1,col="#AC4900")
  
  #circos.segments(sector.index="0",x0=0,y0=max(dataset1$Prop_Win_Total)*2.105,x1=52084813,y1=max(dataset1$Prop_Win_Total)*2.105,lwd=0.3,lty=1,col="#BC5B00")
  circos.segments(sector.index=1,x0=0,y0=max(dataset1$Prop_Win_Total)*2.105,x1=64743805,y1=max(dataset1$Prop_Win_Total)*2.105,lwd=0.3,lty=1,col="#BC5B00")
  circos.segments(sector.index=2,x0=0,y0=max(dataset1$Prop_Win_Total)*2.105,x1=75091153,y1=max(dataset1$Prop_Win_Total)*2.105,lwd=0.3,lty=1,col="#BC5B00")
  circos.segments(sector.index=3,x0=0,y0=max(dataset1$Prop_Win_Total)*2.105,x1=46750189,y1=max(dataset1$Prop_Win_Total)*2.105,lwd=0.3,lty=1,col="#BC5B00")
  circos.segments(sector.index=4,x0=0,y0=max(dataset1$Prop_Win_Total)*2.105,x1=58981928,y1=max(dataset1$Prop_Win_Total)*2.105,lwd=0.3,lty=1,col="#BC5B00")
  circos.segments(sector.index=5,x0=0,y0=max(dataset1$Prop_Win_Total)*2.105,x1=85844054,y1=max(dataset1$Prop_Win_Total)*2.105,lwd=0.3,lty=1,col="#BC5B00")
  circos.segments(sector.index=6,x0=0,y0=max(dataset1$Prop_Win_Total)*2.105,x1=67387004,y1=max(dataset1$Prop_Win_Total)*2.105,lwd=0.3,lty=1,col="#BC5B00")
  circos.segments(sector.index=7,x0=0,y0=max(dataset1$Prop_Win_Total)*2.105,x1=67077172,y1=max(dataset1$Prop_Win_Total)*2.105,lwd=0.3,lty=1,col="#BC5B00")

  #circos.segments(sector.index="0",x0=0,y0=max(dataset1$Prop_Win_Total)*1.905,x1=52084813,y1=max(dataset1$Prop_Win_Total)*1.905,lwd=0.3,lty=1,col="#CD6E00")
  circos.segments(sector.index=1,x0=0,y0=max(dataset1$Prop_Win_Total)*1.905,x1=64743805,y1=max(dataset1$Prop_Win_Total)*1.905,lwd=0.3,lty=1,col="#CD6E00")
  circos.segments(sector.index=2,x0=0,y0=max(dataset1$Prop_Win_Total)*1.905,x1=75091153,y1=max(dataset1$Prop_Win_Total)*1.905,lwd=0.3,lty=1,col="#CD6E00")
  circos.segments(sector.index=3,x0=0,y0=max(dataset1$Prop_Win_Total)*1.905,x1=46750189,y1=max(dataset1$Prop_Win_Total)*1.905,lwd=0.3,lty=1,col="#CD6E00")
  circos.segments(sector.index=4,x0=0,y0=max(dataset1$Prop_Win_Total)*1.905,x1=58981928,y1=max(dataset1$Prop_Win_Total)*1.905,lwd=0.3,lty=1,col="#CD6E00")
  circos.segments(sector.index=5,x0=0,y0=max(dataset1$Prop_Win_Total)*1.905,x1=85844054,y1=max(dataset1$Prop_Win_Total)*1.905,lwd=0.3,lty=1,col="#CD6E00")
  circos.segments(sector.index=6,x0=0,y0=max(dataset1$Prop_Win_Total)*1.905,x1=67387004,y1=max(dataset1$Prop_Win_Total)*1.905,lwd=0.3,lty=1,col="#CD6E00")
  circos.segments(sector.index=7,x0=0,y0=max(dataset1$Prop_Win_Total)*1.905,x1=67077172,y1=max(dataset1$Prop_Win_Total)*1.905,lwd=0.3,lty=1,col="#CD6E00")

  #circos.segments(sector.index="0",x0=0,y0=max(dataset1$Prop_Win_Total)*1.705,x1=52084813,y1=max(dataset1$Prop_Win_Total)*1.705,lwd=0.3,lty=1,col="#DD8000")
  circos.segments(sector.index=1,x0=0,y0=max(dataset1$Prop_Win_Total)*1.705,x1=64743805,y1=max(dataset1$Prop_Win_Total)*1.705,lwd=0.3,lty=1,col="#DD8000")
  circos.segments(sector.index=2,x0=0,y0=max(dataset1$Prop_Win_Total)*1.705,x1=75091153,y1=max(dataset1$Prop_Win_Total)*1.705,lwd=0.3,lty=1,col="#DD8000")
  circos.segments(sector.index=3,x0=0,y0=max(dataset1$Prop_Win_Total)*1.705,x1=46750189,y1=max(dataset1$Prop_Win_Total)*1.705,lwd=0.3,lty=1,col="#DD8000")
  circos.segments(sector.index=4,x0=0,y0=max(dataset1$Prop_Win_Total)*1.705,x1=58981928,y1=max(dataset1$Prop_Win_Total)*1.705,lwd=0.3,lty=1,col="#DD8000")
  circos.segments(sector.index=5,x0=0,y0=max(dataset1$Prop_Win_Total)*1.705,x1=85844054,y1=max(dataset1$Prop_Win_Total)*1.705,lwd=0.3,lty=1,col="#DD8000")
  circos.segments(sector.index=6,x0=0,y0=max(dataset1$Prop_Win_Total)*1.705,x1=67387004,y1=max(dataset1$Prop_Win_Total)*1.705,lwd=0.3,lty=1,col="#DD8000")
  circos.segments(sector.index=7,x0=0,y0=max(dataset1$Prop_Win_Total)*1.705,x1=67077172,y1=max(dataset1$Prop_Win_Total)*1.705,lwd=0.3,lty=1,col="#DD8000")
  
  #circos.segments(sector.index="0",x0=0,y0=max(dataset1$Prop_Win_Total)*1.505,x1=52084813,y1=max(dataset1$Prop_Win_Total)*1.505,lwd=0.3,lty=1,col="#EE9200")
  circos.segments(sector.index=1,x0=0,y0=max(dataset1$Prop_Win_Total)*1.505,x1=64743805,y1=max(dataset1$Prop_Win_Total)*1.505,lwd=0.3,lty=1,col="#EE9200")
  circos.segments(sector.index=2,x0=0,y0=max(dataset1$Prop_Win_Total)*1.505,x1=75091153,y1=max(dataset1$Prop_Win_Total)*1.505,lwd=0.3,lty=1,col="#EE9200")
  circos.segments(sector.index=3,x0=0,y0=max(dataset1$Prop_Win_Total)*1.505,x1=46750189,y1=max(dataset1$Prop_Win_Total)*1.505,lwd=0.3,lty=1,col="#EE9200")
  circos.segments(sector.index=4,x0=0,y0=max(dataset1$Prop_Win_Total)*1.505,x1=58981928,y1=max(dataset1$Prop_Win_Total)*1.505,lwd=0.3,lty=1,col="#EE9200")
  circos.segments(sector.index=5,x0=0,y0=max(dataset1$Prop_Win_Total)*1.505,x1=85844054,y1=max(dataset1$Prop_Win_Total)*1.505,lwd=0.3,lty=1,col="#EE9200")
  circos.segments(sector.index=6,x0=0,y0=max(dataset1$Prop_Win_Total)*1.505,x1=67387004,y1=max(dataset1$Prop_Win_Total)*1.505,lwd=0.3,lty=1,col="#EE9200")
  circos.segments(sector.index=7,x0=0,y0=max(dataset1$Prop_Win_Total)*1.505,x1=67077172,y1=max(dataset1$Prop_Win_Total)*1.505,lwd=0.3,lty=1,col="#EE9200")
  
  #circos.segments(sector.index="0",x0=0,y0=max(dataset1$Prop_Win_Total)*1.305,x1=52084813,y1=max(dataset1$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=1,x0=0,y0=max(dataset1$Prop_Win_Total)*1.305,x1=64743805,y1=max(dataset1$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=2,x0=0,y0=max(dataset1$Prop_Win_Total)*1.305,x1=75091153,y1=max(dataset1$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=3,x0=0,y0=max(dataset1$Prop_Win_Total)*1.305,x1=46750189,y1=max(dataset1$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=4,x0=0,y0=max(dataset1$Prop_Win_Total)*1.305,x1=58981928,y1=max(dataset1$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=5,x0=0,y0=max(dataset1$Prop_Win_Total)*1.305,x1=85844054,y1=max(dataset1$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=6,x0=0,y0=max(dataset1$Prop_Win_Total)*1.305,x1=67387004,y1=max(dataset1$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  circos.segments(sector.index=7,x0=0,y0=max(dataset1$Prop_Win_Total)*1.305,x1=67077172,y1=max(dataset1$Prop_Win_Total)*1.305,lwd=0.3,lty=1,col="#FFA500")
  
  # QTL OW
  #fun_color_range <- colorRampPalette(c("orangered4","orange"))
  #fun_color_range(8)
  #[1] [1] "#8B2500" "#9B3700" "#AC4900" "#BC5B00" "#CD6E00" "#DD8000" "#EE9200" "#FFA500"
  circos.rect(sector.index=1,xleft=29315450,ybottom=max(dataset1$Prop_Win_Total)*2.01,xright=47822809,ytop=max(dataset1$Prop_Win_Total)*2.2,lwd=0.01,lty=1,col="#BC5B00") # lightcoral = 2017
  
  circos.rect(sector.index=2,xleft=155667,ybottom=max(dataset1$Prop_Win_Total)*1.61,xright=2520115,ytop=max(dataset1$Prop_Win_Total)*1.8,lwd=0.01,lty=1,col="#DD8000") #  = 2019
  
  circos.rect(sector.index=3,xleft=8902898,ybottom=max(dataset1$Prop_Win_Total)*2.61,xright=21944360,ytop=max(dataset1$Prop_Win_Total)*2.8,lwd=0.01,lty=1,col="#8B2500") # blue = 2014
  circos.rect(sector.index=3,xleft=21944360,ybottom=max(dataset1$Prop_Win_Total)*2.41,xright=33894541,ytop=max(dataset1$Prop_Win_Total)*2.6,lwd=0.01,lty=1,col="#9B3700") # darkolivegreene4 = 2015
  circos.rect(sector.index=3,xleft=33894541,ybottom=max(dataset1$Prop_Win_Total)*2.21,xright=36727828,ytop=max(dataset1$Prop_Win_Total)*2.4,lwd=0.01,lty=1,col="#AC4900") # darkorchid4 = 2016
  circos.rect(sector.index=3,xleft=21528835,ybottom=max(dataset1$Prop_Win_Total)*2.01,xright=42602509,ytop=max(dataset1$Prop_Win_Total)*2.2,lwd=0.01,lty=1,col="#BC5B00") # lightcoral = 2017
  circos.rect(sector.index=3,xleft=16767733,ybottom=max(dataset1$Prop_Win_Total)*1.81,xright=39576078,ytop=max(dataset1$Prop_Win_Total)*2,lwd=0.01,lty=1,col="#CD6E00") # chocolate4 = 2018
  circos.rect(sector.index=3,xleft=21944395,ybottom=max(dataset1$Prop_Win_Total)*1.61,xright=42602544,ytop=max(dataset1$Prop_Win_Total)*1.8,lwd=0.01,lty=1,col="#DD8000") #  = 2019
  circos.rect(sector.index=3,xleft=21670446,ybottom=max(dataset1$Prop_Win_Total)*1.41,xright=34372538,ytop=max(dataset1$Prop_Win_Total)*1.6,lwd=0.01,lty=1,col="#EE9200") #  = 2020  
  circos.rect(sector.index=3,xleft=24503054,ybottom=max(dataset1$Prop_Win_Total)*1.21,xright=42602544,ytop=max(dataset1$Prop_Win_Total)*1.4,lwd=0.01,lty=1,col="#FFA500") #  = 2021  
  
  circos.rect(sector.index=4,xleft=1013129,ybottom=max(dataset1$Prop_Win_Total)*2.01,xright=58379605,ytop=max(dataset1$Prop_Win_Total)*2.2,lwd=0.01,lty=1,col="#BC5B00") # lightcoral = 2017
  #circos.rect(sector.index=4,xleft=10793638,ybottom=max(dataset1$Prop_Win_Total)*1.21,xright=39431628,ytop=max(dataset1$Prop_Win_Total)*1.4,lwd=0.01,lty=1,col="gray35") # gray35 = mean
  circos.rect(sector.index=4,xleft=712841,ybottom=max(dataset1$Prop_Win_Total)*1.41,xright=50636847,ytop=max(dataset1$Prop_Win_Total)*1.6,lwd=0.01,lty=1,col="#EE9200") #  = 2020  
  
  circos.rect(sector.index=5,xleft=2088736,ybottom=max(dataset1$Prop_Win_Total)*2.21,xright=59728329,ytop=max(dataset1$Prop_Win_Total)*2.4,lwd=0.01,lty=1,col="#AC4900") # darkorchid4 = 2016
  circos.rect(sector.index=5,xleft=1386243,ybottom=max(dataset1$Prop_Win_Total)*2.01,xright=5286233,ytop=max(dataset1$Prop_Win_Total)*2.2,lwd=0.01,lty=1,col="#BC5B00") # lightcoral = 2017
  circos.rect(sector.index=5,xleft=1386243,ybottom=max(dataset1$Prop_Win_Total)*1.81,xright=5286233,ytop=max(dataset1$Prop_Win_Total)*2,lwd=0.01,lty=1,col="#CD6E00") # chocolate4 = 2018
  #circos.rect(sector.index=5,xleft=1386243,ybottom=max(dataset1$Prop_Win_Total)*1.21,xright=3969432,ytop=max(dataset1$Prop_Win_Total)*1.4,lwd=0.01,lty=1,col="gray35") # gray35 = mean
  circos.rect(sector.index=5,xleft=2341457,ybottom=max(dataset1$Prop_Win_Total)*1.61,xright=48861807,ytop=max(dataset1$Prop_Win_Total)*1.8,lwd=0.01,lty=1,col="#DD8000") #  = 2019
  circos.rect(sector.index=5,xleft=2341457,ybottom=max(dataset1$Prop_Win_Total)*1.41,xright=6845198,ytop=max(dataset1$Prop_Win_Total)*1.6,lwd=0.01,lty=1,col="#EE9200") #  = 2020  
  circos.rect(sector.index=5,xleft=146753,ybottom=max(dataset1$Prop_Win_Total)*1.21,xright=32475844,ytop=max(dataset1$Prop_Win_Total)*1.4,lwd=0.01,lty=1,col="#FFA500") #  = 2021  
  
 # position oeak
  #circos.rect(sector.index=1,xleft=NA,ybottom=max(dataset1$Prop_Win_Total)*2.01,xright=NA,ytop=max(dataset1$Prop_Win_Total)*2.2,lwd=0.01,lty=1,col="black") # lightcoral = 2017
  
  circos.rect(sector.index=2,xleft=155667,ybottom=max(dataset1$Prop_Win_Total)*1.61,xright=155667,ytop=max(dataset1$Prop_Win_Total)*1.8,lwd=0.01,lty=1,col="black") #  = 2019
  
  circos.rect(sector.index=3,xleft=16767733,ybottom=max(dataset1$Prop_Win_Total)*2.61,xright=16767733,ytop=max(dataset1$Prop_Win_Total)*2.8,lwd=2,lty=1,col="black") # blue = 2014
  circos.rect(sector.index=3,xleft=25600783,ybottom=max(dataset1$Prop_Win_Total)*2.41,xright=25600783,ytop=max(dataset1$Prop_Win_Total)*2.6,lwd=2,lty=1,col="black") # darkolivegreene4 = 2015
  circos.rect(sector.index=3,xleft=34874101,ybottom=max(dataset1$Prop_Win_Total)*2.21,xright=34874101,ytop=max(dataset1$Prop_Win_Total)*2.4,lwd=2,lty=1,col="black") # darkorchid4 = 2016
  circos.rect(sector.index=3,xleft=33894541,ybottom=max(dataset1$Prop_Win_Total)*2.01,xright=33894541,ytop=max(dataset1$Prop_Win_Total)*2.2,lwd=2,lty=1,col="black") # lightcoral = 2017
  circos.rect(sector.index=3,xleft=33281238,ybottom=max(dataset1$Prop_Win_Total)*1.81,xright=33281238,ytop=max(dataset1$Prop_Win_Total)*2,lwd=2,lty=1,col="black") # chocolate4 = 2018
  circos.rect(sector.index=3,xleft=34874136,ybottom=max(dataset1$Prop_Win_Total)*1.61,xright=34874136,ytop=max(dataset1$Prop_Win_Total)*1.8,lwd=2,lty=1,col="black") #  = 2019
  #circos.rect(sector.index=3,xleft=NA,ybottom=max(dataset1$Prop_Win_Total)*1.41,xright=NA,ytop=max(dataset1$Prop_Win_Total)*1.6,lwd=2,lty=1,col="black") #  = 2020  
  circos.rect(sector.index=3,xleft=32557626,ybottom=max(dataset1$Prop_Win_Total)*1.21,xright=32557626,ytop=max(dataset1$Prop_Win_Total)*1.4,lwd=2,lty=1,col="black") #  = 2021  
  
  circos.rect(sector.index=4,xleft=10793638,ybottom=max(dataset1$Prop_Win_Total)*2.01,xright=10793638,ytop=max(dataset1$Prop_Win_Total)*2.2,lwd=2,lty=1,col="black") # lightcoral = 2017
  ##circos.rect(sector.index=4,xleft=10793638,ybottom=max(dataset1$Prop_Win_Total)*1.21,xright=39431628,ytop=max(dataset1$Prop_Win_Total)*1.4,lwd=0.01,lty=1,col="gray35") # gray35 = mean
  #circos.rect(sector.index=4,xleft=NA,ybottom=max(dataset1$Prop_Win_Total)*1.41,xright=NA,ytop=max(dataset1$Prop_Win_Total)*1.6,lwd=2,lty=1,col="black") #  = 2020  
  
  circos.rect(sector.index=5,xleft=41155983,ybottom=max(dataset1$Prop_Win_Total)*2.21,xright=41155983,ytop=max(dataset1$Prop_Win_Total)*2.4,lwd=2,lty=1,col="black") # darkorchid4 = 2016
  circos.rect(sector.index=5,xleft=2502713,ybottom=max(dataset1$Prop_Win_Total)*2.01,xright=2502713,ytop=max(dataset1$Prop_Win_Total)*2.2,lwd=2,lty=1,col="black") # lightcoral = 2017
  circos.rect(sector.index=5,xleft=2502713,ybottom=max(dataset1$Prop_Win_Total)*1.81,xright=2502713,ytop=max(dataset1$Prop_Win_Total)*2,lwd=2,lty=1,col="black") # chocolate4 = 2018
  #circos.rect(sector.index=5,xleft=1386243,ybottom=max(dataset1$Prop_Win_Total)*1.21,xright=3969432,ytop=max(dataset1$Prop_Win_Total)*1.4,lwd=0.01,lty=1,col="gray35") # gray35 = mean
  circos.rect(sector.index=5,xleft=6845198,ybottom=max(dataset1$Prop_Win_Total)*1.61,xright=6845198,ytop=max(dataset1$Prop_Win_Total)*1.8,lwd=2,lty=1,col="black") #  = 2019
  circos.rect(sector.index=5,xleft=4536400,ybottom=max(dataset1$Prop_Win_Total)*1.41,xright=4536400,ytop=max(dataset1$Prop_Win_Total)*1.6,lwd=2,lty=1,col="black") #  = 2020  
  circos.rect(sector.index=5,xleft=2088771,ybottom=max(dataset1$Prop_Win_Total)*1.21,xright=2088771,ytop=max(dataset1$Prop_Win_Total)*1.4,lwd=2,lty=1,col="black") #  = 2021  
  
  
  for (gene in 1:nrow(candidateRgene)){
    #print(candidateRgene[gene,1])
    circos.segments(sector.index=as.factor(candidateRgene[gene,2]),x0=candidateRgene[gene,3],y0=0,x1=candidateRgene[gene,3],y1=max(dataset1$Prop_Win_Total)*0.85,lwd=0.4,lty=1,col="darkgrey")
  }

  # density gene

  circos.trackLines(densitygene$chr, densitygene$startwin+(densitygene$endwin-densitygene$startwin), densitygene$nb_genes/20000, col="forestgreen", lwd=1.5, type="l")

 # density SNP
  circos.trackLines(dataset1$chr, as.numeric(dataset1$startwin)+window*0.5, dataset1$Prop_Win_Total,lwd=2)


  columntoplot=paste0(trait1,".",model,"_P_value")
  columntofilter=paste0(trait1,".",model,"_P_value_FDR")
  # 2nd track = GWAS general 2014
  circos.trackPlotRegion(factors=  gwasres1$Chrom, y=-log10(gwasres1[[columntoplot]]), track.height=0.11, track.margin=c(0.03,0.01), bg.border = NA, ylim = c(0,max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)+0.3), panel.fun=function(Position,columntoplot){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.8)
  })

  # metaQTL Diana
  circos.rect(sector.index=3,xleft=21605249,ybottom=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*-0.1,xright=24567362,ytop=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=3,xleft=34220246,ybottom=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*-0.1,xright=37772912,ytop=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=2414969,ybottom=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*-0.1,xright=4219224,ytop=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=18827666,ybottom=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*-0.1,xright=24889549,ytop=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  for (gene in 1:nrow(candidateRgene)){
    #print(candidateRgene[gene,1])
    circos.segments(sector.index=as.factor(candidateRgene[gene,2]),x0=candidateRgene[gene,3],y0=0,x1=candidateRgene[gene,3],y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE),lwd=0.4,lty=1,col="darkgrey")
  }
  
  
  # values
  circos.trackPoints(gwasres1$Chrom, gwasres1$Position, -log10(gwasres1[[columntoplot]]),pch=20,col=ifelse(-log10(gwasres1[[columntofilter]])>=-log10(0.05),"firebrick",ifelse(-log10(gwasres1[[columntofilter]])>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres1[[columntofilter]])>=-log10(0.20),"goldenrod2",ifelse(gwasres1$Chrom%%2==0,"navyblue","dodgerblue3")))),cex=ifelse(-log10(gwasres1[[columntofilter]])>=-log10(0.05),2,ifelse(-log10(gwasres1[[columntofilter]])>=-log10(0.10),1.5,ifelse(-log10(gwasres1[[columntofilter]])>=-log10(0.20),1,0.75))))

  # 3rd track = GWAS general 2015
  columntoplot=paste0(trait2,".",model,"_P_value")
  columntofilter=paste0(trait2,".",model,"_P_value_FDR")
  circos.trackPlotRegion(factors=  gwasres2$Chrom, y=-log10(gwasres2[[columntoplot]]), track.height=0.11, track.margin=c(0.03,0.01), bg.border = NA, ylim = c(0,max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE)+0.3), panel.fun=function(Position,columntoplot){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.8)
  })
  
  # metaQTL Diana
  circos.rect(sector.index=3,xleft=21605249,ybottom=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE)*-0.1,xright=24567362,ytop=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=3,xleft=34220246,ybottom=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE)*-0.1,xright=37772912,ytop=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=2414969,ybottom=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE)*-0.1,xright=4219224,ytop=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=18827666,ybottom=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE)*-0.1,xright=24889549,ytop=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  for (gene in 1:nrow(candidateRgene)){
    #print(candidateRgene[gene,1])
    circos.segments(sector.index=as.factor(candidateRgene[gene,2]),x0=candidateRgene[gene,3],y0=0,x1=candidateRgene[gene,3],y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=0.4,lty=1,col="darkgrey")
  }
  
  
  # values
  circos.trackPoints(gwasres2$Chrom, gwasres2$Position, -log10(gwasres2[[columntoplot]]),pch=20,col=ifelse(-log10(gwasres2[[columntofilter]])>=-log10(0.05),"firebrick",ifelse(-log10(gwasres2[[columntofilter]])>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres2[[columntofilter]])>=-log10(0.20),"goldenrod2",ifelse(gwasres2$Chrom%%2==0,"darkolivegreen","darkolivegreen4")))),cex=ifelse(-log10(gwasres2[[columntofilter]])>=-log10(0.05),2,ifelse(-log10(gwasres2[[columntofilter]])>=-log10(0.10),1.5,ifelse(-log10(gwasres2[[columntofilter]])>=-log10(0.20),1,0.75))))
  

  
  # 4rd track = GWAS general 2016
  columntoplot=paste0(trait3,".",model,"_P_value")
  columntofilter=paste0(trait3,".",model,"_P_value_FDR")
  circos.trackPlotRegion(factors=  gwasres3$Chrom, y=-log10(gwasres3[[columntoplot]]), track.height=0.11, track.margin=c(0.03,0.01), bg.border = NA, ylim = c(0,max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)+0.3), panel.fun=function(Position,columntoplot){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    #circos.axis(h=-0.5, direction="inside",labels.facing="inside", lwd=2,labels.cex=0.6)
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.8)
  })
  
  # metaQTL Diana
  circos.rect(sector.index=3,xleft=21605249,ybottom=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)*-0.1,xright=24567362,ytop=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=3,xleft=34220246,ybottom=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)*-0.1,xright=37772912,ytop=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=2414969,ybottom=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)*-0.1,xright=4219224,ytop=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=18827666,ybottom=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)*-0.1,xright=24889549,ytop=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  for (gene in 1:nrow(candidateRgene)){
    #print(candidateRgene[gene,1])
    circos.segments(sector.index=as.factor(candidateRgene[gene,2]),x0=candidateRgene[gene,3],y0=0,x1=candidateRgene[gene,3],y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=0.4,lty=1,col="darkgrey")
  }
  
  
  # values
  circos.trackPoints(gwasres3$Chrom, gwasres3$Position, -log10(gwasres3[[columntoplot]]),pch=20,col=ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.05),"firebrick",ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.20),"goldenrod2",ifelse(gwasres3$Chrom%%2==0,"darkmagenta","darkorchid1")))),cex=ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.05),2,ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.10),1.5,ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.20),1,0.75))))
  
  
  
  # p < 0.01
  circos.trackPlotRegion(factors= dataset3$chr, y=dataset3$PropSNPWinPvalueInf01, track.height=0.11, track.margin=c(0.03,0.01), bg.border = NA, panel.fun=function(startwin,PropSNPWinPvalueInf01){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.6)
  })
  #
  # metaQTL Diana
  circos.rect(sector.index=3,xleft=21605249,ybottom=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*-0.1,xright=24567362,ytop=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=3,xleft=34220246,ybottom=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*-0.1,xright=37772912,ytop=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=2414969,ybottom=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*-0.1,xright=4219224,ytop=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  circos.rect(sector.index=5,xleft=18827666,ybottom=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*-0.1,xright=24889549,ytop=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  for (gene in 1:nrow(candidateRgene)){
    #print(candidateRgene[gene,1])
    circos.segments(sector.index=as.factor(candidateRgene[gene,2]),x0=candidateRgene[gene,3],y0=0,x1=candidateRgene[gene,3],y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=0.4,lty=1,col="darkgrey")
  }
  
  circos.trackLines(dataset1$chr, dataset1$startwin+window*0.5, dataset1$PropSNPWinPvalueInf01,lwd=1,col=ifelse(dataset1$chr%%2==0,"navyblue","dodgerblue3"))
  circos.trackLines(dataset2$chr, dataset2$startwin+window*0.5, dataset2$PropSNPWinPvalueInf01,lwd=1,col=ifelse(dataset2$chr%%2==0,"darkolivegreen","darkolivegreen4"))
  circos.trackLines(dataset3$chr, dataset3$startwin+window*0.5, dataset3$PropSNPWinPvalueInf01,lwd=1,col=ifelse(dataset3$chr%%2==0,"darkmagenta","darkorchid1"))
  
  # track numbers
  points(x=-0.015, y = 0.98 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.98, labels="1",cex=1.1,col="white")
  
  points(x=-0.015, y = 0.88 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.88, labels="2",cex=1.1,col="white")
  
  points(x=-0.015, y = 0.725 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.725, labels="3",cex=1.1,col="white") 
  
  points(x=-0.015, y = 0.575 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.575, labels="4",cex=1.1,col="white") 
  
  points(x=-0.015, y = 0.430 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.430, labels="5",cex=1.1,col="white") 
  
  # internal information
  
  text(x=0, y = 0.26, labels="KEY",cex=1.2,font=2)
  
  segments(x0=-0.18,y0=0.22,x1=-0.18,y1=0.18,col="darkgrey",lwd=2)
  segments(x0=-0.165,y0=0.21,x1=-0.135,y1=0.21,col="black",lwd=2)
  segments(x0=-0.165,y0=0.19,x1=-0.135,y1=0.19,col="forestgreen",lwd=1.5)
  text(x=-0.13, y = 0.208, labels=expression("SNP density (MAF" >= "0.05)"),cex=0.7,col="black",pos=4)
  text(x=-0.13, y = 0.188, labels="Gene density",cex=0.7,col="forestgreen",pos=4)
  points(x=-0.20, y = 0.20 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = 0.20, labels="1",cex=1,col="white")
  
  segments(x0=-0.18,y0=0.18,x1=-0.18,y1=0.14,col="darkgrey",lwd=2)
  points(x=-0.15, y = 0.17 , type = "p", cex=2,col="navyblue",pch=20)
  points(x=-0.15, y = 0.15 , type = "p", cex=2,col="dodgerblue3",pch=20)
  text(x=-0.13, y = 0.16, labels="GWAS Black spot disease severity 2014",cex=0.7,col="blue",pos=4)
  points(x=-0.20, y = 0.16 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = 0.16, labels="2",cex=1,col="white")
  
  
  segments(x0=-0.18,y0=0.13,x1=-0.18,y1=0.09,col="darkgrey",lwd=2)
  points(x=-0.15, y = 0.12 , type = "p", cex=2,col="darkolivegreen",pch=20)
  points(x=-0.15, y = 0.10 , type = "p", cex=2,col="darkolivegreen4",pch=20)
  text(x=-0.13, y = 0.11, labels="GWAS Black spot disease severity 2015",cex=0.7,col="darkolivegreen4",pos=4)
  points(x=-0.20, y = 0.11 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = 0.11, labels="3",cex=1,col="white")
  
  
  segments(x0=-0.18,y0=0.08,x1=-0.18,y1=0.04,col="darkgrey",lwd=2)
  points(x=-0.15, y = 0.07 , type = "p", cex=2,col="darkmagenta",pch=20)
  points(x=-0.15, y = 0.05 , type = "p", cex=2,col="darkorchid1",pch=20)
  text(x=-0.13, y = 0.06, labels="GWAS Black spot disease severity 2016",cex=0.7,col="darkorchid4",pos=4)
  points(x=-0.20, y = 0.06 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = 0.06, labels="4",cex=1,col="white")
  
  segments(x0=-0.18,y0=0.03,x1=-0.18,y1=-0.03,col="darkgrey",lwd=2)
  points(x=-0.15, y = 0.02 , type = "p", cex=2,col="goldenrod2",pch=20)
  text(x=-0.13, y = 0.015, labels=expression("SNP with q-value"<= "0.20"),cex=0.7,col="goldenrod2",pos=4)
  points(x=-0.15, y = 0 , type = "p", cex=2,col="darkorange2",pch=20)
  text(x=-0.13, y = -0.005, labels=expression("SNP with q-value"<= "0.10"),cex=0.7,col="darkorange2",pos=4)
  points(x=-0.15, y = -0.02 , type = "p", cex=2,col="firebrick",pch=20)
  text(x=-0.13, y = -0.025, labels=expression("SNP with q-value"<= "0.05"),cex=0.7,col="firebrick",pos=4)
  points(x=-0.26, y = 0 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.26, y = 0, labels="2",cex=1,col="white")
  points(x=-0.23, y = 0, type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.23, y = 0, labels="3",cex=1,col="white")
  points(x=-0.20, y = 0, type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = 0, labels="4",cex=1,col="white")
  
  segments(x0=-0.18,y0=-0.04,x1=-0.18,y1=-0.13,col="darkgrey",lwd=2)
  segments(x0=-0.165,y0=-0.05,x1=-0.135,y1=-0.05,col="navyblue",lwd=1)
  segments(x0=-0.165,y0=-0.06,x1=-0.135,y1=-0.06,col="dodgerblue3",lwd=1)
  text(x=-0.13, y = -0.058, labels=expression("Density of SNPs with p-value"<="0.01 (2014)"),cex=0.7,col="blue",pos=4)
  
  segments(x0=-0.165,y0=-0.08,x1=-0.135,y1=-0.08,col="darkolivegreen",lwd=1)
  segments(x0=-0.165,y0=-0.09,x1=-0.135,y1=-0.09,col="darkolivegreen4",lwd=1)
  text(x=-0.13, y = -0.088, labels=expression("Density of SNP with p-value"<="0.01 (2015)"),cex=0.7,col="darkolivegreen4",pos=4)
  
  segments(x0=-0.165,y0=-0.11,x1=-0.135,y1=-0.11,col="darkmagenta",lwd=1)
  segments(x0=-0.165,y0=-0.12,x1=-0.135,y1=-0.12,col="darkorchid1",lwd=1)
  text(x=-0.13, y = -0.118, labels=expression("Density of SNP with p-value"<="0.01 (2016)"),cex=0.7,col="darkorchid4",pos=4)

  points(x=-0.20, y = -0.085, type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = -0.085, labels="5",cex=1,col="white")
  
  segments(x0=-0.165,y0=-0.142,x1=-0.135,y1=-0.142,col="#8B2500",lwd=1)
  segments(x0=-0.165,y0=-0.146,x1=-0.135,y1=-0.146,col="#9B3700",lwd=1)
  segments(x0=-0.165,y0=-0.150,x1=-0.135,y1=-0.150,col="#AC4900",lwd=1)
  segments(x0=-0.165,y0=-0.154,x1=-0.135,y1=-0.154,col="#BC5B00",lwd=1)
  segments(x0=-0.165,y0=-0.158,x1=-0.135,y1=-0.158,col="#CD6E00",lwd=1)
  segments(x0=-0.165,y0=-0.162,x1=-0.135,y1=-0.162,col="#DD8000",lwd=1)
  segments(x0=-0.165,y0=-0.166,x1=-0.135,y1=-0.166,col="#EE9200",lwd=1)
  segments(x0=-0.165,y0=-0.170,x1=-0.135,y1=-0.170,col="#FFA500",lwd=1)
  rect(xleft=-0.158,ybottom=-0.144,xright=-0.142, ytop=-0.14,col="#8B2500",border="black",lwd=0.01)
  rect(xleft=-0.158,ybottom=-0.148,xright=-0.142, ytop=-0.144,col="#9B3700",border="black",lwd=0.01)
  rect(xleft=-0.158,ybottom=-0.152,xright=-0.142, ytop=-0.148,col="#AC4900",border="black",lwd=0.01)
  rect(xleft=-0.158,ybottom=-0.156,xright=-0.142, ytop=-0.152,col="#BC5B00",border="black",lwd=0.01)
  rect(xleft=-0.158,ybottom=-0.160,xright=-0.142, ytop=-0.156,col="#CD6E00",border="black",lwd=0.01)
  rect(xleft=-0.158,ybottom=-0.164,xright=-0.142, ytop=-0.160,col="#DD8000",border="black",lwd=0.01)
  rect(xleft=-0.158,ybottom=-0.168,xright=-0.142, ytop=-0.164,col="#EE9200",border="black",lwd=0.01)
  rect(xleft=-0.158,ybottom=-0.172,xright=-0.142, ytop=-0.168,col="#FFA500",border="black",lwd=0.01)
  text(x=-0.205, y = -0.143, labels=expression("2014"),cex=0.45,col="#8B2500",pos=4)
  text(x=-0.205, y = -0.171, labels=expression("2021"),cex=0.45,col="#FFA500",pos=4)
  text(x=-0.13, y = -0.155, labels=expression("QTL (Old Blush * Rosa x wichuriana only)"),cex=0.7,col="black",pos=4)
  
  rect(xleft=-0.165,ybottom=-0.18,xright=-0.135, ytop=-0.195,col="pink",border="black")
  text(x=-0.13, y = -0.191, labels=expression("Meta-QTL from Lopez Arias et al. 2020"),cex=0.7,col="black",pos=4)
  
  segments(x0=-0.165,y0=-0.215,x1=-0.135,y1=-0.215,col="darkgrey",lwd=1)
  text(x=-0.13, y = -0.22, labels=expression("Candidate R genes"),cex=0.7,col="black",pos=4)
  
 dev.off() 
}


