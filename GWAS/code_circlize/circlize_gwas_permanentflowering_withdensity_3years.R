## TL - 161121

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/")
densitygene=read.table("densitygene_slidwin500kb.txt",h=T)

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/permanentflowering/")
#gwasres=read.table("TN2014_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="TN2014"

#trait="Nb_Clust"
#gwasres1=read.table("./A1/Nb_Clust_A1_gwas_results/Nb_Clust_A1_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait1=paste0(trait,"_A1")
#gwasres2=read.table("./A2/Nb_Clust_A2_gwas_results/Nb_Clust_A2_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait2=paste0(trait,"_A2")
#gwasres3=read.table("./A3/Nb_Clust_A3_gwas_results/Nb_Clust_A3_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait3=paste0(trait,"_A3")

#trait="Reb_Mag"
#gwasres1=read.table("./A1/Reb_Mag_A1_gwas_results/Reb_Mag_A1_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait1=paste0(trait,"_A1")
#gwasres2=read.table("./A2/Reb_Mag_A2_gwas_results/Reb_Mag_A2_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait2=paste0(trait,"_A2")
#gwasres3=read.table("./A3/Reb_Mag_A3_gwas_results/Reb_Mag_A3_scores.txt.withoutchrNA.sed",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait3=paste0(trait,"_A3")

#trait="First_Reb"
#trait="Lmax_P"
#trait="Rat_P"
#trait="Min_M"
#trait="DifMax_M"
#trait="Min_V"
#trait="Area_Peak_Inflo"
#trait="Max_Peak"
#trait="Max_Reb"
#trait="Max_Reb"
#trait="Rat_Peak"
#trait="First_Flo"
trait="Cont_Flo"
trait1=paste0(trait,"_A1")
gwasres1=read.table(paste0("./A1/",trait1,"_gwas_results/",trait1,"_scores.txt.withoutchrNA.sed"),h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait2=paste0(trait,"_A2")
gwasres2=read.table(paste0("./A2/",trait2,"_gwas_results/",trait2,"_scores.txt.withoutchrNA.sed"),h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait3=paste0(trait,"_A3")
gwasres3=read.table(paste0("./A3/",trait3,"_gwas_results/",trait3,"_scores.txt.withoutchrNA.sed"),h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE



# for all models gwasres1
for (model in c("additive","diplo.additive","general","diplo.general","1.dom.alt","1.dom.ref","2.dom.alt","2.dom.ref")){
  print(paste("Currently working on model: ",model))
  # year 1 (A1)
  ### slidingwindows R
  ## fist loop by chromosome
  dataset1 <- data.frame(matrix(ncol = 12, nrow = 0))
  colnames(dataset1)=c("chr","startwin","endwin","Nb_SNPs_Win","Prop_Win_Total","NbSNPWinPvalueInf05","NbSNPWinPvalueInf01","NbSNPWinPvalueInf001","SNPswithoutpvalues","PropSNPWinPvalueInf05","PropSNPWinPvalueInf01","PropSNPWinPvalueInf001")
  window=1000000
  for (chr in c(0:7)){
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
  # year2 (A2)
    dataset2 <- data.frame(matrix(ncol = 12, nrow = 0))
    colnames(dataset2)=c("chr","startwin","endwin","Nb_SNPs_Win","Prop_Win_Total","NbSNPWinPvalueInf05","NbSNPWinPvalueInf01","NbSNPWinPvalueInf001","SNPswithoutpvalues","PropSNPWinPvalueInf05","PropSNPWinPvalueInf01","PropSNPWinPvalueInf001")
    window=1000000
    for (chr in c(0:7)){
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
    # year3 (A3)
    dataset3 <- data.frame(matrix(ncol = 12, nrow = 0))
    colnames(dataset3)=c("chr","startwin","endwin","Nb_SNPs_Win","Prop_Win_Total","NbSNPWinPvalueInf05","NbSNPWinPvalueInf01","NbSNPWinPvalueInf001","SNPswithoutpvalues","PropSNPWinPvalueInf05","PropSNPWinPvalueInf01","PropSNPWinPvalueInf001")
    window=1000000
    for (chr in c(0:7)){
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
  
  pdf(file=paste0("circlize_PF_",trait,"_",model,"slidwin","_",format(window,scientific=FALSE),".pdf"),width=12,height=12)
  library("circlize")
  circos.clear()
  par(mar=c(1,1,1,1), lwd=0.1, cex =0.7)
  circos.par("track.height"=0.18,start.degree=90,gap.degree=4,points.overflow.warning=FALSE,cell.padding=c(0, 0, 0, 0),track.margin=c(0.005,0.005))
  circos.initialize(factors = gwasres1$Chrom, x = gwasres1$Position)

  # 1st track density SNP
  circos.trackPlotRegion(factors= dataset1$chr, y=dataset1$Prop_Win_Total, track.height=0.06, track.margin=c(0.03,0.01), bg.border = NA, panel.fun=function(startwin,Prop_Win_Total){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.6)
  })

  # roQSN region
  circos.rect(sector.index=3,xleft=28500000,ybottom=max(dataset1$Prop_Win_Total)*-0.1,xright=33000000,ytop=max(dataset1$Prop_Win_Total)*0.9,lwd=0.01,lty=1,col="pink")

  #ABCDE
  circos.segments(sector.index=1,x0=41514796,y0=0,x1=41514796,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52096899,y0=0,x1=52096899,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52554849,y0=0,x1=52554849,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11833234,y0=0,x1=11833234,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=17461203,y0=0,x1=17461203,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57253197,y0=0,x1=57253197,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57254457,y0=0,x1=57254457,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67016526,y0=0,x1=67016526,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67036489,y0=0,x1=67036489,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33238136,y0=0,x1=33238136,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46474063,y0=0,x1=46474063,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8425314,y0=0,x1=8425314,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8430126,y0=0,x1=8430126,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=61446161,y0=0,x1=61446161,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=58611899,y0=0,x1=58611899,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=60224075,y0=0,x1=60224075,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4414048,y0=0,x1=4414048,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4468385,y0=0,x1=4468385,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=50555186,y0=0,x1=50555186,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=52218507,y0=0,x1=52218507,y1=max(dataset1$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.rect(sector.index=3,xleft=41280000,ybottom=0,xright=41990000,ytop=max(dataset1$Prop_Win_Total,na.rm=TRUE),lwd=0.01,lty=1,col="purple",border=0)
  
  
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

  # roQSN region
  circos.rect(sector.index=3,xleft=28500000,ybottom=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*-0.1,xright=33000000,ytop=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  
  # ABCDE
  circos.segments(sector.index=1,x0=41514796,y0=0,x1=41514796,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52096899,y0=0,x1=52096899,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52554849,y0=0,x1=52554849,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11833234,y0=0,x1=11833234,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=17461203,y0=0,x1=17461203,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57253197,y0=0,x1=57253197,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57254457,y0=0,x1=57254457,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67016526,y0=0,x1=67016526,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67036489,y0=0,x1=67036489,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33238136,y0=0,x1=33238136,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46474063,y0=0,x1=46474063,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8425314,y0=0,x1=8425314,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8430126,y0=0,x1=8430126,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=61446161,y0=0,x1=61446161,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=58611899,y0=0,x1=58611899,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=60224075,y0=0,x1=60224075,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4414048,y0=0,x1=4414048,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4468385,y0=0,x1=4468385,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=50555186,y0=0,x1=50555186,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=52218507,y0=0,x1=52218507,y1=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.7,lwd=2,lty=1,col="purple")
  circos.rect(sector.index=3,xleft=41280000,ybottom=0,xright=41990000,ytop=max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.8,lwd=0.01,lty=1,col="purple",border=0)
  
  circos.text(sector.index=1, 41514796, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "AGL6",cex=0.75,col="purple")
  circos.text(sector.index=1, 52096899, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.9, "TOE2",cex=0.75,col="purple")
  circos.text(sector.index=1, 52554849, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "euAP1",cex=0.75,col="purple")  
  circos.text(sector.index=2, 11833234, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "TM6",cex=0.75,col="purple")  
  circos.text(sector.index=2, 17461203, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.9, "AP2",cex=0.75,col="purple")
  circos.text(sector.index=2, 57253197, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "PLE",cex=0.75,col="purple") 
  circos.text(sector.index=2, 57254457, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.9, "Agamous D",cex=0.75,col="purple")
  circos.text(sector.index=2, 67016526, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "euFUL",cex=0.75,col="purple") 
  circos.text(sector.index=2, 67036489, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.9, "SEP1/2",cex=0.75,col="purple") 
  circos.text(sector.index=3, 33238136, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "AP2/TOE (NP)",cex=0.75,col="purple")
  circos.text(sector.index=4, 46474063, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "SEP3",cex=0.75,col="purple")  
  circos.text(sector.index=5, 8425314, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "AG",cex=0.75,col="purple") 
  circos.text(sector.index=5, 8430126, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.9, "Agamous C",cex=0.75,col="purple")
  circos.text(sector.index=5, 61446161, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "TOE1",cex=0.75,col="purple")
  circos.text(sector.index=6, 58611899, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "PI/GLO",cex=0.75,col="purple")  
  circos.text(sector.index=6, 60224075, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.9, "AP3",cex=0.75,col="purple")  
  circos.text(sector.index=7, 4414048, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "SEP4",cex=0.75,col="purple") 
  circos.text(sector.index=7, 4468385, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.9, "euFUL",cex=0.75,col="purple")  
  circos.text(sector.index=7, 50555186, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE), "D",cex=0.75,col="purple")  
  circos.text(sector.index=7, 52218507, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.9, "AGL12",cex=0.75,col="purple")   
  circos.text(sector.index=3, 41500000, max(-log10(gwasres1[[columntoplot]]),na.rm=TRUE)*0.9, "S locus",cex=0.75,col="purple") 
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
  
  # roQSN region
  circos.rect(sector.index=3,xleft=28500000,ybottom=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE)*-0.1,xright=33000000,ytop=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  # ABCDE
  circos.segments(sector.index=1,x0=41514796,y0=0,x1=41514796,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52096899,y0=0,x1=52096899,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52554849,y0=0,x1=52554849,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11833234,y0=0,x1=11833234,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=17461203,y0=0,x1=17461203,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57253197,y0=0,x1=57253197,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57254457,y0=0,x1=57254457,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67016526,y0=0,x1=67016526,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67036489,y0=0,x1=67036489,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33238136,y0=0,x1=33238136,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46474063,y0=0,x1=46474063,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8425314,y0=0,x1=8425314,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8430126,y0=0,x1=8430126,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=61446161,y0=0,x1=61446161,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=58611899,y0=0,x1=58611899,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=60224075,y0=0,x1=60224075,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4414048,y0=0,x1=4414048,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4468385,y0=0,x1=4468385,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=50555186,y0=0,x1=50555186,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=52218507,y0=0,x1=52218507,y1=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.rect(sector.index=3,xleft=41280000,ybottom=0,xright=41990000,ytop=max(-log10(gwasres2[[columntoplot]]),na.rm=TRUE),lwd=0.01,lty=1,col="purple",border=0)
  
  
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
  
  # roQSN region
  circos.rect(sector.index=3,xleft=28500000,ybottom=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)*-0.1,xright=33000000,ytop=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  # ABCDE
  circos.segments(sector.index=1,x0=41514796,y0=0,x1=41514796,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52096899,y0=0,x1=52096899,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52554849,y0=0,x1=52554849,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11833234,y0=0,x1=11833234,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=17461203,y0=0,x1=17461203,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57253197,y0=0,x1=57253197,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57254457,y0=0,x1=57254457,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67016526,y0=0,x1=67016526,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67036489,y0=0,x1=67036489,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33238136,y0=0,x1=33238136,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46474063,y0=0,x1=46474063,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8425314,y0=0,x1=8425314,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8430126,y0=0,x1=8430126,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=61446161,y0=0,x1=61446161,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=58611899,y0=0,x1=58611899,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=60224075,y0=0,x1=60224075,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4414048,y0=0,x1=4414048,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4468385,y0=0,x1=4468385,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=50555186,y0=0,x1=50555186,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=52218507,y0=0,x1=52218507,y1=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.rect(sector.index=3,xleft=41280000,ybottom=0,xright=41990000,ytop=max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE),lwd=0.01,lty=1,col="purple",border=0)
  
  
  # values
  circos.trackPoints(gwasres3$Chrom, gwasres3$Position, -log10(gwasres3[[columntoplot]]),pch=20,col=ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.05),"firebrick",ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.10),"darkorange2",ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.20),"goldenrod2",ifelse(gwasres3$Chrom%%2==0,"darkmagenta","darkorchid1")))),cex=ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.05),2,ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.10),1.5,ifelse(-log10(gwasres3[[columntofilter]])>=-log10(0.20),1,0.75))))
  
  
  
  # p < 0.01
  # select the track exhiting the highest value (to use it for the y scale)
  
  # then plot
  circos.trackPlotRegion(factors= dataset3$chr, y=dataset3$PropSNPWinPvalueInf01, ylim = c(0,max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01)), track.height=0.11, track.margin=c(0.03,0.01), bg.border = NA, panel.fun=function(startwin,PropSNPWinPvalueInf01){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.6)
  })
  #
  # roQSN region
  circos.rect(sector.index=3,xleft=28500000,ybottom=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*-0.1,xright=33000000,ytop=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  # ABCDE
  circos.segments(sector.index=1,x0=41514796,y0=0,x1=41514796,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52096899,y0=0,x1=52096899,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52554849,y0=0,x1=52554849,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11833234,y0=0,x1=11833234,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=17461203,y0=0,x1=17461203,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57253197,y0=0,x1=57253197,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57254457,y0=0,x1=57254457,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67016526,y0=0,x1=67016526,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67036489,y0=0,x1=67036489,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33238136,y0=0,x1=33238136,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46474063,y0=0,x1=46474063,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8425314,y0=0,x1=8425314,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8430126,y0=0,x1=8430126,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=61446161,y0=0,x1=61446161,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=58611899,y0=0,x1=58611899,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=60224075,y0=0,x1=60224075,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4414048,y0=0,x1=4414048,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4468385,y0=0,x1=4468385,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=50555186,y0=0,x1=50555186,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=52218507,y0=0,x1=52218507,y1=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.rect(sector.index=3,xleft=41280000,ybottom=0,xright=41990000,ytop=max(dataset1$PropSNPWinPvalueInf01,dataset2$PropSNPWinPvalueInf01,dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=0.01,lty=1,col="purple",border=0)
  
  circos.trackLines(dataset1$chr, dataset1$startwin+window*0.5, dataset1$PropSNPWinPvalueInf01,lwd=1,col=ifelse(dataset1$chr%%2==0,"navyblue","dodgerblue3"))
  circos.trackLines(dataset2$chr, dataset2$startwin+window*0.5, dataset2$PropSNPWinPvalueInf01,lwd=1,col=ifelse(dataset2$chr%%2==0,"darkolivegreen","darkolivegreen4"))
  circos.trackLines(dataset3$chr, dataset3$startwin+window*0.5, dataset3$PropSNPWinPvalueInf01,lwd=1,col=ifelse(dataset3$chr%%2==0,"darkmagenta","darkorchid1"))
  
  # track numbers
  points(x=-0.015, y = 1.01 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 1.01, labels="1",cex=1.1,col="white")
  
  points(x=-0.015, y = 0.91 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.91, labels="2",cex=1.1,col="white")
  
  points(x=-0.015, y = 0.755 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.755, labels="3",cex=1.1,col="white") 
  
  points(x=-0.015, y = 0.605 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.605, labels="4",cex=1.1,col="white") 
  
  points(x=-0.015, y = 0.46 , type = "p", cex=4,col="black",pch=20)
  text(x=-0.015, y = 0.46, labels="5",cex=1.1,col="white") 
  
  # internal information
  
  text(x=0, y = 0.275, labels="KEY",cex=1.2,font=2)
  
  segments(x0=-0.18,y0=0.23,x1=-0.18,y1=0.19,col="darkgrey",lwd=2)
  segments(x0=-0.165,y0=0.22,x1=-0.135,y1=0.22,col="black",lwd=2)
  segments(x0=-0.165,y0=0.20,x1=-0.135,y1=0.20,col="forestgreen",lwd=1.5)
  text(x=-0.13, y = 0.218, labels=expression("SNP density (MAF" >= "0.05)"),cex=0.7,col="black",pos=4)
  text(x=-0.13, y = 0.198, labels="Gene density",cex=0.7,col="forestgreen",pos=4)
  points(x=-0.20, y = 0.21 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = 0.21, labels="1",cex=1,col="white")
  
  segments(x0=-0.18,y0=0.17,x1=-0.18,y1=0.13,col="darkgrey",lwd=2)
  points(x=-0.15, y = 0.16 , type = "p", cex=2,col="navyblue",pch=20)
  points(x=-0.15, y = 0.14 , type = "p", cex=2,col="dodgerblue3",pch=20)
  text(x=-0.13, y = 0.15, labels=paste0("GWAS ",trait," Year 1 (A1)"),cex=0.7,col="blue",pos=4)
  points(x=-0.20, y = 0.15 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = 0.15, labels="2",cex=1,col="white")
  
  
  segments(x0=-0.18,y0=0.11,x1=-0.18,y1=0.07,col="darkgrey",lwd=2)
  points(x=-0.15, y = 0.10 , type = "p", cex=2,col="darkolivegreen",pch=20)
  points(x=-0.15, y = 0.08 , type = "p", cex=2,col="darkolivegreen4",pch=20)
  text(x=-0.13, y = 0.09, labels=paste0("GWAS ",trait," Year 2 (A2)"),cex=0.7,col="darkolivegreen4",pos=4)
  points(x=-0.20, y = 0.09 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = 0.09, labels="3",cex=1,col="white")
  
  
  segments(x0=-0.18,y0=0.05,x1=-0.18,y1=0.01,col="darkgrey",lwd=2)
  points(x=-0.15, y = 0.04 , type = "p", cex=2,col="darkmagenta",pch=20)
  points(x=-0.15, y = 0.02 , type = "p", cex=2,col="darkorchid1",pch=20)
  text(x=-0.13, y = 0.03, labels=paste0("GWAS ",trait," Year 3 (A3)"),cex=0.7,col="darkorchid4",pos=4)
  points(x=-0.20, y = 0.03 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = 0.03, labels="4",cex=1,col="white")
  
  segments(x0=-0.18,y0=-0.01,x1=-0.18,y1=-0.07,col="darkgrey",lwd=2)
  points(x=-0.15, y = -0.02 , type = "p", cex=2,col="goldenrod2",pch=20)
  text(x=-0.13, y = -0.025, labels=expression("SNP with q-value"<= "0.20"),cex=0.7,col="goldenrod2",pos=4)
  points(x=-0.15, y = -0.04 , type = "p", cex=2,col="darkorange2",pch=20)
  text(x=-0.13, y = -0.045, labels=expression("SNP with q-value"<= "0.10"),cex=0.7,col="darkorange2",pos=4)
  points(x=-0.15, y = -0.06 , type = "p", cex=2,col="firebrick",pch=20)
  text(x=-0.13, y = -0.065, labels=expression("SNP with q-value"<= "0.05"),cex=0.7,col="firebrick",pos=4)
  points(x=-0.26, y = -0.04 , type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.26, y = -0.04, labels="2",cex=1,col="white")
  points(x=-0.23, y = -0.04, type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.23, y = -0.04, labels="3",cex=1,col="white")
  points(x=-0.20, y = -0.04, type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = -0.04, labels="4",cex=1,col="white")
  
  segments(x0=-0.18,y0=-0.09,x1=-0.18,y1=-0.18,col="darkgrey",lwd=2)
  segments(x0=-0.165,y0=-0.10,x1=-0.135,y1=-0.10,col="navyblue",lwd=1)
  segments(x0=-0.165,y0=-0.11,x1=-0.135,y1=-0.11,col="dodgerblue3",lwd=1)
  text(x=-0.13, y = -0.108, labels=expression("Density of SNPs with p-value"<="0.01 (Year 1)"),cex=0.7,col="blue",pos=4)
  
  segments(x0=-0.165,y0=-0.13,x1=-0.135,y1=-0.13,col="darkolivegreen",lwd=1)
  segments(x0=-0.165,y0=-0.14,x1=-0.135,y1=-0.14,col="darkolivegreen4",lwd=1)
  text(x=-0.13, y = -0.138, labels=expression("Density of SNP with p-value"<="0.01 (Year 2)"),cex=0.7,col="darkolivegreen4",pos=4)
  
  segments(x0=-0.165,y0=-0.16,x1=-0.135,y1=-0.16,col="darkmagenta",lwd=1)
  segments(x0=-0.165,y0=-0.17,x1=-0.135,y1=-0.17,col="darkorchid1",lwd=1)
  text(x=-0.13, y = -0.168, labels=expression("Density of SNP with p-value"<="0.01 (Year 3)"),cex=0.7,col="darkorchid4",pos=4)

  points(x=-0.20, y = -0.135, type = "p", cex=3.5,col="black",pch=20)
  text(x=-0.20, y = -0.135, labels="5",cex=1,col="white")
  
  rect(xleft=-0.165,ybottom=-0.19,xright=-0.135, ytop=-0.205,col="pink",border="black")
  text(x=-0.13, y = -0.201, labels=expression("roQSN region (28.5 - 33 Mb)"),cex=0.7,col="black",pos=4)

  
 dev.off() 
}


