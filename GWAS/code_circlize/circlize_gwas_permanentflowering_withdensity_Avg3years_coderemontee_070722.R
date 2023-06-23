## TL - 161121

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/")
densitygene=read.table("densitygene_slidwin500kb.txt",h=T)

#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/permanentflowering/")
#gwasres3=read.table("./codenumremonteeAvg3years_gwas_results/codenumremonteeAvg3years_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="codenumremonteeAvg3years"
#gwasres3=read.table("./codenumremonteeHarmonise3years_gwas_results/codenumremonteeHarmonise3years_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait="codenumremonteeHarmonise3years"
#trait3=paste0(trait,"")
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/")
gwasres3=read.table("./year_obtained_gwas_results/year_obtained_scores.sed.withoutchrNA",h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
trait="year_obtained"
trait3=paste0(trait,"")

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
#trait="Cont_Flo"
#trait1=paste0(trait,"_A1")
#gwasres1=read.table(paste0("./A1/",trait1,"_gwas_results/",trait1,"_scores.txt.withoutchrNA.sed"),h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait2=paste0(trait,"_A2")
#gwasres2=read.table(paste0("./A2/",trait2,"_gwas_results/",trait2,"_scores.txt.withoutchrNA.sed"),h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE
#trait3=paste0(trait,"_A3")
#gwasres3=read.table(paste0("./A3/",trait3,"_gwas_results/",trait3,"_scores.txt.withoutchrNA.sed"),h=T) #awk '$2 != "NA" {print $0}' INFILE > OUTFILE



# for all models gwasres1
for (model in c("additive","diplo.additive","general","diplo.general","1.dom.alt","1.dom.ref","2.dom.alt","2.dom.ref")){
  print(paste("Currently working on model: ",model))
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
  circos.initialize(factors = gwasres3$Chrom, x = gwasres3$Position)

  # 1st track density SNP
  circos.trackPlotRegion(factors= dataset3$chr, y=dataset3$Prop_Win_Total, track.height=0.06, track.margin=c(0.03,0.01), bg.border = NA, panel.fun=function(startwin,Prop_Win_Total){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.6)
  })

  # roQSN region
  circos.rect(sector.index=3,xleft=28500000,ybottom=max(dataset3$Prop_Win_Total)*-0.1,xright=33000000,ytop=max(dataset3$Prop_Win_Total)*0.9,lwd=0.01,lty=1,col="pink")

  #ABCDE
  circos.segments(sector.index=1,x0=41514796,y0=0,x1=41514796,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52096899,y0=0,x1=52096899,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52554849,y0=0,x1=52554849,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11833234,y0=0,x1=11833234,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=17461203,y0=0,x1=17461203,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57253197,y0=0,x1=57253197,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57254457,y0=0,x1=57254457,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67016526,y0=0,x1=67016526,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67036489,y0=0,x1=67036489,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33238136,y0=0,x1=33238136,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46474063,y0=0,x1=46474063,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8425314,y0=0,x1=8425314,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8430126,y0=0,x1=8430126,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=61446161,y0=0,x1=61446161,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=58611899,y0=0,x1=58611899,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=60224075,y0=0,x1=60224075,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4414048,y0=0,x1=4414048,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4468385,y0=0,x1=4468385,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=50555186,y0=0,x1=50555186,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=52218507,y0=0,x1=52218507,y1=max(dataset3$Prop_Win_Total),lwd=2,lty=1,col="purple")
  circos.rect(sector.index=3,xleft=41280000,ybottom=0,xright=41990000,ytop=max(dataset3$Prop_Win_Total,na.rm=TRUE),lwd=0.01,lty=1,col="purple",border=0)
  
  
  # density gene

  circos.trackLines(densitygene$chr, densitygene$startwin+(densitygene$endwin-densitygene$startwin), densitygene$nb_genes/20000, col="forestgreen", lwd=1.5, type="l")

 # density SNP
  circos.trackLines(dataset3$chr, as.numeric(dataset3$startwin)+window*0.5, dataset3$Prop_Win_Total,lwd=2)

  
  # 4rd track = GWAS general 2016
  columntoplot=paste0(trait,".",model,"_P_value")
  columntofilter=paste0(trait,".",model,"_P_value_FDR")
  circos.trackPlotRegion(factors=  gwasres3$Chrom, y=-log10(gwasres3[[columntoplot]]), track.height=0.18, track.margin=c(0.03,0.01), bg.border = NA, ylim = c(0,max(-log10(gwasres3[[columntoplot]]),na.rm=TRUE)+0.3), panel.fun=function(Position,columntoplot){
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
  circos.trackPlotRegion(factors= dataset3$chr, y=dataset3$PropSNPWinPvalueInf01, ylim = c(0,max(dataset3$PropSNPWinPvalueInf01)), track.height=0.11, track.margin=c(0.03,0.01), bg.border = NA, panel.fun=function(startwin,PropSNPWinPvalueInf01){
    circos.axis(h=0, major.at=NULL, labels=NULL, major.tick=FALSE, lwd=1,labels.cex=0.6,col="grey")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    #circos.text(mean(xlim), max(ylim), paste("Chr_",sector.index), facing = "inside", cex = 1,font=2) # nom des secteurs seulement pour la première
    circos.yaxis("left", lwd=1,labels.cex=0.6)
  })
  #
  # roQSN region
  circos.rect(sector.index=3,xleft=28500000,ybottom=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*-0.1,xright=33000000,ytop=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE)*1.1,lwd=0.01,lty=1,col="pink")
  
  # ABCDE
  circos.segments(sector.index=1,x0=41514796,y0=0,x1=41514796,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52096899,y0=0,x1=52096899,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=1,x0=52554849,y0=0,x1=52554849,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=11833234,y0=0,x1=11833234,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=17461203,y0=0,x1=17461203,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57253197,y0=0,x1=57253197,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=57254457,y0=0,x1=57254457,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67016526,y0=0,x1=67016526,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=2,x0=67036489,y0=0,x1=67036489,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=3,x0=33238136,y0=0,x1=33238136,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=4,x0=46474063,y0=0,x1=46474063,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8425314,y0=0,x1=8425314,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=8430126,y0=0,x1=8430126,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=5,x0=61446161,y0=0,x1=61446161,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=58611899,y0=0,x1=58611899,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=6,x0=60224075,y0=0,x1=60224075,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4414048,y0=0,x1=4414048,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=4468385,y0=0,x1=4468385,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=50555186,y0=0,x1=50555186,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.segments(sector.index=7,x0=52218507,y0=0,x1=52218507,y1=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=2,lty=1,col="purple")
  circos.rect(sector.index=3,xleft=41280000,ybottom=0,xright=41990000,ytop=max(dataset3$PropSNPWinPvalueInf01,na.rm=TRUE),lwd=0.01,lty=1,col="purple",border=0)
  
  circos.trackLines(dataset3$chr, dataset3$startwin+window*0.5, dataset3$PropSNPWinPvalueInf01,lwd=1,col=ifelse(dataset3$chr%%2==0,"darkmagenta","darkorchid1"))
 
 dev.off() 
}


