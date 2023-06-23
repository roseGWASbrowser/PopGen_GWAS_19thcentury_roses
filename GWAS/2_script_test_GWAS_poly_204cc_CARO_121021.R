## October 2021 from Elise April 21- GWAS using genotyping results from fitPoly 
## 205 samples considered as tetraploids only CARO

#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### Load GWASpoly functions and needed packages
#### The tar version of the GWASpoly package makes trouble. I was not able to make a clean install of the package.
#### Thus, here I am loading functions one by one.
#### Attention ! Order matters (functions are calling each others).
#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### 

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/")
source("qqunif.plot.R") # https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R

# note: change path to your path!
#source("~/Software/GWASpoly-master/R/design.score.R")
#source("~/Software/GWASpoly-master/R/get.QTL.R")
#source("~/Software/GWASpoly-master/R/get_x.R")
#source("~/Software/GWASpoly-master/R/GWASpoly.data.R")
#source("~/Software/GWASpoly-master/R/GWASpoly.K.R")
#source("~/Software/GWASpoly-master/R/GWASpoly.R")
#source("~/Software/GWASpoly-master/R/make.full.R")
#source("~/Software/GWASpoly-master/R/makeLOCO.R")
#source("~/Software/GWASpoly-master/R/manhattan.plot.R")
#source("~/Software/GWASpoly-master/R/qq.plot.R")
#source("~/Software/GWASpoly-master/R/qvalue.R")
#source("~/Software/GWASpoly-master/R/read.GWASpoly.R")
#source("~/Software/GWASpoly-master/R/score.calc.R")
#source("~/Software/GWASpoly-master/R/set.K.R")
#source("~/Software/GWASpoly-master/R/set.params.R")
#source("~/Software/GWASpoly-master/R/set.threshold.R")
#source("~/Software/GWASpoly-master/R/write.GWASpoly.R")
#source("~/Software/GWASpoly-master/R/GWASpoly.fitted.R") 
#source("~/Software/GWASpoly-master/R/GWASpoly.thresh.R") 

source("~/Software/GWASpoly-master-elise/design.score.R")
source("~/Software/GWASpoly-master-elise/get.QTL.R")
source("~/Software/GWASpoly-master-elise/get_x.R")
source("~/Software/GWASpoly-master-elise/GWASpoly.data.R")
source("~/Software/GWASpoly-master-elise/GWASpoly.K.R")
source("~/Software/GWASpoly-master-elise/GWASpoly.R")
source("~/Software/GWASpoly-master-elise/make.full.R")
source("~/Software/GWASpoly-master-elise/makeLOCO.R")
source("~/Software/GWASpoly-master-elise/manhattan.plot.R")
source("~/Software/GWASpoly-master-elise/qq.plot.R")
source("~/Software/GWASpoly-master-elise/qvalue.R")
source("~/Software/GWASpoly-master-elise/read.GWASpoly.R")
source("~/Software/GWASpoly-master-elise/score.calc.R")
source("~/Software/GWASpoly-master-elise/set.K.R")
source("~/Software/GWASpoly-master-elise/set.params.R")
source("~/Software/GWASpoly-master-elise/set.threshold.R")
source("~/Software/GWASpoly-master-elise/write.GWASpoly.R")
source("~/Software/GWASpoly-master-elise/GWASpoly.fitted.R") 
source("~/Software/GWASpoly-master-elise/GWASpoly.thresh.R") 
# note: install packages if needed
library(rrBLUP) # requested to run GWASpoly functions !!!!!
library(tidyr)
library(ggplot2)
library(plyr)
library(gplots)
library(data.table)
library(stats)
library(doParallel)

#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### STEP 1/ Prepare inputs files for GWASpoly
#### #### #### #### #### #### #### #### #### #### #### #### ####  
#### #### #### #### #### #### #### #### #### #### #### #### #### 

################################################################
#### SNP data, read snp file and reformate for GWASpoly
################################################################

### read snp file 
#marker <-  fread("matrix_fitPoly_270_scores_formated_only_passing_probes.204samplescc_withoutmono_MAF5pc_AX_ID.txt") # genotype matrix
# version tronque du 120722
marker <-  fread("matrix_fitPoly_270_scores_formated_only_passing_probes.204samplescc_withoutmono_MAF5pc_AX_ID_tronq.txt.txt") # genotype matrix
dim(marker) # 64813 snps vs 96 genotypes
marker[1:5,1:10]
colnames(marker)[1] <- "Marker"

### read snp positions
map <- fread("map_postions_markers_AXIOM_probes.txt") # load position of markers on genome needed by GWASpoly
head(map)

### group snp positions and markers
TableS1 <- merge(map,marker)
TableS1[1:5,1:10]

### create recode function to convert 0 to AAAA 1 to BAAA, 2 to BBAA, 3 to BBBA and 4 to BBBB
recode <-colwise(function(x){return(ifelse(x == 0, "AAAA",
                                    ifelse(x == 1, "BAAA",
                                    ifelse(x == 2, "BBAA",
                                    ifelse(x == 3, "BBBA",
                                    ifelse(x == 4, "BBBB",x))))))})
temp = recode(TableS1[,-c(1:3)])
temp[1:5,1:10]
TableS1=cbind(TableS1[,c(1:3)],temp)
rm(temp)
TableS1[1:5,1:10]
 
### order by Chrom and Position
TableS1 <- TableS1[order(TableS1$Chrom,TableS1$Position),]
TableS1[1:5,1:5]

dim(TableS1) # 46,691 probes x 207 
list_geno = colnames(TableS1)[-(1:3)]
length(list_geno) #204 genotypes
list_geno<- noquote(list_geno) # remove quote if pheno have noquote


#write.csv(TableS1,"TableS1_recodeddataset_204cc_MAF5pc.csv",row.names = F,quote = F) # write final version in csv for import with GWASpoly
write.csv(TableS1,"TableS1_recodeddataset_204cc_MAF5pc_tronq2.csv",row.names = F,quote = F) # write final version in csv for import with GWASpoly

################################################################
#### phenotypic data, read file and reformate for GWASpoly
################################################################

### read phenotype file ready
pheno <-  read.table("phenotype_petals",header=T) # phenotypes
head(pheno)
dim(pheno)

pheno <- pheno[which(pheno$ID %in% list_geno),] # keep only individuals that are in the list, here no change both files with the same 204cc
dim(pheno)

pdf(file="Pheno_petals_histo_204cc_121021.pdf",height=7,width=7)
ggplot(pheno, aes(x=mean_petals))+ 
  geom_histogram(color="darkgrey",fill="grey",binwidth=5)+
  xlab("Averaged number of petals counted per cultivar")+ylab("Number of rose cultivars")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
dev.off()

pheno = pheno[,c(1,10)]
colnames(pheno)[1] <- "Name"


### read structure 
structure <-  read.table("204cc_structure_PCoA_DAPC_MAF5.txt",header=T) # phenotypes
head(structure)
colnames(structure)[1] <- "Name"

### join pheno et structure
TableS2 <- merge(pheno,structure[,c(1:6)]) # use only the columns of the PCoA
TableS2[1:5,]

write.csv(TableS2,"TableS2_petals_PCoA.csv",row.names = F,quote = F) # write final version in csv for import with GWASpoly


#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### 
### STEP 2/ Import data into GWASpoly
#### #### #### #### #### #### #### #### #### #### #### #### ####   format="ACGT"
#### #### #### #### #### #### #### #### #### #### #### #### #### 

#data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_petals_PCoA.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc.csv",format="AB",n.traits=3,delim=",")
#data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_TN_PCoA.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc.csv",format="AB",n.traits=3,delim=",")
#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/scent2017")
#data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_dataPF2017.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc.csv",format="AB",n.traits=150,delim=",")
#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/scent2018")
#data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_dataPF2018.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc.csv",format="AB",n.traits=153,delim=",")
#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/scentAvg2years")
#data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_dataPFavg2years.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc.csv",format="AB",n.traits=149,delim=",")
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/")
#data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_PermanentFlowering_A3_v2.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc.csv",format="AB",n.traits=28,delim=",")
# version tronque du 12/07/22
data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_PermanentFlowering_A3_v2_tronq.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc_tronq.csv",format="AB",n.traits=28,delim=",")
#data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_PermanentFlowering_A3.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc.csv",format="AB",n.traits=14,delim=",")
#data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_PermanentFlowering_A2.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc.csv",format="AB",n.traits=14,delim=",")
#data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_PermanentFlowering_A1.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc.csv",format="AB",n.traits=14,delim=",")
#setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/phenoflorhige")
#data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_phenoflorhigequalitative.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc.csv",format="AB",n.traits=32,delim=",")
#data <-read.GWASpoly(ploidy=4,pheno.file="TableS2_phenoflorhigequantitative.csv",geno.file="TableS1_recodeddataset_204cc_MAF5pc.csv",format="AB",n.traits=16,delim=",")

# take a few minutes

str(data) # list of objects contained in data 
head(data@map) # use @ to look at objects
head(data@fixed)
head(data@pheno)
data@geno[1:5,1:5]


#Number of polymorphic markers: 46691 
#Missing marker data imputed with population mode 
#N = 204 individuals with phenotypic and genotypic information 
#Detected following fixed effects:
#PCo1
#PCo2
#PCo3
#PCo4
#PCo5
#Detected following traits:
#mean_petals



#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### 
### STEP 3/ Computing K matrix and setting structure cofactors
#### #### #### #### #### #### #### #### #### #### #### #### ####  
#### #### #### #### #### #### #### #### #### #### #### #### #### 

### compute K matrix
#Mscale <- scale(data@geno,center=T,scale=F)
#Kval <- tcrossprod(Mscale)
#Kmatrix <- Kval/mean(diag(Kval))

#gid.K <- rownames(Kmatrix)
#gid <- rownames(data@geno)
#missing <- setdiff(gid,gid.K)
#if (length(missing)==0){
#  print("ok")
#  data <- new("GWASpoly.K",map=data@map,pheno=data@pheno,fixed=data@fixed,geno=data@geno,ploidy=data@ploidy,K=K)
#}
#K <- list(all=Kmatrix[gid,gid])

data <- set.K(data) ## this is computing VanRaden kinship from marker data / you could provide your own matrix if you like
str(data) 
dim(data@K)
#View(data@K)
### 3 next lines, mandatory on my session, but not needed for Elise

####
#dim(Kmatrix) # or dim(data@K) without the 3 lines above
#View(data@K)
#View(as.list(data@K))

#data@K <- list(Kmatrix)


### plotting K matrix
hr <- hclust(dist(K), method="complete")
#hr <- hclust(dist(data@K), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.2))

clusterCols1 <- rainbow(length(unique(structure$kmeans))) 
myClusterSideBar1 <- clusterCols1[as.numeric(structure$kmeans)]

pdf(file="Kinship_204cc_121021.pdf",height=18,width=18)
#tiff("Kinship_vanraden.tiff",width = 1000, height = 1000,res=200)
#heatmap(data@K, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hr), scale="none",
heatmap(K, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hr), scale="none",
        RowSideColors= myClusterSideBar1,cexRow=0.65,margins = c(12, 12))
rect(xleft=0.01,xright=0.025,ytop=0.925,ybottom=0.95,col=clusterCols1[1])
rect(xleft=0.01,xright=0.025,ytop=0.9,ybottom=0.925,col=clusterCols1[2])
rect(xleft=0.01,xright=0.025,ytop=0.875,ybottom=0.90,col=clusterCols1[3])
rect(xleft=0.01,xright=0.025,ytop=0.85,ybottom=0.875,col=clusterCols1[4])
rect(xleft=0.01,xright=0.025,ytop=0.825,ybottom=0.85,col=clusterCols1[5])
rect(xleft=0.01,xright=0.025,ytop=0.8,ybottom=0.825,col=clusterCols1[6])
rect(xleft=0.01,xright=0.025,ytop=0.775,ybottom=0.8,col=clusterCols1[7])
rect(xleft=0.01,xright=0.025,ytop=0.75,ybottom=0.775,col=clusterCols1[8])
rect(xleft=0.01,xright=0.025,ytop=0.725,ybottom=0.750,col=clusterCols1[9])
text(x=0.021,y=0.985,"DAPC groups\n(SNP, K=9)",cex=0.7,font=2)
text(x=0.018,y=0.939,"1",cex=0.7,font=2)
text(x=0.018,y=0.914,"2",cex=0.7,font=2)
text(x=0.018,y=0.889,"3",cex=0.7,font=2)
text(x=0.018,y=0.864,"4",cex=0.7,font=2)
text(x=0.018,y=0.839,"5",cex=0.7,font=2)
text(x=0.018,y=0.814,"6",cex=0.7,font=2)
text(x=0.018,y=0.789,"7",cex=0.7,font=2)
text(x=0.018,y=0.764,"8",cex=0.7,font=2)
text(x=0.018,y=0.739,"9",cex=0.7,font=2)
dev.off()

### set params for running GWAS
params <- set.params(fixed=c("PCo1","PCo2","PCo3","PCo4","PCo5"),fixed.type=rep("numeric",5),n.PC=0,MAF=0.05,geno.freq=0.95,P3D=T)
## do not link to data, just create a object that will be used when running GWAS to indicate where/how to find structure
## you can removed fixed effect here to run model with correction for population structure

#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### 
### STEP 4/ LETS GO GWAS ! (rrBLUP need to be loaded)
#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### 
tested_models <- c("additive","diplo-additive","general","diplo-general","1-dom","2-dom") # indicate which models to test
traits <- colnames(data@pheno)[-1] # indicate name of phenotype to analyze

traits <- colnames(data@pheno)[2] # year_obtained
#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### BELOW RUN GWAS ! TAKE LONG TIME IF SEVERAL MODELS ARE TESTED AND SEVERAL TRAITS ANALYZED !
#### BETTER TO RUN ON CLUSTER !


list_GWAS_scores <- list() # list to store scores for several traits
list_GWAS_effects <- list() # list to store effects for several traits


for (i in 1:length(traits[1:1])) {
#for (i in 1:length(traits[1:8])) {  # here run for the 1st 3 compounds

print(paste("running GWAS for",traits[i]))
print("######################")
  
#### run GWAS
data_res <- GWASpoly(data,models=tested_models,traits=traits[i],params=params,n.core=1,quiet=F)
#str(data_res) ## now you have @scores (-log10 pvalue) and @effects added for the various traits analyzed

print(paste("extract results for",traits[i]))
print("######################")

#### extract effects
effects = data_res@effects
effects=as.data.frame(effects)
effects$Marker = row.names(effects)
row.names(effects) <- 1:length(row.names(effects))
effects = effects[,c(length(colnames(effects)),1:(length(colnames(effects))-1))]
colnames(effects)[-1] <-  paste(colnames(effects)[-1],"effect",sep="_")

effects <-  merge(data@map,effects)

list_GWAS_effects[[i]] <-  effects # store effects


#### extract scores
scores = data_res@scores # -log10 pvalues from rrBLUP
scores=as.data.frame(scores)
scores$Marker = row.names(scores)
row.names(scores) <- 1:length(row.names(scores))
scores = scores[,c(length(colnames(scores)),1:(length(colnames(scores))-1))]
colnames(scores)[-1] <-  paste(colnames(scores)[-1],"score",sep="_")

#### get P-values from scores
temp = as.matrix(as.data.frame(data_res@scores))
temp = 10^-temp 
colnames(temp) <- paste(colnames(temp),"P_value",sep="_")
scores = cbind(scores,temp)

#### get FDR corrected P-values from scores
temp_BH = sapply(as.data.frame(temp),p.adjust,method="BH") # need to have as.data.frame otherwise adjust is done across columns
temp_BH=as.data.frame(temp_BH)
colnames(temp_BH) <- paste(colnames(temp_BH),"FDR",sep="_")
scores = cbind(scores,temp_BH)
rm(temp,temp_BH) # clean space

#### store scores into list
scores_pos <- merge(data@map,scores)
list_GWAS_scores[[i]] <-   scores_pos # store effects

####  change directory / one directory per traits
subDir <- paste(traits[i],"_gwas_results",sep="")
mainDir <- getwd()

if (file.exists(subDir)){
  setwd(file.path(mainDir, subDir))
} else {
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir, subDir))
  
}

####  write result files (scores and effects) in trait directory
write.table(scores_pos,paste(traits[i],"_scores.txt",sep=""),quote=F,row.names = F)
write.table(effects,paste(traits[i],"_effects.txt",sep=""),quote=F,row.names = F)

print(paste("plotting results for",traits[i]))
print("######################")


####  qqplot and  manhattan plots (for each model tested)

# bring positions and chrom information
scores_pos = merge(scores, data@map[,1:3])
scores_pos = scores_pos[,c(1,(length(colnames(scores_pos))-1):length(colnames(scores_pos)), 2:length(colnames(scores)))] 
scores_pos = scores_pos[order(scores_pos$Chrom,scores_pos$Position),] # chrom ordered from 0 to NA
scores_pos$Chrom <- paste("Chr",scores_pos$Chrom,sep="_")

# get the number and names of the columns to plot 
tested_models_to_plot <- as.character(sapply(data_res@scores, names))
index <- (length(tested_models_to_plot)+4) : (length(tested_models_to_plot)*2+3) # from data.frame architecture

# plots
for (j in colnames(scores_pos)[index]) { 
 
# for a given model select non NA values  
temp = scores_pos[which(!is.na(scores_pos[,j])),]
rownames(temp)<-1:length(rownames(temp))
head(temp)
temp$pos <- 1:length(rownames(temp))   

## qqplot
tiff(paste(j,"qqplot.tiff",sep="_"),width=800,height=800,res=150)
p<-qqunif.plot(temp[,j],main=j) 
plot(p)
dev.off()

## manhattan
max = tapply(temp$pos, temp$Chrom, max, na.rm=T)
min = tapply(temp$pos, temp$Chrom, min, na.rm=T)
middle = (as.numeric(max)+as.numeric(min))/2

temp$Chrom = as.factor(temp$Chrom)
levels(temp$Chrom)

colors <- rainbow(length(unique(temp$Chrom))) # 
coloring <- colors[as.numeric(temp$Chrom)]

tiff(paste(j,"manhattan.tiff",sep="_"),width=1900,height=900,res=150)
plot(rownames(temp),-log10(temp[,j]),
     pch=19,col=coloring,
     ylab="-log10(Pvalue)",
     xlab="Positions",
     main=j,xaxt='n')
axis(side = 1, at = middle  ,  labels = c("Ch00","Ch01","Ch02","Ch03","Ch04","Ch05","Ch06","Ch07","ChNA"))
dev.off()

}

### back to initial directory
setwd(file.path(mainDir))

### clean space
rm(temp,scores,scores_pos,effects,min,max,middle,index,tested_models_to_plot,data_res,scores_pos)


}




#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### 
### STEP 5/ ALL RESULS ARE STORED IN LIST WITH DATA.FRAME PER TRAITS / CAN USE THESE LIST TO EXTRACT SIGNFICANT SNPS AND DO FURTHER PLOTTING
#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### ####

### save the results as .Rdata for later
save(list_GWAS_scores, file = paste("GWAS_scores",Sys.Date(),"Rdata",sep="."))
save(list_GWAS_effects, file = paste("GWAS_effects",Sys.Date(),"Rdata",sep="."))

### if you want to import back from .Rdata
#load("GWAS_scores.2021-01-18.Rdata")
#load("GWAS_effects.2021-01-18.Rdata")


### plots manhattan all models for a given trait
#####################################################################
########################################################################
str(list_GWAS_scores) # 5 data.frame corresponding to the 5 traits tested
as.character(sapply(list_GWAS_scores[1], names)) # structure different models for a given traits

significant_markers <- list() # list to store significant markers 


for (i in 1:length(traits)) { 

print(traits[i])  
  
#### fdr significant markers
temp_scores_fdr <- as.data.frame(list_GWAS_scores[i])[,c(1,2,3,6:13,22:29)]
temp_scores_fdr <-  temp_scores_fdr[order(temp_scores_fdr$Chrom,TableS1$Position),]
rownames(temp_scores_fdr) <- 1 : length(rownames(temp_scores_fdr))  
colnames(temp_scores_fdr)[4:11] <- paste(c("additive","diplo.additive","general","diplo.general","1.dom.alt","1.dom.ref","2.dom.alt","2.dom.ref"),"_score",sep="")
colnames(temp_scores_fdr)[12:19] <- paste(c("additive","diplo.additive","general","diplo.general","1.dom.alt","1.dom.ref","2.dom.alt","2.dom.ref"),"_FDR",sep="")
temp_scores_fdr$pos <- 1 : length(rownames(temp_scores_fdr))  
temp_scores_fdr$Chrom <- paste("Chr",temp_scores_fdr$Chrom,sep="")
temp_scores_fdr$Chrom <- as.factor(temp_scores_fdr$Chrom)

temp_scores_fdr_5 <- temp_scores_fdr[apply(temp_scores_fdr[,12:19]<0.05,1,any,na.rm=TRUE),]
print(temp_scores_fdr_5)  

significant_markers[[i]] <-   temp_scores_fdr_5 # store effects (with map positions)



#### plots
temp_scores <- as.data.frame(list_GWAS_scores[i])[,c(1:13)]
temp_scores <-  temp_scores[order(temp_scores$Chrom,TableS1$Position),]
rownames(temp_scores) <- 1 : length(rownames(temp_scores))  
colnames(temp_scores)[6:13] <- c("additive","diplo.additive","general","diplo.general","1.dom.alt","1.dom.ref","2.dom.alt","2.dom.ref")
temp_scores$pos <- 1 : length(rownames(temp_scores))  
temp_scores$Chrom <- paste("Chr",temp_scores$Chrom,sep="")
temp_scores$Chrom <- as.factor(temp_scores$Chrom)

max = tapply(temp_scores$pos, temp_scores$Chrom, max, na.rm=T)
min = tapply(temp_scores$pos, temp_scores$Chrom, min, na.rm=T)
middle = (as.numeric(max)+as.numeric(min))/2

temp_scores_long <- temp_scores %>% gather(model, scores, "additive" : "2.dom.ref")
temp_scores_long$Marker <- as.factor(temp_scores_long$Marker)
temp_scores_long$Ref <- as.factor(temp_scores_long$Ref)
temp_scores_long$Alt <- as.factor(temp_scores_long$Alt)
temp_scores_long$model <- as.factor(temp_scores_long$model)

tiff(paste(traits[i],"_all_models.tiff",sep=""),res=150,width = 1200, height = 480)

p<- ggplot(temp_scores_long, aes(x=pos, y=scores,shape=model,color=Chrom)) +
  geom_point(size=1)+theme_bw()+ggtitle(traits[i])+guides( color = FALSE, shape=guide_legend(title="models"))+
  scale_x_continuous(name ="Chromosomes",breaks=middle,labels=c("Ch00","Ch01","Ch02","Ch03","Ch04","Ch05","Ch06","Ch07","ChNA"))+
  scale_shape_manual(values=c(0,1,2,3,4,15,19,11))+ylab("-log10(P-value)")

print(p)  

dev.off()
  
}

#p.adjust(0.00000107085,n=46691,method="BH")
#[1] 0.04999906             
#[1] 5.970271 (-log10 value FDR=0.05)
# p.adjust(0.0000021418,n=46691,method="BH")
#[1] 0.1000028     
# 5.669221 (-log10 value FDR=0.1)

# p.adjust(0.0000042835,n=46691,method="BH")
# [1] 0.2000009
# 5.368151
