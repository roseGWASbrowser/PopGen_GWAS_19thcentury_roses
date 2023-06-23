## 15/10/2020 Elise - analyzing output from fitPoly
## 188 active tetraploids

## Note: I run this script on MSU cluster. Change the paths to your directory if needed.
## I have "fitPoly_188_scores.dat" file in the directory where I am running. This file contains scores I got after fitPoly run.
## Tool files "probes_SNPs_positions_genome_K5_P_Q_names.txt" and "samples_ploidy_original.txt" needed to get positions on genome and expected ploidy levels.

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/fitpoly270/")

### Load packages
#############################################
library(ggplot2)
library(data.table)
library(tidyr) 


### Load results fitPoly
#############################################
df_score <- fread("fitPoly_270CARO_scores.dat", header=T)  # use fread to save time in reading large files 
dim(df_score) #  31512510 by 12
df_score[1:5,1:12] 

str(df_score)
df_score$marker = as.factor(df_score$marker)

### Modify results fitPoly to go into fitTetra_Output() fonction from SNPolisher (if needed)
#############################################
str(df_score) # need to modify slightly so it fit what fitTetra_Ouput()  is expecting

df_score$nsamp <- 270
df_score$select <- 1
df_score$model <- "model" # just as a place older to format the file
colnames(df_score)[2] <- "markername"
colnames(df_score)[3] <- "sample"

colnames(df_score)
df_score = df_score[,c(1,2,3,15,13,14,4:12)]

head(df_score)

df_score$markername=as.factor(as.character(df_score$markername)) ## will transform into factor / was not by default
df_score$sample=as.factor(as.character(df_score$sample)) ## will transform into factor / was not by default


### Stat file at the probe level
#############################################
length(levels(df_score$markername)) # 117,240 look like all markers are there same the one rejected
df_score[1:5,1:15]

df_stat = as.data.frame.matrix(table(df_score$markername,df_score$geno,useNA="ifany")) # make a contingency table into a df
dim(df_stat) # 117,240 x 6

df_stat$markername = rownames(df_stat)
rownames(df_stat) = seq(1 : length(rownames(df_stat)))
df_stat = df_stat[,c(7,1,2,3,4,5,6)]
colnames(df_stat) [2:7] <- c("nulliplex","simplex","duplex","triplex","quadruplex","no_call")
df_stat$nsamp <- rowSums(df_stat[,c("nulliplex","simplex","duplex","triplex","quadruplex","no_call")]) # add nb samples per markers / should be all 188
head(df_stat)
summary(df_stat)

### Nb of probes which passed fitPoly filters/scoring ?
#############################################
# produce ps file for no call probes
list_no_call = subset(df_stat$markername,df_stat$no_call == 270)
length(list_no_call) # 24,706 probes without call

# produce ps file for scored probes
list_scored = subset(df_stat$markername,df_stat$no_call < 270)
length(list_scored) # 92,007 probes with call

##########################################################
##########################################################
# generate data for plotting
##########################################################
##########################################################

# compute the weighted average of the 5 dosage scores, where the weights are the probabilities assigned by fitPoly
df_score$mean_score <- df_score$P0 * 0 + df_score$P1 * 1 +df_score$P2 * 2 + +df_score$P3 * 3 + df_score$P4 * 4

# scored probes indicator
df_score$scored_probes <- NA
df_score$scored_probes[which(df_score$markername %in% list_scored)] <- "yes"
df_score$scored_probes[which(!(df_score$markername %in% list_scored))] <- "no"

# list of genotypes
list_geno <- levels(df_score$sample)

# add ploidy info to df_score
df_ploidy = read.table("samples_ploidy_expected.txt",sep="\t",header=T)
df_ploidy$ploidy=as.factor(df_ploidy$ploidy)
summary(df_ploidy)

df_score = merge(df_score,df_ploidy) # merge on sample column
head(df_score)
dim(df_score)


# get SNP positions to order markers & compute cumulative positions
df_positions = fread("probes_SNPs_positions_genome_K5_P_Q_names.txt",header=T,stringsAsFactors=T,sep="\t")
head(df_positions)

colnames(df_positions) [3] <- "markername" # rename such as in df_score

df_score = merge(df_score,df_positions,by="markername")
head(df_score)
dim(df_score)

# compute cumulative position on OB genome
df_score_order_OB <- df_score[order(df_score$Chromosome_OB,df_score$Position_OB_genome_bp),]
df_score_order_OB$Position_OB_genome_bp_cumulative <- NA
rownames(df_score_order_OB) <- seq(1:length(rownames(df_score_order_OB)))

max_ch1 = max(df_score_order_OB[df_score_order_OB$Chromosome_OB == "Chr01","Position_OB_genome_bp"],na.rm=T)
max_ch2 = max(df_score_order_OB[df_score_order_OB$Chromosome_OB == "Chr02","Position_OB_genome_bp"],na.rm=T)
max_ch3 = max(df_score_order_OB[df_score_order_OB$Chromosome_OB == "Chr03","Position_OB_genome_bp"],na.rm=T)
max_ch4 = max(df_score_order_OB[df_score_order_OB$Chromosome_OB == "Chr04","Position_OB_genome_bp"],na.rm=T)
max_ch5 = max(df_score_order_OB[df_score_order_OB$Chromosome_OB == "Chr05","Position_OB_genome_bp"],na.rm=T)
max_ch6 = max(df_score_order_OB[df_score_order_OB$Chromosome_OB == "Chr06","Position_OB_genome_bp"],na.rm=T)
max_ch7 = max(df_score_order_OB[df_score_order_OB$Chromosome_OB == "Chr07","Position_OB_genome_bp"],na.rm=T)


df_score_order_OB$Position_OB_genome_bp_cumulative <- NA
df_score_order_OB$Position_OB_genome_bp_cumulative[which(df_score_order_OB$Chromosome_OB == "Chr01")] <- df_score_order_OB$Position_OB_genome_bp[which(df_score_order_OB$Chromosome_OB == "Chr01")] 
df_score_order_OB$Position_OB_genome_bp_cumulative[which(df_score_order_OB$Chromosome_OB == "Chr02")] <- df_score_order_OB$Position_OB_genome_bp[which(df_score_order_OB$Chromosome_OB == "Chr02")]  + max_ch1
df_score_order_OB$Position_OB_genome_bp_cumulative[which(df_score_order_OB$Chromosome_OB == "Chr03")] <- df_score_order_OB$Position_OB_genome_bp[which(df_score_order_OB$Chromosome_OB == "Chr03")]  + max_ch1+max_ch2
df_score_order_OB$Position_OB_genome_bp_cumulative[which(df_score_order_OB$Chromosome_OB == "Chr04")] <- df_score_order_OB$Position_OB_genome_bp[which(df_score_order_OB$Chromosome_OB == "Chr04")]  + max_ch1+max_ch2+max_ch3
df_score_order_OB$Position_OB_genome_bp_cumulative[which(df_score_order_OB$Chromosome_OB == "Chr05")] <- df_score_order_OB$Position_OB_genome_bp[which(df_score_order_OB$Chromosome_OB == "Chr05")]  + max_ch1+max_ch2+max_ch3+max_ch4
df_score_order_OB$Position_OB_genome_bp_cumulative[which(df_score_order_OB$Chromosome_OB == "Chr06")] <- df_score_order_OB$Position_OB_genome_bp[which(df_score_order_OB$Chromosome_OB == "Chr06")]  + max_ch1+max_ch2+max_ch3+max_ch4+max_ch5
df_score_order_OB$Position_OB_genome_bp_cumulative[which(df_score_order_OB$Chromosome_OB == "Chr07")] <- df_score_order_OB$Position_OB_genome_bp[which(df_score_order_OB$Chromosome_OB == "Chr07")]  + max_ch1+max_ch2+max_ch3+max_ch4+max_ch5+max_ch6
df_score_order_OB$Position_OB_genome_bp_cumulative[which(df_score_order_OB$Chromosome_OB == "Chr00")] <- df_score_order_OB$Position_OB_genome_bp[which(df_score_order_OB$Chromosome_OB == "Chr00")]  + max_ch1+max_ch2+max_ch3+max_ch4+max_ch5+max_ch6+max_ch7



##########################################################
##########################################################
# plotting / take some time, better to do on server
##########################################################
##########################################################
dir.create("ploidy_diagnostic_plots_perK")

adelaide= subset(df_score,df_score$sample=="adelaide_d_orleans_E10" & df_score$scored_probes == "yes")
aglaia=subset(df_score,df_score$sample=="aglaia_H10" & df_score$scored_probes == "yes")
lafrance=subset(df_score,df_score$sample=="la_france _D7" & df_score$scored_probes == "yes")
lafavorite=subset(df_score,df_score$sample=="la_favorite_lyon_A10" & df_score$scored_probes == "yes")
alfredcolomb=subset(df_score,df_score$sample=="alfred_colomb_C3" & df_score$scored_probes == "yes")
robustabourbon=subset(df_score,df_score$sample=="robusta_bourbon_C11" & df_score$scored_probes == "yes")


for (k in unique(df_score$Chromosome_OB)){
  print(k)
# to be used for the last plot

list_geno = levels(df_score$sample)

for (i in list_geno) {
  list_geno = levels(df_score$sample)

  tiff(paste("ploidy_diagnostic_plots_perK/",i,"_",k,".tiff",sep=""),width = 800, height = 600,res=50, type="cairo")
  
  par(mfrow=c(2,2))   
  
#### plot on OB / AVERAGE SCORE ##############################################
##############################################################################
  
#### hist score ##############################################
##############################################################################
  
temp = subset(df_score,df_score$sample == i & df_score$scored_probes =="yes"  & df_score$Chromosome_OB == k)  
dim(temp) 
head(temp)
 
 
hist(temp$mean_score,col="lightblue",xlim=c(0,4),breaks=seq(0,4,0.01),
       xlab="Weighted average score (allelic dose)",main=paste("probes: ",nrow(temp)),ylab="Probe counts")
  
mtext(side=3, cex=1, paste("expected ploidy :",unique(temp$ploidy)),col="black")
  
  
rm(temp)
  
  
  

#### hist ratio ##############################################
##############################################################################

temp = subset(df_score,df_score$sample == i & df_score$scored_probes =="yes"  & df_score$Chromosome_OB == k)
temp_total = subset(df_score,df_score$sample == i & df_score$scored_probes =="yes")
dim(temp) 
   
hist(temp$ratio,col="lightblue",xlim=c(0,1),breaks=seq(0,1,0.01),
     xlab="Ratio",main=paste("probes: ",nrow(temp)),ylab="Probe counts")
mtext(side=3, cex=1, paste("expected ploidy :",unique(temp$ploidy)),col="black")


### plot with density curved for cultivars with clear 2n, 3n and 4n patterns

plot(density(adelaide$ratio,na.rm=TRUE), col="darkred", lwd=1.5,lty=2,bty='l',main="")
lines(density(aglaia$ratio,na.rm=TRUE),col="red",lwd=1.5,lty=2)
lines(density(lafrance$ratio,na.rm=TRUE), col="goldenrod1", lwd=1.5,lty=2)
lines(density(lafavorite$ratio,na.rm=TRUE),col="darkorange",lwd=1.5,lty=2)
lines(density(alfredcolomb$ratio,na.rm=TRUE), col="green", lwd=1.5,lty=2)
lines(density(robustabourbon$ratio,na.rm=TRUE),col="forestgreen",lwd=1.5,lty=2)
mtext(side=3, cex=1, "red=2n;orange=3n;green=4n",col="black")


### plot with observed density for the focal individual
plot(density(adelaide$ratio,na.rm=TRUE),col="darkred", lwd=1.2,lty=2,bty="l",main="")
lines(density(aglaia$ratio,na.rm=TRUE),col="red",lwd=1.2,lty=2)
lines(density(lafrance$ratio,na.rm=TRUE), col="goldenrod1", lwd=1.2,lty=2)
lines(density(lafavorite$ratio,na.rm=TRUE),col="darkorange",lwd=1.2,lty=2)
lines(density(alfredcolomb$ratio,na.rm=TRUE), col="green", lwd=1.2,lty=2)
lines(density(robustabourbon$ratio,na.rm=TRUE),col="forestgreen",lwd=1.2,lty=2)
lines(density(temp_total$ratio,na.rm=TRUE),col="blue", lwd=1.5,bty='l',lty=2)
lines(density(temp$ratio,na.rm=TRUE),col="blue", lwd=2.5,bty='l')

rm(temp)
rm(temp_total)
dev.off()

cat(i," ")

}

}

  
  
  
  















