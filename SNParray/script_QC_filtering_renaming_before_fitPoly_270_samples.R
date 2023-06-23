#### Loading fitPolyTools as function (tar package archive make troubles !)
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
setwd("/bigvol/benoit/Analyses/Temp_Tibo/Roses")
source("fitPolyTools_20200312.R") ## you can see function are all imported correctly 
## need to have fitPolyTools_20200312.R in same directory or indicate path

library(ggplot2)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
####  Read and convert the array data +  print convert file (take a few minutes)
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

# path to Axiom array data with 188 samples on MSU cluster (change to your path if needed)
#/mnt/home/alberte2/ROSE_SNP_ARRAY/fitPoly_FILTERING/Axiom_output_Rose_188_samples

datAX <- readAxiomSummary(AXdata="/bigvol/benoit/Analyses/Temp_Tibo/Roses/AxiomGT1.270.summary.txt",out="") # take few minutes
# not saving re-formated data, will do later after filtering, thus use out=""

head(datAX)
dim(datAX) # 37202220  / 137786 probes for 270 samples
length(levels(datAX$MarkerName))
length(levels(datAX$SampleName))

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
####  Rename probes so it will be appropriate later for filtering probes
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
library(data.table) # for fast reading of input files

mrktable <- fread("names_P_Q_format.txt")
head(mrktable)

colnames(mrktable)[2] <- "MarkerName"

datAX <- merge(mrktable[,c(2,5)],datAX) ## merge by marker name and add column with P / Q names
dim(datAX) # 37202220  / 137786 probes for 270 samples
head(datAX)

datAX = datAX[,-1]
colnames(datAX)[1] <- "MarkerName" 

head(datAX)
datAX$MarkerName <- as.factor(datAX$MarkerName)
length(levels(datAX$MarkerName)) # still 137786 / ok / just renamed using Voorrips convention


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### Cleaning the data based on total signal intensity
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
head(datAX)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
### 1a : Selection of samples based on total signal intensity R
#### #### #### #### #### #### #### #### #### #### #### ####
library(ggplot2)

sampRmean <- tapply(datAX$R, datAX$SampleName, mean)
temp =as.data.frame(sampRmean)
head(temp)

tiff("270_meanRsignalpersample.tiff",width = 2300, height = 2300,res=300)
ggplot(temp, aes(x=sampRmean)) + 
  geom_histogram(binwidth=0.5)+xlab("Mean R per sample")+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()



#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
### 2a : Selection of markers based on total signal intensity R
#### #### #### #### #### #### #### #### #### #### #### ####

Rstats <- calcRstats(datAX, out=NA) 
head(Rstats) # stats at the marker level

quantile(Rstats$q95, probs=c(0, 0.5, 1)) # min, mean and max of R 95% quantiles for each markers
# large variation depending on the markers

# look at a few markers for various mean R levels
### R levels to look at
Rlevels <- seq(100, 6200, by=400) 
length(Rlevels) # 16 levels

### select random markers corresponding to the various R levels
sel.Rstats <- selMarkers_byR(Rstats, Rlevels, mrkperlevel=3) # select markers
dim(sel.Rstats) # 45 markers selected 

### draw plots
subDir <- "diagnostic_plot_probes_R_95_quantiles_270_samples"
mainDir <- getwd()

if (file.exists(subDir)){
  print("directory already created")} else {
  dir.create(file.path(mainDir, subDir))
 }

drawXYplots(dat=datAX, markers=sel.Rstats,
            out="diagnostic_plot_probes_R_95_quantiles_270_samples/R_95_quantile_probe_level", drawRthresholds=TRUE)

# hard to decide / will cut after marker 12 Rq95 > 1400
tiff("270_Rsignal_probes_quantile95_R.tiff",width = 2300, height = 2300,res=300)
ggplot(Rstats, aes(x=q95)) + 
  geom_histogram(binwidth=5)+xlab("4 signal intensity")+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
  geom_vline(xintercept = 1400,col="red",size=1)
dev.off()

keep <- Rstats$MarkerName[Rstats$q95 >= 1400]
dat <- datAX[datAX$MarkerName %in% keep,]
dim(dat) #  33373890 / 270 = 123607 probes
dim(datAX) # 37202220 /270 = 137786 probes 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
### 3c : Selection based on R values of samples within markers & print file input for fitPoly (large file >> take some time)
#### #### #### #### #### #### #### #### #### #### #### ####
# The following command selects all markers where q95 (the 95% quantile of R) is > 1400, 
# and within these markers it assigns a missing value (NA) to all markers where R < 0.5 * q95.

dat.sel <- makeFitPolyFiles(dat, out="AxiomGT1.270.summary.fitPoly.format",filetype="dat", Rquantile=0.95, marker.threshold=1400,
                            Rthreshold.param=c(0.95, 0.5, 0)) # will produce file to use in fitPoly, take a few minutes 

#dat.sel is in this case a list with only one element; simplify:
dat.sel <- dat.sel[[1]]
head(dat.sel)

sum(is.na(dat.sel$ratio))/length(dat.sel$ratio)*100
# 1475530 / 33373890  * 100 = 4.42 % of the data were set to missing

# missing data per sample
missing_sample_level = as.data.frame.matrix(table(dat.sel$SampleName,is.na(dat.sel$ratio)) )
head(missing_sample_level)
colnames(missing_sample_level) <- c("data","missing")
missing_sample_level$percent_missing = missing_sample_level$missing/(missing_sample_level$data+missing_sample_level$missing)*100
mean(missing_sample_level$percent_missing)
quantile(missing_sample_level$percent_missing, probs=c(0, 0.5,0.75, 0.90,0.95, 1))

tiff("270samples_missing_sample_level.tiff",res=300,width = 3000, height = 3000,type="cairo")
ggplot(missing_sample_level, aes(x=percent_missing)) + 
  geom_histogram(binwidth=0.5)+xlim(0, 10)+
  xlab("Missing data per sample (%)")+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()

# missing data per probe
missing_probe_level = as.data.frame.matrix(table(dat.sel$MarkerName,is.na(dat.sel$ratio)) )
head(missing_probe_level)
colnames(missing_probe_level) <- c("data","missing")
missing_probe_level$percent_missing = missing_probe_level$missing/(missing_probe_level$data+missing_probe_level$missing)*100
mean(missing_probe_level$percent_miss)
quantile(missing_probe_level$percent_miss, probs=c(0, 0.5,0.75, 0.90,0.95, 1)) # 90% of probes with less than 13% missing data

tiff("270samples_missing_probe_level.tiff",res=300,width = 3000, height = 3000,type="cairo")
ggplot(missing_probe_level, aes(x=percent_missing)) + 
  geom_histogram(binwidth=1)+
  xlab("Missing data per probe (%)")+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()


