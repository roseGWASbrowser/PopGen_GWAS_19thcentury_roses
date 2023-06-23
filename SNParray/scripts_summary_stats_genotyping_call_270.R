## 15/10/2020 Elise - analyzing output from fitPoly
## 188 active tetraploids

## Note: I run this script on MSU cluster. Change the paths to your directory if needed. Install packages if needed.
## I have "fitPoly_188_scores.dat" file in the directory where I am running. This is scores I got after fitPoly run.
## I am also using "AxiomGT1.summary_188_samples.fitPoly.format_poly.dat" to get marker stats before fitPoly. This is output from fitPolyTools filtering step.
## Tool files "probes_SNPs_positions_genome_K5_P_Q_names.txt" and "samples_ploidy_original.txt" needed to have positions on genome and expected ploidy levels.

### Load packages
#############################################
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(zoo)


setwd("/bigvol/benoit/Analyses/Temp_Tibo/Roses/output_fitpoly_270samples_CARO/")

### Load results fitPoly
#############################################
df_score <- fread("fitPoly_270CARO_scores.dat", header=T)  # use fread to save time in reading large files 
dim(df_score) #  33103324 by 12 = 33103324/284 => 116,561 probes
df_score[1:5,1:12] 

str(df_score)
df_score$marker = as.factor(df_score$marker)
length(levels(df_score$marker )) # 116,713 probes

### compare to original marker set after filtering
df_original_marker <- fread("/bigvol/benoit/Analyses/Temp_Tibo/Roses/AxiomGT1.270.summary.fitPoly.format_poly.dat")
head(df_original_marker )
length(levels(as.factor(df_original_marker$MarkerName))) # 123,607 probes were kept after fitPolyTools filtering

df_original_marker_stat <- as.data.frame.matrix(table(df_original_marker$MarkerName,is.na(df_original_marker$ratio)) )
head(df_original_marker_stat)
colnames(df_original_marker_stat) <- c("data","missing")
df_original_marker_stat$percent_missing = df_original_marker_stat$missing/(df_original_marker_stat$data+df_original_marker_stat$missing)*100
mean(df_original_marker_stat$percent_missing)  # 4.42
df_original_marker_stat$MarkerName <- rownames(df_original_marker_stat)    
rownames(df_original_marker_stat) <- 1: length(rownames(df_original_marker_stat))
df_original_marker_stat <- df_original_marker_stat[,c(4,1,2,3)]    

`%notin%` <- Negate(`%in%`)
lost_marker = df_original_marker_stat [which(df_original_marker_stat$MarkerName %notin% df_score$MarkerName),]
dim(lost_marker ) # 6894
subset(lost_marker,lost_marker$missing <= 25) # all these 6894 with more than 25% missing, excepted Affx-86783002P Affx-86783051Q Affx-86824872P and Affx-86825785P not sure why. 
# Note: The parameter call.threshold in fitPoly run was set to 0.75, meaning minimum 75% of genotypes to be scored to keep a probe. 
# In log file (fitPoly_270.log), these four probes are noted as : optmodel1=NA - rejected: optmodel1=NA. This suggests no optimal model was found.
df_original_marker[which(df_original_marker$MarkerName=="Affx-86783002P"),] 
df_original_marker[which(df_original_marker$MarkerName=="Affx-86783051Q"),] 

### Modify results fitPoly format to go into fitTetra_Output() fonction from SNPolisher (if needed)
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

fwrite(df_score,"fitPoly_270_scores_formatedSNPolisher.dat",quote=F,row.names = F,sep = "\t") # fwrite for faster writing

### STATS AT THE PROBE LEVEL
#############################################
#############################################
#############################################

### Stat file at the probe level
#############################################
length(levels(df_score$markername)) #  116,713 look like all markers are there same the one rejected
df_score[1:5,1:15]

df_stat = as.data.frame.matrix(table(df_score$markername,df_score$geno,useNA="ifany")) # make a contingency table into a df
dim(df_stat) # 116,713 x 6

df_stat$markername = rownames(df_stat)
rownames(df_stat) = seq(1 : length(rownames(df_stat)))
df_stat = df_stat[,c(7,1,2,3,4,5,6)]
colnames(df_stat) [2:7] <- c("nulliplex","simplex","duplex","triplex","quadruplex","no_call")
df_stat$nsamp <- rowSums(df_stat[,c("nulliplex","simplex","duplex","triplex","quadruplex","no_call")]) # add nb samples per markers / should be all 284
head(df_stat)
summary(df_stat)

# write stat file per probe
fwrite(df_stat,"fitpoly_270_stat_probes.txt",row.names=FALSE,quote=FALSE,sep = "\t") # fwrite for faster writing

### Nb of probes which passed fitPoly filters/scoring ?
#############################################
# produce ps file for no call probes
list_no_call = subset(df_stat$markername,df_stat$no_call == 270)
length(list_no_call) # 25,007 probes without call
write.table(list_no_call,file="fitpoly_270_list_no_call_probes.ps",row.names=FALSE,col.names="probeset_id",quote=FALSE,sep = "\t")

# produce ps file for scored probes
list_scored = subset(df_stat$markername,df_stat$no_call < 270)
length(list_scored) # 91,554 probes with call
write.table(list_scored,file="list_scored_probes.ps",sep="\t",row.names=FALSE,col.names="probeset_id",quote=FALSE)

### Proportion of no call genotypes for passing markers at the probe level
#############################################
df_stat_called = subset(df_stat,df_stat$no_call < 270)
head(df_stat_called)

tiff("fitpoly_270_no_Call_in_passing_probes_ggplot.tiff",width = 2300, height = 2300,res=300)
ggplot(df_stat_called, aes(x=no_call /nsamp*100)) + xlim(0,100)+
        geom_histogram(binwidth=1)+xlab("%NA per probe")+
        geom_vline(xintercept = 25,color="red",size = 1.5)+
        annotate(geom="text", x=60, y=0, label="call.threshold = 0.75", color="red")+
          theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()

tiff("fitpoly_270_no_Call_in_passing_probes_ggplot_2.tiff",width = 2300, height = 2300,res=300)
ggplot(df_stat_called, aes(x=no_call /nsamp*100)) +
    geom_histogram(binwidth=1)+xlab("%NA per probe")+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()

quantile(df_stat_called$no_call/df_stat_called$nsamp*100)
min(df_stat_called$no_call/df_stat_called$nsamp*100) # 0%
max(df_stat_called$no_call/df_stat_called$nsamp*100) # 24.8%
median(df_stat_called$no_call/df_stat_called$nsamp*100) #5.2%

### Are all passing SNP polymorphs ? are all SNP with heteroz calls ?
#############################################

# all nuliplex / none
n = length(subset(df_stat_called$markername,df_stat_called$nulliplex > 0 & df_stat_called$simplex == 0 & df_stat_called$duplex ==0 & df_stat_called$triplex ==0 & df_stat_called$quadruplex ==0))
n # 337
# all simplex / none
s = length(subset(df_stat_called$markername,df_stat_called$nulliplex == 0 & df_stat_called$simplex > 0 & df_stat_called$duplex ==0 & df_stat_called$triplex ==0 & df_stat_called$quadruplex ==0))
s # 32
# all duplex / none
d = length(subset(df_stat_called$markername,df_stat_called$nulliplex == 0 & df_stat_called$simplex == 0 & df_stat_called$duplex > 0 & df_stat_called$triplex ==0 & df_stat_called$quadruplex ==0))
d # 35
# all triplex / none
t = length(subset(df_stat_called$markername,df_stat_called$nulliplex == 0 & df_stat_called$simplex == 0 & df_stat_called$duplex == 0 & df_stat_called$triplex > 0 & df_stat_called$quadruplex ==0))
t # 28
# all quadruplex / none
q = length(subset(df_stat_called$markername,df_stat_called$nulliplex == 0 & df_stat_called$simplex == 0 & df_stat_called$duplex == 0 & df_stat_called$triplex == 0 & df_stat_called$quadruplex > 0))
q # 220

# SNP with no heteroz call / none
hom = length(subset(df_stat_called$markername,  df_stat_called$simplex == 0 & df_stat_called$duplex ==0 & df_stat_called$triplex ==0 ))
hom # 559
# SNP with  heteroz call / all SNP have heteroz calls
het = length(subset(df_stat_called$markername, df_stat_called$simplex > 0 | df_stat_called$duplex > 0 | df_stat_called$triplex > 0 ))
het # 91448

### Export summary of stat at probe level
#############################################
df_sum_stat_probes =data.frame(nulliplex_only = n,
                               simplex_only = s,
                               duplex_only = d,
                               triplex_only = t,
                               quadruplex_only = q,
                               all_homoz = hom,
                               heteroz = het,
                               scored_probes = length(list_scored),
                               no_call_probes = length(list_no_call),
                               scored_probes_less_5_miss = length(subset(df_stat_called$markername,df_stat_called$no_call/df_stat_called$nsamp < 5/100 )),
                               scored_probes_less_10_miss = length(subset(df_stat_called$markername,df_stat_called$no_call/df_stat_called$nsamp < 10/100 )),
                               scored_probes_less_20_miss = length(subset(df_stat_called$markername,df_stat_called$no_call/df_stat_called$nsamp < 20/100 )),
                               scored_probes_less_30_miss = length(subset(df_stat_called$markername,df_stat_called$no_call/df_stat_called$nsamp < 30/100 )),
                               scored_probes_less_40_miss = length(subset(df_stat_called$markername,df_stat_called$no_call/df_stat_called$nsamp < 40/100 ))               )

write.table(df_sum_stat_probes,"fitpoly_270_sum_stat_probes.txt",quote=F,row.names = F)

rm(n,s,d,t,q,hom,het)

### % heterozygocy vs genome position
#############################################

# bring position information
df_stat_called$percent_het_call <- (df_stat_called$simplex + df_stat_called$duplex + df_stat_called$triplex)/(df_stat_called$nsamp - df_stat_called$no_call) * 100
head(df_stat_called)

dim(subset(df_stat_called,df_stat_called$percent_het_call==100)) # 10,195 probes with only heterozygote probes

tiff("fitpoly_270_percent_heteroz_call.tiff",width = 2300, height = 2300,res=300)
ggplot(df_stat_called, aes(x=percent_het_call)) + 
    geom_histogram(binwidth=1)+xlab("%heterozygous calls per probe (i.e. simplex, duplex, triplex)")+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()

df_positions = fread("probes_SNPs_positions_genome_K5_P_Q_names.txt",header=T,stringsAsFactors=T,sep="\t")
head(df_positions) # need to use "unique_name" for merging both

colnames(df_positions)
colnames(df_stat_called)[1] <- "unique_name"

df_stat_called_positions = merge(df_stat_called,df_positions)
df_stat_called_positions_OB = df_stat_called_positions[which(!is.na(df_stat_called_positions$Chromosome_OB )),]

# cumulative position on OB genome

df_stat_called_positions_OB <- df_stat_called_positions_OB[order(df_stat_called_positions_OB$Chromosome_OB,df_stat_called_positions_OB$Position_OB_genome_bp),]
df_stat_called_positions_OB$Position_OB_genome_bp_cumulative <- NA
rownames(df_stat_called_positions_OB) <- seq(1:length(rownames(df_stat_called_positions_OB)))

max_ch1 = max(df_stat_called_positions_OB[df_stat_called_positions_OB$Chromosome_OB == "Chr01","Position_OB_genome_bp"],na.rm=T)
max_ch2 = max(df_stat_called_positions_OB[df_stat_called_positions_OB$Chromosome_OB == "Chr02","Position_OB_genome_bp"],na.rm=T)
max_ch3 = max(df_stat_called_positions_OB[df_stat_called_positions_OB$Chromosome_OB == "Chr03","Position_OB_genome_bp"],na.rm=T)
max_ch4 = max(df_stat_called_positions_OB[df_stat_called_positions_OB$Chromosome_OB == "Chr04","Position_OB_genome_bp"],na.rm=T)
max_ch5 = max(df_stat_called_positions_OB[df_stat_called_positions_OB$Chromosome_OB == "Chr05","Position_OB_genome_bp"],na.rm=T)
max_ch6 = max(df_stat_called_positions_OB[df_stat_called_positions_OB$Chromosome_OB == "Chr06","Position_OB_genome_bp"],na.rm=T)
max_ch7 = max(df_stat_called_positions_OB[df_stat_called_positions_OB$Chromosome_OB == "Chr07","Position_OB_genome_bp"],na.rm=T)


df_stat_called_positions_OB$Position_OB_genome_bp_cumulative <- NA
df_stat_called_positions_OB$Position_OB_genome_bp_cumulative[which(df_stat_called_positions_OB$Chromosome_OB == "Chr01")] <- df_stat_called_positions_OB$Position_OB_genome_bp[which(df_stat_called_positions_OB$Chromosome_OB == "Chr01")] 
df_stat_called_positions_OB$Position_OB_genome_bp_cumulative[which(df_stat_called_positions_OB$Chromosome_OB == "Chr02")] <- df_stat_called_positions_OB$Position_OB_genome_bp[which(df_stat_called_positions_OB$Chromosome_OB == "Chr02")]  + max_ch1
df_stat_called_positions_OB$Position_OB_genome_bp_cumulative[which(df_stat_called_positions_OB$Chromosome_OB == "Chr03")] <- df_stat_called_positions_OB$Position_OB_genome_bp[which(df_stat_called_positions_OB$Chromosome_OB == "Chr03")]  + max_ch1+max_ch2
df_stat_called_positions_OB$Position_OB_genome_bp_cumulative[which(df_stat_called_positions_OB$Chromosome_OB == "Chr04")] <- df_stat_called_positions_OB$Position_OB_genome_bp[which(df_stat_called_positions_OB$Chromosome_OB == "Chr04")]  + max_ch1+max_ch2+max_ch3
df_stat_called_positions_OB$Position_OB_genome_bp_cumulative[which(df_stat_called_positions_OB$Chromosome_OB == "Chr05")] <- df_stat_called_positions_OB$Position_OB_genome_bp[which(df_stat_called_positions_OB$Chromosome_OB == "Chr05")]  + max_ch1+max_ch2+max_ch3+max_ch4
df_stat_called_positions_OB$Position_OB_genome_bp_cumulative[which(df_stat_called_positions_OB$Chromosome_OB == "Chr06")] <- df_stat_called_positions_OB$Position_OB_genome_bp[which(df_stat_called_positions_OB$Chromosome_OB == "Chr06")]  + max_ch1+max_ch2+max_ch3+max_ch4+max_ch5
df_stat_called_positions_OB$Position_OB_genome_bp_cumulative[which(df_stat_called_positions_OB$Chromosome_OB == "Chr07")] <- df_stat_called_positions_OB$Position_OB_genome_bp[which(df_stat_called_positions_OB$Chromosome_OB == "Chr07")]  + max_ch1+max_ch2+max_ch3+max_ch4+max_ch5+max_ch6
df_stat_called_positions_OB$Position_OB_genome_bp_cumulative[which(df_stat_called_positions_OB$Chromosome_OB == "Chr00")] <- df_stat_called_positions_OB$Position_OB_genome_bp[which(df_stat_called_positions_OB$Chromosome_OB == "Chr00")]  + max_ch1+max_ch2+max_ch3+max_ch4+max_ch5+max_ch6+max_ch7

head(df_stat_called_positions_OB)

# get middle cumulative positions
min_cum_ch1 = min(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr01"),"Position_OB_genome_bp_cumulative"])
max_cum_ch1 = max(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr01"),"Position_OB_genome_bp_cumulative"])
center_cum_ch1 = (min_cum_ch1+max_cum_ch1)/2

min_cum_ch2 = min(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr02"),"Position_OB_genome_bp_cumulative"])
max_cum_ch2 = max(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr02"),"Position_OB_genome_bp_cumulative"])
center_cum_ch2 = (min_cum_ch2+max_cum_ch2)/2

min_cum_ch3 = min(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr03"),"Position_OB_genome_bp_cumulative"])
max_cum_ch3 = max(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr03"),"Position_OB_genome_bp_cumulative"])
center_cum_ch3 = (min_cum_ch3+max_cum_ch3)/2

min_cum_ch4 = min(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr04"),"Position_OB_genome_bp_cumulative"])
max_cum_ch4 = max(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr04"),"Position_OB_genome_bp_cumulative"])
center_cum_ch4 = (min_cum_ch4+max_cum_ch4)/2

min_cum_ch5 = min(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr05"),"Position_OB_genome_bp_cumulative"])
max_cum_ch5 = max(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr05"),"Position_OB_genome_bp_cumulative"])
center_cum_ch5 = (min_cum_ch5+max_cum_ch5)/2

min_cum_ch6 = min(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr06"),"Position_OB_genome_bp_cumulative"])
max_cum_ch6 = max(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr06"),"Position_OB_genome_bp_cumulative"])
center_cum_ch6 = (min_cum_ch6+max_cum_ch6)/2

min_cum_ch7 = min(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr07"),"Position_OB_genome_bp_cumulative"])
max_cum_ch7 = max(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr07"),"Position_OB_genome_bp_cumulative"])
center_cum_ch7 = (min_cum_ch7+max_cum_ch7)/2

min_cum_ch0 = min(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr00"),"Position_OB_genome_bp_cumulative"])
max_cum_ch0 = max(df_stat_called_positions_OB[which(df_stat_called_positions_OB$Chromosome_OB=="Chr00"),"Position_OB_genome_bp_cumulative"])
center_cum_ch0 = (min_cum_ch0+max_cum_ch0)/2


# cumulative heterozygote percentage / compute rolling means with different window sizes
df_stat_called_positions_OB_cumul <- df_stat_called_positions_OB %>%
    select(Position_OB_genome_bp_cumulative, percent_het_call) %>%
    mutate(percent_het_call_w10 = rollmean(percent_het_call, k = 10, fill = NA),
           percent_het_call_w50 = rollmean(percent_het_call, k = 50, fill = NA),
           percent_het_call_w100 = rollmean(percent_het_call, k = 100, fill = NA),
           percent_het_call_w200 = rollmean(percent_het_call, k = 200, fill = NA),
           percent_het_call_w1000 = rollmean(percent_het_call, k = 1000, fill = NA),
           percent_het_call_w10000 = rollmean(percent_het_call, k = 10000, fill = NA))

head(df_stat_called_positions_OB_cumul)
summary(df_stat_called_positions_OB_cumul)

# plot window 1000
tiff("fitpoly_270_heterozygote_call_on_genome_sliding_window.tiff",width = 4000, height = 2000,res=300)
ggplot(data=df_stat_called_positions_OB_cumul,aes(Position_OB_genome_bp_cumulative/10^6, percent_het_call_w1000)) +
    geom_line()+ylim(0,100)+ylab("% heterozygous call (window 1000 probes)")+
        geom_vline(xintercept = max_ch1/10^6, color = "black", size=0.5)+
        geom_vline(xintercept = (max_ch1+max_ch2)/10^6, color = "black", size=0.5)+
        geom_vline(xintercept = (max_ch1+max_ch2+max_ch3)/10^6, color = "black", size=0.5)+
        geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4)/10^6, color = "black", size=0.5)+
        geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4+max_ch5)/10^6, color = "black", size=0.5)+
        geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4+max_ch5+max_ch6)/10^6, color = "black", size=0.5)+
        geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4+max_ch5+max_ch6+max_ch7)/10^6, color = "black", size=0.5)+
        scale_x_continuous("positions OB genome",breaks= c(center_cum_ch1,center_cum_ch2,center_cum_ch3,center_cum_ch4,center_cum_ch5,center_cum_ch6,center_cum_ch7,center_cum_ch0)/10^6, labels=c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr00"))+
        theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()    
    
    
### % NA vs genome position
#############################################   
    
# cumulative heterozygote percentage / rolling means
df_stat_called_positions_OB$percent_NA <- df_stat_called_positions_OB$no_call/df_stat_called_positions_OB$nsamp*100 

df_stat_called_positions_OB_cumul <- df_stat_called_positions_OB %>%
    select(Position_OB_genome_bp_cumulative, percent_NA) %>%
    mutate(percent_NA_w10 = rollmean(percent_NA, k = 10, fill = NA),
           percent_NA_w50 = rollmean(percent_NA, k = 50, fill = NA),
           percent_NA_w100 = rollmean(percent_NA, k = 100, fill = NA),
           percent_NA_w200 = rollmean(percent_NA, k = 200, fill = NA),
           percent_NA_w1000 = rollmean(percent_NA, k = 1000, fill = NA),
           percent_NA_w10000 = rollmean(percent_NA, k = 10000, fill = NA))

head(df_stat_called_positions_OB_cumul)
summary(df_stat_called_positions_OB_cumul)

# plot window 50
tiff("NA_call_on_genome_sliding_window.tiff",width = 4000, height = 2000,res=300)
ggplot(data=df_stat_called_positions_OB_cumul,aes(Position_OB_genome_bp_cumulative/10^6, percent_NA_w50)) +
    geom_line(col="red")+theme_bw()+ylim(0,100)+ylab("% NA call (window 50 probes)")+
    geom_vline(xintercept = max_ch1/10^6, color = "black", size=0.5)+ggtitle("% NA call per probe")+
    geom_vline(xintercept = (max_ch1+max_ch2)/10^6, color = "black", size=0.5)+
    geom_vline(xintercept = (max_ch1+max_ch2+max_ch3)/10^6, color = "black", size=0.5)+
    geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4)/10^6, color = "black", size=0.5)+
    geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4+max_ch5)/10^6, color = "black", size=0.5)+
    geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4+max_ch5+max_ch6)/10^6, color = "black", size=0.5)+
    geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4+max_ch5+max_ch6+max_ch7)/10^6, color = "black", size=0.5)+
    scale_x_continuous("positions OB genome",breaks= c(center_cum_ch1,center_cum_ch2,center_cum_ch3,center_cum_ch4,center_cum_ch5,center_cum_ch6,center_cum_ch7,center_cum_ch0)/10^6, labels=c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr00"))+
        theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()       


### compute MAF 
#############################################   
head(df_stat_called_positions_OB)

df_stat_called_positions_OB$nb_A_allele = df_stat_called_positions_OB$nulliplex * 4 +
    df_stat_called_positions_OB$simplex * 3 +
    df_stat_called_positions_OB$duplex * 2 +
    df_stat_called_positions_OB$triplex * 1

df_stat_called_positions_OB$nb_B_allele = df_stat_called_positions_OB$simplex * 1 +
    df_stat_called_positions_OB$duplex * 2 +
    df_stat_called_positions_OB$triplex * 3 +
    df_stat_called_positions_OB$quadruplex * 4

df_stat_called_positions_OB$sum_alleles = df_stat_called_positions_OB$nb_A_allele+df_stat_called_positions_OB$nb_B_allele

head(df_stat_called_positions_OB)

df_stat_called_positions_OB$freq_A_allele = df_stat_called_positions_OB$nb_A_allele/(df_stat_called_positions_OB$nb_A_allele+df_stat_called_positions_OB$nb_B_allele)
df_stat_called_positions_OB$freq_B_allele = df_stat_called_positions_OB$nb_B_allele/(df_stat_called_positions_OB$nb_A_allele+df_stat_called_positions_OB$nb_B_allele)

df_stat_called_positions_OB$MAF <- df_stat_called_positions_OB$freq_A_allele
df_stat_called_positions_OB$MAF[df_stat_called_positions_OB$MAF>0.50]<-df_stat_called_positions_OB$freq_B_allele[df_stat_called_positions_OB$MAF>0.50]

summary(df_stat_called_positions_OB$MAF)

tiff("fitpoly_270_MAF_distribution.tiff",width = 2300, height = 2300,res=300)
ggplot(df_stat_called_positions_OB, aes(x=MAF)) + 
    geom_histogram(binwidth=0.005)+xlab("Minor Allele Frequency (MAF)")+
    geom_vline(xintercept =0.05, color = "red", size=1)+
       theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()


# cumulative MAF percentage / rolling means
df_stat_called_positions_OB_cumul <- df_stat_called_positions_OB %>%
    select(Position_OB_genome_bp_cumulative, MAF) %>%
    mutate(MAF_w10 = rollmean(MAF, k = 10, fill = NA),
           MAF_w50 = rollmean(MAF, k = 50, fill = NA),
           MAF_w100 = rollmean(MAF, k = 100, fill = NA),
           MAF_w200 = rollmean(MAF, k = 200, fill = NA),
           MAF_w1000 = rollmean(MAF, k = 1000, fill = NA),
           MAF_w10000 = rollmean(MAF, k = 10000, fill = NA))

head(df_stat_called_positions_OB_cumul)
summary(df_stat_called_positions_OB_cumul)

# plot window 50
tiff("fitpoly_270_MAF_on_genome_sliding_window.tiff",width = 1800, height = 800,res=200)
ggplot(data=df_stat_called_positions_OB_cumul,aes(Position_OB_genome_bp_cumulative/10^6, MAF_w50)) +
    geom_line(col="blue")+theme_bw()+ylim(0,1)+ylab("MAF (window 50 probes)")+
    geom_vline(xintercept = max_ch1/10^6, color = "black", size=0.5)+
    geom_vline(xintercept = (max_ch1+max_ch2)/10^6, color = "black", size=0.5)+
    geom_vline(xintercept = (max_ch1+max_ch2+max_ch3)/10^6, color = "black", size=0.5)+
    geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4)/10^6, color = "black", size=0.5)+
    geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4+max_ch5)/10^6, color = "black", size=0.5)+
    geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4+max_ch5+max_ch6)/10^6, color = "black", size=0.5)+
    geom_vline(xintercept = (max_ch1+max_ch2+max_ch3+max_ch4+max_ch5+max_ch6+max_ch7)/10^6, color = "black", size=0.5)+
    scale_x_continuous("positions OB genome",breaks= c(center_cum_ch1,center_cum_ch2,center_cum_ch3,center_cum_ch4,center_cum_ch5,center_cum_ch6,center_cum_ch7,center_cum_ch0)/10^6, labels=c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr00"))+
    theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()    


dim(subset(df_stat_called_positions_OB_cumul,df_stat_called_positions_OB_cumul$MAF>=0.05))
# 68283

### Nb of passing probes per chromosome
#############################################

dim(subset(df_stat_called_positions_OB,df_stat_called_positions_OB$Chromosome_OB =="Chr00")) # 3919   
dim(subset(df_stat_called_positions_OB,df_stat_called_positions_OB$Chromosome_OB =="Chr01")) # 9969
dim(subset(df_stat_called_positions_OB,df_stat_called_positions_OB$Chromosome_OB =="Chr02")) # 15758
dim(subset(df_stat_called_positions_OB,df_stat_called_positions_OB$Chromosome_OB =="Chr03")) # 8939 
dim(subset(df_stat_called_positions_OB,df_stat_called_positions_OB$Chromosome_OB =="Chr04")) # 10778    
dim(subset(df_stat_called_positions_OB,df_stat_called_positions_OB$Chromosome_OB =="Chr05")) # 12706    
dim(subset(df_stat_called_positions_OB,df_stat_called_positions_OB$Chromosome_OB =="Chr06")) # 13807    
dim(subset(df_stat_called_positions_OB,df_stat_called_positions_OB$Chromosome_OB =="Chr07")) # 12801    

table(df_stat_called_positions$Chromosome_OB,useNA = "ifany")

### Filter df_score to keep only scored probes
#############################################

length(levels(df_score$markername)) # 116,713
length(unique(list_scored)) # 92,007

df_score_pass = df_score[which(df_score$markername %in% list_scored),]
dim(df_score_pass) # 24841890    15

head(df_score_pass)
str(df_score_pass)


# need to reset levels for markers
df_score_pass[,c(1,2)] <- lapply(df_score_pass[,c(1,2)], as.character) # 
df_score_pass[,c(1,2)] <- lapply(df_score_pass[,c(1,2)], factor) # 

length(levels(df_score_pass$marker))
length(levels(df_score_pass$markername))

# need to reset rownames
rownames(df_score_pass) <- seq(1:length(rownames(df_score_pass)))

# exporting scores for passing probes only in format adapted to SNPolisher
fwrite(df_score_pass,"fitPoly_270_scores_formated_only_passing.dat",quote=F,row.names = F,sep = "\t") # fwrite for faster writing

### STATS AT THE GENOTYPE LEVEL
#############################################
#############################################
#############################################

### Stat file at genotype level 
#############################################
df_indiv = as.data.frame.matrix(table(df_score_pass$sample,df_score_pass$geno,useNA="ifany")) # make a contingency table into a df
dim(df_indiv)

head(df_indiv)

df_indiv$sample = rownames(df_indiv)
rownames(df_indiv) = seq(1 : length(rownames(df_indiv)))
df_indiv = df_indiv[,c(7,1,2,3,4,5,6)]
colnames(df_indiv) [2:7] <- c("nulliplex","simplex","duplex","triplex","quadruplex","no_call")
df_indiv$nprob <- rowSums(df_indiv[,c("nulliplex","simplex","duplex","triplex","quadruplex","no_call")]) # add nb probes per sample / should be all 4529
df_indiv$het_call = rowSums(df_indiv[,c("simplex","duplex","triplex")])

head(df_indiv)
summary(df_indiv)

# write stat file per genotype
fwrite(df_indiv,"fitpoly_270_stat_genotypes.txt",row.names=FALSE,quote=FALSE,sep="\t") # fwrite for fast writing

### Percentage of het loci per genotype
#############################################

tiff("fitpoly_270_percent_het_loci_per_indiv_distribution.tiff",width = 2300, height = 2300,res=300)
ggplot(df_indiv, aes(x=(het_call/nprob)*100)) + 
    geom_histogram(binwidth=1)+xlab("% heterozygous calls per sample")+xlim(0,100)+
    theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()

summary(df_indiv$het_call/df_indiv$nprob*100)
min(df_indiv$het_call/df_indiv$nprob)
mean(df_indiv$het_call/df_indiv$nprob)
max(df_indiv$het_call/df_indiv$nprob)
median(df_indiv$het_call/df_indiv$nprob)


df_sum_stat_genotype_heteroz =data.frame(min = min(df_indiv$het_call/df_indiv$nprob),
                               mean = mean(df_indiv$het_call/df_indiv$nprob),
                               max = max(df_indiv$het_call/df_indiv$nprob),
				sd = sd(df_indiv$het_call/df_indiv$nprob),
                               median = median(df_indiv$het_call/df_indiv$nprob))
                              
# write stat file per genotype
fwrite(df_sum_stat_genotype_heteroz,"fitpoly_270_sum_stat_genotype_heteroz.txt",row.names=FALSE,quote=FALSE,sep="\t") # fwrite for fast writing

### Percentage NA per genotype
#############################################

tiff("fitpoly_270_NA_per_indiv_distribution.tiff",width = 2300, height = 2300,res=300)
ggplot(df_indiv, aes(x=(no_call /nprob)*100)) + 
    geom_histogram(binwidth=1)+theme_bw()+xlab("% NA calls per sample")+xlim(0,25)+
    theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.line=element_line(colour="black",size=1.25),axis.ticks=element_line(colour="black",size=1.25),
        axis.text.x=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y=element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"))
dev.off()

summary(df_indiv$no_call/df_indiv$nprob*100)


###################################################################################################################################################
# END SCRIPT 31/08/21
####################################################################################################################################################


### Genotyping call per genotypes 
#############################################
## proportion nulliplex, simplex, duplex, triplex, quadruplex, no call across genotyped probes per sample


## showing as a barplot, but need to tidy data before
df_indiv_tidy <- df_indiv %>% pivot_longer(c("nulliplex", "simplex","duplex","triplex","quadruplex","no_call"), names_to = "score", values_to = "nb_probes")
df_indiv_tidy=as.data.frame(df_indiv_tidy)
head(df_indiv_tidy)
str(df_indiv_tidy)

df_indiv_tidy[,c(1,4)] <- lapply(df_indiv_tidy[,c(1,4)], factor) # 
df_indiv_tidy$score <- factor(df_indiv_tidy$score,levels(df_indiv_tidy$score)[c(3,5,1,6,4,2)])

# load ploidy levels
df_ploidy = read.table("/mnt/ufs18/home-057/alberte2/ROSE_SNP_ARRAY/fitPoly_FILTERING/Axiom_output_Rose_188_samples/samples_ploidy_original.txt",header=T)
df_ploidy$ploidy=as.factor(df_ploidy$ploidy)
summary(df_ploidy)

df_indiv_tidy = merge(df_indiv_tidy,df_ploidy)
head(df_indiv_tidy)


# order by no call
temp=subset(df_indiv_tidy,df_indiv_tidy$score=="no_call")
head(temp)
temp <- temp[order(temp$ploidy,temp$nb_probes),]

df_indiv_tidy$sample <- factor(df_indiv_tidy$sample, levels=temp$sample)

colors <- temp$ploidy

tiff("scored_probes_per_class_5.tiff", width = 3200, height = 1200,res=200, type="cairo")
ggplot(data=df_indiv_tidy, aes(x=sample , y=nb_probes, fill=score)) +
        geom_bar(stat="identity",width = 1)+theme_minimal()+
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=as.numeric(colors)))+
                ggtitle("Scored probes per allele dosage class") +xlab("188 genotypes")
dev.off()


# order by nulli
temp=subset(df_indiv_tidy,df_indiv_tidy$score=="nulliplex")
head(temp)
temp <- temp[order(temp$ploidy,temp$nb_probes),]

df_indiv_tidy$sample <- factor(df_indiv_tidy$sample, levels=temp$sample)

colors <- temp$ploidy

tiff("scored_probes_per_class_0.tiff", width = 3200, height = 1200,res=200, type="cairo")

ggplot(data=df_indiv_tidy, aes(x=sample , y=nb_probes, fill=score)) +
        geom_bar(stat="identity",width = 1)+theme_minimal()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=as.numeric(colors)))+
        ggtitle("Scored probes per allele dosage class") +xlab("188 genotypes")
dev.off()



# order by simplex
temp=subset(df_indiv_tidy,df_indiv_tidy$score=="simplex")
head(temp)
temp <- temp[order(temp$ploidy,temp$nb_probes),]

df_indiv_tidy$sample <- factor(df_indiv_tidy$sample, levels=temp$sample)

colors <- temp$ploidy

tiff("scored_probes_per_class_1.tiff", width = 3200, height = 1200,res=200, type="cairo")

ggplot(data=df_indiv_tidy, aes(x=sample , y=nb_probes, fill=score)) +
        geom_bar(stat="identity",width = 1)+theme_minimal()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=as.numeric(colors)))+
        ggtitle("Scored probes per allele dosage class") +xlab("188 genotypes")
dev.off()


# order by duplex
temp=subset(df_indiv_tidy,df_indiv_tidy$score=="duplex")
head(temp)
temp <- temp[order(temp$ploidy,temp$nb_probes),]

df_indiv_tidy$sample <- factor(df_indiv_tidy$sample, levels=temp$sample)

colors <- temp$ploidy

tiff("scored_probes_per_class_2.tiff", width = 3200, height = 1200,res=200, type="cairo")

ggplot(data=df_indiv_tidy, aes(x=sample , y=nb_probes, fill=score)) +
        geom_bar(stat="identity",width = 1)+theme_minimal()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=as.numeric(colors)))+
        ggtitle("Scored probes per allele dosage class") +xlab("188 genotypes")
dev.off()


# order by triplex
temp=subset(df_indiv_tidy,df_indiv_tidy$score=="triplex")
head(temp)
temp <- temp[order(temp$ploidy,temp$nb_probes),]

df_indiv_tidy$sample <- factor(df_indiv_tidy$sample, levels=temp$sample)

colors <- temp$ploidy

tiff("scored_probes_per_class_3.tiff", width = 3200, height = 1200,res=200, type="cairo")

ggplot(data=df_indiv_tidy, aes(x=sample , y=nb_probes, fill=score)) +
        geom_bar(stat="identity",width = 1)+theme_minimal()+
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=as.numeric(colors)))+
        ggtitle("Scored probes per allele dosage class") +xlab("188 genotypes")
dev.off()

# order by quadruplex
temp=subset(df_indiv_tidy,df_indiv_tidy$score=="quadruplex")
head(temp)
temp <- temp[order(temp$ploidy,temp$nb_probes),]

df_indiv_tidy$sample <- factor(df_indiv_tidy$sample, levels=temp$sample)

colors <- temp$ploidy

tiff("scored_probes_per_class_4.tiff", width = 3200, height = 1200,res=200, type="cairo")

ggplot(data=df_indiv_tidy, aes(x=sample , y=nb_probes, fill=score)) +
        geom_bar(stat="identity",width = 1)+theme_minimal()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=as.numeric(colors)))+
        ggtitle("Scored probes per allele dosage class") +xlab("188 genotypes")
dev.off()


# average and sd of each call 
df_sum_stat_genotypes <- data.frame(sample=as.character(levels(df_indiv_tidy$score)),row.names=c(1 :length(levels(df_indiv_tidy$score))))
df_sum_stat_genotypes$mean <- tapply(df_indiv_tidy$nb_probes,df_indiv_tidy$score,mean,na.rm=T)
df_sum_stat_genotypes$median <-tapply(df_indiv_tidy$nb_probes,df_indiv_tidy$score,median,na.rm=T)
df_sum_stat_genotypes$sd <-tapply(df_indiv_tidy$nb_probes,df_indiv_tidy$score,sd,na.rm=T)
df_sum_stat_genotypes$min <-tapply(df_indiv_tidy$nb_probes,df_indiv_tidy$score,min,na.rm=T)
df_sum_stat_genotypes$max <-tapply(df_indiv_tidy$nb_probes,df_indiv_tidy$score,max,na.rm=T)

fwrite(df_sum_stat_genotypes,"df_sum_stat_genotypes.txt",quote=F,row.names = F,sep="\t")

### Genotyping call per genotypes / with expected ploidy levels
#############################################

# load ploidy levels
summary(df_ploidy)

list_diploids = subset(df_ploidy$sample,df_ploidy$ploidy==2)
list_diploids=as.factor(as.character(list_diploids))
levels(list_diploids)

list_tetraploids = subset(df_ploidy$sample,df_ploidy$ploidy==4)
list_tetraploids=as.factor(as.character(list_tetraploids))
levels(list_tetraploids)

# subset only diploids
############################
df_indiv_tidy_diploids = subset(df_indiv_tidy,df_indiv_tidy$sample %in% list_diploids)
dim(df_indiv_tidy_diploids)
df_indiv_tidy_diploids$sample=as.factor(as.character(df_indiv_tidy_diploids$sample))
       
temp=subset(df_indiv_tidy_diploids,df_indiv_tidy_diploids$score=="duplex")
head(temp)
temp <- temp[order(temp$nb_probes),]

df_indiv_tidy_diploids$sample <- factor(df_indiv_tidy_diploids$sample, levels=temp$sample)

colors <- temp$ploidy

tiff("scored_probes_per_class_diplo.tiff", width = 2200, height = 1200,res=200, type="cairo")

ggplot(data=df_indiv_tidy_diploids, aes(x=sample , y=nb_probes, fill=score)) +
        geom_bar(stat="identity",width = 1)+theme_minimal()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=colors))+
        ggtitle("Scored probes per allele dosage class") +xlab("58 genotypes supposed to be diploids")
dev.off()


# average and sd of each call in the diploids
df_sum_stat_genotypes_diplo <- data.frame(sample=as.character(levels(df_indiv_tidy_diploids$score)),row.names=c(1 :length(levels(df_indiv_tidy_diploids$score))))
df_sum_stat_genotypes_diplo$mean <- tapply(df_indiv_tidy_diploids$nb_probes,df_indiv_tidy_diploids$score,mean,na.rm=T)
df_sum_stat_genotypes_diplo$median <-tapply(df_indiv_tidy_diploids$nb_probes,df_indiv_tidy_diploids$score,median,na.rm=T)
df_sum_stat_genotypes_diplo$sd <-tapply(df_indiv_tidy_diploids$nb_probes,df_indiv_tidy_diploids$score,sd,na.rm=T)
df_sum_stat_genotypes_diplo$min <-tapply(df_indiv_tidy_diploids$nb_probes,df_indiv_tidy_diploids$score,min,na.rm=T)
df_sum_stat_genotypes_diplo$max <-tapply(df_indiv_tidy_diploids$nb_probes,df_indiv_tidy_diploids$score,max,na.rm=T)

fwrite(df_sum_stat_genotypes_diplo,"df_sum_stat_genotypes_diplo.txt",quote=F,row.names = F,sep="\t")

# subset only tetraploids
#############################
df_indiv_tidy_tetraploids = subset(df_indiv_tidy,df_indiv_tidy$sample %in% list_tetraploids)
dim(df_indiv_tidy_tetraploids)
df_indiv_tidy_tetraploids$sample=as.factor(as.character(df_indiv_tidy_tetraploids$sample))

temp=subset(df_indiv_tidy_tetraploids,df_indiv_tidy_tetraploids$score=="duplex")
head(temp)
temp <- temp[order(temp$nb_probes),]

colors <- temp$ploidy

df_indiv_tidy_tetraploids$sample <- factor(df_indiv_tidy_tetraploids$sample, levels=temp$sample)

tiff("scored_probes_per_class_tetra.tiff", width = 2200, height = 1200,res=200, type="cairo")

ggplot(data=df_indiv_tidy_tetraploids, aes(x=sample , y=nb_probes, fill=score)) +
        geom_bar(stat="identity",width = 1)+theme_minimal()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=colors))+
        ggtitle("Scored probes per allele dosage class") +xlab("134 genotypes supposed to be tetraploids")

dev.off()

# average and sd of each call in the diploids
df_sum_stat_genotypes_tetra <- data.frame(sample=as.character(levels(df_indiv_tidy_tetraploids$score)),row.names=c(1 :length(levels(df_indiv_tidy_tetraploids$score))))
df_sum_stat_genotypes_tetra$mean <- tapply(df_indiv_tidy_tetraploids$nb_probes,df_indiv_tidy_tetraploids$score,mean,na.rm=T)
df_sum_stat_genotypes_tetra$median <-tapply(df_indiv_tidy_tetraploids$nb_probes,df_indiv_tidy_tetraploids$score,median,na.rm=T)
df_sum_stat_genotypes_tetra$sd <-tapply(df_indiv_tidy_tetraploids$nb_probes,df_indiv_tidy_tetraploids$score,sd,na.rm=T)
df_sum_stat_genotypes_tetra$min <-tapply(df_indiv_tidy_tetraploids$nb_probes,df_indiv_tidy_tetraploids$score,min,na.rm=T)
df_sum_stat_genotypes_tetra$max <-tapply(df_indiv_tidy_tetraploids$nb_probes,df_indiv_tidy_tetraploids$score,max,na.rm=T)

fwrite(df_sum_stat_genotypes_tetra,"df_sum_stat_genotypes_tetra.txt",quote=F,row.names = F,sep="\t")
