## Modified by TL - 290621 from April 2021 Elise - PCoA and DAPC on genotyping results from fitPoly with all 96 Delbard genotypes as active samples
## 185 tetraploids 

### "matrix_185_tetra_samples_scores_passing_probes.txt" is the output from script "script_generating_output_as_matrix_185_tetra.R".
### It corresponds to fitPoly results formated as a matrix. 

##### Load packages
#########################
library(plyr)  # data manipulation
library(tidyr) # data manipulation
library(adegenet) # genetic analysis
library(ade4) # PCA, DAPC and others multivariate analysis
library(ggplot2) # plots
library(ape) # tree and phylogeny
library(pegas) # population genetic
library(seqinr) # work with sequence data
library(BBmisc) # tools to normalize data
library(hierfstat)
library(devtools)
library(corrplot)

setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/GWAS_204cc/")

##### Load matrix of probes/snps
###########################################################################
###########################################################################
df_snp = read.table("matrix_fitPoly_270_scores_formated_only_passing_probes.205samplescc",h=T,sep="\t")
str(df_snp)
dim(df_snp) # [1] 92007   205 => 92 007 probes by 204 genotypes 

df_snp=as.data.frame(df_snp)
df_snp[1:5,1:5]
df_snp[1:5,95:97]
dim(df_snp)

##### STEP 1 / Filter data at individual level and SNP level 
##### Individuals with too much missing data need to be filtered out.
##### Monomorphic probes, low MAF probes and probes with too much missing data will make troubles in PCA and scaling [sd = 0 for monomorphic probes] and need to be filtered out. 
###########################################################################
###########################################################################

# Change data to long format to produce stat files
df_snp_long = gather(df_snp, Genotype, Score, adam_F4:wodan_F5, factor_key=TRUE)
head(df_snp_long)


# A ##########################################################################
# Summary stats to identify individuals with excess missing (> 20%)
df_stat_indiv = as.data.frame.matrix(table(df_snp_long$Genotype,df_snp_long$Score,useNA="ifany"))
head(df_stat_indiv)
colnames(df_stat_indiv) <- c("A/A/A/A","B/A/A/A","B/B/A/A","B/B/B/A","B/B/B/B","missing")
df_stat_indiv$Total_SNP = df_stat_indiv$'A/A/A/A' + df_stat_indiv$'B/A/A/A' + df_stat_indiv$'B/B/A/A' + df_stat_indiv$'B/B/B/A'+df_stat_indiv$'B/B/B/B'+df_stat_indiv$missing

# Filter individuals with more than 20% missing #########################
df_stat_indiv$percent_NA = df_stat_indiv$missing/df_stat_indiv$Total_SNP * 100
max(df_stat_indiv$percent_NA) # 16.88024

hist(df_stat_indiv$percent_NA,main="",xlab="NA% per individual",xlim=c(0,100),nclass=200,col="darkblue")
abline(v=20,col="red",lwd=2) 
# not any samples with 20% or more missing data, thus we can keep all samples

pdf(file="missingness_205cc_111021.pdf",height=8,width=8)
ggplot()+geom_density(aes(df_stat_indiv$percent_NA),fill="grey",adjust=0.5)+
  xlab("NA% per individual")+
  xlim(0,100)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_vline(xintercept=20,col="red",lwd=1)
dev.off()

keep_genotypes = subset(rownames(df_stat_indiv),df_stat_indiv$percent_NA <= 20) # list of genotypes to keep
length(keep_genotypes) # 204 genotypes stay in the analysis / we keep all genotypes here

head(df_snp_long)
df_snp_long_filter = df_snp_long[which(df_snp_long$Genotype %in% keep_genotypes),]
dim(df_snp_long) #  18769428        3
dim(df_snp_long_filter) #  18769428        3
# no difference between df_snp_long and df_snp_long_filter here because we kept all genotypes !

# B ##########################################################################
# Summary stats to identify monomorphic probes ######################### start from data with individuals filtered out
df_stat_snp = as.data.frame.matrix(table(df_snp_long_filter$markername,df_snp_long_filter$Score,useNA="ifany"))
head(df_stat_snp)
colnames(df_stat_snp) <- c("A/A/A/A","B/A/A/A","B/B/A/A","B/B/B/A","B/B/B/B","missing")

df_stat_snp$status <- "polymorphic"

# Filter monomorphic probes 
# mono BBBB / 239
df_stat_snp[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'> 0),]
dim(df_stat_snp[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'> 0),])
df_stat_snp$status[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'> 0)] <- "monomorphic"

# mono BBBA / 24
df_stat_snp[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'>0 & df_stat_snp$'B/B/B/B'== 0),]
dim(df_stat_snp[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'>0 & df_stat_snp$'B/B/B/B'== 0),])
df_stat_snp$status[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'>0 & df_stat_snp$'B/B/B/B'== 0)] <- "monomorphic"

# mono BBAA / 22
df_stat_snp[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'>0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'== 0),]
dim(df_stat_snp[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'>0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'== 0),])
df_stat_snp$status[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'>0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'== 0)] <- "monomorphic"

# mono BAAA / 31
df_stat_snp[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'>0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'== 0),]
dim(df_stat_snp[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'>0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'== 0),])
df_stat_snp$status[which(df_stat_snp$'A/A/A/A'==0 & df_stat_snp$'B/A/A/A'>0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'== 0)] <- "monomorphic"

# mono AAAA / 429
df_stat_snp[which(df_stat_snp$'A/A/A/A'>0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'== 0),]
dim(df_stat_snp[which(df_stat_snp$'A/A/A/A'>0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'== 0),])
df_stat_snp$status[which(df_stat_snp$'A/A/A/A'>0 & df_stat_snp$'B/A/A/A'==0 & df_stat_snp$'B/B/A/A'==0 & df_stat_snp$'B/B/B/A'==0 & df_stat_snp$'B/B/B/B'== 0)] <- "monomorphic"

head(df_stat_snp)
table(df_stat_snp$status) # 734 monomorphic

# Filter probes with > 10% missing
df_stat_snp$Total_indiv = df_stat_snp$'A/A/A/A' + df_stat_snp$'B/A/A/A' + df_stat_snp$'B/B/A/A' + df_stat_snp$'B/B/B/A'+df_stat_snp$'B/B/B/B'+df_stat_snp$missing
head(df_stat_snp)

df_stat_snp$percent_NA = df_stat_snp$missing/df_stat_snp$Total_indiv * 100 
max(df_stat_snp$percent_NA) # 28.43137

pdf(file="missingnesslocus_205cc_111021.pdf",height=8,width=8)
ggplot()+geom_density(aes(df_stat_snp$percent_NA),fill="grey",adjust=0.5)+
  xlab("NA% per probe")+ylab("density")+
  xlim(0,30)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
geom_vline(xintercept=10,col="red",lwd=1)+
  annotate("text",x=6,y=0.10,label="64,649 probes",col="red",cex=4)+
  annotate("text",x=18,y=0.10,label="27,358 probes",col="red",cex=4)  
dev.off()

dim(subset(df_stat_snp,df_stat_snp$percent_NA > 10)) # 27,358 probes with more than 10% missing data
dim(subset(df_stat_snp,df_stat_snp$percent_NA <= 10)) # 64,649 probes with no more than 10% missing data

keep_snps = subset(rownames(df_stat_snp),df_stat_snp$status =="polymorphic" & df_stat_snp$percent_NA <= 10)
length(keep_snps) # 64,062 probes

# C ##########################################################################
# filtered matrix / missing data and monomorphic probes ######################### 

length(keep_genotypes) # 204 genotypes stay in the analysis
length(keep_snps) # 64,062 probes stay in the analysis

df_snp[1:5,1:5]

df_snp_filter = df_snp[which(df_snp$markername %in% keep_snps),]
df_snp_filter = df_snp_filter[,which(colnames(df_snp_filter)%in%c("markername",keep_genotypes))]

dim(df_snp_filter) # 64062   204
df_snp_filter[1:5,1:5]

# D ##########################################################################
# Compute minor allele frequency (MAF) and extract probes with MAF >= 5%
df_snp_filter[1:5,1:5]
str(df_snp_filter)

# Transform into matrix
matrix_snp_filter = t(df_snp_filter) # transpose to have individuals as rows and transform to matrix
matrix_snp_filter[1:5,1:5]
dim(matrix_snp_filter) # if transposed nb_ind nb_probes => 204 64062
rownames(matrix_snp_filter) # are individuals?
colnames(matrix_snp_filter)<- matrix_snp_filter[1,] # make SNP names as column names
matrix_snp_filter=matrix_snp_filter[-1,] # remove first row that is just SNP name

dim(matrix_snp_filter) # now the matrix is 204 64062
matrix_snp_filter[1:10,1:10] # only individuals with NA <= 20% and polymorphic SNP with less then 10% missing

## but for some reasons white space introduced during transposition / need to remove them
matrix_snp_filter= matrix(as.numeric(unlist(matrix_snp_filter)),nrow=nrow(matrix_snp_filter))
matrix_snp_filter[1:5,1:5]
dim(matrix_snp_filter)
rownames(matrix_snp_filter) <- colnames(df_snp_filter)[-1]
colnames(matrix_snp_filter) <- df_snp_filter$markername

L = length(keep_genotypes)

df_MAF = data.frame(SNP = colnames(matrix_snp_filter) ,
                    individuals=rep(L,length(colnames(matrix_snp_filter))),
                    NA_value =  colSums(is.na(matrix_snp_filter)),
                    total_allele = (rep(L,length(colnames(matrix_snp_filter)))-colSums(is.na(matrix_snp_filter)))*4,
                    B_allele = colSums (matrix_snp_filter, na.rm = T, dims = 1),
                    row.names=c(1:length(colnames(matrix_snp_filter))))

df_MAF$A_allele = df_MAF$total_allele - df_MAF$B_allele
df_MAF$A_freq = df_MAF$A_allele/df_MAF$total_allele *100
df_MAF$B_freq = df_MAF$B_allele/df_MAF$total_allele *100
head(df_MAF)         

hist(df_MAF$A_freq,main="Allele A frequency",xlab="")
abline(v=50,lwd=2,col="red")

hist(df_MAF$B_freq,main="Allele B frequency",xlab="")
abline(v=50,lwd=2,col="red")

df_MAF$MAF <- df_MAF$A_freq
df_MAF$MAF[df_MAF$MAF > 50 ] <- df_MAF$B_freq[df_MAF$MAF > 50 ]

df_MAF$minor_allele <- "A"
df_MAF$minor_allele[df_MAF$A_freq > 50 ] <- "B"

df_MAF = df_MAF[,c(1,2,3,4,6,5,7,8,9,10)]
head(df_MAF)
summary(df_MAF$MAF)

# plot
pdf(file="MAF_204cc_111021.pdf",height=9,width=9)
ggplot()+geom_density(aes(df_MAF$MAF),fill="grey",adjust=0.5)+
  xlab("Minor Allele Frequency (%)")+
  xlim(0,50)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_vline(xintercept=5,col="red",lwd=1)+
  annotate("text",x=1.2,y=0.10,label="17,438 probes",col="red",cex=4)+
  annotate("text",x=25,y=0.10,label="46,691 probes",col="red",cex=4)  
dev.off()
dim(subset(df_MAF,df_MAF$MAF >= 5)) # 46,691 probes with MAF >= 0.05
dim(subset(df_MAF,df_MAF$MAF <= 5)) # 17,438 probes with a MAF strictly lower than 0.05

keep_snps_MAF_5 <- subset(df_MAF$SNP,df_MAF$MAF >= 5)
length(keep_snps_MAF_5) # keep 46,691 probes at MAF5 for further analysis

df_snp_filter_MAF_5 <- df_snp_filter[which(df_snp_filter$markername %in% keep_snps_MAF_5) ,]
dim(df_snp_filter_MAF_5) # 46,691 probes

write.table(df_MAF,"matrix_fitPoly_270_scores_formated_only_passing_probes.204samplescc.MAF.txt",quote=F,row.names = F)



##### STEP 2 / Recode data to easily import into adegenet for further analysis
###########################################################################
###########################################################################
###########################################################################
df_snp_filter_MAF_5[1:5,1:5]
write.table(df_snp_filter_MAF_5,"matrix_fitPoly_270_scores_formated_only_passing_probes.204samplescc_withoutmono_MAF5pc.txt",quote=F,row.names = F) # save matrix to use in GWAS

# Make recode function to convert 0 to A/A/A/A, 1 to B/A/A/A, 2 to B/B/A/A, 3 to B/B/B/A and 4 to B/B/B/B
recode <-colwise(function(x){return(ifelse(x == 0, "A/A/A/A",
                                    ifelse(x == 1 , "B/A/A/A",
                                    ifelse(x==2 ,"B/B/A/A",
                                    ifelse(x==3,"B/B/B/A",
                                    ifelse(x==4,"B/B/B/B",x))))))})

rownames(df_snp_filter_MAF_5) <- df_snp_filter_MAF_5[,1]
df_snp_filter_MAF_5[1:5,1:5]
# Apply recode function to MAF 5 data
df_snp_filter_MAF_5_recode = recode(df_snp_filter_MAF_5) # take ~ a minute
rownames(df_snp_filter_MAF_5_recode)<- df_snp_filter_MAF_5[,1]


df_snp_filter_MAF_5_recode[1:5,1:5]

# Transpose to have individuals as rows and transform to matrix
matrix_snp_filter_MAF_5_recode = t(df_snp_filter_MAF_5_recode) 
matrix_snp_filter_MAF_5_recode[1:5,1:5]
dim(matrix_snp_filter_MAF_5_recode)

rownames(matrix_snp_filter_MAF_5_recode) # are individuals
colnames(matrix_snp_filter_MAF_5_recode)<- matrix_snp_filter_MAF_5_recode[1,] # make probe names as column names
matrix_snp_filter_MAF_5_recode=matrix_snp_filter_MAF_5_recode[-1,] # remove first row that is just probe names

dim(matrix_snp_filter_MAF_5_recode) # 204 individuals vs 46691
matrix_snp_filter_MAF_5_recode[1:5,1:5] # only individuals with NA <= 20% and polymorphic SNP with MAF >= 5% and no more than 10% missing data
colnames(matrix_snp_filter_MAF_5_recode)<- df_snp_filter_MAF_5[,1]
matrix_snp_filter_MAF_5_recode[1:5,1:5]


#########################
# TL - SCRIPT STOP HERE #
#########################



##### STEP 3 / Read genetic data with adegenet and heterozygosity test
###########################################################################
###########################################################################
###########################################################################

# A ##########################################################################
# Read data with adegenet
genind_CARO_MAF_5 <-df2genind(matrix_snp_filter_MAF_5_recode,ploidy=4,sep="/") # take some time to formate by adegenet
str(genind_CARO_MAF_5)
genind_CARO_MAF_5$tab[1:5,1:5] # genotype table, markers have been split by allele, numbers indicate quantity of a given allele
dim(genind_CARO_MAF_5$tab) #  93382 (*2 because two alleles per SNP)


# B ##########################################################################
# Perform heterozygosity test (see adegenet manual for more details)
stats_CARO_MAF_5 <- summary(genind_CARO_MAF_5) # very long is computing heterozygosity test
str(stats_CARO_MAF_5) #  (heterozygote status is sum ABBB+AABB+AAAB)

barplot(stats_CARO_MAF_5$Hobs[1:50], main="Heterozygosity-observed", # only 50 first probes
        ylab="Hobs",xaxt='n',ylim=c(0,1)) # 

barplot(stats_CARO_MAF_5$Hexp[1:50], main="Heterozygosity-expected", # only 50 first probes
        ylab="Hexp",xaxt='n',ylim=c(0,1)) # 

barplot(stats_CARO_MAF_5$Hexp[1:50]-stats_CARO_MAF_5$Hobs[1:50], main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs",xaxt='n') # excess of heterozygocy !!!! 

bartlett.test(list(stats_CARO_MAF_5$Hexp,stats_CARO_MAF_5$Hobs)) # significant difference between expected and obs

t.test(stats_CARO_MAF_5$Hexp,stats_CARO_MAF_5$Hobs,pair=T,var.equal=TRUE,alter="two.sided") # clearly  Hexp differet from Hobs
t.test(stats_CARO_MAF_5$Hexp,stats_CARO_MAF_5$Hobs,pair=T,var.equal=TRUE,alter="less")# clearly  Hexp < Hobs

# link to clonality and polyploidy of roses? excess of heterozygosity? this population in particular?

#### STEP 4/ perform PCOA (Principal Coordinate Analysis)
# It is equivalent to PCA without scaling data, only centering. 
# PCoA is cleaner than PCA when dealing with missing data, no need to impute missing data a as we do in PCA.
###########################################################################
###########################################################################

# A ##########################################################################
# compute euclidian distance 
genind_CARO_MAF_5[1:5,1:5]

Dist_CARO <- dist(tab(genind_CARO_MAF_5)) # generate euclidian distance matrix
Dist_CARO

Dist_CARO <- cailliez(Dist_CARO, print = FALSE, tol = 1e-07, cor.zero = TRUE) # make euclidian to avoid error


# B ##########################################################################
# pcoa analysis
pco_CARO<- dudi.pco(Dist_CARO, scannf=FALSE,nf=10) # nf=10, only 10 axis conserved for reporting coordinates of individuals

# eigenvalue plot (as %)
pco_CARO$eig[1:10] # 
lambda_values=pco_CARO$eig/sum(pco_CARO$eig)*100 # 
lambda_values[1:89]

pdf(file="204cc_eigen_pcoA_MAF5_121021.pdf",height=9,width=9)
#tiff("204cc_eigen_pcoa_MAF5.tiff",width = 480, height = 480)
barplot(lambda_values[1:50],main="PCoA eigenvalues\n(204cc with MAF>5%)", 
        col=heat.colors(50),ylim=c(0,25),ylab="Variance explained (%)",xlab="PCs")# eigen values
dev.off()

# save PCA indiv coordinates / bring other information to color plot
df.PCoA_MAF5 =as.data.frame(pco_CARO$li)
df.PCoA_MAF5$Accession_name <- row.names(df.PCoA_MAF5)
rownames(df.PCoA_MAF5) <- 1:length(row.names(df.PCoA_MAF5))
df.PCoA_MAF5=df.PCoA_MAF5[,c(11,1:10)]
head(df.PCoA_MAF5)

# import descriptors
df_descriptors <- read.table("df_description_CARO_tetraploids_genotypes.txt",header=T)
head(df_descriptors)

df.PCoA_MAF5 = merge(df_descriptors,df.PCoA_MAF5)

colnames(df.PCoA_MAF5)

# ploting PCOA with various coloring 


# color petal
df.PCoA_MAF5$flower_color_simplify_summary <-  factor(df.PCoA_MAF5$flower_color_simplify_summary,levels(df.PCoA_MAF5$flower_color_simplify_summary)[c(7,5,6,1,4,9,8,2,3)])

plot(df.PCoA_MAF5$A1,df.PCoA_MAF5$A2)
identified <- identify(df.PCoA_MAF5$A1,df.PCoA_MAF5$A2,labels=df.PCoA_MAF5$Accession_name,pos=T)
df.PCoA_MAF5$pos <- NA
df.PCoA_MAF5[identified$ind,]$pos <- identified$pos


tiff("axis_1_axis_2_MAF5_pcoa_experiment_petal_color.tiff",width = 2400, height = 1900,res=200)
ggplot(df.PCoA_MAF5, aes(x=A1, y=A2, group= flower_color_simplify_summary )) +
  geom_point(aes(shape=experiment,color=flower_color_simplify_summary,size=nbreP))+
  ylab(paste("Axis 2 - ",round(lambda_values[2],1),"%",sep="" ))+
  xlab(paste("Axis 1 - ",round(lambda_values[1],1),"%",sep="" ))+
  theme_bw()  +   scale_shape_manual(name="experiment",values=c(19,17))+   
  scale_color_manual(name="petal color",values=c("red","pink","purple","antiquewhite2","orange","yellow","black","green","grey","black"))+
  geom_text(data=subset(df.PCoA_MAF5,pos == 2),aes(label=Accession_name),vjust=1,size=3) +
  geom_text(data=subset(df.PCoA_MAF5,pos == 3),aes(label=Accession_name),vjust=1,size=3) +
  geom_text(data=subset(df.PCoA_MAF5,pos == 4),aes(label=Accession_name),vjust=1,size=3) 

dev.off() 

tiff("axis_1_axis_2_MAF5_pcoa_experiment_petal_color.tiff",width = 2400, height = 1900,res=200)
ggplot(df.PCoA_MAF5, aes(x=A1, y=A2, group= flower_color_simplify_summary )) +
  geom_point(aes(shape=type,color=flower_color_simplify_summary,size=nbreP))+
  ylab(paste("Axis 2 - ",round(lambda_values[2],1),"%",sep="" ))+
  xlab(paste("Axis 1 - ",round(lambda_values[1],1),"%",sep="" ))+
  theme_bw()  +   scale_shape_manual(name="group",values=c(19,17,2))+   
  scale_color_manual(name="petal color",values=c("red","pink","purple","antiquewhite2","orange","yellow","black","green","grey","black"))+
  geom_text(data=subset(df.PCoA_MAF5,pos == 2),aes(label=Accession_name),vjust=1,size=3) +
  geom_text(data=subset(df.PCoA_MAF5,pos == 3),aes(label=Accession_name),vjust=1,size=3) +
  geom_text(data=subset(df.PCoA_MAF5,pos == 4),aes(label=Accession_name),vjust=1,size=3) 

dev.off()

tiff("axis_1_axis_2_MAF5_pcoa_experiment_timeline.tiff",width = 2400, height = 1900,res=200)
ggplot(df.PCoA_MAF5, aes(x=A1, y=A2, group= Classification_JCC2_Simplify )) +
  geom_point(aes(color=Classification_JCC2_Simplify,size=nbreP))+
  ylab(paste("Axis 2 - ",round(lambda_values[2],1),"%",sep="" ))+
  xlab(paste("Axis 1 - ",round(lambda_values[1],1),"%",sep="" ))+
  theme_bw()  +      
  geom_text(data=subset(df.PCoA_MAF5,pos == 2),aes(label=Accession_name),vjust=1,size=3) +
  geom_text(data=subset(df.PCoA_MAF5,pos == 3),aes(label=Accession_name),vjust=1,size=3) +
  geom_text(data=subset(df.PCoA_MAF5,pos == 4),aes(label=Accession_name),vjust=1,size=3) 

dev.off()


#### STEP 5/ perform DAPC 
# See adegenet manual and Jombart et al. 2010 in BMC Genomics
# https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
###########################################################################
###########################################################################

# A ##########################################################################
# find groups using find clusters
grp_MAF_5 <- find.clusters(genind_CARO_MAF_5, max.n.clust=50,n.pca =999) # max 50 clusters # by default center and scale
str(grp_MAF_5) # n.pca = 999 allow to retained all PC
# the BIC graph suggests 5 groups to be a good number, enter 5 !

grp_MAF_5$Kstat
grp_MAF_5$stat # retained 9 groups
grp_MAF_5$grp
grp_MAF_5$size

# B ##########################################################################
# run DAPC first attend with large number of PC, later we will reduce number using a score approach
dapc5 <- dapc(genind_CARO_MAF_5, grp_MAF_5$grp,n.pca=180, n.da=999) # retained 60 PC (n.pca=60) & all 5 discriminant functions (n.da=999)
str(dapc5)
summary(dapc5)

# plots
scatter(dapc5,posi.da="topleft",posi.pca="topright",pch=17:22,cellipse = 1.5,scree.pca=TRUE)
scatter(dapc5,posi.da="topright",posi.pca=NULL,pch=17:22,cellipse = 1.5,scree.pca=TRUE)
# compute a score to select optimal number of PC (see Jombart manual)
a_score <- a.score(dapc5)
a_score

a_score_opt <- optim.a.score(dapc5) # 23 PC is optimal !


# C ##########################################################################
# re-run DAPC with large 23 PC only
dapc5 <- dapc(genind_CARO_MAF_5, grp_MAF_5$grp,n.pca=23, n.da=999) # retained 23 PC (n.pca=23) & all 5 discriminant functions (n.da=999)
str(dapc5)
summary(dapc5)

# plots
pdf(file="204cc_MAF5pc_111021.pdf",height=11,width=11)
scatter(dapc5,posi.da="topleft",posi.pca="bottomleft",pch=20,cellipse = 1.5,scree.pca=TRUE)
dev.off()

# manual plots to link with descriptors
df_DAPC_MAF5=dapc5$ind.coord
df_DAPC_MAF5=as.data.frame(df_DAPC_MAF5)
df_DAPC_MAF5$Accession_name <- as.factor(rownames(df_DAPC_MAF5))
df_DAPC_MAF5=df_DAPC_MAF5[,c(9,1:8)]
rownames(df_DAPC_MAF5) <- 1:length(rownames(df_DAPC_MAF5))
head(df_DAPC_MAF5)
df_DAPC_MAF5$group <- grp_MAF_5$grp

df_DAPC_MAF5 <- merge(df_descriptors , df_DAPC_MAF5)


ggplot(df_DAPC_MAF5, aes(x=LD1, y=LD2)) +
  geom_point(size=2)+   theme_bw() 

df_DAPC_MAF5$flower_color_simplify_summary <-  factor(df_DAPC_MAF5$flower_color_simplify_summary,levels(df_DAPC_MAF5$flower_color_simplify_summary)[c(7,5,6,1,4,9,8,2,3)])


tiff("DAPC_MAF_5_16PC_5grp.tiff",res=200,width=2100,height=1600)
ggplot(df_DAPC_MAF5, aes(x=LD1, y=LD2, group=as.factor(group))) +
  geom_point(aes(shape=type,color=flower_color_simplify_summary,size=nbreP))+ 
  stat_ellipse()+
  theme_bw() +   
  scale_shape_manual("type",values=c(19,17,2))+
  scale_color_manual(name="petal color",values=c("red","pink","purple","antiquewhite2","orange","yellow","black","green","grey"))

dev.off()



#### STEP 6/ export DAPC and PCoA results
###########################################################################
###########################################################################
df_structure_MAF5 = data.frame(Accession_name = df.PCoA_MAF5$Accession_name ,
                               PCo1= df.PCoA_MAF5$A1,
                               PCo2=df.PCoA_MAF5$A2,
                               PCo3=df.PCoA_MAF5$A3,
                               PCo4=df.PCoA_MAF5$A4,
                               PCo5=df.PCoA_MAF5$A5,
                               DA1=dapc5$ind.coord[,1],
                               DA2=dapc5$ind.coord[,2],
                               DA3=dapc5$ind.coord[,3],
                               DA4=dapc5$ind.coord[,4],
                               DA5=dapc5$ind.coord[,5],
                               kmeans=grp_MAF_5$grp)

write.table(df_structure_MAF5,"204cc_structure_PCoA_DAPC_MAF5.txt",quote=F,row.names = F)


M<-cor(df_structure_MAF5[,-c(1,12)])
corrplot(M, method="circle")

pdf(file="204cc_correlation_PCoA_DAPC_axes_121021.pdf",height=9,width=9)
corrplot(M, method="circle")
dev.off()

