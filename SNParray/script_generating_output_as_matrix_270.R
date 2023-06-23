## 22/10/2020 Elise - producing genotyping matrix

### Load packages
#############################################
library(data.table)
library(tidyr)

setwd("/bigvol/benoit/Analyses/Temp_Tibo/Roses/output_fitpoly_270samples_CARO/")

#########################################################################
df = fread("fitPoly_270_scores_formated_only_passing.sed.dat")
passing_probes = fread("list_scored_probes.ps")

#df_pass = df[which(df$MarkerName %in% passing_probes$probeset_id),]
df_pass = df[which(df$markername %in% passing_probes$probeset_id),]
dim(df_pass ) ; dim(df)

matrix_df_pass = spread(df_pass[,c(2,3,15)], key = sample, value = geno)


write.table(matrix_df_pass,"matrix_fitPoly_270_scores_formated_only_passing_probes.txt",quote=F,row.names=F)
