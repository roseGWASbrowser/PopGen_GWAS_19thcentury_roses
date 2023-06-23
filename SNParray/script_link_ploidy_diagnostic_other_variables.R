## 15/12/2020 Elise - linking ploidy diagnostic to design variables
## 188 active tetraploids

#### I created file "df_stat_genotypes_ploidy_diag.txt" wich of "df_stat_genotypes.txt" with design variables and diagnostic ploidy (see SUMMARY_diagnostic_ploidy_Elise.xlsx).

##### load packages
library(ggplot2)

##### load data
df_indiv_ploid_diag = read.table("df_stat_genotypes_ploidy_diag.txt",header=T)
head(df_indiv_ploid_diag)
df_indiv_ploid_diag$percent_NA = df_indiv_ploid_diag$no_call/df_indiv_ploid_diag$nprob*100

##### density plots percent missing genotyping calls
tiff("proportion_missing_per_sample_diag_ploidy.tiff", width = 1500, height = 1200,res=200, type="cairo")
ggplot(df_indiv_ploid_diag, aes(x=percent_NA,fill=diagnostic_ploidy)) + 
  geom_density(alpha=0.6)+theme_bw()+xlab("% NA calls per sample")
dev.off()

##### link to origin of accessions v2
table_JCC2 = table(df_indiv_ploid_diag$diagnostic_ploidy,df_indiv_ploid_diag$Classification_JCC2_Simplify)
table_JCC2=as.data.frame(table_JCC2)
table_JCC2$Var2 <- factor(table_JCC2$Var2,levels(table_JCC2$Var2)[c(3,2,1)])

tiff("diag_ploidy_vs_cultiv_wild.tiff", width = 1500, height = 1200,res=200, type="cairo")
ggplot(table_JCC2, aes(Var2, Var1,  size = Freq)) +
  geom_point(col="cyan") +
  theme(legend.position="bottom") +
  scale_size(range = c(1,15)) +
  guides(color= guide_legend(), size=guide_legend())+theme_bw()+
  xlab("")+ylab("ploidy")+ geom_text(aes(label=Freq))+
theme(axis.text.x = element_text(color="black", size=14),
      axis.text.y = element_text(color="black", 
                                 size=14))+
  theme(axis.text.x = element_text(angle = 90)) + guides(color = FALSE, size = FALSE)
dev.off()

##### link to origin of accessions v1
table_JCC1 = table(df_indiv_ploid_diag$diagnostic_ploidy,df_indiv_ploid_diag$Classification_JCC_Simplify)
table_JCC1=as.data.frame(table_JCC1)

ggplot(table_JCC1, aes(Var2, Var1,  size = Freq)) +
  geom_point(col="cyan") +
  theme(legend.position="bottom") +
  scale_size(range = c(1,15)) +
  guides(color= guide_legend(), size=guide_legend())+theme_bw()+
  xlab("")+ylab("ploidy")+ geom_text(aes(label=Freq))+
theme(axis.text.x = element_text(color="black", size=14),
          axis.text.y = element_text(color="black", 
                                     size=14)) + guides(color = FALSE, size = FALSE)
  
##### link to group genetic Liorzou et al.
table_liorzou = table(df_indiv_ploid_diag$diagnostic_ploidy,df_indiv_ploid_diag$Genetic_group_Liorzou)
table_liorzou=as.data.frame(table_liorzou)

tiff("diag_ploidy_vs_liorzou.tiff", width = 2500, height = 1200,res=200, type="cairo")
ggplot(table_liorzou, aes(Var2, Var1,  size = Freq)) +
  geom_point(col="cyan") +
  theme(legend.position="bottom") +
  scale_size(range = c(1,15)) +
  guides(color= guide_legend(), size=guide_legend())+theme_bw()+
  xlab("")+ylab("ploidy")+ geom_text(aes(label=Freq))+
  theme(axis.text.x = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", 
                                   size=14)) + guides(color = FALSE, size = FALSE)
dev.off()

##### link to horti groups
table_horti = table(df_indiv_ploid_diag$diagnostic_ploidy,df_indiv_ploid_diag$Classification_Horticole_Simplify)
table_horti=as.data.frame(table_horti)

tiff("diag_ploidy_vs_horti.tiff", width = 2500, height = 1200,res=200, type="cairo")
ggplot(table_horti, aes(Var2, Var1,  size = Freq)) +
  geom_point(col="cyan") +
  theme(legend.position="bottom") +
  scale_size(range = c(1,15)) +
  guides(color= guide_legend(), size=guide_legend())+theme_bw()+
  xlab("")+ylab("ploidy")+ geom_text(aes(label=Freq))+
  theme(axis.text.x = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", 
                                   size=14)) + guides(color = FALSE, size = FALSE)+
  theme(axis.text.x = element_text(angle = 90)) 
dev.off()

##### link horti groups to Liorzou groups
table_test = table(df_indiv_ploid_diag$Classification_Horticole_Simplify,df_indiv_ploid_diag$Genetic_group_Liorzou)
table_test=as.data.frame(table_test)

tiff("group_genetic_vs_horti.tiff", width = 2500, height = 1200,res=200, type="cairo")
ggplot(table_test, aes(Var2, Var1,  size = Freq)) +
  geom_point(col="cyan") +
  theme(legend.position="bottom") +
  scale_size(range = c(1,15)) +
  guides(color= guide_legend(), size=guide_legend())+theme_bw()+
  xlab("genetic group")+ylab("horticol group")+ geom_text(aes(label=Freq))+
  theme(axis.text.x = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", 
                                   size=14)) + guides(color = FALSE, size = FALSE)

dev.off()
