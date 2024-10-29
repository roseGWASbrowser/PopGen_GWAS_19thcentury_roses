# TL - 130121 
# Modified 28-29/10/24 to more explicitly show how low coverage impacts the accuracy of the results
library("ggplot2")
library("grid")
library("gridExtra")
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/sequencing/first_round_calling")
coverage=read.table(file="rose32.jointvcf.clean.PASSonly.SNPonly.vcf.2pcsites.vcf.depth.hzposonly.newID2")

ggplot(coverage_to_use, aes(x=V8/V5)) + geom_density()  + 
  ylab("Density")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


# filter per indivdual all sites with cov higher than 20 & compute allelic balance for the minor allele
coverage_to_use=subset(coverage,(coverage$V5>=20))
coverage_to_use$V9<- ifelse((coverage_to_use$V7/coverage_to_use$V5)<=0.5,(coverage_to_use$V7/coverage_to_use$V5),(coverage_to_use$V8/coverage_to_use$V5))

# Compute mode fonction -> MOST FREQUENTLY OBSERVED (EXACT) ALLELE FREQUENCY AMONG ALL VALUES
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}
# loop to compute the most frequent observed value "mode" per individual
mode_val=NULL
for (ind in sort(unique(coverage_to_use$V1))){
  covind=subset(coverage_to_use,coverage_to_use$V1==ind)
  mode_val=rbind(mode_val,cbind(ind,nrow(covind),Mode(covind$V9,na.rm=TRUE)))
}
mode_val<-as.data.frame.matrix(mode_val)
colnames(mode_val)<- c("indiv","number_SNPs_more_20","mode_allelic_balance")

# Most frequent allele balance values per individual:
                          #indiv  number_SNPs_more_20 mode_allelic_balance
#1    01_Duchesse_de_Montebello               22652                 0.25
#2              02_Blush_damask               15115                 0.25
#3      03_Comtesse_de_murinais               11341                 0.25
#4  04_Rosa_gallica_officinalis               85975                 0.25
#5            05_General_kleber                5841                 0.25
#6          06_Rosa_x_damascena               59580                 0.25
#7              07_Fen_zhang_lu               26478                  0.5
#8             08_Jin_ou_fan_lu                6761                 0.25
#9                 09_Old_blush               78655                  0.5
#10 10_Hume_s_blush_tea_scented               49581                  0.5
#11 11_Rosa_chinensis_sanguinea               38358                  0.5
#12 12_Rosa_chinensis_mutabilis               49688                  0.5
#13    13_Rosa_odorata_gigantea               32758                  0.5
#14 14_Rosa_chinensis_spontanea               26617                  0.5
#15      15_Marguerite_de_roman                4798                 0.25
#16 16_Triomphe_de_l_exposition                7789                 0.25
#17           17_Charles_Lawson                5310                 0.25
#18    18_Jacques_Cartier_blanc                5104                 0.25
#19            19_Yellow_Island               73462                 0.25
#20           20_Victor_Verdier                1094                  0.2
#21            21_Lady_Waterlow                3265                  0.2
#22                 22_La_tosca                6320                 0.25
#23              23_La_favorite                4041                 0.25
#24             24_Irene_Bonnet               12767                 0.25
#25                25_La_France               95247    0.333333333333333
#26            26_Rosa_moschata               30814                  0.5
#27  27_Rosa_xanthina_spontanea               18381                  0.5
#28           28_Rosa_pendulina               51475                 0.25
#29         29_Rosa_wichuraiana               30208                  0.5
#30             30_Rosa_majalis               22308                  0.5
#31            31_rosa_arvensis               23897                  0.5
#32         32_Rosa_minutifolia               19486                  0.5


#### Generate plots

#Color code:
cols_restricted<- c("01_Duchesse_de_Montebello" = "#000080",
                    "02_Blush_damask" = "#000080", 
                    "03_Comtesse_de_murinais" = "#000080", 
                    "04_Rosa_gallica_officinalis" = "#000080", 
                    "05_General_kleber" = "#000080", 
                    "06_Rosa_x_damascena" = "#000080", 
                    "07_Fen_zhang_lu" = "#EEB422",
                    "08_Jin_ou_fan_lu" = "#EEB422", 
                    "09_Old_blush" = "#EEB422", 
                    "10_Hume_s_blush_tea_scented" = "#EEB422", 
                    "11_Rosa_chinensis_sanguinea" = "#EEB422", 
                    "12_Rosa_chinensis_mutabilis" = "#EEB422", 
                    "13_Rosa_odorata_gigantea" = "#EEB422",
                    "14_Rosa_chinensis_spontanea" = "#EEB422",
                    "15_Marguerite_de_roman" = "#228B22", 
                    "16_Triomphe_de_l_exposition" = "#228B22", 
                    "17_Charles_Lawson" = "#228B22", 
                    "18_Jacques_Cartier_blanc" = "#228B22", 
                    "19_Yellow_Island" = "#889F22", 
                    "20_Victor_Verdier" = "#889F22", 
                    "21_Lady_Waterlow" = "#889F22", 
                    "22_La_tosca" = "#889F22", 
                    "23_La_favorite" = "#889F22",
                    "24_Irene_Bonnet" = "#889F22",
                    "25_La_France" = "#889F22",
                    "26_Rosa_moschata" = "#FF0000",
                    "27_Rosa_xanthina_spontanea" = "#FF0000",
                    "28_Rosa_pendulina" = "#FF0000",
                    "29_Rosa_wichuraiana" = "#FF0000",
                    "30_Rosa_majalis" = "#FF0000",
                    "31_rosa_arvensis" = "#FF0000",
                    "32_Rosa_minutifolia" = "#FF0000")
  

p1<- ggplot(coverage_to_use, aes(x=V9,color=V1)) + geom_density(adjust=1.5)  + 
  scale_colour_manual(name="Individuals", values = cols_restricted)+  # "#ffe458","#cb745b"))+
  xlab("Allelic balance at heterozygous calls (minor allele)")+
  ylab("Density")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
print(p1)



pdf(file = "allelicbalance_perind_281024_withcolors_v2.pdf",width = 10,height = 8)
for (ind in sort(unique(coverage_to_use$V1))){
  print(ind) 
  p1<- ggplot(coverage_to_use, aes(x=V9,color=V1,linetype=V1,size=V1)) + geom_density(adjust=1.5)  + 
    scale_colour_manual(name="Individuals", values = cols_restricted)+  # "#ffe458","#cb745b"))+
    scale_linetype_manual(name="Individuals", values = ifelse(sort(unique(coverage_to_use$V1))==ind,"solid","longdash"))+  # "#ffe458","#cb745b"))+
    scale_size_manual(name="Individuals", values=ifelse(sort(unique(coverage_to_use$V1))==ind,1.5,0.7))+
    xlab("Allelic balance at heterozygous calls (minor allele)")+
    ylab("Density")+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  print(p1)
  #Sys.sleep(2)
}
dev.off()

### final Extended Fig. S2 with histograms
pdf(file = "allelicbalance_perind_281024_withcolors_withhistograms.pdf",width = 20,height = 9)
for (ind in sort(unique(coverage_to_use$V1))){
  print(ind) 
  p1<- ggplot(coverage_to_use, aes(x=V9,color=V1,linetype=V1,size=V1)) + 
    geom_density(adjust=1.5)  + 
    scale_colour_manual(name="Individuals", values = cols_restricted)+  # "#ffe458","#cb745b"))+
    scale_linetype_manual(name="Individuals", values = ifelse(sort(unique(coverage_to_use$V1))==ind,"solid","longdash"))+  # "#ffe458","#cb745b"))+
    scale_size_manual(name="Individuals", values=ifelse(sort(unique(coverage_to_use$V1))==ind,1.5,0.7))+
    xlab("Allelic balance at heterozygous calls (minor allele)")+
    ylab("Density")+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  #print(p1)
  
  coverage_to_use_ind=subset(coverage_to_use,coverage_to_use$V1==ind)
  coverage_to_use_ind_mode=subset(mode_val,mode_val$indiv==ind)
  p2<- ggplot(coverage_to_use_ind, aes(x=V9)) + 
    geom_histogram(color="darkgrey", fill="lightgrey",bins = 100) +
    ## add line to highlight the most frequently observed value among all allelic balance values
    geom_vline(xintercept = as.numeric(as.character(coverage_to_use_ind_mode$mode_allelic_balance)),color="darkred",lwd=2,lty=2)+ 
    xlab("Allelic balance at heterozygous calls (minor allele)")+
    ylab("#SNPs")+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+ 
    annotate("text", x = as.numeric(as.character(coverage_to_use_ind_mode$mode_allelic_balance))-0.09, y = 100, label = paste0(ind,":\nMost freq. obs. allelic balance value=",round(as.numeric(as.character(coverage_to_use_ind_mode$mode_allelic_balance)),4),"\nn=",coverage_to_use_ind_mode$number_SNPs_more_20," SNPs (DP>=20)"),color="darkred")
  #print(p2)
  grid.arrange(p1, p2, nrow = 1)
  
  #Sys.sleep(2)
}
dev.off()



pdf(file = "allelicbalance_perind_281024.pdf",width = 10,height = 8)
for (ind in sort(unique(coverage_to_use$V1))){
 print(ind) 
 p1<- ggplot(coverage_to_use, aes(x=V9,color=V1)) + geom_density(adjust=1.5)  + 
   scale_colour_manual(values = ifelse(sort(unique(coverage_to_use$V1))==ind,"red2","grey70"))+  # "#ffe458","#cb745b"))+
   xlab("Allelic balance at heterozygous calls (minor allele)")+
    ylab("Density")+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
 print(p1)
 #Sys.sleep(2)
}
dev.off()



pdf(file = "allelicbalance_perind_interchromosomicvar_281024.pdf",width = 10,height = 8)
for (ind in sort(unique(coverage_to_use$V1))){
  coverage_to_use_ind=subset(coverage,(coverage$V1==ind)&(coverage$V5>=20))
  coverage_to_use_ind$V9<- ifelse((coverage_to_use_ind$V7/coverage_to_use_ind$V5)<=0.5,(coverage_to_use_ind$V7/coverage_to_use_ind$V5),(coverage_to_use_ind$V8/coverage_to_use_ind$V5))
  print(ind) 
  p1<- ggplot(coverage_to_use_ind, aes(x=V9,color=V2)) + geom_density(adjust=1.5)  + 
    xlab("Allelic balance at heterozygous calls (minor allele)")+
    ylab("Density")+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
  print(p1+ggtitle(ind) )
  #Sys.sleep(2)
}
dev.off()
