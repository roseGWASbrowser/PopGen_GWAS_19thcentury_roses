# TL - 110921
library(ggplot2)
library("ggrepel")
library("grid")
library("gridExtra")
setwd("/home/thibaultleroy/Rose/documents_transfert_autreordi/PCA_071021/")

data=read.table("matrix_fitPoly_270_scores_formated_only_passing_probes.205samplescc",h=T)

# recode for adegenet
library(plyr)
recode <-colwise(function(x){return(ifelse(x == 0, "A/A/A/A",
                                    ifelse(x == 1 , "B/A/A/A",
                                    ifelse(x==2 ,"B/B/A/A",
                                    ifelse(x==3,"B/B/B/A",
                                    ifelse(x==4,"B/B/B/B",x))))))})

datarec = recode(data) # take ~ a minute
colnames(datarec)

# transpose
datarectransp = t(datarec) 
datarectransp =datarectransp[-1,]
datarectransp[1:5,1:5]

colnames(datarectransp)


# perform PCA
library(adegenet)

dataadegenet <-df2genind(datarectransp,ploidy=4,sep="/") # take some time to formate by adegenet
str(dataadegenet)


x.cows <- tab(dataadegenet, freq=TRUE, NA.method="mean")
pca.cows <- dudi.pca(x.cows, center=TRUE, scale=FALSE)

propvarexp=rbind(pca.cows$eig/sum(pca.cows$eig),cumsum(pca.cows$eig)/sum(pca.cows$eig))

write.table(pca.cows$li,file="res_PCA_071021.txt",sep="\t")

pdf(file="PCA_130921.pdf",width=18,height=18)
par(mfrow=c(2,2))
plot(pca.cows$li[,1],pca.cows$li[,2], xlab=paste("PC 1 - ",round(propvarexp[1,1]*100,2),"%"),ylab=paste("PC 2 - ",round(propvarexp[1,2]*100,2),"%"),pch=20,cex=2,
     col= ifelse(rownames(pca.cows$li)=="duchesse_de_montebello_H7","navyblue",ifelse(rownames(pca.cows$li)=="kean_D12","navyblue",
       ifelse(rownames(pca.cows$li)=="berengere_B9","navyblue",ifelse(rownames(pca.cows$li)=="blush_damask_H1","navyblue",
       ifelse(rownames(pca.cows$li)=="hume_s_blush_tea_scented_china_C7","orange",ifelse(rownames(pca.cows$li)=="old_blush_G12","orange",ifelse(rownames(pca.cows$li)=="rosa_chinensis_mutabilis_G3","orange",
       ifelse(rownames(pca.cows$li)=="mme_lierval_E9","purple",ifelse(rownames(pca.cows$li)=="victor_verdier_chene_H6","purple",
       ifelse(rownames(pca.cows$li)=="marguerite_de_roman_E12","purple",ifelse(rownames(pca.cows$li)=="triomphe_de_l_exposition_A8","purple",ifelse(rownames(pca.cows$li)=="alfred_colomb_C3","purple",                                                                                                                                      
       ifelse(rownames(pca.cows$li)=="la_france_D7","red",ifelse(rownames(pca.cows$li)=="la_favorite_lyon_A10","red","darkgrey")))))))))))))))
text(pca.cows$li[,1],pca.cows$li[,2]+1, labels=rownames(pca.cows$li),cex=0.6,pos=4,
     col= ifelse(rownames(pca.cows$li)=="duchesse_de_montebello_H7","navyblue",ifelse(rownames(pca.cows$li)=="kean_D12","navyblue",
          ifelse(rownames(pca.cows$li)=="berengere_B9","navyblue",ifelse(rownames(pca.cows$li)=="blush_damask_H1","navyblue",
          ifelse(rownames(pca.cows$li)=="hume_s_blush_tea_scented_china_C7","orange",ifelse(rownames(pca.cows$li)=="old_blush_G12","orange",ifelse(rownames(pca.cows$li)=="rosa_chinensis_mutabilis_G3","orange",
          ifelse(rownames(pca.cows$li)=="mme_lierval_E9","purple",ifelse(rownames(pca.cows$li)=="victor_verdier_chene_H6","purple",
          ifelse(rownames(pca.cows$li)=="marguerite_de_roman_E12","purple",ifelse(rownames(pca.cows$li)=="triomphe_de_l_exposition_A8","purple",ifelse(rownames(pca.cows$li)=="alfred_colomb_C3","purple",                                                                                                                                      
          ifelse(rownames(pca.cows$li)=="la_france_D7","red",ifelse(rownames(pca.cows$li)=="la_favorite_lyon_A10","red","darkgrey")))))))))))))))


plot(pca.cows$li[,1],pca.cows$li[,3], xlab=paste("PC 1 - ",round(propvarexp[1,1]*100,2),"%"),ylab=paste("PC 3 - ",round(propvarexp[1,3]*100,2),"%"),pch=20,cex=2,
     col= ifelse(rownames(pca.cows$li)=="duchesse_de_montebello_H7","navyblue",ifelse(rownames(pca.cows$li)=="kean_D12","navyblue",
     ifelse(rownames(pca.cows$li)=="berengere_B9","navyblue",ifelse(rownames(pca.cows$li)=="blush_damask_H1","navyblue",
     ifelse(rownames(pca.cows$li)=="hume_s_blush_tea_scented_china_C7","orange",ifelse(rownames(pca.cows$li)=="old_blush_G12","orange",ifelse(rownames(pca.cows$li)=="rosa_chinensis_mutabilis_G3","orange",
     ifelse(rownames(pca.cows$li)=="mme_lierval_E9","purple",ifelse(rownames(pca.cows$li)=="victor_verdier_chene_H6","purple",
     ifelse(rownames(pca.cows$li)=="marguerite_de_roman_E12","purple",ifelse(rownames(pca.cows$li)=="triomphe_de_l_exposition_A8","purple",ifelse(rownames(pca.cows$li)=="alfred_colomb_C3","purple",                                                                                                                                      
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ifelse(rownames(pca.cows$li)=="la_france_D7","red",ifelse(rownames(pca.cows$li)=="la_favorite_lyon_A10","red","darkgrey")))))))))))))))
text(pca.cows$li[,1],pca.cows$li[,3]+1, labels=rownames(pca.cows$li),cex=0.6,pos=4,
     col= ifelse(rownames(pca.cows$li)=="duchesse_de_montebello_H7","navyblue",ifelse(rownames(pca.cows$li)=="kean_D12","navyblue",
    ifelse(rownames(pca.cows$li)=="berengere_B9","navyblue",ifelse(rownames(pca.cows$li)=="blush_damask_H1","navyblue",
    ifelse(rownames(pca.cows$li)=="hume_s_blush_tea_scented_china_C7","orange",ifelse(rownames(pca.cows$li)=="old_blush_G12","orange",ifelse(rownames(pca.cows$li)=="rosa_chinensis_mutabilis_G3","orange",
    ifelse(rownames(pca.cows$li)=="mme_lierval_E9","purple",ifelse(rownames(pca.cows$li)=="victor_verdier_chene_H6","purple",
    ifelse(rownames(pca.cows$li)=="marguerite_de_roman_E12","purple",ifelse(rownames(pca.cows$li)=="triomphe_de_l_exposition_A8","purple",ifelse(rownames(pca.cows$li)=="alfred_colomb_C3","purple",                                                                                                                                      
    ifelse(rownames(pca.cows$li)=="la_france_D7","red",ifelse(rownames(pca.cows$li)=="la_favorite_lyon_A10","red","darkgrey")))))))))))))))


plot(pca.cows$li[,1],pca.cows$li[,4], xlab=paste("PC 1 - ",round(propvarexp[1,1]*100,2),"%"),ylab=paste("PC 4 - ",round(propvarexp[1,4]*100,2),"%"),pch=20,cex=2,
     col= ifelse(rownames(pca.cows$li)=="duchesse_de_montebello_H7","navyblue",ifelse(rownames(pca.cows$li)=="kean_D12","navyblue",
     ifelse(rownames(pca.cows$li)=="berengere_B9","navyblue",ifelse(rownames(pca.cows$li)=="blush_damask_H1","navyblue",
     ifelse(rownames(pca.cows$li)=="hume_s_blush_tea_scented_china_C7","orange",ifelse(rownames(pca.cows$li)=="old_blush_G12","orange",ifelse(rownames(pca.cows$li)=="rosa_chinensis_mutabilis_G3","orange",
     ifelse(rownames(pca.cows$li)=="mme_lierval_E9","purple",ifelse(rownames(pca.cows$li)=="victor_verdier_chene_H6","purple",
     ifelse(rownames(pca.cows$li)=="marguerite_de_roman_E12","purple",ifelse(rownames(pca.cows$li)=="triomphe_de_l_exposition_A8","purple",ifelse(rownames(pca.cows$li)=="alfred_colomb_C3","purple",                                                                                                                                      
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ifelse(rownames(pca.cows$li)=="la_france_D7","red",ifelse(rownames(pca.cows$li)=="la_favorite_lyon_A10","red","darkgrey")))))))))))))))
text(pca.cows$li[,1],pca.cows$li[,4]+1, labels=rownames(pca.cows$li),cex=0.6,pos=4,
     col= ifelse(rownames(pca.cows$li)=="duchesse_de_montebello_H7","navyblue",ifelse(rownames(pca.cows$li)=="kean_D12","navyblue",
     ifelse(rownames(pca.cows$li)=="berengere_B9","navyblue",ifelse(rownames(pca.cows$li)=="blush_damask_H1","navyblue",
     ifelse(rownames(pca.cows$li)=="hume_s_blush_tea_scented_china_C7","orange",ifelse(rownames(pca.cows$li)=="old_blush_G12","orange",ifelse(rownames(pca.cows$li)=="rosa_chinensis_mutabilis_G3","orange",
     ifelse(rownames(pca.cows$li)=="mme_lierval_E9","purple",ifelse(rownames(pca.cows$li)=="victor_verdier_chene_H6","purple",
     ifelse(rownames(pca.cows$li)=="marguerite_de_roman_E12","purple",ifelse(rownames(pca.cows$li)=="triomphe_de_l_exposition_A8","purple",ifelse(rownames(pca.cows$li)=="alfred_colomb_C3","purple",                                                                                                                                      
     ifelse(rownames(pca.cows$li)=="la_france_D7","red",ifelse(rownames(pca.cows$li)=="la_favorite_lyon_A10","red","darkgrey")))))))))))))))




plot(pca.cows$li[,2],pca.cows$li[,3], xlab=paste("PC 2 - ",round(propvarexp[1,2]*100,2),"%"),ylab=paste("PC 3 - ",round(propvarexp[1,3]*100,2),"%"),pch=20,cex=2,
     col= ifelse(rownames(pca.cows$li)=="duchesse_de_montebello_H7","navyblue",ifelse(rownames(pca.cows$li)=="kean_D12","navyblue",
     ifelse(rownames(pca.cows$li)=="berengere_B9","navyblue",ifelse(rownames(pca.cows$li)=="blush_damask_H1","navyblue",
     ifelse(rownames(pca.cows$li)=="hume_s_blush_tea_scented_china_C7","orange",ifelse(rownames(pca.cows$li)=="old_blush_G12","orange",ifelse(rownames(pca.cows$li)=="rosa_chinensis_mutabilis_G3","orange",
     ifelse(rownames(pca.cows$li)=="mme_lierval_E9","purple",ifelse(rownames(pca.cows$li)=="victor_verdier_chene_H6","purple",
     ifelse(rownames(pca.cows$li)=="marguerite_de_roman_E12","purple",ifelse(rownames(pca.cows$li)=="triomphe_de_l_exposition_A8","purple",ifelse(rownames(pca.cows$li)=="alfred_colomb_C3","purple",                                                                                                                                      
     ifelse(rownames(pca.cows$li)=="la_france_D7","red",ifelse(rownames(pca.cows$li)=="la_favorite_lyon_A10","red","darkgrey")))))))))))))))

text(pca.cows$li[,2],pca.cows$li[,3]+1, labels=rownames(pca.cows$li),cex=0.6,pos=4,
     col= ifelse(rownames(pca.cows$li)=="duchesse_de_montebello_H7","navyblue",ifelse(rownames(pca.cows$li)=="kean_D12","navyblue",
     ifelse(rownames(pca.cows$li)=="berengere_B9","navyblue",ifelse(rownames(pca.cows$li)=="blush_damask_H1","navyblue",
     ifelse(rownames(pca.cows$li)=="hume_s_blush_tea_scented_china_C7","orange",ifelse(rownames(pca.cows$li)=="old_blush_G12","orange",ifelse(rownames(pca.cows$li)=="rosa_chinensis_mutabilis_G3","orange",
     ifelse(rownames(pca.cows$li)=="mme_lierval_E9","purple",ifelse(rownames(pca.cows$li)=="victor_verdier_chene_H6","purple",
     ifelse(rownames(pca.cows$li)=="marguerite_de_roman_E12","purple",ifelse(rownames(pca.cows$li)=="triomphe_de_l_exposition_A8","purple",ifelse(rownames(pca.cows$li)=="alfred_colomb_C3","purple",                                                                                                                                      
     ifelse(rownames(pca.cows$li)=="la_france_D7","red",ifelse(rownames(pca.cows$li)=="la_favorite_lyon_A10","red","darkgrey")))))))))))))))

dev.off()

#write.table(pca.cows$li,file="res_PCA_270genotypes.txt",sep="\t",quote=FALSE)

pca_faststr=read.table("res_PCA_071021_withinfo.txt",h=T)


plot(pca_faststr$Axis1,pca_faststr$Axis2, xlab=paste("PC 1 - ",round(propvarexp[1,1]*100,2),"%"),ylab=paste("PC 2 - ",round(propvarexp[1,2]*100,2),"%"),pch=20,cex=2,
     col= ifelse(pca_faststr$ID_PCA=="duchesse_de_montebello_H7","navyblue",ifelse(pca_faststr$ID_PCA=="kean_D12","navyblue",
     ifelse(pca_faststr$ID_PCA=="berengere_B9","navyblue",ifelse(pca_faststr$ID_PCA=="blush_damask_H1","navyblue",
     ifelse(pca_faststr$ID_PCA=="hume_s_blush_tea_scented_china_C7","orange",ifelse(pca_faststr$ID_PCA=="old_blush_G12","orange",ifelse(pca_faststr$ID_PCA=="rosa_chinensis_mutabilis_G3","orange",
     ifelse(pca_faststr$ID_PCA=="mme_lierval_E9","purple",ifelse(pca_faststr$ID_PCA=="victor_verdier_chene_H6","purple",
     ifelse(pca_faststr$ID_PCA=="marguerite_de_roman_E12","purple",ifelse(pca_faststr$ID_PCA=="triomphe_de_l_exposition_A8","purple",ifelse(pca_faststr$ID_PCA=="alfred_colomb_C3","purple",                                                                                                                                      
     ifelse(pca_faststr$ID_PCA=="la_france_D7","red",ifelse(pca_faststr$ID_PCA=="la_favorite_lyon_A10","red","darkgrey")))))))))))))))

text(pca_faststr$Axis1,pca_faststr$Axis2+1, labels=pca_faststr$ID_PCA,cex=0.6,pos=4,
     col= ifelse(pca_faststr$ID_PCA=="duchesse_de_montebello_H7","navyblue",ifelse(pca_faststr$ID_PCA=="kean_D12","navyblue",
     ifelse(pca_faststr$ID_PCA=="berengere_B9","navyblue",ifelse(pca_faststr$ID_PCA=="blush_damask_H1","navyblue",
     ifelse(pca_faststr$ID_PCA=="hume_s_blush_tea_scented_china_C7","orange",ifelse(pca_faststr$ID_PCA=="old_blush_G12","orange",ifelse(pca_faststr$ID_PCA=="rosa_chinensis_mutabilis_G3","orange",
     ifelse(pca_faststr$ID_PCA=="mme_lierval_E9","purple",ifelse(pca_faststr$ID_PCA=="victor_verdier_chene_H6","purple",
     ifelse(pca_faststr$ID_PCA=="marguerite_de_roman_E12","purple",ifelse(pca_faststr$ID_PCA=="triomphe_de_l_exposition_A8","purple",ifelse(pca_faststr$ID_PCA=="alfred_colomb_C3","purple",                                                                                                                                      
     ifelse(pca_faststr$ID_PCA=="la_france_D7","red",ifelse(pca_faststr$ID_PCA=="la_favorite_lyon_A10","red","darkgrey")))))))))))))))



plot(pca_faststr$Axis1,pca_faststr$Axis2, xlab=paste("PC 1 - ",round(propvarexp[1,1]*100,2),"%"),ylab=paste("PC 2 - ",round(propvarexp[1,2]*100,2),"%"),pch=20,cex=2,
    col= ifelse(pca_faststr$ID_PCA=="duchesse_de_montebello_H7","navyblue",ifelse(pca_faststr$ID_PCA=="kean_D12","navyblue",
    ifelse(pca_faststr$ID_PCA=="berengere_B9","navyblue",ifelse(pca_faststr$ID_PCA=="blush_damask_H1","navyblue",
    ifelse(pca_faststr$ID_PCA=="hume_s_blush_tea_scented_china_C7","orange",ifelse(pca_faststr$ID_PCA=="old_blush_G12","orange",ifelse(pca_faststr$ID_PCA=="rosa_chinensis_mutabilis_G3","orange",
    ifelse(pca_faststr$ID_PCA=="mme_lierval_E9","purple",ifelse(pca_faststr$ID_PCA=="victor_verdier_chene_H6","purple",
    ifelse(pca_faststr$ID_PCA=="marguerite_de_roman_E12","purple",ifelse(pca_faststr$ID_PCA=="triomphe_de_l_exposition_A8","purple",ifelse(pca_faststr$ID_PCA=="alfred_colomb_C3","purple",                                                                                                                                      
     ifelse(pca_faststr$ID_PCA=="la_france_D7","red",ifelse(pca_faststr$ID_PCA=="la_favorite_lyon_A10","red",
    ifelse(pca_faststr$genetic_group=="06","deeppink4", ifelse(pca_faststr$genetic_group=="05","deeppink","darkgrey")))))))))))))))))

text(pca_faststr$Axis1,pca_faststr$Axis2+1, labels=pca_faststr$ID_PCA,cex=0.6,pos=4,
     col= ifelse(pca_faststr$ID_PCA=="duchesse_de_montebello_H7","navyblue",ifelse(pca_faststr$ID_PCA=="kean_D12","navyblue",
          ifelse(pca_faststr$ID_PCA=="berengere_B9","navyblue",ifelse(pca_faststr$ID_PCA=="blush_damask_H1","navyblue",
          ifelse(pca_faststr$ID_PCA=="hume_s_blush_tea_scented_china_C7","orange",ifelse(pca_faststr$ID_PCA=="old_blush_G12","orange",ifelse(pca_faststr$ID_PCA=="rosa_chinensis_mutabilis_G3","orange",
          ifelse(pca_faststr$ID_PCA=="mme_lierval_E9","purple",ifelse(pca_faststr$ID_PCA=="victor_verdier_chene_H6","purple",
          ifelse(pca_faststr$ID_PCA=="marguerite_de_roman_E12","purple",ifelse(pca_faststr$ID_PCA=="triomphe_de_l_exposition_A8","purple",ifelse(pca_faststr$ID_PCA=="alfred_colomb_C3","purple",                                                                                                                                      
          ifelse(pca_faststr$ID_PCA=="la_france_D7","red",ifelse(pca_faststr$ID_PCA=="la_favorite_lyon_A10","red",
          ifelse(pca_faststr$genetic_group=="06","deeppink4", ifelse(pca_faststr$genetic_group=="05","deeppink","darkgrey")))))))))))))))))

# PCA group

ggplot(pca_faststr,aes(x=Axis1,y=Axis2,colour=genetic_group))+geom_point(size=2.5)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = genetic_group, fill = genetic_group), color = 'black', size = 3.5)

# genetic group Liorzou


plotA<-ggplot(pca_faststr,aes(x=Axis1,y=Axis2,colour=as.factor(genetic_group)))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(ID), fill = as.factor(genetic_group)), color = 'black', size = 2)+
  labs(colour="Genetic Groups (SSR Liorzou)",fill="Genetic Groups (SSR Liorzou)")
   
plotB<-ggplot(pca_faststr,aes(x=Axis1,y=Axis3,colour=as.factor(genetic_group)))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 3 (8.0%)")+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))#+
  #geom_label_repel(aes(label = factor(substring(ID, 1, 12)), fill = as.factor(genetic_group)), color = 'black', size = 3.5)

plotC<-ggplot(pca_faststr,aes(x=Axis1,y=Axis4,colour=as.factor(genetic_group)))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 4 (5.6%)")+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))#+
  #geom_label_repel(aes(label = factor(substring(ID, 1, 12)), fill = as.factor(genetic_group)), color = 'black', size = 3.5)

### generate plot with 3 panels

pdf(height=16,width=24,file="PCA_205rosesCC_300622.pdf")
#grid.arrange(plotA, plotB, plotC, widths=c(2,1,1))

pushViewport(viewport(layout=grid.layout(2,3)))
vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col = y)
print(plotA, vp=vplayout(1:2,1:2))
print(plotB, vp=vplayout(1,3))
print(plotC, vp=vplayout(2,3))
dev.off()

# Date d'obtention

plotA<-ggplot(pca_faststr,aes(x=Axis1,y=Axis2,colour=as.factor(period_corjune2017_mathilde)))+geom_point(size=3.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(date_mathilde), fill = as.factor(period_corjune2017_mathilde)), color = 'black', size = 3.2)+
  labs(colour="Temporal breeding pattern",fill="Temporal breeding pattern")


plotB<-ggplot(pca_faststr,aes(x=Axis1,y=Axis3,colour=as.factor(period_corjune2017_mathilde)))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 3 (8.0%)")+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))#+
#geom_label_repel(aes(label = factor(substring(ID, 1, 12)), fill = as.factor(genetic_group)), color = 'black', size = 3.5)

plotC<-ggplot(pca_faststr,aes(x=Axis1,y=Axis4,colour=as.factor(period_corjune2017_mathilde)))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 4 (5.6%)")+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))#+
#geom_label_repel(aes(label = factor(substring(ID, 1, 12)), fill = as.factor(genetic_group)), color = 'black', size = 3.5)

### generate plot with 3 panels

pdf(height=12,width=20,file="PCA_205rosesCC_periodbreeding_300622.pdf")
#grid.arrange(plotA, plotB, plotC, widths=c(2,1,1))

pushViewport(viewport(layout=grid.layout(2,3)))
vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col = y)
print(plotA, vp=vplayout(1:2,1:2))
print(plotB, vp=vplayout(1,3))
print(plotC, vp=vplayout(2,3))
dev.off()

# polynormiale date d'obtention
ggplot(pca_faststr,aes(x=Axis1,y=date_mathilde))+
  stat_smooth(method='lm', formula = y~poly(x,2))+
  geom_point(size=2.5)+
  xlab("Axis 1 (PCA)")+
  ylab("Year")+
  ylim(1750,1910)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))
summary(lm(formula=pca_faststr$date_mathilde~poly(pca_faststr$Axis1,degree=2)))

ggplot(pca_faststr,aes(x=Axis2,y=date_mathilde))+
  stat_smooth(method='lm', formula = y~poly(x,2))+
  geom_point(size=2.5)+
  xlab("Axis 2 (PCA)")+
  ylab("Year")+
  ylim(1750,1910)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

summary(lm(formula=pca_faststr$date_mathilde~poly(pca_faststr$Axis2,degree=2)))


###
plotA <- ggplot(pca_faststr,aes(x=Axis1,y=Axis2,colour=as.factor(genetic_group)))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  ggtitle("'Observed' Ploidy")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(ploidy_obs), fill = as.factor(genetic_group)), color = 'black', size = 2)+
  labs(colour="Genetic Groups (SSR Liorzou)",fill="Genetic Groups (SSR Liorzou)")

plotB <- ggplot(pca_faststr,aes(x=Axis1,y=Axis2,colour=as.factor(genetic_group)))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  ggtitle("Expected Ploidy")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(ploidy_exp), fill = as.factor(genetic_group)), color = 'black', size = 2)+
  labs(colour="Genetic Groups (SSR Liorzou)",fill="Genetic Groups (SSR Liorzou)")

pdf(height=12,width=20,file="PCA_204rosesCC_ploidy_251021.pdf")
#grid.arrange(plotA, plotB, plotC, widths=c(2,1,1))
pushViewport(viewport(layout=grid.layout(2,1)))
vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col = y)
print(plotA, vp=vplayout(1,1))
print(plotB, vp=vplayout(2,1))
dev.off()



#### PCA with info TN
pca_TN=read.table("res_PCA_071021_withinfo_withTN.txt",h=T)

plotA<-ggplot(pca_TN,aes(x=Axis1,y=Axis2,colour=as.factor(genetic_group)))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(ID), fill = as.factor(genetic_group)), color = 'black', size = 2)+
  labs(colour="Genetic Groups (SSR Liorzou)",fill="Genetic Groups (SSR Liorzou)")

plotB<-ggplot(pca_TN,aes(x=Axis1,y=Axis2,colour=TN2015))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = TN2015, fill = TN2015), color = 'black', size = 3.5)+
  scale_fill_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred")+
  scale_colour_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred")


plotC<-ggplot(pca_TN,aes(x=Axis1,y=Axis2,colour=TN2016))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = TN2016, fill = TN2016), color = 'black', size = 3.5)+
  scale_fill_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred")+
  scale_colour_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred")

plotD<-ggplot(pca_TN,aes(x=Axis1,y=Axis2,colour=TN2014))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = TN2014, fill = TN2014), color = 'black', size = 3.5)+
  scale_fill_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred")+
  scale_colour_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred")

plotE<-ggplot(pca_TN,aes(x=Axis1,y=TN2014))+
  stat_smooth(method='lm', formula = y~poly(x,2))+
  geom_point(size=2.5)+
  xlab("Axis 1 (PCA)")+
  ylab("Score Blackspot Disease - Year 2014")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plotF<-ggplot(pca_TN,aes(x=Axis1,y=TN2015))+
  stat_smooth(method='lm', formula = y~poly(x,2))+
  geom_point(size=2.5)+
  xlab("Axis 1 (PCA)")+
  ylab("Score Blackspot Disease - Year 2015")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plotG<-ggplot(pca_TN,aes(x=Axis1,y=TN2016))+
  stat_smooth(method='lm', formula = y~poly(x,2))+
  geom_point(size=2.5)+
  xlab("Axis 1 (PCA)")+
  ylab("Score Blackspot Disease - Year 2016")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))



pdf(height=30,width=20,file="PCA_TN_204rosesCC_181021.pdf")
#grid.arrange(plotA, plotB, plotC, widths=c(2,1,1))

pushViewport(viewport(layout=grid.layout(5,3)))
vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col = y)
print(plotA, vp=vplayout(1,1:3))
print(plotB, vp=vplayout(2,1:3))
print(plotC, vp=vplayout(3,1:3))
print(plotD, vp=vplayout(4,1:3))
print(plotE, vp=vplayout(5,1))
print(plotF, vp=vplayout(5,2))
print(plotG, vp=vplayout(5,3))

dev.off()


summary(lm(formula=pca_TN$TN2014~poly(pca_TN$Axis1,degree=2)))
summary(lm(formula=pca_TN$TN2015~poly(pca_TN$Axis1,degree=2)))
summary(lm(formula=pca_TN$TN2016~poly(pca_TN$Axis1,degree=2)))


### flower

pca_TN_flower=read.table("res_PCA_071021_withinfo_withTN_withflower_withpos.txt",h=T)

ggplot(pca_TN_flower,aes(x=Axis1,y=Axis2,colour=mean_petals))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = substring(mean_petals,1,5), fill = mean_petals), color = 'black', size = 3.5)+
  scale_fill_gradient2(mid = "lightblue", high = "red")+
  scale_colour_gradient2(mid = "lightblue", high = "red")

ggplot(pca_TN_flower,aes(x=Axis1,y=mean_petals))+
  stat_smooth(method='loess')+
  geom_point(size=2.5)+
  xlab("Axis 1 (PCA)")+
  ylab("Averaged petal counted")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))





####
ggplot(pca_faststr,aes(x=Axis2,y=Axis3,colour=as.factor(genetic_group)))+geom_point(size=2.5)+
  xlab("Axis 2 (10.6%)")+
  ylab("Axis 3 (8.0%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(substring(ID, 1, 12)), fill = as.factor(genetic_group)), color = 'black', size = 3.5)

      
ggplot(pca_faststr,aes(x=Axis1,y=Axis2,colour=date_mathilde))+geom_point(size=2.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = date_mathilde, fill = date_mathilde), color = 'black', size = 3.5)+
  scale_fill_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred",midpoint=1750)+
  scale_colour_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred",midpoint=1750)
                                                                                                                                          
# Remontee
pca_remont=read.table("res_PCA_withinforemontee_year3_060722.txt",h=T)

pdf(file="PCA_repeatflowering_060722.pdf",width=16,height=12)
ggplot(pca_remont,aes(x=PCo1,y=PCo2,colour=as.factor(coderemonteeYear3)))+geom_point(size=3.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(Name), fill = as.factor(coderemonteeYear3)), color = 'black', size = 3.2)+
  labs(colour="Repeat flowering\n(NR= - ,R= + ,RF= ++)",fill="Repeat flowering\n(NR= - ,R= + ,RF= ++)")
dev.off()

pca_remont=read.table("res_PCA_withinforemontee_3years_harmonise_060722.txt",h=T)
pdf(file="PCA_repeatflowering_harmonise3years_070722.pdf",width=16,height=12)
ggplot(pca_remont,aes(x=PCo1,y=PCo2,colour=as.factor(CoderemonteeHarmonise3years)))+geom_point(size=3.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(Name), fill = as.factor(CoderemonteeHarmonise3years)), color = 'black', size = 3.2)+
  labs(colour="Repeat flowering\n(NR= - ,R= + ,RF= ++)",fill="Repeat flowering\n(NR= - ,R= + ,RF= ++)")
dev.off()

pdf(file="PCA_repeatflowering_060722_Nb_Clust.pdf",width=16,height=12)
ggplot(pca_remont,aes(x=PCo1,y=PCo2,colour=Nb_Clust_A3))+geom_point(size=3.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(Name), fill = Nb_Clust_A3), color = 'black', size = 3.2)+
  labs(colour="Repeat flowering\n(Nb_clust)",fill="Repeat flowering\n(Max_Reb)")+
  scale_fill_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred")+
  scale_colour_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred")
dev.off()



pdf(file="PCA_repeatflowering_060722_Rat_P.pdf",width=16,height=12)
ggplot(pca_remont,aes(x=PCo1,y=PCo2,colour=Rat_P_A3))+geom_point(size=3.5)+
  xlab("Axis 1 (23.9%)")+
  ylab("Axis 2 (10.6%)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(aes(label = factor(Name), fill = Rat_P_A3), color = 'black', size = 3.2)+
  labs(colour="Repeat flowering\n(Rat_P)",fill="Repeat flowering\n(Rat_P)")+
  scale_fill_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred")+
  scale_colour_gradient2(low = "navyblue", mid="lightgoldenrod1",high = "darkred")
dev.off()