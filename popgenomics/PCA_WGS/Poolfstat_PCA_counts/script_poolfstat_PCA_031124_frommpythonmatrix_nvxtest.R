# TL - 031124
# try to generate a PCA with poolfstat following the suggestion from Mathieu Gautier for PCI Evol Biol

# I start by generating a file with two columns per pop with the AD field from a VCF to import the data 
# see script script_generateallelecounts.py, then I exclude the 4 first columns: awk '$1 = ""; $2= ""; $3 = ""; $4 = ""; {print $0}' rose32.jointvcf.clean.PASSonly.varonly.vcf.allele_counts_table.tsv | sed 's/    //g' | sed 's/ \+/\t/g'
# then a bit of cleaning since it seems there are a few empty lines in the inputfiles "grep -Ev "^ $" rose32.jointvcf.clean.PASSonly.varonly.vcf.allele_counts_table.tsv.withoutSNPinfo > rose32.jointvcf.clean.PASSonly.varonly.vcf.allele_counts_table.tsv.withoutSNPinfo.clean"
# pca is then performed with randomallele

library(poolfstat)
library(vcfR)

#colnames(allele_counts) <- c("Ref_Count_sample1", "Alt_Count_sample1", "Ref_Count_sample2", "Alt_Count_sample2", "Ref_Count_sample3", "Alt_Count_sample3", "Ref_Count_sample4", "Alt_Count_sample4", "Ref_Count_sample5", "Alt_Count_sample5", "Ref_Count_sample6", "Alt_Count_sample6", "Ref_Count_sample7", "Alt_Count_sample7", "Ref_Count_sample8", "Alt_Count_sample8", "Ref_Count_sample9", "Alt_Count_sample9", "Ref_Count_sample10", "Alt_Count_sample10", "Ref_Count_sample11", "Alt_Count_sample11", "Ref_Count_sample12", "Alt_Count_sample12", "Ref_Count_sample13", "Alt_Count_sample13", "Ref_Count_sample14", "Alt_Count_sample14", "Ref_Count_sample15", "Alt_Count_sample15", "Ref_Count_sample16", "Alt_Count_sample16","Ref_Count_sample17", "Alt_Count_sample17", "Ref_Count_sample18", "Alt_Count_sample18", "Ref_Count_sample19", "Alt_Count_sample19", "Ref_Count_sample20", "Alt_Count_sample20", "Ref_Count_sample21", "Alt_Count_sample21", "Ref_Count_sample22", "Alt_Count_sample22", "Ref_Count_sample23", "Alt_Count_sample23", "Ref_Count_sample24", "Alt_Count_sample24", "Ref_Count_sample25", "Alt_Count_sample25", "Ref_Count_sample26", "Alt_Count_sample26", "Ref_Count_sample27", "Alt_Count_sample27", "Ref_Count_sample28", "Alt_Count_sample28", "Ref_Count_sample29", "Alt_Count_sample29", "Ref_Count_sample30", "Alt_Count_sample30", "Ref_Count_sample31", "Alt_Count_sample31", "Ref_Count_sample32", "Alt_Count_sample32")

pool_data <-genotreemix2countdata(genotreemix.file = "rose32.jointvcf.clean.PASSonly.varonly.vcf.allele_counts_table.tsv.withoutSNPinfo.clean.inputtreemix.withheader.68M") #rose32.jointvcf.clean.PASSonly.varonly.vcf.allele_counts_table.tsv.withoutSNPinfo.clean.inputtreemix.withheader")

#pool_data <- pooldata(
#  refallele.count = allele_counts[, c("Ref_Count_sample1", "Ref_Count_sample2", "Ref_Count_sample3","Ref_Count_sample4", "Ref_Count_sample5", "Ref_Count_sample6","Ref_Count_sample7", "Ref_Count_sample8", "Ref_Count_sample9","Ref_Count_sample10", "Ref_Count_sample11", "Ref_Count_sample12","Ref_Count_sample13","Ref_Count_sample14", "Ref_Count_sample15", "Ref_Count_sample16","Ref_Count_sample17", "Ref_Count_sample18", "Ref_Count_sample19","Ref_Count_sample20","Ref_Count_sample21", "Ref_Count_sample22","Ref_Count_sample23","Ref_Count_sample24", "Ref_Count_sample25", "Ref_Count_sample26","Ref_Count_sample27", "Ref_Count_sample28", "Ref_Count_sample29","Ref_Count_sample30","Ref_Count_sample31", "Ref_Count_sample32")],
#  altallele.count = allele_counts[, c("Alt_Count_sample1", "Alt_Count_sample2", "Alt_Count_sample3","Alt_Count_sample4", "Alt_Count_sample5", "Alt_Count_sample6","Alt_Count_sample7", "Alt_Count_sample8", "Alt_Count_sample9","Alt_Count_sample10", "Alt_Count_sample11", "Alt_Count_sample12","Alt_Count_sample13","Alt_Count_sample14", "Alt_Count_sample15", "Alt_Count_sample16","Alt_Count_sample17", "Alt_Count_sample18", "Alt_Count_sample19","Alt_Count_sample20","Alt_Count_sample21", "Alt_Count_sample22","Alt_Count_sample23","Alt_Count_sample24", "Alt_Count_sample25", "Alt_Count_sample26","Alt_Count_sample27", "Alt_Count_sample28", "Alt_Count_sample29","Alt_Count_sample30","Alt_Count_sample31", "Alt_Count_sample32")],
#  poolsize = ploidy_levels
#)

#ploidy_levels <- c(4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,2,3,2,4,2,2,2,4,4,2,2,2,4,2,2,4,4) # inferred ploidy in the exact same order than in the vcf
# 175580	175581	175583	175584	175585	175586	175587	175588	175623	175624	175626	175627	175628	177142	177143	HumeBlush	LaFrance	OB2n	Rarve	Rchimut	Rchisan	Rchispo	Rgaloff	Rmaja	Rminut	Rmosch	Rogig	Rpend	Rwich	Rxanth	Rxdama	YellowIsland

# poolfstat data object
#pool_data <- pooldata_from_counts(allele_counts, ploidy_levels)

# perform PCA from poolfstat following Mathieu Gautier's suggestion
pca_result <- randomallele.pca(pool_data)
write.table(pca_result$pop.loadings,file="rose32.jointvcf.clean.PASSonly.varonly.vcf.allele_counts_table.tsv.withoutSNPinfo.68M.PCAres.randomallele.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(pca_result$perc.var,file="rose32.jointvcf.clean.PASSonly.varonly.vcf.allele_counts_table.tsv.withoutSNPinfo.68M.PCAres.explainedvar_randomallele.txt",sep="\t",quote=FALSE,row.names=FALSE)

# generate basic PCA plots
# plot PCA for PC1 and PC2
pdf(width=10,height=10,file="rose32.jointvcf.clean.PASSonly.varonly.vcf.allele_counts_table.tsv.withoutSNPinfo.68M.PCA_poolfstat_randomallele_ploidy_level.pdf")
plot(
  pca_result$pop.loadings[,1], pca_result$pop.loadings[,2],
  xlab = "PC1", ylab = "PC2",
  #main = "PCA of Samples with Variable Ploidy",
  col = "blue", pch = 20)
#text(pca_result$pop.loadings[,1], pca_result$pop.loadings[,2], labels = rownames(pca_result$scores), pos = 3, cex = 0.8)

# plot PCA for PC1 and PC3
plot(
  pca_result$pop.loadings[,1], pca_result$pop.loadings[,3],
  xlab = "PC1", ylab = "PC3",
  #main = "PCA of Samples with variable ploidy",
  col = "blue", pch = 20)
#text(pca_result$scores[,1], pca_result$scores[,3], labels = rownames(pca_result$scores), pos = 3, cex = 0.8)

# plot PCA for PC1 and PC4
#plot(
#  pca_result$scores[,1], pca_result$scores[,4],
#  xlab = "PC1", ylab = "PC4",
#  #main = "PCA of Samples with variable ploidy",
#  col = "blue", pch = 20)
#text(pca_result$scores[,1], pca_result$scores[,4], labels = rownames(pca_result$scores), pos = 3, cex = 0.8)

# plot PCA for PC2 and PC3
#plot(
#  pca_result$scores[,2], pca_result$scores[,3],
#  xlab = "PC2", ylab = "PC3",
#  #main = "PCA of Samples with variable ploidy",
#  col = "blue", pch = 20)
#text(pca_result$scores[,2], pca_result$scores[,3], labels = rownames(pca_result$scores), pos = 3, cex = 0.8)


dev.off()
