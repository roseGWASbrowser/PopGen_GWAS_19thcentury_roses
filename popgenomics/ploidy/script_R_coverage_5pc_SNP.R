 # TL - 281024
library(ggplot2)
setwd("/home/tleroy/Papiers/Rose_19eme/PCIEvolBiol/First_round_PCI")
datacov=read.table("rose32.jointvcf.clean.PASSonly.SNPonly.vcf.5pcsites.vcf.bcftools.depth.withIDs",h=T)

apply(datacov,2,median,na.rm=TRUE) # here 2 meas column-wise, while 1 would mean row-wise
# doesn't work



# get mode of the distribution
# Create the function.
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


for (element in 3:length(colnames(datacov))){
  colID=colnames(datacov[element])
  mediancov=median(as.numeric(as.character(datacov[,element])),na.rm=TRUE)
  meancov=mean(as.numeric(as.character(datacov[,element])),na.rm=TRUE)
  modecov=getmode(as.numeric(as.character(datacov[,element])))
  print(c(colID,mediancov,meancov,modecov))
}


result <- getmode(v)
print(result)
