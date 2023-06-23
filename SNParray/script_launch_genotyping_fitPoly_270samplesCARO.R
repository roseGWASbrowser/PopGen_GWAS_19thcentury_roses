# TL - 230421
# fitpoly used with ploidy = 4, so is equivalent to fitTetra.

library(fitPoly)
library(data.table)
setwd("/bigvol/benoit/Analyses/Temp_Tibo/Roses/")

# load data (data is generated from the script_QC_filtering_renaming...R)
markerData <- fread("AxiomGT1.270.summary.fitPoly.format_poly.dat")

markerData[,c(1,2)] <- lapply(markerData[,c(1,2)], factor) # need to have SampleName and MarkerName as factors
str(markerData)

### Create output directory

subDir <- "output_fitpoly_270samples_CARO"
mainDir <- getwd()

if (file.exists(subDir)){
  setwd(file.path(mainDir, subDir))
}else{
  dir.create(file.path(mainDir,subDir))
  setwd(file.path(mainDir,subDir))
}

# Run fitpoly
saveMarkerModels(ploidy=4,data=markerData,maxiter=500,ncores=30,try.HW=T,p.threshold=0.90,call.threshold=0.75,dip.filter=1,peak.threshold=1,
                 filePrefix="fitPoly_270CARO",rdaFiles=F,plot="fitted",allModelsFile=T)
