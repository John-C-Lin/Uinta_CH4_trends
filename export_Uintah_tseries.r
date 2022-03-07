# Exports Uintah Basin CH4/CO2 time series 
# December 14th, 2020 by John C. Lin (John.Lin@utah.edu)

####################
YEARs <- 2018:2020
species <- c("CH4","CO2")
SITEs <- c("FRU","ROO","CSP","HPL")   # sites to output
obsdir <- "/uufs/chpc.utah.edu/common/home/lin-group4/jcl/SimCity/"
####################

for(ss in species){
  specie <- ss
for(i in 1:length(YEARs)){
  objname <- paste0("SimCity_",specie,"_allsites_hrly_",YEARs[i])
  print(paste("Reading in.....",objname))
  tmp <- getr(objname,path=obsdir)
  colnms <- colnames(tmp)
  sel <- colnms %in% paste0(specie,"_",SITEs)
  dat <- data.frame(Time_UTC=tmp$Time,tmp[,sel])

  outputfilenm <- paste0("Uintah_",specie,"_",YEARs[i],".csv")
  write.csv(dat,outputfilenm,row.names=FALSE)
  print(paste(outputfilenm,"generated"))
  gc()
} #for(i in 1:length(YEARs)){
} #for(ss in species){



