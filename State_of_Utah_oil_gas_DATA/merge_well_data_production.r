#Merges well data (includes well location information) with well-level detailed production information
#  to create detailed database 
#March 21st, 2021 by John C. Lin (John.Lin@utah.edu)

#####################
counties<-c("Duchesne","Uintah")
#Years<-1984:2018
#Years<-2000:2020
Years<-2000:2021
obsYrs<-2015:2020   #yrs when have surface CH4 obs
resultname<-"well_data_gas.oil.production_monthly.rds"
#####################

dat.all<-NULL
for(i in 1:length(counties)){
  county<-counties[i]
  xfiles<-list.files(pattern=paste0(county,"_monthly"))
  dat.all.sub<-NULL
  
  # load well history data
  wdfile <- paste0("well_data_",county,"_all.csv")
  wdat <- read.csv(wdfile)
for(ff in 1:length(xfiles)){
  #datfile<-"prod_Uintah_monthly_2011_2013.csv"
  datfile <- xfiles[ff]
  print(paste("Reading in:",datfile))
  wprod <- read.csv(datfile,stringsAsFactors=FALSE)
  
  dat <- merge(x=wprod,y=wdat,by="API.Well.Number",all.x=TRUE)
  # write.csv(dat,file="welldat.csv",row.names=FALSE)
  
  dat.all.sub<-rbind(dat.all.sub,dat)
  gc()
} # for(ff in 1:length(xfiles)){
  dat.all<-rbind(dat.all,dat.all.sub)
  gc()
} # for(i in 1:length(counties)){

# LOTS of duplicated data especially in 2010;  remove these first
label <- paste(dat.all$Report.Period,dat.all$API.Well.Number,dat.all$Bbls.Oil,dat.all$Mcf.Gas,dat.all$Bbls.Water)
sel <- duplicated(label)
#  figure out which years have the most duplicated instances  
tapply(substring(label[sel],7,10),substring(label[sel],7,10),length) # mostly in 2010
#  figure out whether replicated wells are mostly oil or gas wells
tapply(dat.all$Well.Type.x[sel],dat.all$Well.Type.x[sel],length)     # mostly oil wells
#sel2 <- label == "07/01/2010 4301315111 225 168 358"   # look out at some examples of repeated times
#sel2 <- label == "06/01/2020 4304756292 0 0 0"
#dat.all[sel2,]
dat.all <- dat.all[!sel,]   # remove duplicated rows

YYYY <- substring(dat.all$Report.Period,7,10)
MM <- substring(dat.all$Report.Period,1,2)
DD <- substring(dat.all$Report.Period,4,5)
YYYYMMDD <- paste0(YYYY,MM,DD)
dat.all <- data.frame(YYYYMMDD,dat.all)
saveRDS(dat.all,resultname)
print(paste(resultname,"generated"))

