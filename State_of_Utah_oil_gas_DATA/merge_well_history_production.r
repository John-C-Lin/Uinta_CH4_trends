#Merges well history (but does NOT include location information) with well-level detailed production information
#  to create detailed database 
#March 21st, 2021 by John C. Lin (John.Lin@utah.edu)

#####################
counties<-c("Duchesne","Uintah")
#Years<-1984:2018
Years<-2000:2020
obsYrs<-2015:2020   #yrs when have surface CH4 obs
resultname<-"gas.oil.production_monthly.rds"
#####################

dat.all<-NULL
for(i in 1:length(counties)){
  county<-counties[i]
  xfiles<-list.files(pattern=paste0(county,"_monthly"))
  dat.all.sub<-NULL
  
  # load well history data
  whfile <- paste0("well_history_",county,"_all.csv")
  whist <- read.csv(whfile)
for(ff in 1:length(xfiles)){
  #datfile<-"prod_Uintah_monthly_2011_2013.csv"
  datfile <- xfiles[ff]
  print(paste("Reading in:",datfile))
  wprod <- read.csv(datfile,stringsAsFactors=FALSE)
  
  dat <- merge(x=wprod,y=whist,by="API.Well.Number",all.x=TRUE)
  
  dat.all.sub<-rbind(dat.all.sub,dat)
} # for(ff in 1:length(xfiles)){
  dat.all.sub<-dat.all.sub[order(dat.all.sub$Time),]
  dat.all<-rbind(dat.all,dat.all.sub)
} # for(i in 1:length(counties)){

