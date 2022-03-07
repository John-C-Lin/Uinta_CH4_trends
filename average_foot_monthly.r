# Generate MONTHLY-averaged footprint MAPS
# V2(210428): introduce 'mettype' to distinguish between met files driving STILT (to assess transport error)--e.g., "HRRR", "NARR", "NAM12", "WRF27" 
# April 18th, 2021 by John C. Lin (John.Lin@utah.edu)


################
YEARs <- 2015:2020
#SITE<-"HPL"
#SITE<-"ROO"
SITE<-"CSP"
if(SITE=="CSP")YEARs <- c(2016,2019,2020)
if(SITE=="ROO")YEARs <- 2015:2019

# meteorology driving STILT
mettype <- "HRRR"
#mettype <- "NAM12"
#mettype <- "WRF27"
if(mettype=="HRRR")STILTdir<-"/uufs/chpc.utah.edu/common/home/lin-group7/jcl/CH4_inversion/CH4_inversion_Uintah/CH4_inversion_Uintah_HYSPLIT-STILT_HRRR/out/"
if(mettype=="NARR")STILTdir<-"/uufs/chpc.utah.edu/common/home/lin-group7/jcl/CH4_inversion/CH4_inversion_Uintah/CH4_inversion_Uintah_HYSPLIT-STILT_NARR/out/"
if(mettype=="NAM12")STILTdir<-"/uufs/chpc.utah.edu/common/home/lin-group7/jcl/CH4_inversion/CH4_inversion_Uintah/CH4_inversion_Uintah_HYSPLIT-STILT_NAM12/out/"
if(mettype=="WRF27")STILTdir<-"/uufs/chpc.utah.edu/common/home/lin-group7/jcl/CH4_inversion/CH4_inversion_Uintah/CH4_inversion_Uintah_HYSPLIT-STILT_WRF27km/out/"

obsdir<-"/uufs/chpc.utah.edu/common/home/lin-group4/jcl/SimCity/"  #where CH4 observations reside (to filter out times when SITE measured NAs)

# subset of hours [UTC]; if NULL, then all hours
# UThrs <- NULL   
UThrs <- c(20,21,22,23)  #[UTC]  only analyze afternoon hours, following Foster papers

finalresultname <- paste0(SITE,'_foot_monave_',mettype,'.rds')
################

outputdir<-paste0(STILTdir,"/",SITE,"/by-id/")

require(ncdf4);require(ggplot2)
# Footprint grid settings from 'run_stilt_hysplit-stiltV1.r'
xmn <- -112.5; xmx <- -108
ymn <- 38.5; ymx <- 41.5
xres <- 0.01; yres <- xres
FLAT <- seq(ymn+yres/2,ymx-yres/2,yres)
FLON <- seq(xmn+xres/2,xmx-xres/2,xres)

XFILES <- list.files(outputdir)

# Months and Years to grid
MMs <- formatC(1:12,flag="0",width=2)
YYYYMMs <- paste0(rep(YEARs,each=length(MMs)),rep(MMs,length(YEARs)))
RESULT.MON <- array(data=0, dim=c(length(FLON),length(FLAT),length(YYYYMMs)))
dimnames(RESULT.MON) <- list(FLON,FLAT,YYYYMMs)

for(yy in 1:length(YEARs)){
  YEAR <- YEARs[yy]
  # Simulation timing, yyyy-mm-dd HH:MM:SS (UTC)
  t_start <- paste0(YEAR,'-01-01 00:00:00')
  t_end   <- paste0(YEAR,'-12-31 23:00:00')
  run_times <- seq(from = as.POSIXct(t_start, tz = 'UTC'), to   = as.POSIXct(t_end, tz = 'UTC'), by   = 'hour')
  if(!is.null(UThrs)){
    sel <- substring(run_times,12,13)%in%UThrs  #only analyze afternoon hours, following Foster papers
    run_times <- run_times[sel]
  } # if(!is.null(UThrs)){

  # read in observations
  objname <- paste0("SimCity_CH4_allsites_hrly_",YEARs[yy])
  print(paste("Reading in.....",objname))
  obs <- getr(objname,path=obsdir)[,c("Time",paste0("CH4_",SITE))]

  # skip footprint if CH4_SITE is NA 
  run_times.sub <- subset(run_times,run_times%in%obs$Time)
  Times2rm <- obs[is.na(obs[,paste0("CH4_",SITE)]),'Time']
  sel <- run_times.sub%in%Times2rm
  print(paste('Skipping',sum(sel),'hrs out of total of',length(run_times.sub),'hrs'))
  run_times.sub <- run_times.sub[!sel]

  #!!! skip footprint if transport error is too large !!!

for(mm in 1:12){
  YYYYMM <- paste0(YEAR,formatC(mm,flag="0",width=2))
  run_times.tt <- run_times.sub[format(run_times.sub,"%Y%m")==YYYYMM]
  foot.sum <- 0
  N <- 0   #start counter
for(tt in 1:length(run_times.tt)){
  YYYYMMDDHH <- format(run_times.tt[tt],"%Y%m%d%H")
  sel <- substring(XFILES,1,10)==YYYYMMDDHH
  if(sum(is.na(XFILES[sel])>0)){print(paste("XFILES[sel] is NA for:",YYYYMMDDHH));next}
  ncfilenm <- paste0(outputdir,"/",XFILES[sel],"/",XFILES[sel],"_foot.nc")
  if(length(ncfilenm)>1){print(paste("more than 1 ncfilenm:",paste(ncfilenm,collapse=";  ")))}
  if(!file.exists(ncfilenm)){print(paste("no footprint output for:",YYYYMMDDHH));next}
  # retrieve footprint
  footfile <- nc_open(ncfilenm)
  # print(paste(N,"Reading in:",ncfilenm))
  foot.all <- ncvar_get(footfile,"foot")
  foot <- apply(foot.all,c(1,2),sum)   #sum over back times
  foot.sum <- foot.sum + foot
  N <- N+1
  print(paste(N,"Processed",YYYYMMDDHH))
  nc_close(footfile)
} # for(tt in 1:length(run_times.tt)){
  RESULT.MON[,,YYYYMM] <- foot.sum/N
  saveRDS(RESULT.MON,finalresultname)
  print(paste(finalresultname,"saved"))
  gc()
} # for(mm in 1:12){
  gc()
} # for(yy in 1:length(YEARs)){


if(FALSE){
  require(fields)

  SITE <- "CSP"
  if(SITE=="HPL"){LAT <- 40.1434; LON <- -109.4680; zagl <- 4.06}
  if(SITE=="ROO"){LAT <- 40.2941; LON <- -110.0090; zagl <- 4.06}
  if(SITE=="CSP"){LAT <- 40.0509; LON <- -110.0194; zagl <- 4.06}

  #YYYYMM <- "201604"
  YYYYMM <- "202004"
  resultnm <- paste0(SITE,'_foot_monave.rds')
  dat.all <- readRDS(resultnm)
  dat <- dat.all[,,YYYYMM]
  dat[dat==0] <- NA

  #a) log-10 plot
  dat.log10 <- log10(dat)
  zlims <- c(-7,-1)
  #zlims <- NULL
  flon <- as.numeric(rownames(dat))
  flat <- as.numeric(colnames(dat))
  if(!is.null(zlims)){dat.log10[dat.log10<zlims[1]]<-zlims[1];  dat.log10[dat.log10>zlims[2]]<-zlims[2]}
  dev.new();par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
  image.plot(flon,flat,dat.log10,main=paste("Site=",SITE,"\n",YYYYMM),zlim=zlims)
  points(LON,LAT,pch=17,col="black")

  #b) linear plot
  #zlims <- NULL
  zlims <- c(0.00001,0.0001)
  if(!is.null(zlims)){dat[dat<zlims[1]]<-zlims[1];  dat[dat>zlims[2]]<-zlims[2]}
  dev.new();par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
  image.plot(flon,flat,dat,main=paste("Site=",SITE,"\n",YYYYMM),zlim=zlims)
  points(LON,LAT,pch=17,col="black")

} # if(FALSE){
