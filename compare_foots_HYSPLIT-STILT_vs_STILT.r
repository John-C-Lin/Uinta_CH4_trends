#Compare footprints from newly merged HYSPLIT-STILT executable versus STILT executable
#Based on "plot_foot_CH4_inversion_UintahV1.r"
#September 6th, 2020 by John C. Lin (John.Lin@utah.edu)

#########################
#SITE<-"INDY"
#SITE<-"HPL"
#SITE<-"ROO"
SITE<-"CSP"
#SITE<-"SLC"

# HYSPLIT-STILT output directory
outputdir<-paste0("CH4_inversion_Uintah_HYSPLIT-STILT/out/",SITE,"/")
# STILT output directory
outputdir.s<-paste0("CH4_inversion_Uintah_STILT/out/",SITE,"/")

#outputdir<-paste0("out/",SITE,"/by-id/")
pngfileTF<-FALSE; overwritepngfileTF<-FALSE
GoogleMapTF<-FALSE

#lat/lon range of plot
XLIMS<-c(-90,-80);YLIMS<-c(38,42)
#range of footprint plot (log10 limits)
ZLIMS<-c(-7,-1)

# Simulation timing, yyyy-mm-dd HH:MM:SS (UTC)
# t_start <- '2015-07-01 00:00:00'
# t_end   <- '2020-08-31 23:00:00'
t_start <- '2016-01-01 00:00:00'
t_end   <- '2016-12-31 23:00:00'
run_times <- seq(from = as.POSIXct(t_start, tz = 'UTC'),
                 to   = as.POSIXct(t_end, tz = 'UTC'),
                 by   = 'hour')
sel<-substring(run_times,12,13)%in%c("20","21","22","23")  #only analyze afternoon hours, following Foster papers
run_times<-run_times[sel]

calc.footTOT.TF <- FALSE
#########################


require(fields);require(ncdf4)
require(maps);require(RgoogleMaps)

#if(!SITE%in%c("WBB","BOS","INDY","JES","SLC","PORT","LA","SF"))stop(paste("Not valid site:",SITE))
if(SITE=="HPL"){lati <- 40.1434; long <- -109.4680; zagl <- 4.06}
if(SITE=="ROO"){lati <- 40.2941; long <- -110.0090; zagl <- 4.06}
if(SITE=="CSP"){lati <- 40.0509; long <- -110.0194; zagl <- 4.06}
if(SITE=="WBB"){lati <- 40.766 ; long <- -111.849; zagl <- 35}
if(SITE=="BOS"){lati <- 42.3470; long <- -71.0840; zagl <- 20}
if(SITE=="INDY"){lati <- 39.7833; long <- -86.1652; zagl <- 20}
if(SITE=="JES"){lati <- 39.1723; long <- -76.7765; zagl <- 20}
if(SITE=="SLC"){lati <- 40.6539; long <- -111.8878; zagl <- 20}
if(SITE=="LA"){lati <- 33.8700; long <- -118.2800; zagl <- 20}
if(SITE=="PORT"){lati <- 45.5100; long <- -122.6000; zagl <- 20}
if(SITE=="SF"){lati <- 37.7956; long <- -122.2794; zagl <- 20}

xmn <- long-2; xmx <- long+2
ymn <- lati-2; ymx <- lati+2
xres <- 0.01
yres <- xres

if(calc.footTOT.TF){
#loop over all receptor times 
YYYYMMDDHHMMs <- format(run_times,format="%Y%m%d%H%M")
RESULT <- data.frame(Time=run_times,HYSPLIT.STILT=NA,STILT=NA)
for(i in 1:length(run_times)){
  YYYYMMDDHHMM<-YYYYMMDDHHMMs[i]
  print(paste(" Date:",YYYYMMDDHHMM))
  ###############################
  #I.  HYSPLIT-STILT results
  receptor <- paste0(YYYYMMDDHHMM,"_",long,"_",lati,"_",zagl)
  ncfilenm<-paste0(outputdir,"/by-id/",receptor,"/",receptor,"_foot.nc") # e.g., 201802042100_-109.468_40.1434_4.06_foot.nc
  if(file.exists(ncfilenm)){
  #read in footprint and calculate total
    footfile<-nc_open(ncfilenm)
    flat<-ncvar_get(footfile,"lat");flon<-ncvar_get(footfile,"lon")
    foot.all<-ncvar_get(footfile,"foot")
    footsum<-apply(foot.all,c(1,2),sum)   #sum over backward time, but preserve spatial structure
    footTOT<-sum(footsum)                 #TOTAL footprint
    #  foot.log<-log10(footsum);foot.log[footsum==0]<-NA
    RESULT[i,"HYSPLIT.STILT"] <- footTOT
  } else { 
    print(paste("no HYSPLIT-STILT footprint output for:",receptor))
  } # if(file.exists(ncfilenm)){

  ###############################
  #II.  STILT results
  receptor <- paste0(substring(YYYYMMDDHHMM,1,10),"_",long,"_",lati,"_",zagl)
  ncfilenm<-paste0(outputdir.s,"/by-id/",receptor,"/",receptor,"_foot.nc") # e.g., 2018020421_-109.468_40.1434_4.06_foot.nc
  if(file.exists(ncfilenm)){
  #read in footprint and calculate total
    footfile<-nc_open(ncfilenm)
    flat<-ncvar_get(footfile,"lat");flon<-ncvar_get(footfile,"lon")
    foot.all<-ncvar_get(footfile,"foot")
    footsum<-apply(foot.all,c(1,2),sum)   #sum over backward time, but preserve spatial structure
    footTOT<-sum(footsum)                 #TOTAL footprint
    #  foot.log<-log10(footsum);foot.log[footsum==0]<-NA
    RESULT[i,"STILT"] <- footTOT
  } else { 
    print(paste("no STILT footprint output for:",receptor))
  } # if(file.exists(ncfilenm)){

  saveRDS(RESULT,"footTOT_HYSPLIT-STILT_STILT.RDS")
  gc()
} # for(i in 1:length(run_times)){
} # if(calc.footTOT.TF){

dat <- readRDS("footTOT_HYSPLIT-STILT_STILT.RDS")
# plot time series
dev.new(); par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
plot(dat[,c("Time","HYSPLIT.STILT")],pch=16,cex=0.5,xlab="Time",ylab="Total Footprint [ppm/(umole/m2/s)]",main=paste("Site:",SITE))
points(dat[,c("Time","STILT")],pch=16,cex=0.5,col="darkgray")
legend("topright",c("HYSPLIT-STILT","STILT"),text.col=c("black","darkgray"),col=c("black","darkgray"),pch=16)

# plot differences
dfootTOT <- dat$HYSPLIT.STILT-dat$STILT
dev.new(); par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
plot(x=dat$Time,y=dfootTOT,pch=16,cex=0.5,xlab="Time",ylab="Total Footprint Difference (HYSPLIT.STILT minus STILT)",
     main=paste("Site:",SITE))

# plot differences--subset of months
MONs <- c(4:10)
sel <- as.numeric(format(dat$Time,format="%m"))%in%MONs
dev.new(); par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
xmain <- paste("Site:",SITE,"\nMonths:",paste(MONs,collapse=" "))
plot(x=dat$Time[sel],y=dfootTOT[sel],pch=16,cex=0.5,xlab="Time",
     ylab="Total Footprint Diff (HYSPLIT.STILT minus STILT)",main=xmain)
# plot distribution of differences
dev.new(); par(cex.axis=1.3,cex.lab=1.2,cex.main=1.3)
hist(dfootTOT[sel],xlab="Total Footprint Diff (HYSPLIT.STILT minus STILT)\n[ppm/(umole/m2/s)]",main=xmain)
xbias <- mean(dfootTOT[sel],na.rm=T)
xstdev <- stdev(dfootTOT[sel],na.rm=T)
abline(v=xbias,col="red")
legend("topleft",c(paste("mean(HYSPLIT.STILT - STILT):",signif(xbias,3)),
                   paste("stdev(HYSPLIT.STILT - STILT):",signif(xstdev,3))),cex=1.1,bty="n")
# plot distribution of PERCENTAGE differences
dfootTOT.perc <- 100*(dat$HYSPLIT.STILT-dat$STILT)/dat$HYSPLIT.STILT
dev.new(); par(cex.axis=1.3,cex.lab=1.2,cex.main=1.3)
hist(dfootTOT.perc[sel],xlab="% Footprint Diff: 100*(HYSPLIT.STILT - STILT)/HYSPLIT.STILT",main=xmain)
xbias <- mean(dfootTOT.perc[sel],na.rm=T)
xstdev <- stdev(dfootTOT.perc[sel],na.rm=T)
abline(v=xbias,col="red")
legend("topleft",c(paste("mean:",signif(xbias,3),"%"),
                   paste("stdev:",signif(xstdev,3),"%")),cex=1.2,bty="n")






