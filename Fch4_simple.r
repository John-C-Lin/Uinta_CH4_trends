# Simple STILT analyses of Uintah Basin methane emissions
# Idea is to calculate AFTERNOON enhancements (CH4_HPL - CH4_FRU), & then divide by total footprint strength:
#   F_CH4 = (CH4_HPL - CH4_FRU)/sum(foot_i,j)
#
#  Then see how stable F_CH4 is over different seasons and years.
# V2(200911): use new HYSPLIT-STILT executable instead of STILT executable
# V3(200920): instead of annual production data, use monthly production data
# V4(210320): add CSP as receptor, save results that are site-specific
# V5(210412): add option to fill in gaps in FRU time series and merge in wind obs + sim to filter out times when HRRR-STILT was off using quadrant method
# V6(210416): filter out days based on within 45o of observed wind direction (better than quadrant method above, since not throw way data close to quadrant boundaries)
# V7(210422): calculate Energy-Normalized Methane Average (ENMA) emissions, make CH4 volume content explicit variable, combine monthly and annual emissions (for paper figure), create function to convert from [MCF NatGas] => [kg CH4]
# V8(210423): introduce 'mettype' to distinguish between met files driving STILT (to assess transport error)--e.g., "HRRR", "NARR", "NAM12"
# V9(211004): sensitivity analysis to change model domain
# May 29th, 2019 by John C. Lin (John.Lin@utah.edu)

require(ncdf4); require(fields); require(geosphere)
####################
YEARs <- 2015:2020
MONsub <- 4:9 #subset of months to calculate fluxes and leak rates
#MONsub <- 6:8 #subset of months to calculate fluxes and leak rates
HRs<-20:23  #[UTC]  only analyze afternoon hours, following Foster papers
#HRs<-17:23  #[UTC]  only analyze afternoon hours, following Foster papers
SITE<-"HPL" #HPL is site to focus on for calculating long-term F_CH4
#SITE<-"ROO" 
#SITE<-"CSP" 
if(SITE=="CSP")YEARs <- c(2016,2019,2020)
if(SITE=="ROO")YEARs <- 2015:2019

obsdir <- "/uufs/chpc.utah.edu/common/home/lin-group4/jcl/SimCity/"
winddir <- "/uufs/chpc.utah.edu/common/home/lin-group7/jcl/Transporterr_stiltread/Output"  #where wind obs and sim values are found

mettype <- "HRRR"
#mettype <- "HRRR_2000particles"
#mettype <- "NAM12"
#mettype <- "WRF27"
# Where STILT output is found
# a) STILT executable
#  STILTdir<-"/uufs/chpc.utah.edu/common/home/lin-group7/jcl/CH4_inversion/CH4_inversion_Uintah/CH4_inversion_Uintah_STILT/out/"
# b) HYSPLIT-STILT executable
#mettype <- toupper(mettype)
if(mettype=="HRRR")STILTdir<-"/uufs/chpc.utah.edu/common/home/lin-group7/jcl/CH4_inversion/CH4_inversion_Uintah/CH4_inversion_Uintah_HYSPLIT-STILT_HRRR/out/"
if(mettype=="HRRR_2000particles")STILTdir<-"/uufs/chpc.utah.edu/common/home/lin-group7/jcl/CH4_inversion/CH4_inversion_Uintah/CH4_inversion_Uintah_HYSPLIT-STILT_HRRR/out/"
if(mettype=="NAM12")STILTdir<-"/uufs/chpc.utah.edu/common/home/lin-group7/jcl/CH4_inversion/CH4_inversion_Uintah/CH4_inversion_Uintah_HYSPLIT-STILT_NAM12/out/"
if(mettype=="WRF27")STILTdir<-"/uufs/chpc.utah.edu/common/home/lin-group7/jcl/CH4_inversion/CH4_inversion_Uintah/CH4_inversion_Uintah_HYSPLIT-STILT_WRF27km/out/"
print(paste("STILTdir=",STILTdir))

resultname<-paste("Fch4_",SITE,"_daily_",mettype,".rds",sep="")

Nday.min <- 10          # minimum days necessary for a monthly average to be retained (otherwise assigned NA)

preparedata.TF <- TRUE  # run lines to prepare data?
fillFRU.TF <- TRUE      # fill in gaps in background (FRU) time series?
filterUV.TF <- TRUE     # filter times based on U/V (filter out times when HRRR is off--i.e., large transport errors) ?
domainsens.TF <- FALSE  # sensitivity analysis varyiing the domain size?
CH4.VOLFRAC <- 0.89     # volume fraction of CH4 in natural gas of 0.89 [Karion et al. 2013]
XLIMS<-c(-110.6,-109) # lon range of Uintah Basin--default
YLIMS<-c(39.9,40.5)   # lat range of Uintah Basin--default

if(domainsens.TF){
  print("vary domain from default; remember to set preparedata.TF to TRUE!")
  #  preparedata.TF <- TRUE  
  XLIMS<-c(-110,-109)    # lon range of Uintah Basin--change for sensitivity analysis (in response to Reviewer #1)
} # if(domainsens.TF){

####################

if(filterUV.TF&substring(mettype,1,4)!="HRRR")stop(paste0("Only extracted HRRR wind vectors; filterUV.TF needs to be set to FALSE for mettype=",mettype))

if(SITE=="HPL"){LAT <- 40.1434; LON <- -109.4680; zagl <- 4.06; WINDsite <- "UBHSP"; WINDsite.2015 <- "A1633"}  # UBHSP site not available in 2015, so use A1633 (RedWash) site instead in 2015
if(SITE=="ROO"){LAT <- 40.2941; LON <- -110.0090; zagl <- 4.06; WINDsite <- "QRS"}
if(SITE=="CSP"){LAT <- 40.0509; LON <- -110.0194; zagl <- 4.06; WINDsite <- "UBCSP"}

# Function that converts from natural gas volume in thousands of [MCF]  to [kg of CH4]
MCF2kg <- function(MCF,CH4.VOLFRAC=0.89){
  NatGas<-MCF*CH4.VOLFRAC   #[MCF natural gas] => [MCF CH4], using volume fraction of CH4 in natural gas of 0.89 [Karion et al. 2013]
  NatGas<-NatGas*1000   #[MCF CH4] => [ft^3 CH4]
  NatGas<-NatGas*0.0283 #[ft^3 CH4] => [m^3 CH4]
  R.CH4<-8.3143*1000/16.043   #Ideal Gas constant of CH4 [J/K/kg]: (Rg/mCH4), where Rg is universal gas constant and mCH4 is molar mass of CH4
  rho.CH4<-(101300)/(R.CH4*288.7)   #density of CH4 [kg/m3], using P=101.3kPa and T=288.7K [Karion et al. 2013]
  NatGas<-NatGas*rho.CH4 #[m^3 CH4] => [kg CH4]
  return(NatGas)
} # MCF2kg <- function(MCF,CH4.VOLFRAC=0.89){

if(preparedata.TF){
#########################################################
#I.  Calculate CH4 enhancements over background site (FRU)
dat.all<-NULL
for(i in 1:length(YEARs)){
  objname<-paste0("SimCity_CH4_allsites_hrly_",YEARs[i])
  print(paste("Reading in.....",objname))
  tmp<-getr(objname,path=obsdir)[,c("Time",paste0("CH4_",c(SITE,"FRU")))]
  print(colnames(tmp))
  dat.all<-rbind(dat.all,tmp)
  gc()
} #for(i in 1:length(YEARs)){

# filter for specified hrs [UTC]
YYYYMMDDHH <- format(dat.all$Time,"%Y%m%d%H")
dat.all <- data.frame(YYYYMMDDHH,dat.all)
sel <- as.numeric(substring(dat.all$YYYYMMDDHH,9,10))%in%HRs
dat.all <- dat.all[sel,]

# fill in gaps in background (FRU) time series
if(fillFRU.TF){
  dat <- dat.all
  YYYYMM <- substring(dat$YYYYMMDDHH,1,6)
  dat <- data.frame(YYYYMM,dat)
  # calculate % of NAs in each month
  N <- tapply(dat$YYYYMM,dat$YYYYMM,length)
  sumNA <- function(x)return(sum(is.na(x)))
  N.NA <- tapply(dat$CH4_FRU,dat$YYYYMM,sumNA)
  print("% of data with  NAs at FRU in each month:")
  print(round(100*N.NA/N,2))

  # calculate CV (coefficient of variation)--what % is variability in FRU compared to magnitude of CH4 enhancement?
  xsigma <- sqrt(tapply(dat$CH4_FRU,dat$YYYYMM,var,na.rm=T))
  dCH4.ave <- tapply(dat[,paste0('CH4_',SITE)]-dat$CH4_FRU,dat$YYYYMM,mean,na.rm=T)
  print('% stdev(FRU.hrly)/mean(SITE-FRU):')
  print(round(100*xsigma/dCH4.ave,2))
  #  calculate CV for DAILY-AVERAGED CH4 (for selected UThrs)
  FRU.day <- tapply(dat$CH4_FRU,substring(dat$YYYYMMDDHH,1,8),mean,na.rm=T)
  dCH4.day <- tapply(dat[,paste0('CH4_',SITE)]-dat$CH4_FRU,substring(dat$YYYYMMDDHH,1,8),mean,na.rm=T)
  YYYYMMDD <- names(FRU.day)
  xsigma <- sqrt(tapply(FRU.day,substring(YYYYMMDD,1,6),var,na.rm=T))
  dCH4.ave <- tapply(dCH4.day,substring(YYYYMMDD,1,6),mean,na.rm=T)
  print('% stdev(FRU.day)/mean(SITE-FRU):')
  print(round(100*xsigma/dCH4.ave,2))

  # fill in gaps in FRU time series with monthly average
  print("fill in gaps in FRU time series with monthly average")
  f <- function(x){
    x.ave <- mean(x,na.rm=T)
    x[is.na(x)] <- x.ave
    return(x)
  } # f <- function(x){
  FRU.gapfilled <- unlist(tapply(dat$CH4_FRU,dat$YYYYMM,f))
  if(TRUE){
    FRU.monave <- tapply(dat$CH4_FRU,dat$YYYYMM,mean,na.rm=T)
    tt <- as.numeric(substring(names(FRU.monave),1,4))+as.numeric(substring(names(FRU.monave),5,6))/12
    dev.new(); plot(tt,FRU.monave,pch=16,type="o")
    dev.new(); plot(dat[,c("Time","CH4_FRU")],pch=16)
    isNA <- is.na(dat$CH4_FRU)
    points(dat$Time[isNA],FRU.gapfilled[isNA],pch=16,col="orange")
    legend("topright",c("all","gap-filled"),pch=16,col=c("black","orange"))
  } #if(FALSE){
  dat$CH4_FRU <- FRU.gapfilled
  dat.all <- dat
} # if(fillFRU.TF){

# calculate CH4 enhancement over FRU [ppm]
dCH4 <- dat.all[,paste0('CH4_',SITE)]-dat.all$CH4_FRU
dat.all <- data.frame(dCH4,dat.all)

#calculate DAILY dCH4
dCH4 <- data.frame(Time=dat.all[,"Time"],dCH4)   #[ppm]
Time.day<-format(dCH4[,"Time"],format="%Y%m%d")
dCH4.ave<-tapply(dCH4[,"dCH4"],Time.day,mean,na.rm=T)
Time<-strptime(names(dCH4.ave),"%Y%m%d",tz="GMT")
dCH4.ave<-data.frame(YYYYMMDD=names(dCH4.ave),Time,dCH4.ave,stringsAsFactors=FALSE)
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
plot(dCH4.ave[,c("Time","dCH4.ave")],main=paste(SITE,"\nUThrs:",paste(HRs,collapse=" ")),pch=16,
     xlab="Time",ylab="Daily Averge Enhancement over Background [ppm]")
dCH4.ave<-data.frame(dCH4.ave,footsum=NA,foot.basin=NA)
sel<-substring(dCH4.ave[,"YYYYMMDD"],1,6)%in%c("201501","201502","201503","201504","201505","201506")   #lack of obs--remove
dCH4.ave<-dCH4.ave[!sel,]


dat.all <- data.frame(dat.all,footsum=NA,foot.basin=NA)
#########################################################
#II.  Calculate total footprint strengths and subset of footprint strength that just fall within Uintah Basin
footdir<-paste0(STILTdir,"/",SITE,"/footprints/")
if(mettype=="HRRR_2000particles") footdir<-paste0(STILTdir,"/",SITE,"_2000particles/footprints/")  #sensitive analysis with 10x base simulation particle number
for(i in 1:nrow(dat.all)){
  YYYYMMDDHH<-dat.all[i,"YYYYMMDDHH"]
  footfile<-list.files(footdir,pattern=YYYYMMDDHH)
  if(length(footfile)==0){print(paste("Cannot find any STILT footprints for",YYYYMMDDHH,"; skip"));next}
  if(length(footfile)> 1){print(paste("More than ONE STILT footprint file found for",YYYYMMDDHH,"; skip"));next}
  xfile <- paste0(footdir,footfile)
  print(paste("Reading in....",footfile))
  xnc <- nc_open(xfile)
  xlon <- ncvar_get(xnc,varid="lon"); xlat<-ncvar_get(xnc,varid="lat")
  xfoot <- apply(ncvar_get(xnc,varid="foot"),c(1,2),sum)
  # total footprint ALL doomain
  footsum <- sum(xfoot)  #[ppm/(umole/m2/s)] 
  dat.all[i,"footsum"] <- footsum
  # total footprint within Uintah Basin 
  sel.x <- xlon>=XLIMS[1]&xlon<=XLIMS[2]
  sel.y <- xlat>=YLIMS[1]&xlat<=YLIMS[2]
  foot.basin <- sum(xfoot[sel.x,sel.y])
  dat.all[i,"foot.basin"]<-foot.basin

  saveRDS(dat.all,file=paste0("Fch4_",SITE,"_",mettype,".rds"))
  nc_close(xnc)
  gc()
} #for(i in 1:nrow(dCH4)){


#########################################################
#III.  Merge with wind data (obs and HRRR simulated winds) to enable times when HRRR are off to be filtered out
dat.all <- readRDS(paste0("Fch4_",SITE,"_",mettype,".rds"))
wind.all <- NULL
for(yy in 1:length(YEARs)){
  # read in HRRR winds + observed values
  windname <- paste0("HRRR_obs_",YEARs[yy],"01010000to",YEARs[yy],"12312300.RDS")
  if(YEARs[yy]=="2018")windname <- paste0("HRRR_obs_",YEARs[yy],"01010000to",YEARs[yy],"12310000.RDS")  # typo during 2018
  colnms <- c("Time",paste0(c("Usim.","Vsim.","Tsim.","Uobs.","Vobs.","Tobs."),WINDsite))
  # UBHSP site not available in 2015, so use A1633 (RedWash) site instead in 2015
  if(SITE=="HPL")colnms <- c("Time",paste0(c("Usim.","Vsim.","Tsim.","Uobs.","Vobs.","Tobs."),WINDsite),
                                    paste0(c("Usim.","Vsim.","Tsim.","Uobs.","Vobs.","Tobs."),WINDsite.2015))
  if(file.exists(paste0(winddir,"/",windname))) {
    print(paste("Reading in.....",windname))
    tmp0 <- readRDS(paste0(winddir,"/",windname))$dat
    tmp <- matrix(NA,nrow=nrow(tmp0),ncol=length(colnms))
    colnames(tmp) <- colnms
    tmp <- data.frame(tmp)
    colnms.sel <- colnames(tmp0)%in%colnms
    tmp[,colnames(tmp0)[colnms.sel]] <- tmp0[,colnms.sel]
  } else {
    print(paste("Doesn't exists:",windname,"; replaces with NAs"))
    tmp.dat <- dat.all[substring(dat.all$YYYYMM,1,4)==YEARs[yy],]
    tmp <- matrix(NA,nrow=nrow(tmp.dat),ncol=length(colnms))
    colnames(tmp) <- colnms
    tmp <- data.frame(tmp)
    tmp$Time <- tmp.dat$Time
  } # if(file.exists(paste0(winddir,"/",windname)){
  wind.all <- rbind(wind.all,tmp)

  gc()
} # for(yy in 1:length(YEARs)){

  # merge dat.all with wind data
  result <- merge(x=dat.all,y=wind.all,by='Time',all.x=TRUE,all.y=FALSE)
  dat.all <- result
  saveRDS(dat.all,file=paste0("Fch4_",SITE,"_",mettype,".rds"))
  print(paste0("Fch4_",SITE,"_",mettype,".rds"," written out"))
  write.csv(dat.all,file=paste0("Fch4_",SITE,"_",mettype,".csv"),row.names=FALSE)
  print(paste0("Fch4_",SITE,"_",mettype,".csv"," written out"))

} # if(preparedata.TF){

#########################################################
#IV.   Filter times based on U/V (filter out times when HRRR is off--i.e., large transport errors) ?
dat.all <- readRDS(paste0("Fch4_",SITE,"_",mettype,".rds"))
if(filterUV.TF){ 
  dat <- dat.all
  Usim <- dat[,paste0("Usim.",WINDsite)]; Vsim <- dat[,paste0("Vsim.",WINDsite)]
  Uobs <- dat[,paste0("Uobs.",WINDsite)]; Vobs <- dat[,paste0("Vobs.",WINDsite)]
  if(SITE=="HPL"){
    # UBHSP site not available in 2015, so use A1633 (RedWash) site instead in 2015
    sel <- substring(dat$YYYYMM,1,4)=="2015"
    Usim[sel] <- dat[sel,paste0("Usim.",WINDsite.2015)]
    Vsim[sel] <- dat[sel,paste0("Vsim.",WINDsite.2015)]
    Uobs[sel] <- dat[sel,paste0("Uobs.",WINDsite.2015)]
    Vobs[sel] <- dat[sel,paste0("Vobs.",WINDsite.2015)]
  } # if(SITE=="HPL"){
  isNA <- is.na(Usim)|is.na(Vsim)|is.na(Uobs)|is.na(Vobs)
  # a) U/V filter:  both simulated U & V have to have the same sign as the observed (i.e., wind vector has to be same QUADRANT as observed)
  #  SEL <- sign(Usim)!=sign(Uobs)
  #  SEL <- sel|(sign(Vsim)!=sign(Vobs))

  # b) U/V filter:  filter out days based on +/-45o of observed (better than quadrant method above, since not throw way data close to quadrant boundaries
  WSPD.sim <- sqrt(Usim^2 + Vsim^2); WSPD.obs <- sqrt(Uobs^2 + Vobs^2)
  print("Observed windspeeds (quantile):")
  print(signif(quantile(WSPD.obs,c(0,0.05,0.1,0.25,0.5,0.75,1.0),na.rm=T),3))
  #  see:  http://tornado.sfsu.edu/geosciences/classes/m430/Wind/WindDirection.html
  WDIR.sim <- (180/pi)*atan2(y=-1*Vsim,x=-1*Usim)   # wind direction (-180o to +180o); direction wind is blowing FROM
  WDIR.obs <- (180/pi)*atan2(y=-1*Vobs,x=-1*Uobs)   # wind direction (-180o to +180o); direction wind is blowing FROM
  dWDIR <- WDIR.sim - WDIR.obs
  sel <- abs(dWDIR) > 180
  dWDIR[sel&!isNA] <- sign(dWDIR[sel&!isNA])*(360 - abs(dWDIR[sel&!isNA]))
  SEL <- abs(dWDIR) > 45   # simulated wind direction has to be within 45o of observed
  SEL[WSPD.obs < 1.0] <- FALSE   # not remove data if windspeed is too low (wind anemometer not reliable)

  print(paste("Out of",sum(!isNA),"non-NA U/V times,",sum(SEL[!isNA]),"failed to pass U/V filter: ",signif(100*sum(SEL[!isNA])/sum(!isNA),3),"%"))
  # assign NAs to footsum and foot.basin when failed to pass U/V filter
  print(paste("Before filtering, footsum has",sum(!is.na(dat$footsum)),"non-NA values"))
  dat$footsum[SEL] <- NA
  print(paste(" After filtering, footsum has",sum(!is.na(dat$footsum)),"non-NA values"))
  dat$foot.basin[SEL] <- NA
  dat.all <- dat
} # if(filterUV.TF){ 

#########################################################
#V.   Calculate daily averages
dat <- dat.all
Time.day<-format(dat$Time,format="%Y%m%d")
dCH4.ave<-tapply(dat$dCH4,Time.day,mean,na.rm=T)
Time<-strptime(names(dCH4.ave),"%Y%m%d",tz="GMT")
footsum<-tapply(dat$footsum,Time.day,mean,na.rm=T)
foot.basin<-tapply(dat$foot.basin,Time.day,mean,na.rm=T)

result<-data.frame(YYYYMMDD=names(dCH4.ave),Time,dCH4.ave,footsum,foot.basin,stringsAsFactors=FALSE)

#########################################################
#VI.   Divide CH4 enhancement by footprint strength to calculate CH4 fluxes
#      e.g., F_CH4 = (CH4_HPL - CH4_FRU)/sum(foot_i,j)
Fch4.sum<-result[,"dCH4.ave"]/result[,"footsum"]  #Fch4 calculated with total footprint
sel<-result[,"footsum"]<quantile(result[,"footsum"],prob=0.10,na.rm=T) #filter out cases when footprint is too weak...
Fch4.sum[Fch4.sum<0|sel]<-NA  
Fch4.basin<-result[,"dCH4.ave"]/result[,"foot.basin"] #Fch4 calculated with footprint only within Uintah Basin
sel<-result[,"foot.basin"]<quantile(result[,"foot.basin"],prob=0.10,na.rm=T) #filter out cases when footprint is too weak...
Fch4.basin[Fch4.basin<0|sel]<-NA
result<-data.frame(result,Fch4.sum,Fch4.basin)

saveRDS(result,resultname)
print(paste(resultname,"written out"))

dat<-readRDS(resultname)
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
plot(dat[,c("Time","Fch4.sum")],main=SITE,pch=16)
points(dat[,c("Time","Fch4.basin")],pch=16,col="blue")
legend("topright",c("sum(foot)","sum(foot.basin)"),col=c("black","blue"),text.col=c("black","blue"),pch=16,cex=1.3)

#filter out winter months
sel<-substring(dat[,"YYYYMMDD"],5,6)%in%c("12","01","02")
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
XMAIN <- paste0(SITE,"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
plot(dat[!sel,c("Time","Fch4.sum")],main=XMAIN,ylim=c(0,1),pch=16)
points(dat[!sel,c("Time","Fch4.basin")],pch=16,col="blue")
legend("topright",c("sum(foot)","sum(foot.basin)"),col=c("black","blue"),text.col=c("black","blue"),pch=16,cex=1.3)

#########################################################
#VII.   Calculate monthly averages, plot seasonal cycle and monthly averages
MM<-substring(dat[,"YYYYMMDD"],5,6)
Fch4.sum.seas<-tapply(dat[,"Fch4.sum"],MM,mean,na.rm=T)
Fch4.basin.seas<-tapply(dat[,"Fch4.basin"],MM,mean,na.rm=T)
NN<-tapply(dat[!is.na(dat[,"Fch4.sum"]),"Fch4.sum"],MM[!is.na(dat[,"Fch4.sum"])],length)
XMAIN <- paste0(SITE,"  Years: ",paste(YEARs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
plot(1:12,Fch4.sum.seas,main=XMAIN,pch=16,type="o")
lines(1:12,Fch4.basin.seas,pch=16,type="o",col="blue")
legend("topright",c("sum(foot)","sum(foot.basin)"),col=c("black","blue"),text.col=c("black","blue"),pch=16,cex=1.3,lty=1)

#plot time series of monthly averages
Time.mon<-format(dat[,"Time"],format="%Y%m")
Fch4.sum<-tapply(dat[,"Fch4.sum"],Time.mon,mean,na.rm=T)
Fch4.basin<-tapply(dat[,"Fch4.basin"],Time.mon,mean,na.rm=T)
Fch4.sum[is.nan(Fch4.sum)] <- NA; Fch4.basin[is.nan(Fch4.basin)] <- NA
#  calculate how many valid days there are in each month
Ntmp <- tapply(dat[!is.na(dat$Fch4.basin),"Fch4.basin"],Time.mon[!is.na(dat$Fch4.basin)],length)
N <- Fch4.basin; N[1:length(N)] <- NA;  N[names(Ntmp)] <- Ntmp
Fch4.basin[N<Nday.min] <- NA   # assign NA if # of valid days less than minimum
Fch4.sum[N<Nday.min] <- NA   # assign NA if # of valid days less than minimum
Fch4.basin.sd<-tapply(dat[,"Fch4.basin"],Time.mon,sd,na.rm=T)
Fch4.basin.sd[N<Nday.min] <- NA   # assign NA if # of valid days less than minimum

YYYYMM<-names(Fch4.sum)
frYr<-as.numeric(substring(YYYYMM,1,4))+(as.numeric(substring(YYYYMM,5,6))-0.5)/12
XMAIN <- paste0(SITE,"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
XSUB <- paste0("mettype=",mettype,"; Nday.min=",Nday.min)
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
plot(frYr,Fch4.sum,main=XMAIN,pch=16,type="o",ylim=c(0,0.2),xlab="Year",ylab="Fch4 [umole/m2/s]",sub=XSUB)
lines(frYr,Fch4.basin,type="o",pch=16,col="blue")
legend("topright",c("sum(foot)","sum(foot.basin)"),col=c("black","blue"),text.col=c("black","blue"),pch=16,cex=1.3,lty=1)

#########################################################
#VIII.   Convert to total emissions at Uintah Basin scale
dx1<-distHaversine(c(XLIMS[1],YLIMS[1]),c(XLIMS[2],YLIMS[1]))   #[m]
dx2<-distHaversine(c(XLIMS[1],YLIMS[2]),c(XLIMS[2],YLIMS[2]))   #[m]
dy1<-distHaversine(c(XLIMS[1],YLIMS[1]),c(XLIMS[1],YLIMS[2]))   #[m]
dy2<-distHaversine(c(XLIMS[2],YLIMS[1]),c(XLIMS[2],YLIMS[2]))   #[m]
AREA<-mean(c(dx1,dx2))*mean(c(dy1,dy2))   #[m^2]
Ech4.sum<-Fch4.sum*AREA*3600  #[umole/m2/s]=>[umole/hr]
Ech4.sum<-Ech4.sum*(12.0111+4*1.008)*1E-6/1000/(1000)  #[umole/hr]=>[10^3 kg/hr]
Ech4.basin<-Fch4.basin*AREA*3600  #[umole/m2/s]=>[umole/hr]
Ech4.basin<-Ech4.basin*(12.0111+4*1.008)*1E-6/1000/(1000)  #[umole/hr]=>[10^3 kg/hr]
Ech4.basin.sd<-Fch4.basin.sd*AREA*3600  #[umole/m2/s]=>[umole/hr]
Ech4.basin.sd<-Ech4.basin.sd*(12.0111+4*1.008)*1E-6/1000/(1000)  #[umole/hr]=>[10^3 kg/hr]

#########################################################
#IX.   Plot all months
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
XMAIN <- paste0(SITE,"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
plot(frYr,Ech4.sum,main=XMAIN,pch=16,type="o",xlab="Year",ylab="Emissions of CH4 from Uintah Basin [10^3 kg/hr]",ylim=c(0,100),sub=XSUB)
lines(frYr,Ech4.basin,type="o",pch=16,col="blue")
#abline(h=55,col="red",lty=2) #Uintah Basin-wide CH4 emissions [10^3 kg/hr], as reported by Karion et al. [2013]
legend("topright",c("sum(foot)","sum(foot.basin)"),col=c("black","blue"),text.col=c("black","blue"),pch=16,cex=1.3,lty=1)
dev.copy(png,"Fch4_simple_0.png");dev.off()

#Generate plots used for NCGG-8 presentation
Ech4.basin.seas<-Fch4.basin.seas*AREA*3600   #[umole/m2/s]=>[umole/hr]
Ech4.basin.seas<-Ech4.basin.seas*(12.0111+4*1.008)*1E-6/1000/(1000)  #[umole/hr]=>[10^3 kg/hr]
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3,mar=c(5,5,4,5))
XMAIN <- paste0(SITE,";  UThrs: ",paste(HRs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
plot(as.numeric(names(Ech4.basin.seas)),Ech4.basin.seas,main=XMAIN,sub=XSUB,
      pch=16,type="o",xlab="Month",ylab="Emissions of CH4 from Uintah Basin [10^3 kg/hr]")
#abline(h=55,col="red",lty=2,lwd=2) #Uintah Basin-wide CH4 emissions [10^3 kg/hr], as reported by Karion et al. [2013]
par(new=T)
plot(as.numeric(names(Fch4.basin.seas)),Fch4.basin.seas,xlab="",ylab="",pch=16,type="o",axes=F)
axis(4)
mtext("Average CH4 Flux [umole/m2/s]",side=4,cex=1.3,line=2)

#########################################################
#X.   Look at time series of Basin-Wide Emissions only for a subset of months
YYYYMM<-names(Ech4.basin)
frYr<-as.numeric(substring(YYYYMM,1,4))+(as.numeric(substring(YYYYMM,5,6))-0.5)/12
SEL<-as.numeric(substring(YYYYMM,5,6))%in%MONsub
#calculate error bars, sample numbers
dat<-readRDS(resultname)
isNA<-is.na(dat[,"Fch4.basin"]); dat2<-dat[!isNA,]
YYYYMM<-substring(dat2[,"YYYYMMDD"],1,6)
NN<-tapply(YYYYMM,YYYYMM,length)
Ech4.basin.stderr<-Ech4.basin.sd
Ech4.basin.stderr[names(Ech4.basin.stderr)]<-Ech4.basin.sd[names(Ech4.basin.stderr)]/sqrt(NN[names(Ech4.basin.stderr)])
Ech4.basin[names(NN)][NN<Nday.min]<-NA   #too small a sample size, so assign NA
XMAIN <- paste0(SITE,";  UThrs: ",paste(HRs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
XSUB <- paste0("mettype=",mettype,";  Mons: ",paste(MONsub,collapse=","),";  Nday.min=",Nday.min)
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3,mar=c(5,5,4,5))
plot(frYr[SEL],Ech4.basin[SEL],main=XMAIN,sub=XSUB,pch=16,lwd=2,
     xlab="Year",ylab="Emissions of CH4 from Uintah Basin [10^3 kg/hr]",type="p")
#abline(h=55,col="red",lty=2,lwd=2) #Uintah Basin-wide CH4 emissions [10^3 kg/hr], as reported by Karion et al. [2013]
segments(x0=frYr[SEL],y0=Ech4.basin[SEL]-Ech4.basin.stderr[SEL],x1=frYr[SEL],y1=Ech4.basin[SEL]+Ech4.basin.stderr[SEL])
devCUR <- dev.cur()
dev.copy(png,"Fch4_simple_1.png");dev.off()
#superimpose oil production to see whether there is any relationship
dev.set(devCUR)
frYr1 <- frYr[SEL]
tmp <- readRDS("State_of_Utah_oil_gas_DATA/gas.oil.production_monthly.rds")
counties <- unique(tmp$County)
Oil <- tapply(tmp[,"Oil..BBLs."],list(tmp[,"Time"]),sum,na.rm=TRUE)  #sum monthly values over counties
ylab <- "Oil [BBLs]"
YYYYMM <- paste0(substring(names(Oil),1,4),substring(names(Oil),6,7))
names(Oil) <- YYYYMM
frYr<-as.numeric(substring(YYYYMM,1,4))+(as.numeric(substring(YYYYMM,5,6))-0.5)/12
par(new=T)
SEL2 <- frYr>=min(frYr1) & frYr <=max(frYr1)
plot(frYr[SEL2],Oil[SEL2],pch=16,type="o",xlab="",ylab="",col="orange",axes=F,lwd=2)
axis(4,col="orange",col.axis="orange")
mtext(ylab,side=4,cex=1.3,line=2,col="orange")
dev.copy(png,"Fch4_simple_1.5.png");dev.off()


#Look at time series ONLY for a subset of months--plot DEVIATION from annual average emissinos
YYYYMM<-names(Ech4.basin)
frYr<-as.numeric(substring(YYYYMM,1,4))+(as.numeric(substring(YYYYMM,5,6))-0.5)/12
SEL<-as.numeric(substring(YYYYMM,5,6))%in%MONsub
#calculate error bars, sample numbers
dat<-readRDS(resultname)
isNA<-is.na(dat[,"Fch4.basin"]); dat2<-dat[!isNA,]
YYYYMM<-substring(dat2[,"YYYYMMDD"],1,6)
NN<-tapply(YYYYMM,YYYYMM,length)
Ech4.basin.stderr<-Ech4.basin.sd
Ech4.basin.stderr[names(Ech4.basin.stderr)]<-Ech4.basin.sd[names(Ech4.basin.stderr)]/sqrt(NN[names(Ech4.basin.stderr)])
Ech4.basin[names(NN)][NN<Nday.min]<-NA   #too small a sample size, so assign NA
XMAIN <- paste0(SITE,";  UThrs: ",paste(HRs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
XSUB <- paste0("mettype=",mettype,";  Mons: ",paste(MONsub,collapse=","),";  Nday.min=",Nday.min)
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3,mar=c(5,7,4,5))
Ech4.yrave <- tapply(Ech4.basin[SEL],substring(names(Ech4.basin[SEL]),1,4),mean,na.rm=T)
dEch4.basin <- Ech4.basin - Ech4.yrave[substring(names(Ech4.basin),1,4)]
names(dEch4.basin) <- names(Ech4.basin)
plot(frYr[SEL],dEch4.basin[SEL],main=XMAIN,sub=XSUB,pch=16,lwd=2,
     xlab="Year",ylab="DEVIAION from Annual Ave CH4 Emissions \nfrom Uintah Basin [10^3 kg/hr]",type="o")
#abline(h=55,col="red",lty=2,lwd=2) #Uintah Basin-wide CH4 emissions [10^3 kg/hr], as reported by Karion et al. [2013]
#superimpose oil production to see whether there is any relationship
frYr1 <- frYr[SEL]
tmp <- readRDS("State_of_Utah_oil_gas_DATA/gas.oil.production_monthly.rds")
counties <- unique(tmp$County)
Oil <- tapply(tmp[,"Oil..BBLs."],list(tmp[,"Time"]),sum,na.rm=TRUE)  #sum monthly values over counties
ylab <- "Oil [BBLs]"
YYYYMM <- paste0(substring(names(Oil),1,4),substring(names(Oil),6,7))
names(Oil) <- YYYYMM
frYr<-as.numeric(substring(YYYYMM,1,4))+(as.numeric(substring(YYYYMM,5,6))-0.5)/12
par(new=T)
SEL2 <- frYr>=min(frYr1) & frYr <=max(frYr1)
plot(frYr[SEL2],Oil[SEL2],pch=16,type="o",xlab="",ylab="",col="orange",axes=F,lwd=2)
axis(4,col="orange",col.axis="orange")
mtext(ylab,side=4,cex=1.3,line=2,col="orange")
dev.copy(png,"Fch4_simple_1.6.png");dev.off()
  R <- cor(dEch4.basin[SEL], Oil[names(dEch4.basin[SEL])], use="na.or.complete")



#########################################################
#XI.  Plot emissions as % of CH4 from Natural Gas PRODUCTION (monthly)
#convert units of natural gas & oil production
tmp<-readRDS("State_of_Utah_oil_gas_DATA/gas.oil.production_monthly.rds")
counties<-unique(tmp$County)
NatGas <- tapply(tmp[,"Natural.Gas..MCF."],list(tmp[,"Time"]),sum,na.rm=TRUE)  #sum monthly values over counties
NatGas <- MCF2kg(MCF=NatGas,CH4.VOLFRAC=CH4.VOLFRAC)  #[MCF Natural Gas per Month] => [kg CH4 per Month]
NatGas <- NatGas/1000  #[kg CH4 per Month] => [10^3 kg CH4 per Month]
NatGas <- NatGas/(24*(365/12)) # [10^3 kg CH4 per Month] => [10^3 kg CH4 per Hour]
ylab<-"Natural Gas [10^3 kg/hr]"
YYYYMM<-paste0(substring(names(NatGas),1,4),substring(names(NatGas),6,7))
names(NatGas)<-YYYYMM
production<-NatGas[names(Ech4.basin)]
YYYYMM<-names(Ech4.basin)
Year<-as.numeric(substring(YYYYMM,1,4))
MM<-as.numeric(substring(YYYYMM,5,6))
percEprod<-100*Ech4.basin/production
SEL<-MM%in%MONsub
frYr<-Year+(MM-0.5)/12
XMAIN <- paste0(SITE,";  UThrs: ",paste(HRs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
XSUB <- paste0("mettype=",mettype,";  Mons: ",paste(MONsub,collapse=","),";  Nday.min=",Nday.min,";  CH4.VOLFRAC=",CH4.VOLFRAC)
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3,mar=c(5,5,4,5))
plot(frYr[SEL],percEprod[SEL],main=XMAIN,sub=XSUB,pch=16,
     xlab="Year",ylab="Emissions of CH4 from Uintah Basin [% Production]")
dev.copy(png,"Fch4_simple_2.png");dev.off()

#Look at time series for a subset of years, and months
dat<-readRDS(resultname)
isNA<-is.na(dat[,"Fch4.basin"]); dat2<-dat[!isNA,]
dat<-dat2[as.numeric(substring(dat2[,"YYYYMMDD"],5,6))%in%MONsub,]
#plot time series of monthly averages
Time.mon<-format(dat[,"Time"],format="%Y%m")
Fch4.basin<-tapply(dat[,"Fch4.basin"],Time.mon,mean,na.rm=T)
Fch4.basin[is.nan(Fch4.basin)] <- NA
#  calculate how many valid days there are in each month
Ntmp <- tapply(dat[!is.na(dat$Fch4.basin),"Fch4.basin"],Time.mon[!is.na(dat$Fch4.basin)],length)
N <- Fch4.basin; N[1:length(N)] <- NA;  N[names(Ntmp)] <- Ntmp
Fch4.basin[N<Nday.min] <- NA   # assign NA if # of valid days less than minimum
YYYY<-substring(names(Fch4.basin),1,4)
#calculate error bars, sample numbers
NN<-tapply(YYYY,YYYY,length)
CONV<-AREA*3600*(12.0111+4*1.008)*1E-6/1000/(1000)  #[umole/m2/s]=>[10^3 kg/hr]
Fch4.yr<-tapply(Fch4.basin,YYYY,mean,na.rm=T)
Fch4.yr.sd<-tapply(Fch4.basin,YYYY,sd,na.rm=T)
Ech4.yr<-Fch4.yr*CONV
Ech4.yr.sd<-Fch4.yr.sd*CONV
Ech4.yr.stderr<-Ech4.yr.sd/sqrt(NN)
YYYY<-names(Ech4.yr)
frYr<-as.numeric(substring(YYYY,1,4))+0.5
XMAIN <- paste0(SITE,";  UThrs: ",paste(HRs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
XSUB <- paste0("mettype=",mettype,";  Mons: ",paste(MONsub,collapse=","),";  Nday.min=",Nday.min)
XSUB <- paste0(XSUB,"\n XLIMS=",round(XLIMS[1],2),",",round(XLIMS[2],2),"; YLIMS=",round(YLIMS[1],2),",",round(YLIMS[2],2))
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3,mar=c(5,5,4,5))
ylims <- c(18,60)
plot(frYr,Ech4.yr,main=XMAIN,sub=XSUB,pch=16,type="o",ylim=ylims,
     xlab="Year",ylab="Emissions of CH4 from Uintah Basin [10^3 kg/hr]")
#abline(h=55,col="red",lty=2,lwd=2) #Uintah Basin-wide CH4 emissions [10^3 kg/hr], as reported by Karion et al. [2013]
segments(x0=frYr,y0=Ech4.yr-Ech4.yr.stderr,x1=frYr,y1=Ech4.yr+Ech4.yr.stderr)
dev.copy(png,"Fch4_simple_4.png");dev.off()

#Plot emissions as % of CH4 from Natural Gas PRODUCTION (annual)
production<-readRDS("State_of_Utah_oil_gas_DATA/gas.oil.production_annual.rds")
rownames(production)<-production$Year
NatGas <- production[as.character(production$Year) %in% rownames(Ech4.yr),"Natural.Gas..MCF."]   # [MCF NatGas/Year]
NatGas <- MCF2kg(MCF=NatGas,CH4.VOLFRAC=CH4.VOLFRAC)  #[MCF NatGas/Year] => [kg CH4/Year]
NatGas <- NatGas/1000  #[kg CH4/Year] => [10^3 kg CH4/Year]
NatGas <- NatGas/(24*365) # [10^3 kg CH4/Year] => [10^3 kg CH4/Hour]
Ech4.yr.perc<-100*Ech4.yr/NatGas   #emission as % of CH4 in natural gas production
Ech4.yr.perc.stderr<-100*Ech4.yr.stderr/NatGas   #stderr as % of CH4 in natural gas production
XMAIN <- paste0(SITE,";  UThrs: ",paste(HRs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
XSUB <- paste0("mettype=",mettype,"; Mons: ",paste(MONsub,collapse=","),";  Nday.min=",Nday.min,";  CH4.VOLFRAC=",CH4.VOLFRAC)
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3,mar=c(5,5,4,5))
ylims<-c(0,10)
plot(frYr,Ech4.yr.perc,main=XMAIN,sub=XSUB,pch=16,type="o",
     xlab="Year",ylab="Emissions of CH4 from Uintah Basin [% Production]",ylim=ylims)
segments(x0=frYr,y0=Ech4.yr.perc-Ech4.yr.perc.stderr,x1=frYr,y1=Ech4.yr.perc+Ech4.yr.perc.stderr)
dev.copy(png,"Fch4_simple_5.png");dev.off()


#########################################################
#XII.  Plot emissions in ENERGY UNITS from Basin (as % of energy in CH4 emissions from Natural Gas + Oil PRODUCTION (annual))
#      This is called the Energy-Normalized Methane Average (ENMA) emissions in Robertson et al. (2017) ES&T paper
production <- readRDS("State_of_Utah_oil_gas_DATA/gas.oil.production_annual.rds")
Oil <- production[as.character(production$Year) %in% rownames(Ech4.yr),"Oil..BBLs."]
names(Oil) <- rownames(Ech4.yr)
Oil.BTU <- Oil * 5.8E6  # oil production in energy units [Barrels/Yr] => [BTU/Yr]; from Robertson et al. (2017)
NatGas <- production[as.character(production$Year) %in% rownames(Ech4.yr),"Natural.Gas..MCF."]   # [MCF NatGas/Year]
NatGas.BTU <- NatGas * 1E6   # gas production in energy units [MCF Natural Gas/Year] => [BTU/Year] from Robertson et al. (2017)
#  convert methane emissions back from [10^3 kg CH4 per Hour] => [MCF (thousands of cubic feet) natural gas/Year]
Egas <- data.frame(Ech4.yr,Ech4.yr.stderr)  # [10^3 kg CH4/hr]
Egas <- Egas*1000*24*365                 # [10^3 kg CH4/hr] => [kg CH4/Year]
CONV <- MCF2kg(MCF=1,CH4.VOLFRAC=CH4.VOLFRAC)  #[MCF Natural Gas] => [kg CH4] conversion factor
Egas <- Egas/CONV                        # [kg CH4/Year] => [MCF Natural Gas/Year]
Egas.BTU <- Egas * 1E6   # gas emissions in energy units [MCF Natural Gas/Year] => [BTU/Year] from Robertson et al. (2017)
ENBA <- 100*Egas.BTU/(Oil.BTU+NatGas.BTU)
dev.new();par(cex.axis=1.3,cex.lab=1.2,cex.main=1.3,mar=c(5,5,4,5))
ylims <- c(0,10)
#ylims <- NULL
plot(frYr,ENBA$Ech4.yr,main=XMAIN,sub=XSUB,pch=16,type="o",
     xlab="Year",ylab="Energy-Normalized Methane Average (ENMA) emissions\n [% Production]",ylim=ylims)
segments(x0=frYr,y0=ENBA$Ech4.yr-ENBA$Ech4.yr.stderr,x1=frYr,y1=ENBA$Ech4.yr+ENBA$Ech4.yr.stderr)
dev.copy(png,"Fch4_simple_6.png");dev.off()

#########################################################
#XIII.  Combine ENMA emissions as % of total NatGas+Oil production of energy with CH4 emissions as % of NatGas production [both %] in same plot
YYYY<-names(Ech4.yr.perc)
frYr<-as.numeric(substring(YYYY,1,4))+0.5
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3,mar=c(5,5,4,5))
ylims<-c(0,10)
#ylims<-c(0,18)
#ylims<-NULL
XMAIN <- paste0(SITE,";  UThrs: ",paste(HRs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
XSUB <- paste0("mettype=",mettype,";  Mons: ",paste(MONsub,collapse=","),";  Nday.min=",Nday.min,";  CH4.VOLFRAC=",CH4.VOLFRAC)
plot(frYr,Ech4.yr.perc,main=XMAIN,sub=XSUB,pch=16,type="o",
     xlab="Year",ylab="CH4 Emissions [% NatGas Production or % Energy Produced]",ylim=ylims,lwd=2)
segments(x0=frYr,y0=Ech4.yr.perc-Ech4.yr.perc.stderr,x1=frYr,y1=Ech4.yr.perc+Ech4.yr.perc.stderr)
lines(frYr,ENBA$Ech4.yr,main=XMAIN,sub=XSUB,pch=16,type="o",col="darkgray",lwd=2)
segments(x0=frYr,y0=ENBA$Ech4.yr-ENBA$Ech4.yr.stderr,x1=frYr,y1=ENBA$Ech4.yr+ENBA$Ech4.yr.stderr,col="darkgray")
legend("bottomleft",c("% NatGas Production","% Energy (NatGas+Oil)"),col=c("black","darkgray"),pch=16,lwd=2,text.col=c("black","darkgray"),cex=1.2,bty="n")
dev.copy(png,"Fch4_simple_7.png");dev.off()


#########################################################
#XIV.  Superimpose time series of monthly Basin-Wide Emissions with ANNUAL values
YYYYMM<-names(Ech4.basin)
frYr<-as.numeric(substring(YYYYMM,1,4))+(as.numeric(substring(YYYYMM,5,6))-0.5)/12
SEL<-as.numeric(substring(YYYYMM,5,6))%in%MONsub
# calculate error bars, sample numbers
dat<-readRDS(resultname)
isNA<-is.na(dat[,"Fch4.basin"]); dat2<-dat[!isNA,]
YYYYMM<-substring(dat2[,"YYYYMMDD"],1,6)
NN<-tapply(YYYYMM,YYYYMM,length)
Ech4.basin.stderr<-Ech4.basin.sd
Ech4.basin.stderr[names(Ech4.basin.stderr)]<-Ech4.basin.sd[names(Ech4.basin.stderr)]/sqrt(NN[names(Ech4.basin.stderr)])
Ech4.basin[names(NN)][NN<Nday.min]<-NA   #too small a sample size, so assign NA
XMAIN <- paste0(SITE,":  UThrs=",paste(HRs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
XSUB <- paste0("mettype=",mettype,";  Mons: ",paste(MONsub,collapse=","),";  Nday.min=",Nday.min,"; XLIMS=",
                round(XLIMS[1],2),",",round(XLIMS[2],2),"; YLIMS=",round(YLIMS[1],2),",",round(YLIMS[2],2))
dev.new();par(cex.axis=1.5,cex.lab=1.5,cex.main=1.5,cex.sub=1.0,mar=c(5,5,4,5))
ylims <- c(18,60)
xlims <- c(2015,2021)
plot(frYr[SEL],Ech4.basin[SEL],main=XMAIN,sub=XSUB,pch=16,lwd=2,ylim=ylims,xlim=xlims,
     #xlab="Year",ylab="Emissions of CH4 from Uintah Basin [10^3 kg/hr]",type="p",col="gray")
     xlab="Year",ylab=expression(paste("CH"[4]," Emissions [10"^3," kg hr"^-1,"]")),type="p",col="gray")
#abline(h=55,col="red",lty=2,lwd=2) #Uintah Basin-wide CH4 emissions [10^3 kg/hr], as reported by Karion et al. [2013]
segments(x0=frYr[SEL],y0=Ech4.basin[SEL]-Ech4.basin.stderr[SEL],x1=frYr[SEL],y1=Ech4.basin[SEL]+Ech4.basin.stderr[SEL],col="gray")
# add in ANNUAL means (with errors)
YYYY<-substring(names(Ech4.basin),1,4)
#calculate error bars, sample numbers
isNA <- is.na(Ech4.basin)
NN <- tapply(YYYY[SEL&!isNA],YYYY[SEL&!isNA],length)
Ech4.basin.yr <- tapply(Ech4.basin[SEL&!isNA],YYYY[SEL&!isNA],mean,na.rm=T)
Ech4.basin.yr.var <- tapply((Ech4.basin.stderr[SEL&!isNA])^2,YYYY[SEL&!isNA],sum)
Ech4.basin.yr.stderr <- sqrt(Ech4.basin.yr.var/NN)
YYYY<-names(Ech4.basin.yr)
frYr<-as.numeric(substring(YYYY,1,4))+0.5
lines(frYr,Ech4.basin.yr,pch=16,type="o",lwd=3)
segments(x0=frYr,y0=Ech4.basin.yr-Ech4.basin.yr.stderr,x1=frYr,y1=Ech4.basin.yr+Ech4.basin.yr.stderr)
dev.copy(png,"Fch4_simple_8.png");dev.off()


# sensitivity analysis varyiing the domain size?
if(domainsens.TF){
#########################################################
#XIV.5  Superimpose time series of monthly Basin-Wide Emissions with ANNUAL values
YYYYMM<-names(Ech4.basin)
frYr<-as.numeric(substring(YYYYMM,1,4))+(as.numeric(substring(YYYYMM,5,6))-0.5)/12
SEL<-as.numeric(substring(YYYYMM,5,6))%in%MONsub
# calculate error bars, sample numbers
dat<-readRDS(resultname)
isNA<-is.na(dat[,"Fch4.basin"]); dat2<-dat[!isNA,]
YYYYMM<-substring(dat2[,"YYYYMMDD"],1,6)
NN<-tapply(YYYYMM,YYYYMM,length)
Ech4.basin.stderr<-Ech4.basin.sd
Ech4.basin.stderr[names(Ech4.basin.stderr)]<-Ech4.basin.sd[names(Ech4.basin.stderr)]/sqrt(NN[names(Ech4.basin.stderr)])
Ech4.basin[names(NN)][NN<Nday.min]<-NA   #too small a sample size, so assign NA
XMAIN <- paste0(SITE,":  UThrs=",paste(HRs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
XSUB <- paste0("mettype=",mettype,";  Mons: ",paste(MONsub,collapse=","),";  Nday.min=",Nday.min,"; XLIMS=",
                round(XLIMS[1],2),",",round(XLIMS[2],2),"; YLIMS=",round(YLIMS[1],2),",",round(YLIMS[2],2))
dev.new();par(cex.axis=1.5,cex.lab=1.5,cex.main=1.5,cex.sub=1.0,mar=c(5,5,4,5))
ylims <- c(18,60)
xlims <- c(2015,2021)
plot(frYr[SEL],Ech4.basin[SEL],main=XMAIN,sub=XSUB,pch=16,lwd=2,ylim=ylims,xlim=xlims,
     xlab="Year",ylab=expression(paste("CH"[4]," Emissions [10"^3," kg hr"^-1,"]")),type="n",col="gray")
#segments(x0=frYr[SEL],y0=Ech4.basin[SEL]-Ech4.basin.stderr[SEL],x1=frYr[SEL],y1=Ech4.basin[SEL]+Ech4.basin.stderr[SEL],col="gray")
# add in ANNUAL means (with errors)
YYYY<-substring(names(Ech4.basin),1,4)
#calculate error bars, sample numbers
isNA <- is.na(Ech4.basin)
NN <- tapply(YYYY[SEL&!isNA],YYYY[SEL&!isNA],length)
Ech4.basin.yr <- tapply(Ech4.basin[SEL&!isNA],YYYY[SEL&!isNA],mean,na.rm=T)
Ech4.basin.yr.var <- tapply((Ech4.basin.stderr[SEL&!isNA])^2,YYYY[SEL&!isNA],sum)
Ech4.basin.yr.stderr <- sqrt(Ech4.basin.yr.var/NN)
YYYY<-names(Ech4.basin.yr)
frYr<-as.numeric(substring(YYYY,1,4))+0.5
  #############################################################################################
  # Correction due to smaller model domain--i.e., need to add in emissions from Uintah County #
  #############################################################################################
  #1) since new model domain used in sensitivity analysis limited to Uintah County, dominated by gas wells, so need to scale up to whole Basin by assuming that gas wells only emit ~85% of total, following Rella et al. [2015] stable isotope measurements
  #Ech4.basin.yr.domainsens <- Ech4.basin.yr/0.85  

  # prepare natural gas production data by county
  tmp<-readRDS("State_of_Utah_oil_gas_DATA/gas.oil.production_monthly.rds")
  counties<-unique(tmp$County)  # DUCHESNE or UINTAH
  NatGas.yr.bycounty <- NULL
  for(cc in 1:length(counties)){
    sel <- tmp$County == counties[cc]
    NatGas <- tapply(tmp[sel,"Natural.Gas..MCF."],list(tmp[sel,"Time"]),sum,na.rm=TRUE)  #sum monthly values over counties
    NatGas <- MCF2kg(MCF=NatGas,CH4.VOLFRAC=CH4.VOLFRAC)  #[MCF Natural Gas per Month] => [kg CH4 per Month]
    NatGas <- NatGas/1000  #[kg CH4 per Month] => [10^3 kg CH4 per Month]
    NatGas <- NatGas/(24*(365/12)) # [10^3 kg CH4 per Month] => [10^3 kg CH4 per Hour]
    YYYYMM <- paste0(substring(names(NatGas),1,4),substring(names(NatGas),6,7))
    names(NatGas) <- YYYYMM
    NatGas <- NatGas[names(Ech4.basin)]
    YYYYMMsel <- paste0(rep(YEARs,each=length(MONsub)),rep(formatC(MONsub,width=2,flag="0"),length(YEARs)))
    NatGas <- NatGas[YYYYMMsel]
    #  average production during selected months each year
    NatGas.yr <- tapply(NatGas,substring(names(NatGas),1,4),mean,na.rm=T)   # [10^3 kg CH4/Hr]
    NatGas.yr.bycounty[[cc]] <- NatGas.yr
  } # for(cc in 1:length(counties)){
  names(NatGas.yr.bycounty) <- counties

  #2) scale gas production in Duchesne county based on same emission characteristics as Uintah County using ratio in gas production
  # Ech4.basin.yr.domainsens <- Ech4.basin.yr*(1+NatGas.yr.bycounty$DUCHESNE/NatGas.yr.bycounty$UINTAH)

  #3) calculate emissions from Duchesne county as solely due to leakage from gas production from oil wells (using leakage rate from oil wells of 14.86%)
  Ech4.basin.yr.domainsens <- Ech4.basin.yr+NatGas.yr.bycounty$DUCHESNE*0.1486
  #############################################################################################
lines(frYr,Ech4.basin.yr,pch=16,type="o",lwd=2,col="gray",lty=2)
lines(frYr,Ech4.basin.yr.domainsens,pch=16,type="o",lwd=3,col="gray")
segments(x0=frYr,y0=Ech4.basin.yr.domainsens-Ech4.basin.yr.stderr,x1=frYr,y1=Ech4.basin.yr.domainsens+Ech4.basin.yr.stderr,col="gray")
# add in default emissions (without modifying model domain--output from "Fch4_simple_multimettype.r"
default <- readRDS(paste0(mettype,".rds"))
Ech4.yr <- default$Ech4.yr; YYYY<-names(Ech4.yr)
frYr<-as.numeric(substring(YYYY,1,4))+0.5
lines(frYr,Ech4.yr,lwd=2,type="o",pch=16)
legend("topright",c("Base case(Uintah&Duchesne)","Domain(Uintah)","Domain(Uintah) + Duchesne"),
       lty=c(1,2,1),lwd=c(3,2,2),col=c("black","gray","gray"),text.col=c("black","gray","gray"),cex=1.1)
dev.copy(png,"Fch4_simple_8.5.png");dev.off()
} # if(domainsens.TF){

#########################################################
#XV.  Plot emissions as % of CH4 from Natural Gas PRODUCTION during SELECTED MONTHS, instead of using ANNUAL PRODUCTION values
#convert units of natural gas & oil production
tmp<-readRDS("State_of_Utah_oil_gas_DATA/gas.oil.production_monthly.rds")
counties<-unique(tmp$County)
NatGas <- tapply(tmp[,"Natural.Gas..MCF."],list(tmp[,"Time"]),sum,na.rm=TRUE)  #sum monthly values over counties
NatGas <- MCF2kg(MCF=NatGas,CH4.VOLFRAC=CH4.VOLFRAC)  #[MCF Natural Gas per Month] => [kg CH4 per Month]
NatGas <- NatGas/1000  #[kg CH4 per Month] => [10^3 kg CH4 per Month]
NatGas <- NatGas/(24*(365/12)) # [10^3 kg CH4 per Month] => [10^3 kg CH4 per Hour]
YYYYMM <- paste0(substring(names(NatGas),1,4),substring(names(NatGas),6,7))
names(NatGas) <- YYYYMM
NatGas <- NatGas[names(Ech4.basin)]
YYYYMMsel <- paste0(rep(YEARs,each=length(MONsub)),rep(formatC(MONsub,width=2,flag="0"),length(YEARs)))
NatGas <- NatGas[YYYYMMsel]
#  average production during selected months each year
NatGas.yr <- tapply(NatGas,substring(names(NatGas),1,4),mean,na.rm=T)   # [10^3 kg CH4/Hr]
frYr<-as.numeric(names(NatGas.yr))+0.5
# Ech4.yr is in units of [10^3 kg CH4/hr]
Ech4.yr.perc <- Ech4.yr*100/NatGas.yr
XMAIN <- paste0(SITE,";  UThrs: ",paste(HRs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
XSUB <- paste0("mettype=",mettype,"; Mons: ",paste(MONsub,collapse=","),";  Nday.min=",Nday.min,";  CH4.VOLFRAC=",CH4.VOLFRAC)
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3,mar=c(5,5,4,5))
ylims<-c(1,10)
plot(frYr,Ech4.yr.perc,main=XMAIN,sub=XSUB,pch=16,type="o",
     xlab="Year",ylab="Emissions of CH4 from Uintah Basin [% Production]",ylim=ylims)
segments(x0=frYr,y0=Ech4.yr.perc-Ech4.yr.perc.stderr,x1=frYr,y1=Ech4.yr.perc+Ech4.yr.perc.stderr)
dev.copy(png,"Fch4_simple_9.png");dev.off()

#########################################################
#XVI.  Combine ENMA emissions as % of total NatGas+Oil production of energy with CH4 emissions as % of NatGas production [both %] in same plot
tmp<-readRDS("State_of_Utah_oil_gas_DATA/gas.oil.production_monthly.rds")
counties<-unique(tmp$County)
Oil <- tapply(tmp[,"Oil..BBLs."],list(tmp[,"Time"]),sum,na.rm=TRUE)  #sum monthly values over counties
Oil.BTU <- Oil * 5.8E6   # oil production in energy units [Barrels/Month] => [BTU/Month]; from Robertson et al. (2017)
NatGas <- tapply(tmp[,"Natural.Gas..MCF."],list(tmp[,"Time"]),sum,na.rm=TRUE)  #sum monthly values over counties
NatGas.BTU <- NatGas*1E6 # gas productioon in energy units [MCF Natural Gas/Month] => {BTU/Month];  from Robertson et al. (2017)  
#  select subset based on Year/Month
YYYYMMsel <- paste0(rep(YEARs,each=length(MONsub)),"-",rep(formatC(MONsub,width=2,flag="0"),length(YEARs)),"-01")
NatGas.BTU <- NatGas.BTU[YYYYMMsel] ; Oil.BTU <- Oil.BTU[YYYYMMsel] 
# average each year
NatGas.BTU.monave <- tapply(NatGas.BTU,substring(names(NatGas.BTU),1,4),mean,na.rm=T)
Oil.BTU.monave <- tapply(Oil.BTU,substring(names(Oil.BTU),1,4),mean,na.rm=T)
NatGas.BTU.hr <- NatGas.BTU.monave/(24*(365/12))  # [BTU per Month] => [BTU per Hour]
Oil.BTU.hr <- Oil.BTU.monave/(24*(365/12))        # [BTU per Month] => [BTU per Hour]
#  convert methane emissions back from [10^3 kg CH4 per Hour] => [MCF (thousands of cubic feet) natural gas/Year]
Egas <- data.frame(Ech4.yr,Ech4.yr.stderr) # [10^3 kg CH4/hr]
Egas <- Egas*1000                          # [10^3 kg CH4/hr] => [kg CH4/hr]
CONV <- MCF2kg(MCF=1,CH4.VOLFRAC=CH4.VOLFRAC)  #[MCF Natural Gas] => [kg CH4] conversion factor
Egas <- Egas/CONV                          # [kg CH4/hr] => [MCF Natural Gas/hr]
Egas <- Egas*CH4.VOLFRAC                   # [MCF Natural Gas/hr] => [MCF methane/hr]
Egas.BTU.hr <- Egas * 1E6   # gas emissions in energy units [MCF Natural Gas/hr] OR [MCF methane/hr]=> [BTU/hr] from Robertson et al. (2017)
ENBA <- 100*Egas.BTU.hr/(Oil.BTU.hr+NatGas.BTU.hr)

YYYY<-names(Ech4.yr.perc)
frYr<-as.numeric(substring(YYYY,1,4))+0.5
dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3,mar=c(5,5,4,5))
ylims<-c(1,10)
#ylims<-NULL
XMAIN <- paste0(SITE,";  UThrs: ",paste(HRs,collapse=","),"\nfillFRU.TF=",fillFRU.TF,";  filterUV.TF=",filterUV.TF)
XSUB <- paste0("mettype=",mettype,";  Mons: ",paste(MONsub,collapse=","),";  Nday.min=",Nday.min,";  CH4.VOLFRAC=",CH4.VOLFRAC)
plot(frYr,Ech4.yr.perc,main=XMAIN,sub=XSUB,pch=16,type="o",
     xlab="Year",ylab="CH4 Emissions [% NatGas Production or % Energy Produced]",ylim=ylims,lwd=2)
segments(x0=frYr,y0=Ech4.yr.perc-Ech4.yr.perc.stderr,x1=frYr,y1=Ech4.yr.perc+Ech4.yr.perc.stderr)
lines(frYr,ENBA$Ech4.yr,main=XMAIN,sub=XSUB,pch=16,type="o",col="darkgray",lwd=2)
segments(x0=frYr,y0=ENBA$Ech4.yr-ENBA$Ech4.yr.stderr,x1=frYr,y1=ENBA$Ech4.yr+ENBA$Ech4.yr.stderr,col="darkgray")
legend("bottomleft",c("% NatGas Production","% Energy (NatGas+Oil)"),col=c("black","darkgray"),pch=16,lwd=2,text.col=c("black","darkgray"),cex=1.2,bty="n")
dev.copy(png,"Fch4_simple_10.png");dev.off()


