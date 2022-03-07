# Plots time series & diurnal cycle of CH4 in Uintah Basin and contrast 2020 (COVID-shutdown period) against previous years
# Based on "Fch4_simple.r"
# September 2nd, 2020 by John C. Lin (John.Lin@utah.edu)

require(ncdf4); require(fields); require(geosphere)
####################
YEARs<-2015:2020
MONs<-1:12  
HRs<-20:23  #[UTC]  only analyze afternoon hours, following Foster papers

SITE<-"HPL" #site to focus on for calculating F_CH4
#SITE<-"CSP" #site to focus on for calculating F_CH4
obsdir<-"/uufs/chpc.utah.edu/common/home/lin-group4/jcl/SimCity/"
STILTdir<-"/uufs/chpc.utah.edu/common/home/lin-group7/jcl/CH4_inversion/CH4_inversion_Uintah/CH4_inversion_Uintah_STILT/out/"
resultname<-paste("Fch4.",SITE,".rds",sep="")
#XLIMS<-c(-111,-108.6)   #longitude range of Uintah Basin
#YLIMS<-c(39.6,40.6)     #latitude range of Uintah Basin
XLIMS<-c(-110.6,-109)   #longitude range of Uintah Basin
YLIMS<-c(39.9,40.5)     #latitude range of Uintah Basin
####################

#########################################################
#I.  Calculate CH4 enhancements over background site (FRU)
YYYYMM.COVID <- c("202006","202007","202008")[2]   # COVID period, in YYYYMM format
MMsel <- substring(YYYYMM.COVID,5,6)
dat.all<-NULL
for(i in 1:length(YEARs)){
  objname<-paste0("SimCity_CH4_allsites_hrly_",YEARs[i])
  print(paste("Reading in.....",objname))
  tmp<-getr(objname,path=obsdir)
  if(paste0("CH4_",SITE)%in%colnames(tmp)){
    tmp2 <- tmp[,c("Time",paste0("CH4_",c(SITE,"FRU")))]
  } else {
    tmp2 <- data.frame(Time=tmp[,"Time"],NA,CH4_FRU=tmp[,"CH4_FRU"])
    colnames(tmp2)[2] <- paste0("CH4_",SITE)
  } # if(paste0("CH4_",SITE)%in%colnames(tmp)){
  dat.all<-rbind(dat.all,tmp2)
  gc()
} #for(i in 1:length(YEARs)){

dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
plot(dat.all[,c("Time","CH4_FRU")],main=paste("FRU\nUThrs:",paste(HRs,collapse=" ")),pch=16,
     xlab="Time")

dCH4<-dat.all[,paste0("CH4_",SITE)]-dat.all[,"CH4_FRU"]
dCH4<-data.frame(Time=dat.all[,"Time"],dCH4)   #[ppm]
YYYYMMDDHH<-format(dCH4[,"Time"],format="%Y%m%d%H")
sel1<-as.numeric(substring(YYYYMMDDHH,5,6))%in%MONs
sel2<-as.numeric(substring(YYYYMMDDHH,9,10))%in%HRs
dCH4<-data.frame(YYYYMMDDHH,dCH4,stringsAsFactors=FALSE)
dCH4<-dCH4[sel1&sel2,]


#calculate DAILY dCH4
Time.day<-format(dCH4[,"Time"],format="%Y%m%d")
dCH4.ave<-tapply(dCH4[,"dCH4"],Time.day,mean,na.rm=T)
Time<-strptime(names(dCH4.ave),"%Y%m%d",tz="GMT")
dCH4.ave<-data.frame(YYYYMMDD=names(dCH4.ave),Time,dCH4.ave,stringsAsFactors=FALSE)

dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)
plot(dCH4.ave[,c("Time","dCH4.ave")],main=paste(SITE,"\nUThrs:",paste(HRs,collapse=" ")),pch=16,
     xlab="Time",ylab="Daily Averge Enhancement over Background [ppm]")
# plot up months with COVID in 2020 (but also in previous years--want to remove periods with significant enhancements)
sel <- substring(dCH4.ave[,"YYYYMMDD"],5,6)%in%MMsel
points(dCH4.ave[sel,c("Time","dCH4.ave")],pch=16,col="red")
legend("topright",paste("months:",paste(MMsel,collapse=" ")),pch=16,col="red",bty="n",text.col="red",cex=1.3)

dCH4.ave<-data.frame(dCH4.ave,footsum=NA,foot.basin=NA)
sel<-substring(dCH4.ave[,"YYYYMMDD"],1,6)%in%c("201501","201502","201503","201504","201505")   #lack of obs--remove
dCH4.ave<-dCH4.ave[!sel,]


#########################################################
#II.  Compare DIURNAL CYCLES @ CSP from 2020 COVID period versus those from previous years
YYYYMM.COVID <- c("202003","202004","202005")[-1]   # COVID period, in YYYYMM format
MMsel<-substring(YYYYMM.COVID,5,6)
SITE <- "CSP"
dat.all<-NULL
for(i in 1:length(YEARs)){
  objname<-paste0("SimCity_CH4_allsites_hrly_",YEARs[i])
  print(paste("Reading in.....",objname))
  tmp<-getr(objname,path=obsdir)
  if(paste0("CH4_",SITE)%in%colnames(tmp)){
    tmp2 <- tmp[,c("Time",paste0("CH4_",c(SITE,"FRU")))]
  } else {
    tmp2 <- data.frame(Time=tmp[,"Time"],NA,CH4_FRU=tmp[,"CH4_FRU"])
    colnames(tmp2)[2] <- paste0("CH4_",SITE)
  } # if(paste0("CH4_",SITE)%in%colnames(tmp)){
  dat.all<-rbind(dat.all,tmp2)
  gc()
} #for(i in 1:length(YEARs)){
YYYYMMDDHH <- format(dat.all[,"Time"],format="%Y%m%d%H")
sel <- substring(YYYYMMDDHH,1,6)%in%c("201501","201502","201503","201504","201505")   #lack of obs--remove
dat.all <- dat.all[!sel,]

if(SITE=="CSP"){
  YYYYsel<-"2016"   
  YYYYMM.other <- paste0(YYYYsel,MMsel)
} #if(SITE=="CSP"){

sitecolnm<-paste0("CH4_",SITE)

# check which months have most non-NA data
YYYYMMDDHH<-format(dat.all[,"Time"],format="%Y%m%d%H")
YYYYMM<-substring(YYYYMMDDHH,1,6)
sumNotNA<-function(x)return(sum(!is.na(x)))
dum <- tapply(dat.all[,sitecolnm],YYYYMM,sumNotNA)
print(dum)


sel<-YYYYMM%in%YYYYMM.COVID
Cdat<-dat.all[sel,]
sel<-YYYYMM%in%YYYYMM.other
Odat<-dat.all[sel,]

# plot time series at SITE and FRU
YLIMs<-c(1.9,5)
dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3)
plot(Cdat[,c("Time",sitecolnm)],pch=16,col=c("orange"),ylim=YLIMs,ylab="CH4 [ppm]")
points(Cdat[,c("Time",paste0("CH4_FRU"))],pch=16,col=c("darkgray"))
legend("topright",pch=16,col=c("orange","darkgray"),legend=c(sitecolnm,"CH4_FRU"))
title(main=paste(YYYYMM.COVID,collapse="  "))

dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3)
plot(Odat[,c("Time",sitecolnm)],pch=16,col=c("orange"),ylim=YLIMs,ylab="CH4 [ppm]")
points(Odat[,c("Time",paste0("CH4_FRU"))],pch=16,col=c("darkgray"))
legend("topright",pch=16,col=c("orange","darkgray"),legend=c(sitecolnm,"CH4_FRU"))
title(main=paste(YYYYMM.other,collapse="  "))

# dev.new(); plot(Odat[1:48,c("Time","CH4_FRU")],pch=16,type="o",col=c("darkgray")); title(main=paste(YYYYMM.other,collapse="  "))
# dev.new(); plot(Cdat[1:48,c("Time","CH4_FRU")],pch=16,type="o",col=c("darkgray"));title(main=paste(YYYYMM.COVID,collapse="  "))

# calculate diurnal cycle of CH4 
diur.C <- tapply(Cdat[,sitecolnm],format(Cdat[,"Time"],format="%H"),mean,na.rm=T)
diur.O <- tapply(Odat[,sitecolnm],format(Odat[,"Time"],format="%H"),mean,na.rm=T)
UThrs <- as.numeric(names(diur.C))
dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3)
plot(UThrs,diur.C,ylim=range(c(diur.C,diur.O)),pch=16,type="o",col="red",xlab="Hour of Day [UT]",ylab="CH4 [ppm]")
lines(UThrs,diur.O,pch=16,type="o",col="black")
title(main=paste("Average Diurnal Cycle at",SITE))
legend("topright",c(paste(YYYYMM.other,collapse=","),paste(YYYYMM.COVID,collapse=",")),
       pch=16,lty=1,text.col=c("black","red"),col=c("black","red"),bty="n")


# calculate diurnal cycle of CH4 ENHANCEMENTS (CH4_CSP - CH4_FRU)
diurDCH4.C <- tapply(Cdat[,sitecolnm]-Cdat[,"CH4_FRU"],format(Cdat[,"Time"],format="%H"),mean,na.rm=T)
sigma.C <- tapply(Cdat[,sitecolnm]-Cdat[,"CH4_FRU"],format(Cdat[,"Time"],format="%H"),stdev,na.rm=T)
N.C <- tapply(Cdat[,sitecolnm]-Cdat[,"CH4_FRU"],format(Cdat[,"Time"],format="%H"),sumNotNA)
stderr.C <- sigma.C/sqrt(N.C)

diurDCH4.O <- tapply(Odat[,sitecolnm]-Odat[,"CH4_FRU"],format(Odat[,"Time"],format="%H"),mean,na.rm=T)
sigma.O <- tapply(Odat[,sitecolnm]-Odat[,"CH4_FRU"],format(Odat[,"Time"],format="%H"),stdev,na.rm=T)
N.O <- tapply(Odat[,sitecolnm]-Odat[,"CH4_FRU"],format(Odat[,"Time"],format="%H"),sumNotNA)
stderr.O <- sigma.O/sqrt(N.O)

dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3)
UThrs <- as.numeric(names(diurDCH4.C))
plot(UThrs,diurDCH4.C,ylim=range(c(diurDCH4.C,diurDCH4.O)),pch=16,
     type="o",col="red",xlab="Hour of Day [UT]",ylab="CH4 enhancement over FRU [ppm]")
segments(UThrs,diurDCH4.C-stderr.C,UThrs,diurDCH4.C+stderr.C,col="red")
lines(UThrs,diurDCH4.O,pch=16,type="o",col="black")
segments(UThrs,diurDCH4.O-stderr.O,UThrs,diurDCH4.O+stderr.O,col="black")
title(main=paste("Average Diurnal Cycle:",SITE,"- FRU"))
legend("topright",c(paste(YYYYMM.other,collapse=","),paste(YYYYMM.COVID,collapse=",")),
       pch=16,lty=1,text.col=c("black","red"),col=c("black","red"),bty="n")
#  number of non-NA data in different months--make sure there are no systematic differences between 2016 and 2020, since there is difference between March, April, and May
N.C <- tapply(Cdat[,sitecolnm]-Cdat[,"CH4_FRU"],format(Cdat[,"Time"],format="%m"),sumNotNA)
N.O <- tapply(Odat[,sitecolnm]-Odat[,"CH4_FRU"],format(Odat[,"Time"],format="%m"),sumNotNA)
Mons <- as.numeric(names(N.C))
dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3)
plot(Mons,N.C,ylim=range(c(N.C,N.O)),pch=16,type="o",col="red",xlab="Month",ylab="# of non-NA data points")
lines(Mons,N.O,pch=16,type="o",col="black")
legend("topright",c(paste(YYYYMM.other,collapse=","),paste(YYYYMM.COVID,collapse=",")),
       pch=16,lty=1,text.col=c("black","red"),col=c("black","red"),bty="n")
title(main=paste("# of non-NA CH4 obs during different months of:",SITE,"- FRU"))



#########################################################
#III.  Compare DIURNAL CYCLES @ HPL & WBB from 2020 COVID period versus those from previous years
##########
#YYYYMM.COVID <- c("202003","202004","202005")[-1]   # COVID period, in YYYYMM format
YYYYMM.COVID <- c("202006","202007","202008")[2]   # COVID period, in YYYYMM format
MMsel <- substring(YYYYMM.COVID,5,6)
SITEs <- c("WBB","HPL","ROO","FRU")[2]
##########

YEARs <- 2015:2020
dat.all<-NULL
for(i in 1:length(YEARs)){
  objname<-paste0("SimCity_CH4_allsites_hrly_",YEARs[i])
  print(paste("Reading in.....",objname))
  tmp <- getr(objname,path=obsdir)
  sel <- paste0("CH4_",SITEs)%in%colnames(tmp)
  if(sum(sel)==length(SITEs)){
    tmp2 <- tmp[,c("Time",paste0("CH4_",c(SITEs,"FRU")))]
  } else {
    tmp2 <- data.frame(Time=tmp[,"Time"],tmp[,paste0("CH4_",c(SITEs[sel],"FRU"))])
    tmp2 <- data.frame(tmp2,NA)
    colnames(tmp2)[ncol(tmp2)] <- paste0("CH4_",SITE[!sel])
  } # if(paste0("CH4_",SITE)%in%colnames(tmp)){
  dat.all<-rbind(dat.all,tmp2)
  gc()
} #for(i in 1:length(YEARs)){

for(ss in 1:length(SITEs)){
  SITE <- SITEs[ss]
  sitecolnm<-paste0("CH4_",SITE)
  YYYYMMDDHH <- format(dat.all[,"Time"],format="%Y%m%d%H")
  dat <- dat.all[substring(YYYYMMDDHH,5,6)%in%MMsel,]
  YYYYMMDDHH <- format(dat[,"Time"],format="%Y%m%d%H")

  dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3)
  ylims <- c(1.95,2.7)
  if(SITE=="WBB") ylims <- c(1.9,2.0)
  if(SITE=="FRU") ylims <- c(1.9,2.1)
  plot(0,0,type="n",xlab="Hour of Day [UT]",ylab="CH4 [ppm]",ylim=ylims,xlim=c(0,24))
  title(main=paste("Average Diurnal Cycle at",SITE,"\n Months:",paste(MMsel,collapse=",")))
  Years <- unique(substring(YYYYMMDDHH,1,4))
  cols <- c("darkblue","lightblue","darkgreen","lightgreen","orange","red")
for(yy in 1:length(Years)){
  sel <- substring(YYYYMMDDHH,1,4)==Years[yy]
  diur <- tapply(dat[sel,sitecolnm],format(dat[sel,"Time"],format="%H"),mean,na.rm=T)
  UThrs <- as.numeric(names(diur))
  lines(UThrs,diur,pch=16,type="l",col=cols[yy],lwd=2) 
} #for(yy in 1:length(Years)){
  legend("topright",Years,bty="n",lwd=2,text.col=cols,col=cols,lty=1)

  dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3)
  ylims <- c(0,1.0)
  if(SITE=="WBB") ylims <- c(0,0.1)
  plot(0,0,type="n",xlab="Hour of Day [UT]",ylab="CH4 enhancement over FRU [ppm]",ylim=ylims,xlim=c(0,24))
  title(main=paste("Average Diurnal Cycle:",SITE,"- FRU,\n Months:",paste(MMsel,collapse=",")))
  Years <- unique(substring(YYYYMMDDHH,1,4))
  cols <- c("darkblue","lightblue","darkgreen","lightgreen","orange","red")
for(yy in 1:length(Years)){
  sel <- substring(YYYYMMDDHH,1,4)==Years[yy]
  diur <- tapply(dat[sel,sitecolnm]-dat[sel,"CH4_FRU"],format(dat[sel,"Time"],format="%H"),mean,na.rm=T)
  UThrs <- as.numeric(names(diur))
  lines(UThrs,diur,pch=16,type="l",col=cols[yy],lwd=2)
} #for(yy in 1:length(Years)){
  legend("topright",Years,bty="n",lwd=2,text.col=cols,col=cols,lty=1)
} #for(ss in 1:length(SITEs)){


