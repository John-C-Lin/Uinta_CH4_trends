# Plots CH4 time series, diurnal, and seasonal cycle in Uintah Basin 
# V2(210521): prepare figures for publicatoin:  a) change colors  b) use subscripts c) remove UThrs, mons labelling
# April 26th, 2021 by John C. Lin (John.Lin@utah.edu)

####################
UTsel <- 20:23  # [UTC]  only analyze afternoon hours, following Foster papers
#obsdir <- "/uufs/chpc.utah.edu/common/home/lin-group15/jcl/SimCity/"
obsdir <- "/uufs/chpc.utah.edu/common/home/u0791084/PROJECTS/SimCity/GHGsites_scripts/"
col.gas <- "#B91C1C"  # red
col.oil <- "#1D4ED8"  # blue
####################

#########################################################
#I.  Compare CH4 time series 
##########
MONsel <- 4:9 # months to examine
YEARs <- 2015:2023
SITEs <- c("FRU","CSP","HPL")
##########
# combine data from different years
dat.all<-NULL
colnms.sel <- paste0("CH4_",SITEs)
for(i in 1:length(YEARs)){
  objname<-paste0("SimCity_CH4_allsites_hrly_",YEARs[i],".rds")
  print(paste("Reading in.....",objname))
  tmp <- readRDS(paste0(obsdir,"/",objname))
  dat.add <- matrix(NA,nrow=nrow(tmp),ncol=length(colnms.sel))
  colnames(dat.add) <- colnms.sel
  dat.add <- data.frame(dat.add)

  sel <- colnames(tmp)%in%colnms.sel
  dat.add[,colnames(tmp)[sel]] <- tmp[,sel]
  dat.add <- data.frame(Time=tmp$Time,dat.add)
  dat.all<-rbind(dat.all,dat.add)
  gc()
} #for(i in 1:length(YEARs)){

# average over UTsel
YYYYMMDDHH <- format(dat.all[,"Time"],format="%Y%m%d%H")
tmp <- dat.all[substring(YYYYMMDDHH,9,10)%in%formatC(UTsel,width=2,flag="0"),]
YYYYMMDD <- format(tmp[,"Time"],format="%Y%m%d")
CH4_FRU <- tapply(tmp$CH4_FRU,YYYYMMDD,mean,na.rm=T)
CH4_CSP <- tapply(tmp$CH4_CSP,YYYYMMDD,mean,na.rm=T)
CH4_HPL <- tapply(tmp$CH4_HPL,YYYYMMDD,mean,na.rm=T)
Time <- as.POSIXct(strptime(names(CH4_HPL),"%Y%m%d"),tz="GMT")

#ylims <- 1000*range(c(CH4_FRU,CH4_CSP,CH4_HPL),na.rm=T)
#ylims <- c(1900,3000)
ylims <- c(1850,2500)
#ylims <- c(1850,2300)
#ylims <- NULL
dev.new(width=8,height=6);par(cex.axis=1.5,cex.main=1.5,cex.lab=1.5,mar=c(5,5,4,2))
plot(Time,1000*CH4_HPL,pch=16,ylim=ylims,cex=0.6,col=col.gas,xlab="Year",ylab=expression("CH"[4]*" [ppb]"))
points(Time,1000*CH4_CSP,pch=16,cex=0.6,col=col.oil)
points(Time,1000*CH4_FRU,pch=16,cex=0.5,col="black")
title(main=expression("Observed Daily CH"[4]*" (afternoon-averaged)"))  # \n UThrs:",paste(UTsel,collapse=",")))
# draw gray boxes for times not included in analyses
xleft <- as.POSIXct(strptime(x=paste0(c(YEARs),"-",formatC(max(MONsel)+1,flag="0",width=2),"-","01"),
                    format="%Y-%m-%d"),tz="GMT")
xright <- as.POSIXct(strptime(x=paste0(c(YEARs+1),"-",formatC(min(MONsel),flag="0",width=2),"-","01"),
                     format="%Y-%m-%d"),tz="GMT")
rect(xleft=xleft,ybottom=ylims[1]-100,xright=xright,ytop=ylims[2]+100,col=gray(level=0.3,alpha=0.5),border="darkgray")
legend("topleft",c("HPL","CSP","FRU"),col=c(col.gas,col.oil,"black"),
       text.col=c(col.gas,col.oil,"black"),pch=16,bg="white",cex=1.3)
dev.copy(png,filename="CH4_tseries_allsites.png",width=480*8/6,height=480);dev.off()


#########################################################
#II.  Compare DIURNAL CYCLES @ HPL in different years
##########
YEARs <- 2015:2023
COLS <- c("darkblue","lightblue","darkgreen","lightgreen","orange","red","darkred","pink","magenta")
names(COLS) <- YEARs
MONsel <- 4:9 # months to examine
#MONsel <- 4:5 # months to examine
#MONsel <- 7:8 # months to examine
SITE <- "HPL"
#ylims <- c(0,0.7)      # [ppm]
#ylims <- c(0,0.7)*1000  # [ppb]
ylims <- c(0,0.5)*1000  # [ppb]
##########
MONsel.txt <- formatC(MONsel,width=2,flag="0")
# combine data from different years
dat.all<-NULL
for(i in 1:length(YEARs)){
  objname<-paste0("SimCity_CH4_allsites_hrly_",YEARs[i],".rds")
  print(paste("Reading in.....",objname))
  tmp <- readRDS(paste0(obsdir,"/",objname))
  sel <- paste0("CH4_",SITE)%in%colnames(tmp)
  if(sum(sel)==length(SITE)){
    tmp2 <- tmp[,c("Time",paste0("CH4_",c(SITE,"FRU")))]
  } else {
    tmp2 <- data.frame(Time=tmp[,"Time"],tmp[,paste0("CH4_",c(SITE,"FRU"))])
    tmp2 <- data.frame(tmp2,NA)
    colnames(tmp2)[ncol(tmp2)] <- paste0("CH4_",SITE[!sel])
  } # if(paste0("CH4_",SITE)%in%colnames(tmp)){
  dat.all<-rbind(dat.all,tmp2)
  gc()
} #for(i in 1:length(YEARs)){

sitecolnm<-paste0("CH4_",SITE)
YYYYMMDDHH <- format(dat.all[,"Time"],format="%Y%m%d%H")
dat <- dat.all[substring(YYYYMMDDHH,5,6)%in%MONsel.txt,]
YYYYMMDDHH <- format(dat[,"Time"],format="%Y%m%d%H")

dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,mar=c(5,5,4,2))
plot(0,0,type="n",xlab="",ylab="CH4 enhancement over FRU [ppb]",ylim=ylims,xlim=c(0,24),axes=FALSE)
title(main=paste("Average Diurnal Cycle:",SITE,"- FRU,\n Months:",paste(MONsel,collapse=",")))
Years <- unique(substring(YYYYMMDDHH,1,4))
for(yy in 1:length(Years)){
  sel <- substring(YYYYMMDDHH,1,4)==Years[yy]
  diur <- tapply(dat[sel,sitecolnm]-dat[sel,"CH4_FRU"],format(dat[sel,"Time"],format="%H"),mean,na.rm=T)
  diur <- 1000*diur # [ppm] => [ppb]
  UThrs <- as.numeric(names(diur))
  MST <- UThrs - 7; MST[MST<0] <- MST[MST<0] + 24
  #  rearrange order according to MST
  diur <- diur[order(MST)]
  MST <- MST[order(MST)]
  UThrs <- as.numeric(names(diur))
  lines(MST,diur,pch=16,type="l",col=COLS[Years[yy]],lwd=2)
} #for(yy in 1:length(Years)){
box();axis(2)
axis(1,at=MST,labels=MST,line=0)
mtext("Hour of Day [MST]",side=1,line=2.3,cex=1.3)
#axis(1,at=MST,labels=UThrs,line=3.5);mtext("[UT]",side=1,line=4)
#abline(v=range(UTsel)-7,col="darkgray",lwd=2,lty=2)  # selected hours [MST]
rect(xleft=-1,ybottom=ylims[1]-100,xright=UTsel[1]-7,ytop=ylims[2]+100,col=gray(level=0.3,alpha=0.5),border="darkgray")
rect(xleft=max(UTsel)+1-7,ybottom=ylims[1]-100,xright=24+1,ytop=ylims[2]+100,col=gray(level=0.3,alpha=0.5),border="darkgray")
legend("bottomleft",Years,bty="n",lwd=2,text.col=COLS,col=COLS,lty=1,bg="white")


#########################################################
#III.  Compare SEASONAL CYCLES @ HPL in different years
##########
##########
dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,mar=c(5,5,4,2))
ylims <- c(0,0.1)*1000  # [ppb]
plot(0,0,type="n",xlab="Month",ylab="CH4 enhancement over FRU [ppb]",ylim=ylims,xlim=range(MONsel))
title(main=paste("Average Seasonal Cycle:",SITE,"- FRU,\n UThrs:",paste(UTsel,collapse=",")))
Years <- unique(substring(YYYYMMDDHH,1,4))
for(yy in 1:length(Years)){
  sel <- substring(YYYYMMDDHH,1,4)==Years[yy]
  sel <- sel&substring(YYYYMMDDHH,9,10)%in%as.character(UTsel)
  seas <- tapply(dat[sel,sitecolnm]-dat[sel,"CH4_FRU"],format(dat[sel,"Time"],format="%m"),mean,na.rm=T)
  seas <- 1000*seas # [ppm] => [ppb]
  Mons <- as.numeric(names(seas))
  lines(Mons,seas,pch=16,type="l",col=COLS[Years[yy]],lwd=2)
} #for(yy in 1:length(Years)){
legend("topleft",Years,bty="n",lwd=2,text.col=COLS,col=COLS,lty=1)




#########################################################
#IV.  Compare DIURNAL CYCLES @ CSP in different years
##########
YEARs <- (2015:2023)[c(2,6,7,8,9)]
COLS <- c("darkblue","lightblue","darkgreen","lightgreen","orange","red","darkred","pink","magenta")[c(2,6,7,8,9)]
names(COLS) <- YEARs
MONsel <- 4:9 # months to examine
#MONsel <- 4:5 # months to examine
SITE <- "CSP"
ylims <- c(0,0.5)*1000  # [ppb]
##########
MONsel.txt <- formatC(MONsel,width=2,flag="0")
# combine data from different years
dat.all<-NULL
for(i in 1:length(YEARs)){
  objname<-paste0("SimCity_CH4_allsites_hrly_",YEARs[i],".rds")
  print(paste("Reading in.....",objname))
  tmp <- readRDS(paste0(obsdir,"/",objname))
  sel <- paste0("CH4_",SITE)%in%colnames(tmp)
  if(sum(sel)==length(SITE)){
    tmp2 <- tmp[,c("Time",paste0("CH4_",c(SITE,"FRU")))]
  } else {
    tmp2 <- data.frame(Time=tmp[,"Time"],tmp[,paste0("CH4_",c(SITE,"FRU"))])
    tmp2 <- data.frame(tmp2,NA)
    colnames(tmp2)[ncol(tmp2)] <- paste0("CH4_",SITE[!sel])
  } # if(paste0("CH4_",SITE)%in%colnames(tmp)){
  dat.all<-rbind(dat.all,tmp2)
  gc()
} #for(i in 1:length(YEARs)){
sitecolnm<-paste0("CH4_",SITE)
YYYYMMDDHH <- format(dat.all[,"Time"],format="%Y%m%d%H")
dat <- dat.all[substring(YYYYMMDDHH,5,6)%in%MONsel.txt,]
YYYYMMDDHH <- format(dat[,"Time"],format="%Y%m%d%H")

dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,mar=c(5,5,4,2))
plot(0,0,type="n",xlab="",ylab="CH4 enhancement over FRU [ppb]",ylim=ylims,xlim=c(0,24),axes=FALSE)
title(main=paste("Average Diurnal Cycle:",SITE,"- FRU,\n Months:",paste(MONsel,collapse=",")))
Years <- unique(substring(YYYYMMDDHH,1,4))
for(yy in 1:length(Years)){
  sel <- substring(YYYYMMDDHH,1,4)==Years[yy]
  diur <- tapply(dat[sel,sitecolnm]-dat[sel,"CH4_FRU"],format(dat[sel,"Time"],format="%H"),mean,na.rm=T)
  diur <- 1000*diur # [ppm] => [ppb]
  UThrs <- as.numeric(names(diur))
  MST <- UThrs - 7; MST[MST<0] <- MST[MST<0] + 24
  #  rearrange order according to MST
  diur <- diur[order(MST)]
  MST <- MST[order(MST)]
  UThrs <- as.numeric(names(diur))
  lines(MST,diur,pch=16,type="l",col=COLS[Years[yy]],lwd=2)
} #for(yy in 1:length(Years)){
box();axis(2)
axis(1,at=MST,labels=MST,line=0)
mtext("Hour of Day [MST]",side=1,line=2.3,cex=1.3)
#axis(1,at=MST,labels=UThrs,line=3.5);mtext("[UT]",side=1,line=4)
#abline(v=range(UTsel)-7,col="darkgray",lwd=2,lty=2)  # selected hours [MST]
rect(xleft=-1,ybottom=ylims[1]-100,xright=UTsel[1]-7,ytop=ylims[2]+100,col=gray(level=0.3,alpha=0.5),border="darkgray")
rect(xleft=max(UTsel)+1-7,ybottom=ylims[1]-100,xright=24+1,ytop=ylims[2]+100,col=gray(level=0.3,alpha=0.5),border="darkgray")
legend("bottomleft",Years,bty="n",lwd=2,text.col=COLS,col=COLS,lty=1,bg="white")


#########################################################
#V.  Compare SEASONAL CYCLES @ CSP in different years
##########
##########
dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,mar=c(5,5,4,2))
ylims <- c(0,0.1)*1000  # [ppb]
plot(0,0,type="n",xlab="Month",ylab="CH4 enhancement over FRU [ppb]",ylim=ylims,xlim=range(MONsel))
title(main=paste("Average Seasonal Cycle:",SITE,"- FRU,\n UThrs:",paste(UTsel,collapse=",")))
Years <- unique(substring(YYYYMMDDHH,1,4))
for(yy in 1:length(Years)){
  sel <- substring(YYYYMMDDHH,1,4)==Years[yy]
  sel <- sel&substring(YYYYMMDDHH,9,10)%in%as.character(UTsel)
  seas <- tapply(dat[sel,sitecolnm]-dat[sel,"CH4_FRU"],format(dat[sel,"Time"],format="%m"),mean,na.rm=T)
  seas <- 1000*seas # [ppm] => [ppb]
  Mons <- as.numeric(names(seas))
  lines(Mons,seas,pch=16,type="l",col=COLS[Years[yy]],lwd=2)
} #for(yy in 1:length(Years)){
legend("topleft",Years,bty="n",lwd=2,text.col=COLS,col=COLS,lty=1)
