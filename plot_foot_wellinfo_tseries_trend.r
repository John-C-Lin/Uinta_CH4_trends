# Plot time series of footprint-convolved well # or footprint-convolved CH4 production and also determine trends
#
# V2(210417): filter out days based on within 45o of observed wind direction (filterUV.TF)
# April 5th, 2021 by John C. Lin (John.Lin@utah.edu)

################
SITE<-"HPL"
#SITE<-"ROO"
#SITE<-"CSP"
# subset of months
MONsub<-4:9
# subset of hours [UTC]; if NULL, then all hours
UThrs <- c(20,21,22,23)  # [UTC]  only analyze afternoon hours, following Foster papers
#datname <- paste0(SITE,'_CH4_foot.wellinfo_ALL.rds') # output from 'merge_CH4_well_info.r'
datname <- paste0(SITE,'_CH4_UVwind_foot.wellinfo_ALL.rds') # output from 'merge_CH4_UVwind_wellinfo.r'
filterUV.TF <- TRUE     # filter times based on U/V (filter out times when HRRR is off--i.e., large transport errors) ?
################

require(ggplot2); require(reshape2); require(geosphere)

if(SITE=="HPL"){LAT <- 40.1434; LON <- -109.4680; zagl <- 4.06; WINDsite <- "UBHSP"; WINDsite.2015 <- "A1633"}  # UBHSP site not available in 2015, so use A1633 (RedWash) site instead in 2015
if(SITE=="ROO"){LAT <- 40.2941; LON <- -110.0090; zagl <- 4.06; WINDsite <- "QRS"}
if(SITE=="CSP"){LAT <- 40.0509; LON <- -110.0194; zagl <- 4.06; WINDsite <- "UBCSP"}

DAT.all <- readRDS(datname)

# filter subset of dataset based on Months & UThrs
YYYYMMDDHH <- format(DAT.all$Time,'%Y%m%d%H')
DAT.all <- data.frame(YYYYMMDDHH, DAT.all)
rownames(DAT.all) <- YYYYMMDDHH
sel <- substring(DAT.all$YYYYMMDDHH, 5, 6)%in%formatC(MONsub,width=2,flag='0')
sel <- sel&(substring(DAT.all$YYYYMMDDHH, 9, 10)%in%formatC(UThrs,width=2,flag='0'))
DAT.all <- DAT.all[sel,]
# filter out times with missing obs at SITE
sel <- !is.na(DAT.all[,paste0("CH4_",SITE)])
DAT.all <- DAT.all[sel,]

if(filterUV.TF){
  dat <- DAT.all
  Usim <- dat[,paste0("Usim.",WINDsite)]; Vsim <- dat[,paste0("Vsim.",WINDsite)]
  Uobs <- dat[,paste0("Uobs.",WINDsite)]; Vobs <- dat[,paste0("Vobs.",WINDsite)]
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
  # assign NAs to wellinfo when failed to pass U/V filter
  print(paste("Before filtering, Mcf.Gas has",sum(!is.na(dat$Mcf.Gas)),"non-NA values"))
  dat[SEL&!isNA,COLNMS] <- NA
  print(paste(" After filtering, Mcf.Gas has",sum(!is.na(dat$Mcf.Gas)),"non-NA values"))
  DAT.all <- dat
} # if(filterUV.TF){ 

# multiply by AREA of footprint grid to remove [1/m^2] factor and to make the slope values larger 
# footprint grid settings from 'run_stilt_hysplit-stiltV1.r'
xmn <- -112.5; xmx <- -108
ymn <- 38.5; ymx <- 41.5
xres <- 0.01; yres <- xres
flat <- seq(ymn+yres/2,ymx-yres/2,yres)
flon <- seq(xmn+xres/2,xmx-xres/2,xres)
# lat/lons refer to footprint grid CENTERS, so move them to lower left (southwest) corner to facilitate gridding
flat.ll <- flat - yres/2; flon.ll <- flon - xres/2
# determine footprint gridcell areas
flon.vec <- rep(flon,length(flat))
flat.vec <- rep(flat,each=length(flon))
dx <- distHaversine(p1=cbind(flon.vec,flat.vec),p2=cbind(flon.vec+xres,flat.vec))
dy <- mean(distHaversine(p1=cbind(flon.vec,flat.vec),p2=cbind(flon.vec,flat.vec+yres)))   # gridcell dimension in y-direction is almost a constant
AREA <- dx*dy  # gridcell area [m^2]
AREA.mat <- matrix(AREA,nrow=length(flon)) #; image.plot(AREA.mat)
#!!!!! average footprint gridcell area is going to be used to remove [1/m^2] factor !!!!!#
AREA.ave <- mean(AREA.mat)

COLNMS <- c('Mcf.Gas','Bbls.Oil','Gas.Well','Oil.Well','Gas.Well_Producing','Oil.Well_Producing')
DAT.all[,COLNMS] <- DAT.all[,COLNMS]*AREA.ave


DAT0 <- DAT.all[,c("Time",COLNMS)]

#----------------------------------------- All Selected Hours ---------------------------------------#
DAT <- DAT0
wellvars <- COLNMS
Time <- DAT$Time
data.table <- NULL
for(i in 1:length(wellvars)){
  Wellinfo <- wellvars[i]
  well <- DAT[,Wellinfo]
  xlm <- lm(well ~ Time)
  m <- coef(xlm)[2]
  p.value <- summary(xlm)$coefficients[,"Pr(>|t|)"]["Time"]
  data.table <- rbind(data.table,data.frame(Wellinfo,m,p.value))
} # for(i in 1:length(wellvars)){
rownames(data.table) <- NULL
# slope ~ [1/s] due to Time being in seconds, so change to [1/year]
data.table$m <- data.table$m*(60*60*24*365)

dev.new(); theme_set(theme_bw())
dat <- melt(DAT, id.vars='Time')
colnames(dat)[colnames(dat)=="variable"] <- "Wellinfo"
g <- ggplot(dat, aes(x=Time,y=value)) + geom_point(,size=1.5) +
     stat_smooth(method = "lm") + geom_abline() + facet_grid(. ~ Wellinfo)
g <- g + labs(title=paste0('SITE =',SITE,';  Time series of sum(foot*wellinfo)',
                           '\n UThrs=',paste(sort(unique(substring(rownames(DAT),9,10))),collapse=","),';  MONs: ',paste(MONsub,collapse=",")),
                            caption=paste0('filterUV.TF=',filterUV.TF))
g <- g + theme(strip.text.x = element_text(size = 14, colour = "black"),plot.caption=element_text(size=12),
               axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=12), axis.text.y=element_text(size=14))
g <- g + geom_text(data = data.table, hjust = 0, col='blue',
                  aes(x = dat$Time[1], y = Inf, label = paste0('\nSlope=',signif(m,3),';  p-val=',signif(p.value,2))))
g + facet_wrap(.~Wellinfo, ncol=2, scales='free')


#----------------------------------------- Average Selected Hours within a Day ---------------------------------------#
YYYYMMDD <- substring(DAT.all$YYYYMMDDHH,1,8)
dat <- DAT.all[,c(COLNMS)]
dum <- NULL
for(cc in 1:ncol(dat)){
  tmp <- tapply(dat[,cc],YYYYMMDD,mean,na.rm=T)
  dum <- cbind(dum,tmp)
} # for(cc in 1:ncol(dat)){
colnames(dum) <- c(COLNMS)
dum[is.nan(dum)] <- NA

Time <- as.POSIXct(rownames(dum),format="%Y%m%d",tz="GMT")
DAT <- data.frame(Time,dum)
wellvars <- COLNMS
data.table <- NULL
for(i in 1:length(wellvars)){
  Wellinfo <- wellvars[i]
  well <- DAT[,Wellinfo]
  xlm <- lm(well ~ Time)
  m <- coef(xlm)[2]
  p.value <- summary(xlm)$coefficients[,"Pr(>|t|)"]["Time"]
  data.table <- rbind(data.table,data.frame(Wellinfo,m,p.value))
} # for(i in 1:length(wellvars)){
rownames(data.table) <- NULL
# slope ~ [1/s] due to Time being in seconds, so change to [1/year]
data.table$m <- data.table$m*(60*60*24*365)

dev.new(); theme_set(theme_bw())
dat <- melt(DAT, id.vars='Time')
colnames(dat)[colnames(dat)=="variable"] <- "Wellinfo"
g <- ggplot(dat, aes(x=Time,y=value)) + geom_point(,size=1.5) +
     stat_smooth(method = "lm") + geom_abline() + facet_grid(. ~ Wellinfo)
g <- g + labs(title=paste0('SITE =',SITE,';  Time series of sum(foot*wellinfo)',
                           '\n mean(UThrs=',paste(UThrs,collapse=","),')',';  MONs: ',paste(MONsub,collapse=",")),
                            caption=paste0('filterUV.TF=',filterUV.TF))
g <- g + theme(strip.text.x = element_text(size = 14, colour = "black"),plot.caption=element_text(size=12),
               axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=12), axis.text.y=element_text(size=14))
g <- g + geom_text(data = data.table, hjust = 0, col='blue',
                  aes(x = dat$Time[1], y = Inf, label = paste0('\nSlope=',signif(m,3),';  p-val=',signif(p.value,2))))
g + facet_wrap(.~Wellinfo, ncol=2, scales='free')
figfilenm <- paste0(SITE,"_foot_wellinfo_tseries_daily.png")
ggsave(figfilenm);print(paste(figfilenm,"generated"))

#----------------------------------------- Average over an entire MONTH ---------------------------------------#
YYYYMM <- substring(DAT.all$YYYYMMDDHH,1,6)
dat <- DAT.all[,c(COLNMS)]
dum <- NULL
for(cc in 1:ncol(dat)){
  tmp <- tapply(dat[,cc],YYYYMM,mean,na.rm=T)
  dum <- cbind(dum,tmp)
} # for(cc in 1:ncol(dat)){
colnames(dum) <- c(COLNMS)
dum[is.nan(dum)] <- NA

YYYYMMDD <- paste0(rownames(dum),"15") # create this, since as.POSIXct doesn't work with just %Y%m, need DAYS too
Time <- as.POSIXct(YYYYMMDD,format="%Y%m%d",tz="GMT")
DAT <- data.frame(Time,dum)
wellvars <- COLNMS
data.table <- NULL
for(i in 1:length(wellvars)){
  Wellinfo <- wellvars[i]
  well <- DAT[,Wellinfo]
  xlm <- lm(well ~ Time)
  m <- coef(xlm)[2]
  p.value <- summary(xlm)$coefficients[,"Pr(>|t|)"]["Time"]
  data.table <- rbind(data.table,data.frame(Wellinfo,m,p.value))
} # for(i in 1:length(wellvars)){
rownames(data.table) <- NULL
# slope ~ [1/s] due to Time being in seconds, so change to [1/year]
data.table$m <- data.table$m*(60*60*24*365)

dev.new(); theme_set(theme_bw())
dat <- melt(DAT, id.vars='Time')
colnames(dat)[colnames(dat)=="variable"] <- "Wellinfo"
g <- ggplot(dat, aes(x=Time,y=value)) + geom_point(,size=1.5) +
     stat_smooth(method = "lm") + geom_abline() + facet_grid(. ~ Wellinfo)
g <- g + labs(title=paste0('SITE =',SITE,';  Time series of sum(foot*wellinfo)',
                           '\n MONTHLY mean;  MONs: ',paste(MONsub,collapse=",")),
                            caption=paste0('filterUV.TF=',filterUV.TF))
g <- g + theme(strip.text.x = element_text(size = 14, colour = "black"),plot.caption=element_text(size=12),
               axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=12), axis.text.y=element_text(size=14))
g <- g + geom_text(data = data.table, hjust = 0, col='blue',
                  aes(x = dat$Time[1], y = Inf, label = paste0('\nSlope=',signif(m,3),';  p-val=',signif(p.value,2))))
g + facet_wrap(.~Wellinfo, ncol=2, scales='free')
figfilenm <- paste0(SITE,"_foot_wellinfo_tseries_monthly.png")
ggsave(figfilenm);print(paste(figfilenm,"generated"))

#----------------------------------------- Average over an entire YEAR ----------------------------------------#
YYYY <- substring(DAT.all$YYYYMMDDHH,1,4)
dat <- DAT.all[,c(COLNMS)]
dum <- NULL
for(cc in 1:ncol(dat)){
  tmp <- tapply(dat[,cc],YYYY,mean,na.rm=T)
  dum <- cbind(dum,tmp)
} # for(cc in 1:ncol(dat)){
colnames(dum) <- c(COLNMS)
dum[is.nan(dum)] <- NA

YYYYMMDD <- paste0(rownames(dum),"0615") # create this, since as.POSIXct doesn't work with just %Y%m, need DAYS too
Time <- as.POSIXct(YYYYMMDD,format="%Y%m%d",tz="GMT")
DAT <- data.frame(Time,dum)
wellvars <- COLNMS
data.table <- NULL
for(i in 1:length(wellvars)){
  Wellinfo <- wellvars[i]
  well <- DAT[,Wellinfo]
  xlm <- lm(well ~ Time)
  m <- coef(xlm)[2]
  p.value <- summary(xlm)$coefficients[,"Pr(>|t|)"]["Time"]
  data.table <- rbind(data.table,data.frame(Wellinfo,m,p.value))
} # for(i in 1:length(wellvars)){
rownames(data.table) <- NULL
# slope ~ [1/s] due to Time being in seconds, so change to [1/year]
data.table$m <- data.table$m*(60*60*24*365)

dev.new(); theme_set(theme_bw())
dat <- melt(DAT, id.vars='Time')
colnames(dat)[colnames(dat)=="variable"] <- "Wellinfo"
g <- ggplot(dat, aes(x=Time,y=value)) + geom_point(,size=1.5) +
     stat_smooth(method = "lm") + geom_abline() + facet_grid(. ~ Wellinfo)
g <- g + labs(title=paste0('SITE =',SITE,';  Time series of sum(foot*wellinfo)',
                           '\n ANNUAL mean;  MONs: ',paste(MONsub,collapse=",")),
                            caption=paste0('filterUV.TF=',filterUV.TF))
g <- g + theme(strip.text.x = element_text(size = 14, colour = "black"),plot.caption=element_text(size=12),
               axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=12), axis.text.y=element_text(size=14))
g <- g + geom_text(data = data.table, hjust = 0, col='blue',
                  aes(x = dat$Time[1], y = Inf, label = paste0('\nSlope=',signif(m,3),';  p-val=',signif(p.value,2))))
g + facet_wrap(.~Wellinfo, ncol=2, scales='free')
figfilenm <- paste0(SITE,"_foot_wellinfo_tseries_yrly.png")
ggsave(figfilenm);print(paste(figfilenm,"generated"))



#----------------------------------------- Separate out Seasons ---------------------------------------#

