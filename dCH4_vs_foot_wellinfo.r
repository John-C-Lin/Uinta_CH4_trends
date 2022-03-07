# Plot dCH4 as function of footprint-convolved well # or footprint-convolved CH4 production
# y-axis is in [ppm], and x-axis is either 
# a)	[well # x ppm/(umole/m2/s)] or 
# b)	[umole CH4 produced x ppm/(umole/m2/s)]
#
# So the SLOPE from the linear regression has units of 
# a) [(umole/m2s)/(well/m2/day)] or 
# b) [(umole/m2/s)/(umole CH4 produced/m2/day)], which is a direct estimate of leak rate!
#
# This allows one to scale emissions to the Basin level using a) well # or b) CH4 produced
# b) is a direct measure of the overall leak rate!
#
# TO DO:  see whether the slope from b) (leak rate) varies from year to year
# TO DO:  see whether a strong linear relationship that emerges?  
#
# V2(210403): add option to fill in gaps in FRU time series
# V3(210405): add calculation of leak from slope in each plot
# V4(210416): filter out days based on within 45o of observed wind direction (filterUV.TF), use Model-II regression, bootstrap to determine slope errs
# V5(210427): introduce 'mettype' to distinguish between met files driving STILT (to assess transport error)--e.g., "HRRR", "NARR", "NAM12", "WRF27" 
# V6(210523): prepare figures for manuscript submission-e.g., Greek symbols, subscript, trim outliers, report N=?
# March 30th, 2021 by John C. Lin (John.Lin@utah.edu)

################
#SITE<-"HPL"
#SITE<-"ROO"
SITE<-"CSP"
# subset of months
MONsub<-4:9
# subset of hours [UTC]; if NULL, then all hours
UThrs <- c(20,21,22,23)  # [UTC]  only analyze afternoon hours, following Foster papers

# meteorology driving STILT
mettype <- "HRRR"
# mettype <- "NAM12"
# mettype <- "WRF27"

#datname <- paste0(SITE,'_CH4_UVwind_foot.wellinfo_ALL.rds') # output from 'merge_CH4_UVwind_wellinfo.r'
datname <- paste0(SITE,'_CH4_UVwind_foot.wellinfo_ALL_',mettype,'.rds') # output from 'merge_CH4_UVwind_wellinfo.r'
fillFRU.TF <- TRUE      # whether or not to fill in gaps in background (FRU) time series
filterUV.TF <- TRUE     # filter times based on U/V (filter out times when HRRR is off--i.e., large transport errors) ?
CH4.VOLFRAC <- 0.89  # volume fraction of CH4 in natural gas of 0.89 [Karion et al. 2013]

col.gas <- "#B91C1C"  # red
col.oil <- "#1D4ED8"  # blue
################

require(ggplot2); require(reshape2); require(lmodel2)

if(SITE=="HPL"){LAT <- 40.1434; LON <- -109.4680; zagl <- 4.06; WINDsite <- "UBHSP"; WINDsite.2015 <- "A1633"}  # UBHSP site not available in 2015, so use A1633 (RedWash) site instead in 2015
if(SITE=="ROO"){LAT <- 40.2941; LON <- -110.0090; zagl <- 4.06; WINDsite <- "QRS"}
if(SITE=="CSP"){LAT <- 40.0509; LON <- -110.0194; zagl <- 4.06; WINDsite <- "UBCSP"}

DAT.all <- readRDS(datname)

# filter subset of dataset based on Months & UThrs
YYYYMMDDHH <- format(DAT.all$Time,'%Y%m%d%H')
DAT.all <- data.frame(YYYYMMDDHH, DAT.all)
sel <- substring(DAT.all$YYYYMMDDHH, 5, 6)%in%formatC(MONsub,width=2,flag='0')
sel <- sel&(substring(DAT.all$YYYYMMDDHH, 9, 10)%in%formatC(UThrs,width=2,flag='0'))
DAT.all <- DAT.all[sel,]

if(fillFRU.TF){
  dat <- DAT.all
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
  if(FALSE){
    FRU.monave <- tapply(dat$CH4_FRU,dat$YYYYMM,mean,na.rm=T)
    tt <- as.numeric(substring(names(FRU.monave),1,4))+as.numeric(substring(names(FRU.monave),5,6))/12
    dev.new(); plot(tt,FRU.monave,pch=16,type="o")
    dev.new(); plot(dat[,c("Time","CH4_FRU")],pch=16)
    isNA <- is.na(dat$CH4_FRU)
    points(dat$Time[isNA],FRU.gapfilled[isNA],pch=16,col="orange")
  } #if(FALSE){
  dat$CH4_FRU <- FRU.gapfilled
  DAT.all <- dat
} # if(fillFRU.TF){

dCH4 <- DAT.all[,paste0("CH4_",SITE)] - DAT.all[,"CH4_FRU"]
DAT.all <- data.frame(dCH4,DAT.all)
rownames(DAT.all) <- DAT.all$YYYYMMDDHH

if(filterUV.TF){
  dat <- DAT.all
  Usim <- dat[,paste0("Usim.",WINDsite)]; Vsim <- dat[,paste0("Vsim.",WINDsite)]
  Uobs <- dat[,paste0("Uobs.",WINDsite)]; Vobs <- dat[,paste0("Vobs.",WINDsite)]
  if(SITE=="HPL"){
    # UBHSP site not available in 2015, so use A1633 (RedWash) site instead in 2015
    sel <- substring(as.character(dat$Time),1,4)=="2015"
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
  # assign NAs to dCH4 when failed to pass U/V filter
  print(paste("Before filtering, dCH4 has",sum(!is.na(dat$dCH4)),"non-NA values"))
  dat$dCH4[SEL] <- NA
  print(paste(" After filtering, dCH4 has",sum(!is.na(dat$dCH4)),"non-NA values"))
  DAT.all <- dat
} # if(filterUV.TF){ 

COLNMS <- c('Mcf.Gas','Bbls.Oil','Gas.Well','Oil.Well','Gas.Well_Producing','Oil.Well_Producing')[-c(3,4)]
DAT0 <- DAT.all[,c("dCH4",COLNMS)]

slope2leak <- function(dydx,CH4.VOLFRAC=0.89) {
  # calculate leak rate based on regression slope between dCH4 and Mcf.Gas
  dydx <- dydx*(3600*24)  # [umole*day]/[MCFnatgas*s] => [umole]/[MCFnatgas]
  dydx <- dydx/CH4.VOLFRAC       # [umole]/MCFnatgas] => [umole]/[MCF CH4]
  dydx <- dydx/1000       # [umole]/[MCF CH4] => [umole]/[ft^3 CH4]
  dydx <- dydx/0.0283     # [umole]/[ft^3 CH4] => [umole]/[m^3 CH4]
  R.CH4<-8.3143*1000/16.043 #Ideal Gas constant of CH4 [J/K/kg]: (Rg/mCH4), where Rg is universal gas constant & mCH4 is molar mass of CH4
  rho.CH4<-(101300)/(R.CH4*288.7)     #density of CH4 [kg/m3], using P=101.3kPa and T=288.7K [Karion et al. 2013]
  dydx <- dydx/rho.CH4    # [umole]/[m^3 CH4] => [umole]/[kg CH4]
  dydx <- dydx/1000       # [umole]/[kg CH4] => [umole]/[g CH4]
  dydx <- dydx*1E-6       # [umole]/[g CH4] => [mole]/[g CH4]
  dydx <- dydx*16         # [mole]/[g CH4] => [g CH4]/[g CH4]
  dydx <- dydx*100        # [% CH4 emitted per CH4 in Natural Gas production]
  leak <- dydx
  return(leak)
} # slope2leak <- function(dydx) {

bootlmodel2 <- function(x,y,N=1000,method="SMA") {
  if(length(x)!=length(y))stop(paste("x and y need to be the same length!"))
  # bootstrap to calculate errors in regression slope and intercept, based on Model-II regression
  result <- matrix(NA,nrow=N,ncol=2); colnames(result) <- c("m","interc")
  result <- data.frame(result)
  for(i in 1:N){ 
    ind.s <- sample(1:length(x),replace=TRUE)
    x2 <- x[ind.s]
    y2 <- y[ind.s]
    #xlm <- lmodel2(formula=y ~ x)
    xlm <- lmodel2(formula=y2 ~ x2)
    sel <- xlm$regression.results[,"Method"]==method #from Nick Murdoch's test code, appears that the SMA method (standard major axis regression)
    result[i,"m"] <- xlm$regression.results[sel,"Slope"]
    result[i,"interc"] <- xlm$regression.results[sel,"Intercept"]
  } # for(i in 1:N){ 
  return(result)
} # bootlmodel2 <- (x,y,N=1000) {

#----------------------------------------- All Selected Hours ---------------------------------------#
DAT <- DAT0
wellvars <- COLNMS
#  trim outliers in dCH4
sel <- DAT$dCH4>quantile(DAT$dCH4,probs=0.99,na.rm=T)
DAT$dCH4[sel] <- NA
dCH4 <- DAT$dCH4
data.table <- NULL
for(i in 1:length(wellvars)){
  Wellinfo <- wellvars[i]
  #  trim outliers in Wellinfo
  sel <- DAT[,Wellinfo]>quantile(DAT[,Wellinfo],probs=0.99,na.rm=T)
  DAT[sel&!is.na(DAT[,Wellinfo]),Wellinfo] <- NA
  well <- DAT[,Wellinfo]
  R <- cor(dCH4, well, use="na.or.complete")
  # a) OLS regression
  #xlm <- lm(dCH4 ~ well)
  #m <- coef(xlm)[2]; interc <- coef(xlm)[1]
  # b) Model-II regression
  xlm <- lmodel2(formula=dCH4 ~ well)
  sel <- xlm$regression.results[,"Method"]=="SMA" #from Nick Murdoch's test code, appears that the SMA method (standard major axis regression)
  m <- xlm$regression.results[sel,"Slope"]
  interc <- xlm$regression.results[sel,"Intercept"]
  N <- xlm$n
  tmp <- bootlmodel2(x=well,y=dCH4)
  m.sd <- sqrt(var(tmp$m))   # error in regression slope from bootstrap
  data.table <- rbind(data.table,data.frame(Wellinfo,R,m,m.sd,interc,N))
} # for(i in 1:length(wellvars)){
rownames(data.table) <- NULL
data.table$Wellinfo <- factor(data.table$Wellinfo, levels=c('Mcf.Gas', 'Gas.Well_Producing', 'Bbls.Oil', 'Oil.Well_Producing'))
levels(data.table$Wellinfo) <- c("Gas Production","Gas Well Density","Oil Production","Oil Well Density")

dydx <- data.table$m[data.table$Wellinfo=="Gas Production"]  # [ppm]/[MCF natgas day-1 m-2 * (ppm/umole/m2/s)] = [umole*day]/[MCFnatgas*s]
dydx <- slope2leak(dydx,CH4.VOLFRAC=CH4.VOLFRAC)       # [umole*day]/[MCFnatgas*s] => [% CH4 emitted per CH4 in Natural Gas production]
dydx.sd <- data.table$m.sd[data.table$Wellinfo=="Gas Production"]   # [ppm]/[MCF natgas day-1 m-2 * (ppm/umole/m2/s)] = [umole*day]/[MCFnatgas*s]
dydx.sd <- slope2leak(dydx.sd,CH4.VOLFRAC=CH4.VOLFRAC) # [umole*day]/[MCFnatgas*s] => [% CH4 emitted per CH4 in Natural Gas production]

dat <- melt(DAT, id.vars='dCH4')
colnames(dat)[colnames(dat)=="variable"] <- "Wellinfo"
#  adjust order of well info variable
dat$Wellinfo <- factor(dat$Wellinfo, levels=c('Mcf.Gas', 'Gas.Well_Producing', 'Bbls.Oil', 'Oil.Well_Producing'))
levels(dat$Wellinfo) <- c("Gas Production","Gas Well Density","Oil Production","Oil Well Density")

YYYYMMDDHH <- rownames(DAT)
dev.new()
theme_set(theme_bw())
g <- ggplot(dat, aes(x=value,y=dCH4)) + geom_point(,size=1.5) +
     stat_smooth(method = "lm") + geom_abline() + facet_grid(. ~ Wellinfo) 
g <- g + labs(title=paste0('SITE =',SITE,';  dCH4 vs sum(foot*wellinfo);  mons: ',paste(unique(substring(YYYYMMDDHH,5,6)),collapse=","),
                           '\n UThrs=',paste(unique(substring(YYYYMMDDHH,9,10)),collapse=","),
                           '; Yrs=',paste(sort(unique(substring(YYYYMMDDHH,1,4))),collapse=",")),
                           caption=paste0('mettype=',mettype,';  fillFRU.TF=',fillFRU.TF,';  filterUV.TF=',filterUV.TF))
g <- g + labs(subtitle=paste("CH4 leak rate from dCH4 vs GasProd slope:",round(dydx,2),"+/-",round(dydx.sd,2),"[%]"))
g <- g + theme(strip.text.x = element_text(size = 14, colour = "black"),plot.caption=element_text(size=12),
               axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=10), axis.text.y=element_text(size=14))
ymax <- max(dat$dCH4,na.rm=T)
g <- g + ylim(c(0,ymax))
g <- g + geom_text(data = data.table, hjust = 0, col='blue',
                  aes(x = 0.0, y = 0.8*ymax, label = paste0('R=', round(R, 2), '\n',
                      'Slope=', signif(m, 3),'+/-',signif(m.sd,3), '\n','N=',N)))
g + facet_wrap(.~Wellinfo, nrow = 2, scales='free')
ggsave(paste0(SITE,'_dCH4_vs_foot_wellinfo_hrly_',mettype,'.png'))


#----------------------------------------- Average Selected Hours within a Day ---------------------------------------#
YYYYMMDD <- substring(DAT.all$YYYYMMDDHH,1,8)
dat <- DAT.all[,c('dCH4',COLNMS)]
dum <- NULL
for(cc in 1:ncol(dat)){
  tmp <- tapply(dat[,cc],YYYYMMDD,mean,na.rm=T)
  dum <- cbind(dum,tmp)
} # for(cc in 1:ncol(dat)){
colnames(dum) <- c('dCH4',COLNMS)
dum[is.nan(dum)] <- NA

DAT <- data.frame(dum)
wellvars <- COLNMS
#  trim outliers in dCH4
sel <- DAT$dCH4>quantile(DAT$dCH4,probs=0.99,na.rm=T)
DAT$dCH4[sel] <- NA
dCH4 <- DAT$dCH4
data.table <- NULL
for(i in 1:length(wellvars)){
  Wellinfo <- wellvars[i]
  #  trim outliers in Wellinfo
  sel <- DAT[,Wellinfo]>quantile(DAT[,Wellinfo],probs=0.99,na.rm=T)
  DAT[sel&!is.na(DAT[,Wellinfo]),Wellinfo] <- NA
  well <- DAT[,Wellinfo]
  R <- cor(dCH4, well, use="na.or.complete")  # Pearson correlation coefficient
  # Model-II regression
  xlm <- lmodel2(formula=dCH4 ~ well)
  sel <- xlm$regression.results[,"Method"]=="SMA" #from Nick Murdoch's test code, appears that the SMA method (standard major axis regression)
  m <- xlm$regression.results[sel,"Slope"]
  interc <- xlm$regression.results[sel,"Intercept"]
  N <- xlm$n
  tmp <- bootlmodel2(x=well,y=dCH4)
  m.sd <- sqrt(var(tmp$m))   # error in regression slope from bootstrap
  data.table <- rbind(data.table,data.frame(Wellinfo,R,m,m.sd,interc,N))
} # for(i in 1:length(wellvars)){
rownames(data.table) <- NULL
data.table$Wellinfo <- factor(data.table$Wellinfo, levels=c('Mcf.Gas', 'Gas.Well_Producing', 'Bbls.Oil', 'Oil.Well_Producing'))
levels(data.table$Wellinfo) <- c("Gas Production","Gas Well Density","Oil Production","Oil Well Density")

dydx <- data.table$m[data.table$Wellinfo=="Gas Production"]  # [ppm]/[MCF natgas day-1 m-2 * (ppm/umole/m2/s)] = [umole*day]/[MCFnatgas*s]
dydx <- slope2leak(dydx,CH4.VOLFRAC=CH4.VOLFRAC)       # [umole*day]/[MCFnatgas*s] => [% CH4 emitted per CH4 in Natural Gas production]
dydx.sd <- data.table$m.sd[data.table$Wellinfo=="Gas Production"]  # [ppm]/[MCF natgas day-1 m-2 * (ppm/umole/m2/s)] = [umole*day]/[MCFnatgas*s]
dydx.sd <- slope2leak(dydx.sd,CH4.VOLFRAC=CH4.VOLFRAC) # [umole*day]/[MCFnatgas*s] => [% CH4 emitted per CH4 in Natural Gas production]

dat <- melt(DAT, id.vars='dCH4')
colnames(dat)[colnames(dat)=="variable"] <- "Wellinfo"
#  adjust order of well info variable
dat$Wellinfo <- factor(dat$Wellinfo, levels=c('Mcf.Gas', 'Gas.Well_Producing', 'Bbls.Oil', 'Oil.Well_Producing'))
levels(dat$Wellinfo) <- c("Gas Production","Gas Well Density","Oil Production","Oil Well Density")
saveRDS(dat,file="dat.rds");print("dat.rds written out")
YYYYMMDD <- rownames(DAT); isNA <- is.na(DAT$dCH4)

dev.new()
#options(scipen=-1)   # tend to use Scientific Notation
theme_set(theme_bw())
g <- ggplot(dat, aes(x=value,y=dCH4)) + geom_point(,size=1.2)  + facet_wrap(. ~ Wellinfo) +
     stat_smooth(method = "lm") + geom_abline() # + facet_wrap(. ~ Wellinfo) 
#g <- g + labs(title=paste0('SITE =',SITE,';  dCH4 vs sum(foot*wellinfo);  mons: ',paste(unique(substring(YYYYMMDD[!isNA],5,6)),collapse=","),
#                           '\n mean(UThrs=',paste(unique(substring(DAT.all$YYYYMMDDHH,9,10)),collapse=","),
#                           '); Yrs=',paste(sort(unique(substring(YYYYMMDD[!isNA],1,4))),collapse=","),
#                           '\n CH4 leak rate from dCH4 vs GasProd slope:',round(dydx,2),"+/-",round(dydx.sd,2),"[%]"))
MONuniq <- paste(unique(substring(YYYYMMDD[!isNA],5,6)),collapse=",")
HHuniq <- paste(unique(substring(DAT.all$YYYYMMDDHH,9,10)),collapse=",")
YYYYuniq <- paste(sort(unique(substring(YYYYMMDD[!isNA],1,4))),collapse=",")
g <- g + labs(title=substitute(paste(SITE," daily; ",Delta,"CH"[4]," vs GasProd slope(leak rate)=",dydx,"+/-",dydx.sd,"[%]"),
                           list(SITE=SITE,HHuniq=HHuniq,dydx=round(dydx,2),dydx.sd=round(dydx.sd,2))))
g <- g + labs(caption=paste0('mettype=',mettype,';  fillFRU.TF=',fillFRU.TF,';  filterUV.TF=',filterUV.TF))
g <- g + scale_x_continuous(labels = scales::scientific)  # forces use of scientific notatiion on x-axis
g <- g + theme(strip.text.x = element_text(size = 14, colour = "black"), axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=11.0), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),plot.title=element_text(size=15.5))
ymax <- max(dat$dCH4,na.rm=T)
g <- g + ylim(c(0,ymax))
g <- g + geom_text(data = data.table, hjust = 0, col='blue',
                  aes(x = 0.0, y = 0.8*ymax, label = paste0('R=', round(R, 2), '\n',
                      'Slope=', signif(m, 3),'+/-',signif(m.sd,3), '\n','N=',N)))
g <- g + labs(y=expression(paste(Delta,"CH4 [ppm]")))
g + facet_wrap(.~Wellinfo, nrow = 2, scales='free')
ggsave(paste0(SITE,'_dCH4_vs_foot_wellinfo_daily_',mettype,'.png'))


# Multiple plot function
# Copied from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}


# Plot each panel separately, in order to more closely control plotting and to enable different x-axis labels
theme_set(theme_bw())
wellvar <- "Gas Production"
g <- ggplot(subset(dat,Wellinfo==wellvar), aes(x=value,y=dCH4)) + geom_point(,size=1.2)
SEL <- data.table$Wellinfo==wellvar
g <- g + geom_abline(slope=data.table$m[SEL],intercept=data.table$interc[SEL],colour="blue")
g <- g + labs(title=paste0(SITE,":  ",wellvar))
if(SITE=="HPL")g <- g + scale_x_continuous(labels = scales::scientific)  # forces use of scientific notatiion on x-axis
g <- g + theme(strip.text.x = element_text(size = 14, colour = "black"), axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=10.0), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),plot.title=element_text(size=14.0))
ymax <- max(dat$dCH4,na.rm=T)
g <- g + ylim(c(0,ymax))
g <- g + geom_text(data = data.table[SEL,], hjust = 0, col='blue', aes(x = 0.0, y = 0.8*ymax, label = paste0('R=', round(R, 2), '\n', 'Slope=', signif(m, 3),'+/-',signif(m.sd,3), '\n','N=',N)),size=4)
xpos <- 1.2E-4
if(SITE=="CSP")xpos <- 5.0E-5
g <- g + geom_text( col=col.gas, aes(x = xpos, y = ymax, label = paste0("Leakage rate: ",round(dydx,2),"+/-",round(dydx.sd,2)," [%]")),size=4.5)
g1 <- g + labs(y=expression(paste(Delta,"CH"[4]," [ppm]")), x=expression(paste("[Mcf day"^-1," ppm/(",mu,"mole s"^-1,")]")))

theme_set(theme_bw())
wellvar <- "Oil Production"
dat.sub <- dat[dat$Wellinfo==wellvar,]
g <- ggplot(subset(dat,Wellinfo==wellvar), aes(x=value,y=dCH4)) + geom_point(,size=1.2)
SEL <- data.table$Wellinfo==wellvar
g <- g + geom_abline(slope=data.table$m[SEL],intercept=data.table$interc[SEL],colour="blue")
g <- g + labs(title=paste0(SITE,":  ",wellvar))
g <- g + scale_x_continuous(labels = scales::scientific)  # forces use of scientific notatiion on x-axis
g <- g + theme(strip.text.x = element_text(size = 14, colour = "black"), axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=10.5), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),plot.title=element_text(size=14.0))
ymax <- max(dat$dCH4,na.rm=T)
g <- g + ylim(c(0,ymax))
g <- g + geom_text(data = data.table[SEL,], hjust = 0, col='blue', aes(x = 0.0, y = 0.8*ymax, label = paste0('R=', round(R, 2), '\n', 'Slope=', signif(m, 3),'+/-',signif(m.sd,3), '\n','N=',N)),size=4)
g2 <- g + labs(y=expression(paste(Delta,"CH"[4]," [ppm]")), x=expression(paste("[Barrels day"^-1," ppm/(",mu,"mole s"^-1,")]")))

theme_set(theme_bw())
wellvar <- "Gas Well Density"
dat.sub <- dat[dat$Wellinfo==wellvar,]
g <- ggplot(subset(dat,Wellinfo==wellvar), aes(x=value,y=dCH4)) + geom_point(,size=1.2)
SEL <- data.table$Wellinfo==wellvar
g <- g + geom_abline(slope=data.table$m[SEL],intercept=data.table$interc[SEL],colour="blue")
g <- g + labs(title=paste0(SITE,":  ",wellvar))
g <- g + scale_x_continuous(labels = scales::scientific)  # forces use of scientific notatiion on x-axis
g <- g + theme(strip.text.x = element_text(size = 14, colour = "black"), axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=10.5), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),plot.title=element_text(size=14.0))
ymax <- max(dat$dCH4,na.rm=T)
g <- g + ylim(c(0,ymax))
g <- g + geom_text(data = data.table[SEL,], hjust = 0, col='blue', aes(x = 0.0, y = 0.8*ymax, label = paste0('R=', round(R, 2), '\n', 'Slope=', signif(m, 3),'+/-',signif(m.sd,3), '\n','N=',N)),size=4)
g3 <- g + labs(y=expression(paste(Delta,"CH"[4]," [ppm]")), x=expression(paste("[Well # day"^-1," ppm/(",mu,"mole s"^-1,")]")))

theme_set(theme_bw())
wellvar <- "Oil Well Density"
dat.sub <- dat[dat$Wellinfo==wellvar,]
g <- ggplot(subset(dat,Wellinfo==wellvar), aes(x=value,y=dCH4)) + geom_point(,size=1.2)
SEL <- data.table$Wellinfo==wellvar
g <- g + geom_abline(slope=data.table$m[SEL],intercept=data.table$interc[SEL],colour="blue")
g <- g + labs(title=paste0(SITE,":  ",wellvar))
g <- g + scale_x_continuous(labels = scales::scientific)  # forces use of scientific notatiion on x-axis
g <- g + theme(strip.text.x = element_text(size = 14, colour = "black"), axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=10.5), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),plot.title=element_text(size=14.0))
ymax <- max(dat$dCH4,na.rm=T)
g <- g + ylim(c(0,ymax))
g <- g + geom_text(data = data.table[SEL,], hjust = 0, col='blue', aes(x = 0.0, y = 0.8*ymax, label = paste0('R=', round(R, 2), '\n', 'Slope=', signif(m, 3),'+/-',signif(m.sd,3), '\n','N=',N)),size=4)
g4 <- g + labs(y=expression(paste(Delta,"CH"[4]," [ppm]")), x=expression(paste("[Well # day"^-1," ppm/(",mu,"mole s"^-1,")]")))

multiplot(g1,g2,g3,g4,cols=2)
#ggsave(paste0(SITE,'_dCH4_vs_foot_wellinfo_daily_',mettype,'_multiplot.png'))  #ggsave only saves one of the panels, so use dev.copy instead...
#dev.copy(png,paste0(SITE,'_dCH4_vs_foot_wellinfo_daily_',mettype,'_multiplot.png'));dev.off()
dev.copy(pdf,paste0(SITE,'_dCH4_vs_foot_wellinfo_daily_',mettype,'_multiplot.pdf'));dev.off()


#----------------------------------------- Analyze Year-by-Year, and Average Selected Hours within a Day ---------------------------------------#
if(FALSE){
YEARs <- unique(substring(DAT.all$YYYYMMDD,1,4))

for(yy in 1:length(YEARs)){

YEAR <- YEARs[yy]
sel <- substring(DAT.all$YYYYMMDD,1,4)==YEAR
YYYYMMDD <- substring(DAT.all$YYYYMMDDHH[sel],1,8)
dat <- DAT.all[sel,c('dCH4',COLNMS)]
dum <- NULL
for(cc in 1:ncol(dat)){
  tmp <- tapply(dat[,cc],YYYYMMDD,mean,na.rm=T)
  dum <- cbind(dum,tmp)
} # for(cc in 1:ncol(dat)){
colnames(dum) <- c('dCH4',COLNMS)
dum[is.nan(dum)] <- NA

DAT <- data.frame(dum)
wellvars <- COLNMS
#  trim outliers in dCH4
sel <- DAT$dCH4>quantile(DAT$dCH4,probs=0.99,na.rm=T)
DAT$dCH4[sel] <- NA
dCH4 <- DAT$dCH4
data.table <- NULL
for(i in 1:length(wellvars)){
  Wellinfo <- wellvars[i]
  #  trim outliers in Wellinfo
  sel <- DAT[,Wellinfo]>quantile(DAT[,Wellinfo],probs=0.99,na.rm=T)
  DAT[sel&!is.na(DAT[,Wellinfo]),Wellinfo] <- NA
  well <- DAT[,Wellinfo]
  R <- cor(dCH4, well, use="na.or.complete")
  # a) OLS regression
  #xlm <- lm(dCH4 ~ well)
  #m <- coef(xlm)[2]; interc <- coef(xlm)[1]
  # b) Model-II regression
  xlm <- lmodel2(formula=dCH4 ~ well)
  sel <- xlm$regression.results[,"Method"]=="SMA" #from Nick Murdoch's test code, appears that the SMA method (standard major axis regression)
  m <- xlm$regression.results[sel,"Slope"]
  interc <- xlm$regression.results[sel,"Intercept"]
  N <- xlm$n
  tmp <- bootlmodel2(x=well,y=dCH4)
  m.sd <- sqrt(var(tmp$m))   # error in regression slope from bootstrap
  data.table <- rbind(data.table,data.frame(Wellinfo,R,m,m.sd,interc,N))
} # for(i in 1:length(wellvars)){
rownames(data.table) <- NULL
data.table$Wellinfo <- factor(data.table$Wellinfo, levels=c('Mcf.Gas', 'Gas.Well_Producing', 'Bbls.Oil', 'Oil.Well_Producing'))
levels(data.table$Wellinfo) <- c("Gas Production","Gas Well Density","Oil Production","Oil Well Density")

dydx <- data.table$m[data.table$Wellinfo=="Gas Production"]  # [ppm]/[MCF natgas day-1 m-2 * (ppm/umole/m2/s)] = [umole*day]/[MCFnatgas*s]
dydx <- slope2leak(dydx,CH4.VOLFRAC=CH4.VOLFRAC)       # [umole*day]/[MCFnatgas*s] => [% CH4 emitted per CH4 in Natural Gas production]
dydx.sd <- data.table$m.sd[data.table$Wellinfo=="Gas Production"]  # [ppm]/[MCF natgas day-1 m-2 * (ppm/umole/m2/s)] = [umole*day]/[MCFnatgas*s]
dydx.sd <- slope2leak(dydx.sd,CH4.VOLFRAC=CH4.VOLFRAC) # [umole*day]/[MCFnatgas*s] => [% CH4 emitted per CH4 in Natural Gas production]

dat <- melt(DAT, id.vars='dCH4')
colnames(dat)[colnames(dat)=="variable"] <- "Wellinfo"
#  adjust order of well info variable
dat$Wellinfo <- factor(dat$Wellinfo, levels=c('Mcf.Gas', 'Gas.Well_Producing', 'Bbls.Oil', 'Oil.Well_Producing'))
levels(dat$Wellinfo) <- c("Gas Production","Gas Well Density","Oil Production","Oil Well Density")

YYYYMMDD <- rownames(DAT); isNA <- is.na(DAT$dCH4)
dev.new()
theme_set(theme_bw())
g <- ggplot(dat, aes(x=value,y=dCH4)) + geom_point(,size=1.5) +
     stat_smooth(method = "lm") + geom_abline() + facet_grid(. ~ Wellinfo) 
g <- g + labs(title=paste0('SITE =',SITE,';  dCH4 vs sum(foot*wellinfo);  mons: ',paste(unique(substring(YYYYMMDD[!isNA],5,6)),collapse=","),
                           '\n mean(UThrs=',paste(unique(substring(DAT.all$YYYYMMDDHH,9,10)),collapse=","),
                           '); Yrs=',paste(sort(unique(substring(YYYYMMDD[!isNA],1,4))),collapse=",")),
                           caption=paste0('mettype=',mettype,';  fillFRU.TF=',fillFRU.TF,';  filterUV.TF=',filterUV.TF))
g <- g + labs(subtitle=paste("CH4 leak rate from dCH4 vs GasProd slope:",round(dydx,2),"+/-",round(dydx.sd,2),"[%]"))
g <- g + theme(strip.text.x = element_text(size = 14, colour = "black"), axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=10), axis.text.y=element_text(size=14),plot.caption=element_text(size=12))
ymax <- max(dat$dCH4,na.rm=T)
g <- g + ylim(c(0,ymax))
g <- g + geom_text(data = data.table, hjust = 0, col='blue',
                  aes(x = 0.0, y = 0.8*ymax, label = paste0('R=', round(R, 2), '\n',
                      'Slope=', signif(m, 3),'+/-',signif(m.sd,3), '\n','N=',N)))
g + facet_wrap(.~Wellinfo, nrow = 2, scales='free')
ggsave(paste0(SITE,'_dCH4_vs_foot_wellinfo_daily_',YEAR,'_',mettype,'.png'))

} # for(yy in 1:length(YEARs)){

} #if(FALSE){

