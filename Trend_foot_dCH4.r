# Look at trends in footprint strength to make sure trends in Fch4 is not due to trend in footprint, as well as trend in dCH4
# Takes output from 'Fch4_simple.r'
# V2(210505):  add SITE, mettype, filterUV.TF
# 3/16/2021 by John C. Lin (John.Lin@utah.edu)

require(ggplot2)

###################
SITE<-"HPL" #HPL is site to focus on for calculating long-term F_CH4
#SITE<-"ROO" 
#SITE<-"CSP" 
MONsub <- 4:9                  # subset of months to examine (filter out wintertime)
HRs<-20:23  #[UTC]  only analyze afternoon hours, following Foster papers

mettype <- "HRRR"
#mettype <- "NAM12"
#mettype <- "WRF27"
filterUV.TF <- TRUE     # filter times based on U/V (filter out times when HRRR is off--i.e., large transport errors) ?
###################


if(filterUV.TF&mettype!="HRRR")stop(paste0("Only extracted HRRR wind vectors; filterUV.TF needs to be set to FALSE for mettype=",mettype))

if(SITE=="HPL"){LAT <- 40.1434; LON <- -109.4680; zagl <- 4.06; WINDsite <- "UBHSP"; WINDsite.2015 <- "A1633"}  # UBHSP site not available in 2015, so use A1633 (RedWash) site instead in 2015
if(SITE=="ROO"){LAT <- 40.2941; LON <- -110.0090; zagl <- 4.06; WINDsite <- "QRS"}
if(SITE=="CSP"){LAT <- 40.0509; LON <- -110.0194; zagl <- 4.06; WINDsite <- "UBCSP"}

datfilenm <- paste("Fch4_",SITE,"_daily_",mettype,".rds",sep="")  # output from 'Fch4_simple.r'
dat.all <- readRDS(datfilenm)
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


dat <- subset(dat.all, as.numeric(substring(dat.all$YYYYMMDD,5,6))%in%MONsub)

#dev.new()
#theme_set(theme_bw())
#g <- ggplot(dat, aes(x=Time,y=foot.basin)) + geom_point(,size=1.5) + geom_smooth(method = "lm") 
#g <- g + theme(plot.title=element_text(size=25), axis.title.y=element_text(size=20),  
#            axis.title.x=element_text(size=20),  axis.text.x=element_text(size=18,angle = 30, vjust=.5),  axis.text.y=element_text(size=18))  # Y axis text
#plot(g)

# Custom function to create regression line from fit (https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/)
ggplotRegression <- function (fit,title='',caption='',ylims=NULL,xlims=NULL) {
  require(ggplot2)
  dev.new()
  theme_set(theme_bw())

  xslope <- signif(fit$coef[[2]], 3)
  xp <- signif(summary(fit)$coef[2,4], 3)
  adjR2 <- signif(summary(fit)$adj.r.squared, 3)
  N <- nobs(fit)
  #xsubtitle <- paste0("slope=",signif(fit$coef[[2]], 3),";  p=",signif(summary(fit)$coef[2,4], 3),";  adjR2=",signif(summary(fit)$adj.r.squared, 3))
  xsubtitle <- substitute(paste("slope=",xslope,";  p=",xp,";  adj-",R^2,"=",adjR2,";  N=",N),list(xslope=xslope,xp=xp,adjR2=adjR2,N=N))

  g <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
      geom_point() +
      stat_smooth(method = "lm", col = "blue") +
      labs(title=title) +
      labs(subtitle = xsubtitle)
      #labs(subtitle = paste0("slope=",signif(fit$coef[[2]], 3),";  p=",signif(summary(fit)$coef[2,4], 3), ";  adjR2=",signif(summary(fit)$adj.r.squared, 3)))
  g <- g + theme(plot.title=element_text(size=18,hjust=0.5), plot.subtitle=element_text(size=17,hjust=0.5,color="blue"), 
                 plot.caption=element_text(size=12),
                 axis.title.y=element_text(size=18), axis.text.y=element_text(size=18),
                 axis.title.x=element_text(size=18), axis.text.x=element_text(size=18,angle = 0, vjust=.5))  
  g <- g + ylim(ylims) + xlim(xlims)
  plot(g)
  return(g)
} # ggplotRegression <- function (fit) {



####################################################
# I.  Trends in Basin-summed footprint             #
####################################################

xtitle <- paste0(SITE,':  Trend in footprint summed over Uinta Basin\nmons: ',paste(MONsub,collapse=","),
                           ';  UThrs: ',paste(HRs,collapse=","), '\n mettype=',mettype)
xcaption <- paste0('filterUV.TF=',filterUV.TF)

#########################
# a) Daily time series
#xfit <- lm(foot.basin ~ Time, data=dat)
#ggplotRegression(xfit)



#########################
# b) Monthly time series
foot.basin.mon <- tapply(dat$foot.basin, dat$YYYYMM, mean, na.rm=T)
Time <- as.POSIXct(paste0(names(foot.basin.mon),"01"), format("%Y%m%d"),tz="GMT")
xfit <- lm(foot.basin.mon ~ Time)
xfit$coef[[2]] <- xfit$coef[[2]]*(365*24*60*60)  # slope is in [1/s]; convert to [1/year]
ylims <- c(0.3,0.9)
xlims <- as.POSIXct(strptime(x=paste0(c(2015,2021),"-01-01"), format="%Y-%m-%d"),tz="GMT")
xbreaks <- as.POSIXct(strptime(x=paste0(2015:2021,"-01-01"), format="%Y-%m-%d"),tz="GMT")
g <- ggplotRegression(xfit,title=xtitle,caption=xcaption,ylim=ylims,xlim=xlims)
g <- g + labs(x="Year",y=expression(paste("Total Footprint in Basin [",mu,"mole m"^-2," s"^-1,"]")))
g <- g + scale_x_continuous(name="Year", breaks=xbreaks, labels=2015:2021)
plot(g)
ggsave(paste0('Trend_foot_month_',mettype,'.png'))



#########################
# c) Annual time series
foot.basin.yr <- tapply(dat$foot.basin, substring(dat$YYYYMM,1,4), mean, na.rm=T)
Year <-  as.numeric(names(foot.basin.yr))
ylims <- c(0.3,0.9)
xlims <- c(2015,2021)
g <- ggplotRegression(lm(foot.basin.yr ~ Year),title=xtitle,caption=xcaption,ylim=ylims,xlim=xlims)
g <- g + labs(x="Year",y=expression(paste("Total Footprint in Basin [",mu,"mole m"^-2," s"^-1,"]")))
g <- g + scale_x_continuous(name="Year", breaks=2015:2021, labels=2015:2021)
plot(g)
ggsave(paste0('Trend_foot_annual_',mettype,'.png'))



####################################################
# II.  Trends in Methane Enhancement (dCH4)        #
####################################################

# for published version, turn off the additional information and caption
#xtitle <- paste0('SITE =',SITE,';  Trend in CH4 enhancement over FRU\nmons: ',paste(MONsub,collapse=","),';  UThrs: ',paste(HRs,collapse=","))
#xcaption <- paste0('filterUV.TF=',filterUV.TF)
xtitle <- substitute(paste(SITE,":  Trend in CH"[4]," enhancement over FRU"),list(SITE=SITE))
xcaption <- paste0('')
#########################
# b) Monthly time series
dCH4 <- 1000*dat$dCH4   # CH4 enhancement [ppb]
dCH4.mon <- tapply(dCH4, dat$YYYYMM, mean, na.rm=T)
Time <- as.POSIXct(paste0(names(dCH4.mon),"01"), format("%Y%m%d"),tz="GMT")
xfit <- lm(dCH4.mon ~ Time)
xfit$coef[[2]] <- xfit$coef[[2]]*(365*24*60*60)  # slope is in [1/s]; convert to [1/year]
ylims <- c(0,100)
xlims <- as.POSIXct(strptime(x=paste0(c(2015,2021),"-01-01"), format="%Y-%m-%d"),tz="GMT")
xbreaks <- as.POSIXct(strptime(x=paste0(2015:2021,"-01-01"), format="%Y-%m-%d"),tz="GMT")
g <- ggplotRegression(xfit,title=xtitle,caption=xcaption,ylim=ylims,xlim=xlims)
g <- g + labs(x="Year",y=expression(paste(Delta,"CH"[4]," [ppb]")))
g <- g + scale_x_continuous(name="Year", breaks=xbreaks, labels=2015:2021)
#g <- g + xlim(xlims)
plot(g)
ggsave(paste0('Trend_dCH4_month_',mettype,'.png'))



#########################
# c) Annual time series
dCH4 <- 1000*dat$dCH4   # CH4 enhancement [ppb]
dCH4.yr <- tapply(dCH4, substring(dat$YYYYMM,1,4), mean, na.rm=T)
Year <-  as.numeric(names(dCH4.yr))
xfit <- lm(dCH4.yr ~ Year)
ylims <- c(0,100)
xlims <- c(2015,2022)
g <- ggplotRegression(xfit,title=xtitle,caption=xcaption,ylim=ylims,xlim=xlims)
g <- g + labs(x="Year",y=expression(paste(Delta,"CH4 [ppb]")))
g <- g + scale_x_continuous(name="Year", breaks=2015:2021, labels=2015:2021)
plot(g)
ggsave(paste0('Trend_dCH4_annual_',mettype,'.png'))



