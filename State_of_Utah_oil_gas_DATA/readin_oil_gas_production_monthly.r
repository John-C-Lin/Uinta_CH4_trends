# Plot time series of oil and natural gas production
# Based on "readin_oil_gas_production_annualV1.R"
# V2(210430): add annually-averaged value
# V3(210521): change annual average to annual mean, and alter colors for gas & oil
# V4(210521): add # of producing wells to plot along with productiion
# June 4th, 2019 by John C. Lin (John.Lin@utah.edu)

#####################
counties <- c("Duchesne","Uintah")
#Years <- 1984:2018
#Years <- 2000:2020
#Years <- 2005:2020
Years <- 2005:2021
obsYrs <- 2015:2021   #yrs when have surface CH4 obs
MONsel <- 4:9 # months to examine

resultname <- "gas.oil.production_monthly.rds"

col.gas <- "#B91C1C"  # red
col.oil <- "#1D4ED8"  # blue

plotwellnumTF <- TRUE # plot number of producing wells?
welldatname <- "well_data_gas.oil.production_monthly.rds"
#####################

dat.all<-NULL
for(i in 1:length(counties)){
  county<-counties[i]
  xfiles<-list.files(pattern=paste0(county,"_monthly"))
  dat.all.sub<-NULL
for(ff in 1:length(xfiles)){
  #datfile<-"prod_Uintah_monthly_2011_2013.csv"
  datfile<-xfiles[ff]
  print(paste("Reading in:",datfile))
  tmp<-read.csv(datfile,stringsAsFactors=FALSE)
  #  average over monthly timescales
  Natural.Gas..MCF.<-tapply(tmp$Mcf.Gas,tmp$Report.Period,sum,na.rm=TRUE)
  Time<-as.POSIXct(x=names(Natural.Gas..MCF.),format="%m/%d/%Y",tz="MST")
  Natural.Gas..MCF.<-Natural.Gas..MCF.[order(Time)]
  Time<-Time[order(Time)]
  Oil..BBLs.<-tapply(tmp$Bbls.Oil,tmp$Report.Period,sum,na.rm=TRUE)
  Time<-as.POSIXct(x=names(Oil..BBLs.),format="%m/%d/%Y",tz="MST")
  Oil..BBLs.<-Oil..BBLs.[order(Time)]
  Time<-Time[order(Time)]
  dat<-data.frame(Time,Natural.Gas..MCF.,Oil..BBLs.,County=toupper(county))
  dat.all.sub<-rbind(dat.all.sub,dat)
} # for(ff in 1:length(xfiles)){
  dat.all.sub<-dat.all.sub[order(dat.all.sub$Time),]
  dat.all<-rbind(dat.all,dat.all.sub)
} # for(i in 1:length(counties)){

Year<-format(dat.all[,"Time"],"%Y")
dat.all<-data.frame(Year,dat.all)
#!!!!!(200920) Aug 2020 data seems incomplete still, so remove !!!!!
#sel<-format(dat.all$Time,"%Y-%m-%d")=="2020-08-01"
#dat.all<-dat.all[!sel,]
#!!!!!(200920) Aug 2020 data seems incomplete still, so remove !!!!!

#Year 2010 data in DUCHESNE were weird (numbers, so remove)
# 04/01/2010  2010 2010-04-01           2801989     850124 DUCHESNE
# 04/01/20101 2010 2010-04-01           2801989     850124 DUCHESNE
sel <- dat.all$County=="DUCHESNE" & dat.all$Year==2010
sel <- sel & nchar(row.names(dat.all))==11
dat.all <- dat.all[!sel,]
saveRDS(dat.all,file=resultname)
print(paste(resultname,"written out"))

dat.all<-dat.all[dat.all[,"Year"]%in%Years,]
vars.1<-c("Natural.Gas..MCF.","Oil..BBLs.")

#convert into MONTHLY production values
#dat.all[-nrow(dat.all),vars.1]<-dat.all[-nrow(dat.all),vars.1]/12  
#  since last year has not been through yet, divide by the number of months that have elapsed when data were downloaded
#dat.all[nrow(dat.all),vars.1]<-dat.all[nrow(dat.all),vars.1]/Mons.in.Last.Year

if(plotwellnumTF) {welldat.all <- readRDS(welldatname); welldat.all <- welldat.all[substring(welldat.all$YYYYMMDD,1,4)%in%Years,]}

#result<-data.frame(Year=Years)
vars<-vars.1
for(j in 1:length(vars)){
  var<-vars[j]
  # first plot TOTAL (sum of both counties)
  sel <- dat.all$County%in%toupper(counties)
  Tot <- tapply(dat.all[sel,var],dat.all$Time[sel],sum)   # generate monthly totals
  ylims <- range(c(Tot,dat.all[sel,var]),na.rm=T)
  Time <- as.POSIXct(names(Tot),format="%Y-%m-%d",tz="MST")
  #  calculate annual mean
  #tot.yr <- tapply(Tot,substring(names(Tot),1,4),sum)     # generate annual totals
  tot.yr <- tapply(Tot,substring(names(Tot),1,4),mean)     # generate annual totals
  Year.x <- as.POSIXct(paste0(names(tot.yr),"-01-01"),tz="MST")
  tt <- paste0(names(tot.yr),"-07-01")   # assign annual total to middle of year
  Tot.yr <- Tot; Tot.yr[1:length(Tot)] <- NA
  Tot.yr[tt] <- tot.yr   # only assign value to middle of year (otherwise =NA)


  #  plot whole time series 
  if(substring(var,1,3)=="Nat"){xmain <- "Natural Gas Production in Uinta Basin"; ylab <- "Production [Mcf]"}
  if(substring(var,1,3)=="Oil"){xmain <- "Oil Production in Uinta Basin"; ylab <- "Production [Barrels]"}
  dev.new(width=12,height=10);par(cex.axis=1.5,cex.lab=1.5,cex.main=1.5)  
  par(mar=c(5,5,4,5))
  plot(Time,Tot.yr,col="black",pch=16,xlab="Year",cex=1.5,ylab=ylab,ylim=ylims)
  lines(Time[!is.na(Tot.yr)],Tot.yr[!is.na(Tot.yr)],lwd=3)
  axis(side=1,at=Year.x,labels=FALSE)
  lines(Time,Tot,lwd=1,col="black")

  title(main=xmain)
  cols <- c(col.oil,col.gas)
  legend("topleft",c("PRODUCTION:","Total (annual mean)","Total",
                     paste(counties,"county")),lwd=c(NA,3,1,1,1),bty="n",pch=c(NA,16,NA,NA,NA),
                     col=c(NA,"black","black",cols),text.col=c(NA,"black","black",cols),cex=1.3)
  # then plot by county
  for(i in 1:length(counties)){
    county<-counties[i]
    dat<-dat.all[toupper(dat.all[,"County"])==toupper(county),]
    lines(dat[,c("Time",var)],lwd=2,col=cols[i])
  } #for(i in 1:length(counties)){

if(plotwellnumTF){
  if(substring(var,1,3)=="Nat"){sel <- welldat.all$Well.Status.x == "Producing" & welldat.all$Well.Type.x == "Gas Well"}
  if(substring(var,1,3)=="Oil"){sel <- welldat.all$Well.Status.x == "Producing" & welldat.all$Well.Type.x == "Oil Well"}
  welldat <- welldat.all[sel,]

  #check whether there's duplicated oil wells in 2010:  looks weird
  #JCL(210521): duplicated rows removed in 'merge_well_data_production.r'
  timelab <- paste(welldat$Report.Period,welldat$API.Well.Number)
  sel <- duplicated(timelab)
  print(paste("Duplicated data in each year for:",var))
  print(tapply(substring(timelab[sel],7,10),substring(timelab[sel],7,10),length))
  #write.csv(welldat[timelab%in%timelab[sel],],"welldat_duplicated.csv",row.names=FALSE)

  tmp <- tapply(welldat$YYYYMMDD,welldat$YYYYMMDD,length)
  #Time.Nwell <- Time
  Time.Nwell <- as.POSIXct(names(tmp),format="%Y%m%d",tz="MST")
  YYYYMMDD <- format(Time.Nwell,"%Y%m%d")
  Nwell <- tmp[YYYYMMDD]

  # dev.new();plot(Time.Nwell,Nwell,pch=16,type="o",main=var)

  # add well number to existing plot
  par(new=TRUE)
  plot(Time.Nwell,Nwell,type="l",lwd=2,col="darkgreen",axes=F,xlab="",ylab="")
  axis(4,cex=1.5,col="darkgreen",col.axis="darkgreen")
  mtext("Number of Producing Wells",side=4,line=3,cex=1.5,col="darkgreen")
} #  if(plotwellnumTF){


  # dashed lines for years included in study
  xleft <- as.POSIXct(strptime(x=paste0(min(obsYrs),"-",formatC(min(MONsel),flag="0",width=2),"-","01"),
                    format="%Y-%m-%d"),tz="GMT")
  xright <- as.POSIXct(strptime(x=paste0(max(obsYrs),"-",formatC(max(MONsel)+1,flag="0",width=2),"-","01"),
                     format="%Y-%m-%d"),tz="GMT")
  abline(v=c(xleft,xright),lwd=2,lty=2)



  figfilenm <- paste0(var,"_monthly_all.png")
  dev.copy(png,figfilenm,width=540*12/10,height=540);dev.off(); print(paste(figfilenm,"generted"))

if(FALSE){
  #  zoom into years when CH4 obs available 
  dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)  
  # first plot TOTAL (sum of both counties)
  sel <- format(Time,"%Y")%in%obsYrs  
  sel2 <- format(dat.all$Time,"%Y")%in%obsYrs  
  ylims <- range(c(tot[sel],dat.all[sel2,var]),na.rm=T)
  plot(Time[sel],tot[sel],type="l",lwd=2,col="black",ylab=var,ylim=ylims,xlab="Year")
  title(main=var)
  legend("topleft",c("Total",counties),lwd=2,col=c("black",cols),text.col=c("black",cols))
  # then plot by county
  cols<-c("brown","darkgreen")
  for(i in 1:length(counties)){
    county<-counties[i]
    dat<-dat.all[toupper(dat.all[,"County"])==toupper(county),]
    lines(dat[,c("Time",var)],lwd=2,col=cols[i])
  } #for(i in 1:length(counties)){
  figfilenm <- paste0(var,"_monthly_subset.png")
  dev.copy(png,figfilenm);dev.off(); print(paste(figfilenm,"generted"))
} # if(FALSE){

} #for(j in 1:length(vars)){



if(TRUE){
# Compare production values between monthly and earlier annual values

#convert units of natural gas & oil production
tmp<-readRDS("gas.oil.production_monthly.rds")
NatGas<-tapply(tmp[,"Natural.Gas..MCF."],list(tmp[,"Year"]),sum,na.rm=TRUE)  #sum monthly values over entire year
#Nmonth<-tapply(tmp[,"Natural.Gas..MCF."],list(tmp[,"Year"]),length)/length(counties)  # number of months in an entire year
#NatGas<-NatGas/Nmonth  #[MCF natural gas per Year] => [MCF natural gas per Month]
#NatGas<-NatGas*0.89   #[MCF natural gas per Month] => [MCF CH4 per Month], using volume fraction of CH4 in natural gas of 0.89 [Karion et al. 2013]
#NatGas<-NatGas/(24*(365/12))   #[MCF CH4 per Month] => [MCF CH4 per Hour]
#NatGas<-NatGas*1000   #[MCF CH4 per Hour] => [ft^3 CH4 per Hour]
#NatGas<-NatGas*0.0283 #[ft^3 CH4 per Hour] => [m^3 CH4 per Hour]
#R.CH4<-8.3143*1000/16.043   #Ideal Gas constant of CH4 [J/K/kg]: (Rg/mCH4), where Rg is universal gas constant and mCH4 is molar mass of CH4
#rho.CH4<-(101300)/(R.CH4*288.7)     #density of CH4 [kg/m3], using P=101.3kPa and T=288.7K [Karion et al. 2013]
#NatGas<-NatGas*rho.CH4 #[m^3 CH4 per Hour] => [kg CH4 per Hour]
#NatGas<-NatGas/1000    #[kg CH4 per Hour] => [10^3 kg CH4 per Hour]
#ylab<-"Natural Gas [10^3 kg/hr]"
ylab <- "Natural.Gas..MCF."
dat.tot<-data.frame(Year=as.numeric(names(NatGas)),NatGas)
colnames(dat.tot)[ncol(dat.tot)] <- ylab
Oil<-tapply(tmp[,"Oil..BBLs."],list(tmp[,"Year"]),sum,na.rm=TRUE)  #sum monthly values over entire year
#Nmonth<-tapply(tmp[,"Oil..BBLs."],list(tmp[,"Year"]),length)/length(counties)  # number of months in an entire year
#Oil<-Oil/Nmonth  #[Barrels] => [Barrels per Month]
#Oil<-Oil/1000   #[Barrels per Month] => [Thousand Barrels per Month]
#Oil<-Oil/(365/12) #[Thousand Barrels per Month] => [Thousand Barrels per Day]
#ylab<-"Oil [Thousand Barrels per Day]"
ylab <- "Oil..BBLs."
dat.tot <- data.frame(dat.tot,Oil)
colnames(dat.tot)[ncol(dat.tot)] <- ylab

# read in annual values
dat.orig <- readRDS("gas.oil.production_annual.rds")  #generated by "readin_oil_gas_production_annualV1.R" based on annual production data

dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)  
plot(dat.orig[,c("Year","Natural.Gas..MCF.")],pch=16,type="o",lwd=4)
lines(dat.tot[,c("Year","Natural.Gas..MCF.")],col="orange")
legend("topleft",c("Annual data","Monthly data"),col=c("black","orange"),lwd=c(4,1),text.col=c("black","orange"))

dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)  
plot(dat.orig[,c("Year","Oil..BBLs.")],pch=16,type="o",lwd=4)
lines(dat.tot[,c("Year","Oil..BBLs.")],col="orange")
legend("topleft",c("Annual data","Monthly data"),col=c("black","orange"),lwd=c(4,1),text.col=c("black","orange"))

} #if(FALSE){
