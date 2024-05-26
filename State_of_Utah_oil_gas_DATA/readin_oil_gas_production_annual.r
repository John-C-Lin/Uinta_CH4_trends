# Reads in time series of oil and natural gas production at ANNUAL time scales
# V2(210422): leave units of natural gas production as MCF (thousands of cubic feet), to preserve consistency with monthly data (see 'readin_oil_gas_production_monthly.r'), and not convert to CH4 (do that conversion in other scripts); also save the amount of WATER produced 
# June 4th, 2019 by John C. Lin (John.Lin@utah.edu)

#####################
counties<-c("Duchesne","Uintah")
#Years<-1984:2018
Years<-2000:2023
obsYrs<-2015:2023   #yrs when have surface CH4 obs
Mons.in.Last.Year <- 12   #number of months in the last year, when data were downloaded (to properly scale last year's production)
resultname<-"gas.oil.production_annual.rds"
#####################

dat.all<-NULL
for(i in 1:length(counties)){
  county<-counties[i]
  datfile<-paste0("prod_",county,"_all_annual.csv")
  tmp<-read.csv(datfile,stringsAsFactors=FALSE)
  tmp<-tmp[order(tmp$Year),]
  dat.all<-rbind(dat.all,tmp)
} #for(i in 1:length(counties)){
dat.all<-dat.all[dat.all[,"Year"]%in%Years,]
vars.1<-c("Natural.Gas..MCF.","Oil..BBLs.","Water..BBLs.")

#downscale from annual to MONTHLY production values
#dat.all[-nrow(dat.all),vars.1]<-dat.all[-nrow(dat.all),vars.1]/12  
#  since last year has not been through yet, divide by the number of months that have elapsed when data were downloaded
#dat.all[nrow(dat.all),vars.1]<-dat.all[nrow(dat.all),vars.1]/Mons.in.Last.Year

result<-data.frame(Year=Years)
vars<-vars.1
for(j in 1:length(vars)){
  var<-vars[j]
  
  dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)  
  cols<-c("brown","darkgreen")
  for(i in 1:length(counties)){
    ylims<-range(dat.all[,var],na.rm=T)
    county<-counties[i]
    dat<-dat.all[toupper(dat.all[,"County"])==toupper(county),]
    if(i==1){
  	  plot(dat[,c("Year",var)],ylim=ylims,type="l",lwd=2,col=cols[i])
  	  dat.tot<-dat[,c("Year",var)]
    }else{
  	  lines(dat[,c("Year",var)],lwd=2,col=cols[i])
  	  dat.tot[,var]<-dat.tot[,var]+dat[,var]
    } #if(i==1){
  } #for(i in 1:length(counties)){
  title(main="Annual Production Values");legend("topleft",counties,lwd=2,col=cols,text.col=cols)

  #convert units
  if(var=="Natural.Gas..MCF."){
    # a) to compare against EIA data
    #dat.tot[,var]<-dat.tot[,var]*1000   #[MCF per Month] => [CF per MOnth];[MCF] stands for THOUSAND Cubic Feet, with the "M" coming from Roman Numeral
    #dat.tot[,var]<-dat.tot[,var]/1E6    #[CF per Month] => [Million Cubic Feet per Month]
    #ylab<-"Natural Gas [Million Cubic Feet per Month]"
    # b) to compare against Karion et al. [2013] results
    #dat.tot[,var]<-dat.tot[,var]*0.89   #[MCF natural gas per Month] => [MCF CH4 per Month], using volume fraction of CH4 in natural gas of 0.89 [Karion et al. 2013]
    #dat.tot[,var]<-dat.tot[,var]/(24*(365/12))   #[MCF CH4 per Month] => [MCF CH4 per Hour]
    #dat.tot[,var]<-dat.tot[,var]*1000   #[MCF CH4 per Hour] => [ft^3 CH4 per Hour]
    #dat.tot[,var]<-dat.tot[,var]*0.0283 #[ft^3 CH4 per Hour] => [m^3 CH4 per Hour]
    #R.CH4<-8.3143*1000/16.043   #Ideal Gas constant of CH4 [J/K/kg]: (Rg/mCH4), where Rg is universal gas constant and mCH4 is molar mass of CH4
    #rho.CH4<-(101300)/(R.CH4*288.7)     #density of CH4 [kg/m3], using P=101.3kPa and T=288.7K [Karion et al. 2013]
    #dat.tot[,var]<-dat.tot[,var]*rho.CH4 #[m^3 CH4 per Hour] => [kg CH4 per Hour]
    #dat.tot[,var]<-dat.tot[,var]/1000    #[kg CH4 per Hour] => [10^3 kg CH4 per Hour]
    #ylab<-"Natural Gas [10^3 kg CH4/hr]"
    # c) 
  } #if(var=="Natural.Gas..MCF.){
  if(var=="Oil..BBLs."){
    # to compare against EIA data
    #dat.tot[,var]<-dat.tot[,var]/1000   #[Barrels per Month] => [Thousand Barrels per Month]
    #dat.tot[,var]<-dat.tot[,var]/(365/12) #[Thousand Barrels per Month] => [Thousand Barrels per Day]
    #ylab<-"Oil [Thousand Barrels per Day]"
  } #if(var==""Natural.Gas..MCF.""){
  
  dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3)  
  plot(dat.tot[,c("Year",var)],pch=16,type="o",lwd=2,ylab=var,
      main=paste("Production Values\nCounties: ",paste(counties,collapse=",")))
  abline(v=range(obsYrs),lty=2)
  
  
  result<-data.frame(result,dat.tot[,var])
  #colnames(result)[ncol(result)]<-ylab
  colnames(result)[ncol(result)] <- var
  saveRDS(result,resultname);print(paste(resultname,"written out"))  
} #for(j in 1:length(vars)){

  
