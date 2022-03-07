# Tally the natural gas production from oil versus gas wells 
# May 23rd, 2021 by John C. Lin (John.Lin@utah.edu)

#####################
counties <- c("Duchesne","Uintah")
Years <- 2005:2020
obsYrs <- 2015:2020   #yrs when have surface CH4 obs
MONsel <- 4:9 # months to examine

col.gas <- "#B91C1C"  # red
col.oil <- "#1D4ED8"  # blue

plotwellnumTF <- TRUE # plot number of producing wells?
welldatname <- "well_data_gas.oil.production_monthly.rds"
#####################

welldat.all <- readRDS(welldatname)
welldat <- welldat.all[substring(welldat.all$YYYYMMDD,1,4)%in%Years,]

GasProd <- tapply(welldat$Mcf.Gas,list(welldat$YYYYMMDD,welldat$Well.Type.x),sum,na.rm=TRUE)
frYr <- as.numeric(substring(rownames(GasProd),1,4))+(as.numeric(substring(rownames(GasProd),5,6))-0.5)/12

ylims <- range(c(GasProd[,"Gas Well"],GasProd[,"Oil Well"]))
dev.new();par(cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
plot(frYr,GasProd[,"Gas Well"],ylim=ylims,col=col.gas,type="l",lwd=2,
     ylab="Gas Production [Mcf]",xlab="Year")
lines(frYr,GasProd[,"Oil Well"],lwd=2,col=col.oil)
axis(1,at=c(Years,2021),labels=FALSE,cex.axis=1.5)
title(main="Natural Gas Production from both Gas & Oil Wells\nin Uinta Basin")
legend("topright",c("Gas Well","Oil Well"),col=c(col.gas,col.oil),
       text.col=c(col.gas,col.oil),lwd=2,cex=1.5,bty="n")

# calculate fraction of natural gas production produced from gas versus oil wells during study period
sel <- frYr>=min(obsYrs)&frYr<=(max(obsYrs)+1)
xsum <- apply(GasProd[sel,],2,sum,na.rm=T)
print(xsum)
print(round((xsum*100/sum(xsum)),2))

# assume leak rate determined from dCH4 vs foot-convolved gas production slopes (from "dCH4_vs_foot_wellinfo.r"), and see % of emissinos from oil versus gas wells
Egaswell <- 0.07*xsum["Gas Well"]
Eoilwell <- 0.15*xsum["Oil Well"]
print(paste0(round(Eoilwell*100/(Egaswell+Eoilwell),2),"[%] of CH4 emissions may be from oil wells"))

# calculate resulting leak rate from BOTH oil and gas wells
leakrate <- 100*(Egaswell+Eoilwell)/(xsum["Gas Well"]+xsum["Oil Well"])
print(paste0("leak rate from both gas and oil wells:  ",round(leakrate,2),"[%]"))





