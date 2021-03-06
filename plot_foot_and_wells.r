# Plots footprint of receptors and adds in well locations
# V2(210418): instead of calculating averaged footprints on the fly, reads in output from 'average_foot_monthly.r'
# V3(210420): add AVERAGE footprint map with well locations and site locations (meant for paper)
# March 21st, 2021 by John C. Lin (John.Lin@utah.edu)

#########################
#SITE<-"HPL"
#SITE<-"ROO"
SITE<-"CSP"

# Uintah Basin domain
#XLIMS <- c(-111.05,-108.9);YLIMS <- c(39.4,41.0)  #set lat/lon limits for limits of map;  set to NULL if dynamically determine
XLIMS <- c(-111.55,-108.9);YLIMS <- c(39.4,41.0)  #set lat/lon limits for limits of map;  set to NULL if dynamically determine
#XLIMS <- c(-111.15,-108.7);YLIMS <- c(39.4,41.0)  #set lat/lon limits for limits of map;  set to NULL if dynamically determine
#XLIMS <- c(-111.10,-108.8);YLIMS <- c(39.4,41.0)  #set lat/lon limits for limits of map;  set to NULL if dynamically determine
#range of footprint plot (log10 limits)

HRs <- 20:23  #[UTC]  only analyze afternoon hours, following Foster papers

welldatname <- "State_of_Utah_oil_gas_DATA/well_data_gas.oil.production_monthly.rds" # generated by 'State_of_Utah_oil_gas_DATA/merge_well_data_production.r'

# meteorology driving STILT
mettype <- "HRRR"
#mettype <- "NAM12"
#mettype <- "WRF27"

#footmonavenm <- paste0(SITE,'_foot_monave.rds')   # output from 'average_foot_monthly.r'
footmonavenm <- paste0(SITE,'_foot_monave_',mettype,'.rds')   # output from 'average_foot_monthly.r'

#YYYYMMs <- c("201604","201605","202004","202005","202008")   # selected Years/Months
MONs <- 4:9; YEARs <- 2015:2020
if(SITE=="CSP") YEARs <- c(2016,2019,2020)
YYYYMMs <- paste0(rep(YEARs,each=length(MONs)),rep(formatC(MONs,width=2,flag="0"),length(YEARs)))
YYYYMMs <- YYYYMMs[!YYYYMMs%in%c("201504","201505","201506")]
#log10TF <- TRUE; ZLIMS <- c(-7,-1)     
#log10TF <- FALSE; ZLIMS <- c(0.00001,0.0001)
log10TF <- FALSE; ZLIMS <- c(0.000005,0.0001)
#########################

require(ncdf4);require(fields)
require(maps);require(RgoogleMaps)
if(SITE=="HPL"){LAT <- 40.1434; LON <- -109.4680; zagl <- 4.06}
if(SITE=="ROO"){LAT <- 40.2941; LON <- -110.0090; zagl <- 4.06}
if(SITE=="CSP"){LAT <- 40.0509; LON <- -110.0194; zagl <- 4.06}

SITES <- data.frame(matrix(NA,nrow=4,ncol=3,dimnames=list(c("HPL","ROO","CSP","FRU"),c("LAT","LON","zagl"))))
for(i in 1:nrow(SITES)){
  if(rownames(SITES)[i]=="HPL"){SITES$LAT[i] <- 40.1434; SITES$LON[i] <- -109.4680; SITES$zagl[i] <- 4.06}
  if(rownames(SITES)[i]=="ROO"){SITES$LAT[i] <- 40.2941; SITES$LON[i] <- -110.0090; SITES$zagl[i] <- 4.06}
  if(rownames(SITES)[i]=="CSP"){SITES$LAT[i] <- 40.0509; SITES$LON[i] <- -110.0194; SITES$zagl[i] <- 4.06}
  if(rownames(SITES)[i]=="FRU"){SITES$LAT[i] <- 40.2087; SITES$LON[i] <- -110.8403; SITES$zagl[i] <- 4.06}
} #for(i in nrow(SITES)){
# not show ROO in the overall plot due to contamination from nearby well...
SITES <- SITES[rownames(SITES)!="ROO",]

DAT <- readRDS(welldatname)

foot.ave.ALL <- readRDS(footmonavenm)

#-------------------------------------------------#
# average of ALL of the selected Year/Months #
foot.AVE <- apply(foot.ave.ALL[,,YYYYMMs],c(1,2),mean,na.rm=T)


#-------------------------------------------------#
#               Plot month-by-month               #
#for(i in 0:length(YYYYMMs)){
for(i in 0){
  if(i==0){
    foot <- foot.AVE
    wdat <- DAT
    #  Well.Type.y and Well.Status.y is the type/status of the well at the time when data downloaded (2021)--March 2021
    #  relevant files are "CH4_inversion_Uintah/State_of_Utah_oil_gas_DATA/well_data_<county>_all.csv"
    oil.producing <- wdat[wdat$Well.Type.y=="Oil Well"&wdat$Well.Status.y=="Producing",]
    gas.producing <- wdat[wdat$Well.Type.y=="Gas Well"&wdat$Well.Status.y=="Producing",]
    YYYYMM <- paste("\nfootprint averaged over selected years/mons")
    saveRDS(oil.producing,file="oil.producing.rds")
    saveRDS(gas.producing,file="gas.producing.rds")
    saveRDS(foot,file=paste0(SITE,"_footAVE.rds"))
    system(paste0("tar cvf ~/public_html/Exchange/rds.tar oil.producing.rds gas.producing.rds ",SITE,"_footAVE.rds"))
  } else { 
    YYYYMM <- YYYYMMs[i]
    foot <- foot.ave.ALL[,,YYYYMM]
    # prepare well data
    sel <- substring(DAT$YYYYMMDD,1,6) == YYYYMM
    wdat <- DAT[sel,]
    oil.producing <- wdat[wdat$Well.Type.x=="Oil Well"&wdat$Well.Status.x=="Producing",]
    gas.producing <- wdat[wdat$Well.Type.x=="Gas Well"&wdat$Well.Status.x=="Producing",]
  } # } else { 

  flon <- as.numeric(rownames(foot)); flat <- as.numeric(colnames(foot))
  foot[foot==0] <- NA
  foot.log10 <- log10(foot)

  xmain<-paste0(SITE,":  ",YYYYMM)
  xsub<-paste0("UThrs: ",paste(HRs,collapse=","),"; mettype=",mettype)
  #plot on Google Maps
  alpha<-0.70  #transparency
  bb <- qbbox(lat=YLIMS, lon=XLIMS, TYPE="quantile")
  par(pty="s");MyMap <- GetMap.bbox(bb$lonR, bb$latR,urlBase = "http://mt1.google.com/vt/lyrs=m", tileDir= "./mapTiles/Google/")
  PlotOnStaticMap(MyMap,lat=LAT,lon=LON,col="red",pch=17,FUN=points,TrueProj = TRUE,mar=c(2,2,3,5))
  PlotOnStaticMap(MyMap,lat=gas.producing$Latitude,lon=gas.producing$Longitude,
                  col="blue",pch=16,cex=0.4,FUN=points,TrueProj = TRUE,add=T)
  PlotOnStaticMap(MyMap,lat=oil.producing$Latitude,lon=oil.producing$Longitude,
                  col="orange",pch=16,cex=0.4,FUN=points,TrueProj = TRUE,add=T)
  title(main=xmain,cex.main=1.5)
  mtext(xsub,side=1,cex=1.1)

  lats <- flat; lons <- flon
  lons.mat<-matrix(rep(lons,times=length(lats)),ncol=length(lats))
  lats.mat<-matrix(rep(lats,times=length(lons)),byrow=T,nrow=length(lons))
  image.coords <- LatLon2XY.centered(MyMap, lat=as.vector(lats.mat), lon=as.vector(lons.mat))
  xnew.mat<-matrix(image.coords$newX,ncol=ncol(lons.mat),nrow=nrow(lons.mat))
  xnew<-xnew.mat[,1]
  ynew.mat<-matrix(image.coords$newY,ncol=ncol(lats.mat),nrow=nrow(lats.mat))
  ynew<-ynew.mat[1,]

  zmat <- foot
  #legend.lab <- "footprint"
  legend.lab <- expression(paste("footprint [",mu,"mole ",m^-2," ",s^-1,"]",sep=""))
  if(log10TF){zmat <- foot.log10;legend.lab <- "log10(footprint)"}
  #figure out color scheme & zlims
  COLS<-c("white","violet","lightblue","blue","yellow","orange","red") #JCL:  just a subset of colors--not to conflict with trajectory color
  # ZLIMS<-range(zmat,na.rm=T)
  #colsc<-designer.colors(n=64,alpha=alpha,col=COLS); colsc.legend<-designer.colors(n=64,alpha=1.0,col=COLS)
  #colsc <- tim.colors(n=64,alpha=alpha); colsc.legend <- tim.colors(n=64,alpha=1.0)
  colsc <- gray(seq(1,0,-0.02),alpha=alpha); colsc.legend <- gray(seq(1,0,-0.02),alpha=1.0)
  zmat2<-zmat
  zmat2[zmat2>ZLIMS[2]]<-ZLIMS[2]   #make sure everything gets plotted on colorscale
  #zmat2[zmat2<ZLIMS[1]]<-ZLIMS[1]
  zmat2[zmat2<ZLIMS[1]]<-NA   # cut off lower end
  image(x=xnew,y=ynew,z=zmat2,add=TRUE,col=colsc,zlim=ZLIMS)
  image.plot(x=xnew,y=ynew,z=zmat2,legend.only=TRUE,add=TRUE,col=colsc.legend,
            legend.lab=legend.lab,zlim=ZLIMS,legend.cex=1.4,legend.line=-2,legend.mar=6)
  PlotOnStaticMap(MyMap,lat=LAT,lon=LON,col="red",pch=17,FUN=points,TrueProj = TRUE,add=T)
  PlotOnStaticMap(MyMap,lat=LAT+0.05,lon=LON,col="red",FUN=text,TrueProj = TRUE,add=TRUE,labels=SITE,cex=1.3)
  legend("bottomleft",c("gas wells","oil wells"),col=c("blue","orange"),text.col=c("blue","orange"),pch=16,cex=1.4,bg="white")
  if(i==0){
    PlotOnStaticMap(MyMap,lat=SITES$LAT,lon=SITES$LON,col="red",pch=17,FUN=points,TrueProj = TRUE,add=T)
    PlotOnStaticMap(MyMap,lat=SITES$LAT+0.05,lon=SITES$LON,col="red",FUN=text,TrueProj = TRUE,add=TRUE,labels=rownames(SITES),cex=1.3)
    # add ALL sites to map
    figfilenm <- paste0(SITE,"_foot_and_wells","_ALL_AVE.png")
  } else { 
    figfilenm <- paste0(SITE,"_foot_and_wells","_",YYYYMM,".png")
  } # if(i==0){
  dev.copy(png,filename=figfilenm);dev.off();print(paste(figfilenm,"generated"))

} # for(i in 1:length(YYYYMMs)){

 # create movie to facilitate visualization of changes in well activity
  system(paste0("convert -dispose Background -loop 0 -delay 40 ",SITE,"_foot_and_wells_??????.png ~/public_html/Exchange/",SITE,"_foot_and_wells_",mettype,".gif"))




