# Plot the gridded output from the convolution between footprint and well data, at monthly timescales generated by "convolve_foot_wellinfo.r"
# April 6th, 2021 by John C. Lin (John.Lin@utah.edu)

#########################
SITE<-"HPL"
#SITE<-"ROO"
#SITE<-"CSP"

# Uintah Basin domain
XLIMS <- c(-111.05,-108.9);YLIMS <- c(39.4,41.0)  #set lat/lon limits for limits of map;  set to NULL if dynamically determine
#XLIMS <- c(-111.15,-108.7);YLIMS <- c(39.4,41.0)  #set lat/lon limits for limits of map;  set to NULL if dynamically determine
#XLIMS <- c(-111.10,-108.8);YLIMS <- c(39.4,41.0)  #set lat/lon limits for limits of map;  set to NULL if dynamically determine

UThrs <- 20:23  #[UTC]  only analyze afternoon hours, following Foster papers

# meteorology driving STILT
mettype <- "HRRR"
#mettype <- "NAM12"
#mettype <- "WRF27"

VAR <- c("Mcf.Gas","Bbls.Oil","Gas.Well","Oil.Well","Gas.Well_Producing","Oil.Well_Producing")[1]
plotwellsTF <- FALSE   #whether or not to plot well locations
welldatname <- "State_of_Utah_oil_gas_DATA/well_data_gas.oil.production_monthly.rds" # generated by 'State_of_Utah_oil_gas_DATA/merge_well_data_production.r'
YYYYMMsel <- NULL  #if set to NULL, then plot ALL the year/months
#YYYYMMsel <- "201605"
#########################

require(fields)
require(maps);require(RgoogleMaps)
require(viridis)
require(geosphere)

if(SITE=="HPL"){LAT <- 40.1434; LON <- -109.4680; zagl <- 4.06}
if(SITE=="ROO"){LAT <- 40.2941; LON <- -110.0090; zagl <- 4.06}
if(SITE=="CSP"){LAT <- 40.0509; LON <- -110.0194; zagl <- 4.06}

if(plotwellsTF)welldat.all <- readRDS(welldatname)

resultnm <- paste0(SITE,"_foot_",VAR,"_gridded_",mettype,".rds")
dat.all <- readRDS(resultnm)
FLON <- as.numeric(rownames(dat.all)); FLAT <- as.numeric(colnames(dat.all))

# determine footprint gridcell areas
xres <- unique(signif(diff(FLON),3)); yres <- unique(signif(diff(FLAT),3))
flon.vec <- rep(FLON,ncol(dat.all))
flat.vec <- rep(FLAT,each=nrow(dat.all))
dx <- distHaversine(p1=cbind(flon.vec,flat.vec),p2=cbind(flon.vec+xres,flat.vec))
dy <- mean(distHaversine(p1=cbind(flon.vec,flat.vec),p2=cbind(flon.vec,flat.vec+yres)))   # gridcell dimension in y-direction is almost a constant
AREA <- dx*dy  # gridcell area [m^2]
AREA.mat <- matrix(AREA,nrow=length(FLON)) #; image.plot(AREA.mat)

YYYYMMs <- YYYYMMsel
if(is.null(YYYYMMsel))YYYYMMs <- dimnames(dat.all)[[3]]

for(i in 1:length(YYYYMMs)){
  YYYYMM <- YYYYMMs[i]
  dum <- unique(as.vector(dat.all[,,YYYYMM]))
  if(length(dum)==1){print(paste("skip blank map for",YYYYMM));next}
  #plot on Google Maps
  #alpha<-0.50  #transparency
  alpha<-1.00  #transparency
  bb <- qbbox(lat=YLIMS, lon=XLIMS, TYPE="quantile")
  par(pty="s");MyMap <- GetMap.bbox(bb$lonR, bb$latR,urlBase = "http://mt1.google.com/vt/lyrs=m", tileDir= "./mapTiles/Google/")
  PlotOnStaticMap(MyMap,lat=LAT,lon=LON,col="black",pch=17,FUN=points,TrueProj = TRUE,mar=c(2,2,3,2),cex=1.5)
  # PlotOnStaticMap(MyMap,lat=LAT+0.05,lon=LON,col="blue",FUN=text,TrueProj = TRUE,add=TRUE,labels=SITE)
  if(plotwellsTF){
    welldat <- welldat.all[substring(welldat.all$YYYYMMDD,1,6)==YYYYMM,]
    welldat <- subset(welldat,welldat[,VAR]>0)
    sel <- welldat$Well.Type.x == "Oil Well"
    PlotOnStaticMap(MyMap,lat=welldat$Latitude[sel],lon=welldat$Longitude[sel],col="darkgray",FUN=points,TrueProj = TRUE,add=T,cex=0.4,pch=16)
    sel <- welldat$Well.Type.x == "Gas Well"
    PlotOnStaticMap(MyMap,lat=welldat$Latitude[sel],lon=welldat$Longitude[sel],col="black",FUN=points,TrueProj = TRUE,add=T,cex=0.4,pch=16)
  } # if(plotwellsTF)

  lats <- FLAT; lons <- FLON
  lons.mat<-matrix(rep(lons,times=length(lats)),ncol=length(lats))
  lats.mat<-matrix(rep(lats,times=length(lons)),byrow=T,nrow=length(lons))
  image.coords <- LatLon2XY.centered(MyMap, lat=as.vector(lats.mat), lon=as.vector(lons.mat))
  xnew.mat<-matrix(image.coords$newX,ncol=ncol(lons.mat),nrow=nrow(lons.mat))
  xnew<-xnew.mat[,1]
  ynew.mat<-matrix(image.coords$newY,ncol=ncol(lats.mat),nrow=nrow(lats.mat))
  ynew<-ynew.mat[1,]

  dat <- dat.all[,,YYYYMM]
  dat <- dat*AREA.mat # remove [m^2] from denominator
  TOT <- sum(dat,na.rm=T)
  dat[dat==0] <- NA 
  zmat <- dat 
  #figure out color scheme & zlims
  #COLS<-c("white","violet","lightblue","blue","yellow","orange","red") #JCL:  just a subset of colors--not to conflict with trajectory color
  #colsc<-designer.colors(n=64,alpha=alpha,col=COLS);colsc.legend<-designer.colors(n=64,alpha=1.0,col=COLS)
  #colsc <- plasma(n=20,alpha=alpha);colsc.legend <- plasma(n=20,alpha=1.0)
  colsc<-tim.colors(n=60,alpha=alpha);colsc.legend<-tim.colors(n=60,alpha=1.0)

  zmat2 <- zmat
  ZLIMS <- range(zmat2,na.rm=T)
  if(VAR=='Mcf.Gas')ZLIMS<-c(0.05,0.5)
  if(VAR=='Bbls.Oil')ZLIMS<-c(0.05,0.5)
  zmat2[zmat2<ZLIMS[1]]<-ZLIMS[1];zmat2[zmat2>ZLIMS[2]]<-ZLIMS[2]   #make sure everything gets plotted on colorscale
  image(x=xnew,y=ynew,z=zmat2,add=TRUE,col=colsc,zlim=ZLIMS)
  image.plot(x=xnew,y=ynew,z=zmat2,legend.only=TRUE,add=TRUE,col=colsc.legend,legend.lab=VAR,zlim=ZLIMS,horizontal=TRUE)
  PlotOnStaticMap(MyMap,lat=LAT,lon=LON,col="black",pch=17,FUN=points,TrueProj = TRUE,add=T,cex=1.5)
  PlotOnStaticMap(MyMap,lat=LAT+0.05,lon=LON,col="black",FUN=text,TrueProj = TRUE,add=TRUE,labels=SITE)
  xmain <- paste0(SITE,"  ",VAR,"  ",YYYYMM,"\nUThrs: ",paste(UThrs,collapse=","),";  total=",signif(TOT,3),";  mettype=",mettype)
  title(main=xmain,cex.main=1.3)
  #mtext(xsub,side=1,cex=1.1)

  if(plotwellsTF)legend("topleft",c("Gas Well","Oil Well"),pch=16,col=c("black","darkgray"))

  figfilenm <- paste0(SITE,"_",VAR,"_foot_wellinfo_convolved_",YYYYMM,"_",mettype,".png")
  dev.copy(png,filename=figfilenm);dev.off();print(paste(figfilenm,"generated"))
  gc()
} # for(i in 1:length(YYYYMMs)){



