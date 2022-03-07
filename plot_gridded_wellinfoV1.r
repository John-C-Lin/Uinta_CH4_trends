# Plot the gridded well info onto Google Map
# March 28th, 2021 by John C. Lin (John.Lin@utah.edu)

################
# output generated by 'grid_well_info_foot.r'
VARS <- c("Mcf.Gas","Bbls.Oil","Gas.Well","Oil.Well","Gas.Well_Producing","Oil.Well_Producing")[c(2,5,6)]

# year/month(s) that want to plot
YYYYMMs <- NULL # if NULL, then ALL of the YYYMMs
#YYYYMMs <- c("201604","202004")
#YYYYMMs <- c("201507","201707","201807","202007")[3]

# Uintah Basin domain
XLIMS <- c(-111.05,-108.9);YLIMS <- c(39.4,41.0)  #set lat/lon limits for limits of map;  set to NULL if dynamically determine
################

require(ncdf4);require(fields)
require(maps);require(RgoogleMaps)
require(viridis);require(geosphere)


SITEs <- c("CSP","HPL","ROO")
LATs <- NULL; LONs <- NULL
for(SITE in SITEs){
  if(SITE=="HPL"){LAT <- 40.1434; LON <- -109.4680; zagl <- 4.06}
  if(SITE=="ROO"){LAT <- 40.2941; LON <- -110.0090; zagl <- 4.06}
  if(SITE=="CSP"){LAT <- 40.0509; LON <- -110.0194; zagl <- 4.06}
  LATs <- c(LATs,LAT); LONs <- c(LONs,LON)
} #for(ss in 1:length(SITEs)){


for(vv in 1:length(VARS)){
  var <- VARS[vv]
  resultname <- paste0(var,"_gridded.rds")
  dat.all <- readRDS(resultname)

  # determine footprint gridcell areas
  FLON <- as.numeric(rownames(dat.all)); FLAT <- as.numeric(colnames(dat.all))
  xres <- unique(signif(diff(FLON),3)); yres <- unique(signif(diff(FLAT),3))
  flon.vec <- rep(FLON,ncol(dat.all))
  flat.vec <- rep(FLAT,each=nrow(dat.all))
  dx <- distHaversine(p1=cbind(flon.vec,flat.vec),p2=cbind(flon.vec+xres,flat.vec))
  dy <- mean(distHaversine(p1=cbind(flon.vec,flat.vec),p2=cbind(flon.vec,flat.vec+yres)))   # gridcell dimension in y-direction is almost a constant
  AREA <- dx*dy  # gridcell area [m^2]
  AREA.mat <- matrix(AREA,nrow=length(FLON)) #; image.plot(AREA.mat)

  YYYYMMsels <- YYYYMMs
  if(is.null(YYYYMMsels))YYYYMMsels <- dimnames(dat.all)[[3]]
for(i in 1:length(YYYYMMsels)){
  YYYYMM <- YYYYMMsels[i]
  dat <- dat.all[,,YYYYMM]
  # remove [m^2] from denominator, resulting in units of [var/day]
  dat <- dat*AREA.mat 

  xmain<-paste0(var,":  ",YYYYMM,"\ntotal=",signif(sum(dat,na.rm=T),3)," [per day]")
  xsub<-""
  #plot on Google Maps
  alpha<-0.50  #transparency
  bb <- qbbox(lat=YLIMS, lon=XLIMS, TYPE="quantile")
  par(pty="s");MyMap <- GetMap.bbox(bb$lonR, bb$latR,urlBase = "http://mt1.google.com/vt/lyrs=m", tileDir= "./mapTiles/Google/")
  PlotOnStaticMap(MyMap,lat=LATs,lon=LONs,col="black",pch=17,FUN=points,TrueProj = TRUE,mar=c(2,2,3,2))
  #PlotOnStaticMap(MyMap,lat=LATs+0.03,lon=LONs,col="black",FUN=text,TrueProj = TRUE,add=TRUE,labels=SITEs)
  title(main=xmain,cex.main=1.5)
  mtext(xsub,side=1,cex=1.1)

  lats <- as.numeric(colnames(dat)); lons <- as.numeric(rownames(dat))
  lons.mat<-matrix(rep(lons,times=length(lats)),ncol=length(lats))
  lats.mat<-matrix(rep(lats,times=length(lons)),byrow=T,nrow=length(lons))
  image.coords <- LatLon2XY.centered(MyMap, lat=as.vector(lats.mat), lon=as.vector(lons.mat))
  xnew.mat<-matrix(image.coords$newX,ncol=ncol(lons.mat),nrow=nrow(lons.mat))
  xnew<-xnew.mat[,1]
  ynew.mat<-matrix(image.coords$newY,ncol=ncol(lats.mat),nrow=nrow(lats.mat))
  ynew<-ynew.mat[1,]

  # set 0s to NA to increase transparency of plot
  zmat <- dat
  zmat[zmat==0] <- NA
  ZLIMS<-range(zmat,na.rm=TRUE)
  if(var=='Mcf.Gas')ZLIMS <- c(100,1E4)
  if(var=='Bbls.Oil')ZLIMS <- c(0,150)
  if(var=='Gas.Well'|var=='Gas.Well_Producing')ZLIMS <- c(0,1.5)  # number of gas wells in each gridcell
  if(var=='Oil.Well'|var=='Oil.Well_Producing')ZLIMS <- c(0,0.5)  # number of oil wells in each gridcell
  #figure out color scheme & zlims
  # COLS<-c("white","violet","lightblue","blue","yellow","orange","red") #JCL:  just a subset of colors--not to conflict with trajectory color
  # colsc<-designer.colors(n=64,alpha=alpha,col=COLS)
  # colsc.legend<-designer.colors(n=64,alpha=1.0,col=COLS)
  #colsc <- tim.colors(n=64,alpha=alpha)
  #colsc.legend <- tim.colors(n=64,alpha=1.0)
  colsc <- plasma(n=20,alpha=alpha)
  colsc.legend <- plasma(n=20,alpha=1.0)

  #make sure everything gets plotted on colorscale
  # zmat[zmat<ZLIMS[1]] <- ZLIMS[1]  #if 0, leave transparent
  zmat[zmat>ZLIMS[2]] <- ZLIMS[2]   
  image(x=xnew,y=ynew,z=zmat,add=TRUE,col=colsc,zlim=ZLIMS)
  image.plot(x=xnew,y=ynew,z=zmat,legend.only=TRUE,add=TRUE,col=colsc.legend,legend.lab=paste(var,"[per day]"),zlim=ZLIMS,horizontal=TRUE)
  PlotOnStaticMap(MyMap,lat=LAT,lon=LON,col="black",pch=16,FUN=points,TrueProj = TRUE,add=T)

  #plot site locations/names
  PlotOnStaticMap(MyMap,lat=LATs,lon=LONs,col="black",pch=17,FUN=points,TrueProj = TRUE,mar=c(2,2,3,2),add=TRUE)
  PlotOnStaticMap(MyMap,lat=LATs+0.03,lon=LONs,col="black",FUN=text,TrueProj = TRUE,add=TRUE,labels=SITEs)

  #plot outliers for OIL production
  #PlotOnStaticMap(MyMap,lat=c(40.255,40.255),lon=c(-110.315,-110.035),col="red",FUN=points,pch=16,TrueProj = TRUE,add=TRUE,cex=1.5)
  
  # figfilenm <- paste0(var,"_",YYYYMM,".jpg"); dev.copy(jpeg,filename=figfilenm)
  figfilenm <- paste0(var,"_",YYYYMM,".png"); dev.copy(png,filename=figfilenm)
  dev.off();print(paste(figfilenm,"generated"))
  gc()
} # for(i in 1:length(YYYYMMs)){
  # create movie to facilitate visualization of changes in well activity
  system(paste0("convert -dispose Background -loop 0 -delay 40 ",VARS[vv],"_??????.png ~/public_html/Exchange/",VARS[vv],".gif"))
  gc()
} # for(vv in 1:length(VARS)){

