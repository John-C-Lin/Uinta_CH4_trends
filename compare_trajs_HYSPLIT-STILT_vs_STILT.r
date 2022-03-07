#Compare trajectories from newly merged HYSPLIT-STILT executable versus STILT executable
#Based on "plot_foot_CH4_inversion_UintahV1.r"
#September 6th, 2020 by John C. Lin (John.Lin@utah.edu)

#########################
#SITE<-"INDY"
SITE<-"HPL"
#SITE<-"ROO"
#SITE<-"CSP"
#SITE<-"SLC"

# HYSPLIT-STILT output directory
outputdir<-paste0("CH4_inversion_Uintah_HYSPLIT-STILT/out/",SITE,"/by-id/")
# STILT output directory
outputdir.s<-paste0("CH4_inversion_Uintah_STILT/out/",SITE,"/by-id/")

#outputdir<-paste0("out/",SITE,"/by-id/")
pngfileTF<-FALSE; overwritepngfileTF<-FALSE
GoogleMapTF<-FALSE

#lat/lon range of plot
#XLIMS<-c(-115,-109);YLIMS<-c(37,42.5)
#XLIMS<-c(-113.5,-111.2);YLIMS<-c(40,41.7)
#XLIMS<-c(-112.3,-111.5);YLIMS<-c(40,41.3)
#XLIMS<-c(-112.7,-111.5);YLIMS<-c(40,41.3)
XLIMS<-c(-90,-80);YLIMS<-c(38,42)
#range of footprint plot (log10 limits)
ZLIMS<-c(-7,-1)
#########################

require(fields);require(ncdf4)
require(maps);require(RgoogleMaps)

#if(!SITE%in%c("WBB","BOS","INDY","JES","SLC","PORT","LA","SF"))stop(paste("Not valid site:",SITE))
if(SITE=="HPL"){lati <- 40.1434; long <- -109.4680; zagl <- 4.06}
if(SITE=="ROO"){lati <- 40.2941; long <- -110.0090; zagl <- 4.06}
if(SITE=="CSP"){lati <- 40.0509; long <- -110.0194; zagl <- 4.06}
if(SITE=="WBB"){lati <- 40.766 ; long <- -111.849; zagl <- 35}
if(SITE=="BOS"){lati <- 42.3470; long <- -71.0840; zagl <- 20}
if(SITE=="INDY"){lati <- 39.7833; long <- -86.1652; zagl <- 20}
if(SITE=="JES"){lati <- 39.1723; long <- -76.7765; zagl <- 20}
if(SITE=="SLC"){lati <- 40.6539; long <- -111.8878; zagl <- 20}
if(SITE=="LA"){lati <- 33.8700; long <- -118.2800; zagl <- 20}
if(SITE=="PORT"){lati <- 45.5100; long <- -122.6000; zagl <- 20}
if(SITE=="SF"){lati <- 37.7956; long <- -122.2794; zagl <- 20}

xmn <- long-2; xmx <- long+2
ymn <- lati-2; ymx <- lati+2
xres <- 0.01
yres <- xres

#loop over the YYYYMMDD folders
YYYYMMDDs<-list.files(outputdir,pattern="2018")     #e.g., "201801012000_-109.468_40.1434_4.06"
YYYYMMDDs.s<-list.files(outputdir.s,pattern="2018") #e.g., "2018010120_-109.468_40.1434_4.06"  NOTE:  does NOT have MINUTE in receptor time
YYYYMMDDs<-YYYYMMDDs[substring(YYYYMMDDs,1,10)%in%substring(YYYYMMDDs.s,1,10)]   #make sure that both STILT & HYSPLIT-STILT results are available
#  YYYYMMDDs<-sample(YYYYMMDDs,3)   #randomly sample a few dates
YYYYMMDDs<-YYYYMMDDs[substring(YYYYMMDDs,1,10)%in%c("2018063020","2018010123","2018053121")]
for(i in 1:length(YYYYMMDDs)){
#for(i in 875:length(YYYYMMDDs)){
#for(i in c(1,50)){
#for(i in 13:length(YYYYMMDDs)){
  ncfilenm<-paste0(outputdir,"/",YYYYMMDDs[i],"/",YYYYMMDDs[i],"_foot.nc")
  if(!file.exists(ncfilenm)){print(paste("no footprint output for:",YYYYMMDDs[i]));next}

  #Generate footprint
  footfile<-nc_open(ncfilenm)
  lat<-ncvar_get(footfile,"lat");lon<-ncvar_get(footfile,"lon")
  foot.all<-ncvar_get(footfile,"foot")
  flat<-lat;flon<-lon
  footsum<-apply(foot.all,c(1,2),sum)   #sum over backward time
  foot.log<-log10(footsum);foot.log[footsum==0]<-NA

  #Read in trajectories 
  trajfile<-paste0(outputdir,"/",YYYYMMDDs[i],"/",YYYYMMDDs[i],"_traj.rds")
  if(!file.exists(trajfile)){print(paste("cannot find HYSPLIT-STILT trajectory file",trajfile));next}
  # change from YYYYMMDDHHmm to YYYYMMDDHH for STILT receptor
  trajname.s <- paste0(substring(YYYYMMDDs[i],1,10),substring(YYYYMMDDs[i],13,nchar(YYYYMMDDs[i])))
  trajfile.s<-paste0(outputdir.s,"/",trajname.s,"/",trajname.s,"_traj.rds")
  if(!file.exists(trajfile.s)){print(paste("cannot find STILT trajectory file",trajfile.s));next}
  # HYSPLIT-STILT trajectories
  dat.all<-readRDS(file=trajfile)
  receptor<-dat.all$receptor
  lon0<-receptor$long;lat0<-receptor$lati;agl0<-receptor$zagl
  receptor$run_time<-strftime(receptor$run_time,"%Y-%m-%d %H:%M",tz="UTC")
  dat<-dat.all$particle
  # STILT trajectories
  dat.all.s<-readRDS(file=trajfile.s)
  dat.s<-dat.all.s$particle

  #Check whether PNG file exists or not
  if(pngfileTF){
    pngfile<-paste0(outputdir,"/",YYYYMMDDs[i],"/",YYYYMMDDs[i],".png")
    if(!overwritepngfileTF){if(file.exists(pngfile)){print(paste("PNG file already exists!; skip:",pngfile));next}}
  } #if(pngfileTF){

  lats<-dat$lati;lons<-dat$long
if(!GoogleMapTF){
  dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.1)
  #image.plot(flon,flat,foot.log,xlab="Lon",ylab="Lat")
  #image.plot(flon,flat,foot.log,xlab="Lon",ylab="Lat",col=terrain.colors(10))
  #image.plot(flon,flat,foot.log,xlab="Lon",ylab="Lat",col=snow.colors(10))
  #image.plot(flon,flat,foot.log,xlab="Lon",ylab="Lat",col=larry.colors())
  image.plot(flon,flat,foot.log,xlab="Lon",ylab="Lat",col=rev(gray(0:10/10)))
  #add trajectories
  # a) STILT-HYSPLIT trajectories
  points(lons,lats,pch=16,col="darkblue",cex=0.4)
  points(lons[dat$foot>0],lats[dat$foot>0],pch=16,col="lightblue",cex=0.2)
  # b) STILT trajectories
  points(dat.s$long,dat.s$lati,pch=16,col="darkgreen",cex=0.2)
  points(dat.s$long[dat.s$foot>0],dat.s$lati[dat.s$foot>0],pch=16,col="lightgreen",cex=0.2)
  #add receptor location
  points(lon0,lat0,pch=17,col="red",cex=1.5)
  map("state",add=TRUE)  #;map("world",add=TRUE)
  map.cities(minpop=100000,cex=1.2)
  title(main=paste("Footprint (STILT model) at receptor (red triangle):\n",receptor$run_time,"UTC   ",
                     round(lon0,digits=3),"LON ",round(lat0,digits=3),"LAT ",round(agl0,digits=3),"m-AGL"),cex.main=1.1)
  # legend("topright",c("HYSPLIT-STILT","STILT"),lty=1,pch=16,lwd=2,col=c("black","darkgray"))
  legend("topright",c("HYSPLIT-STILT","HYSPLIT-STILT(foot>0)","STILT","STILT(foot>0)"),
         lty=1,pch=16,lwd=2,col=c("darkblue","lightblue","darkgreen","lightgreen"),text.col=c("darkblue","lightblue","darkgreen","lightgreen"))
}else{
  #map defined by bounding box
  #bb <- qbbox(lat=lats, lon=lons)   #use trajectories to define map
  #bb <- qbbox(lat=flat, lon=flon)    #use footprint to define map
  bb <- qbbox(lat=YLIMS, lon=XLIMS)    #use footprint to define map
  MyMap <- GetMap.bbox(bb$lonR, bb$latR, maptype="terrain", frame=TRUE,GRAYSCALE=FALSE)
  PlotOnStaticMap(MyMap,mar=c(2,1,3,1))

  #Add trajectories
  site.xy<- LatLon2XY.centered(MyMap, lat=as.numeric(lats), lon=as.numeric(lons))
  points(site.xy$newX,site.xy$newY,pch=16,col="darkgray",cex=0.2)

  #Add receptor location
  site.xy<- LatLon2XY.centered(MyMap, lat=as.numeric(lat0), lon=as.numeric(lon0))
  points(site.xy$newX,site.xy$newY,pch=17,col="black",cex=1.5)

  title(main=paste("CO2-USA Site Footprint (STILT model) at receptor (black triangle):\n",receptor$run_time,"UTC   ",
                     round(lon0,digits=3),"LON ",round(lat0,digits=3),"LAT ",round(agl0,digits=3),"m-AGL"),cex.main=1.1)

  #Add footprint
  alpha<-0.45  #transparency
  lats.m<-matrix(rep(flat,length(flon)),byrow=T,ncol=length(flat))
  lons.m<-matrix(rep(flon,length(flat)),ncol=length(flat))
  xys<- LatLon2XY.centered(MyMap, lat=as.vector(lats.m), lon=as.vector(lons.m)) #converting lat/lon coordinates to map coordinates
  image.coords.x<-matrix(xys$newX,ncol=length(flat))
  image.coords.y<-matrix(xys$newY,ncol=length(flat))
  COLS<-c("white","darkgray","lightblue","RoyalBlue","darkgreen","yellow","orange","red") #colorscale
  colsc<-designer.colors(n=64,alpha=alpha,col=COLS)
  colsc.legend<-designer.colors(n=64,alpha=1.0,col=COLS)

  #generate image
  image(x=image.coords.x[,1],y=image.coords.y[1,],z=foot.log,
          col=colsc,ylim=c(MyMap$BBOX$ll[1],MyMap$BBOX$ur[1]),xlim=c(MyMap$BBOX$ll[2],MyMap$BBOX$ur[2]),add=T,zlim=ZLIMS)

  #Add legend 
  image.plot(x=image.coords.x[,1],y=image.coords.y[1,],z=foot.log,
        col=colsc.legend,ylim=c(MyMap$BBOX$ll[1],MyMap$BBOX$ur[1]),xlim=c(MyMap$BBOX$ll[2],MyMap$BBOX$ur[2]),add=T,
       legend.lab="log10(footprint)",legend.line=2,zlim=ZLIMS,legend.only=T,horizontal=T,legend.width=0.8)
} #if(!GoogleMapTF){

  if(pngfileTF){
    dev.copy(png,filename=pngfile);dev.off()
    print(paste(pngfile,"produced"))
  } #if(pngfileTF){
  nc_close(footfile)
  gc()
} #for(i in 1:length(YYYYMMDDs)){


