# Merges CH4 observations and wind (observed + simulated UV) with results of convolution of footprints with well information
# March 30th, 2021 by John C. Lin (John.Lin@utah.edu)

#################
YEARs <- 2015:2020
#SITE<-"HPL"
#SITE<-"ROO"
SITE<-"CSP"
if(SITE=="CSP")YEARs <- c(2016,2019,2020)
if(SITE=="ROO")YEARs <- 2015:2019

# subset of hours [UTC]; if NULL, then all hours
# UThrs <- NULL   
UThrs <- c(20,21,22,23)  #[UTC]  only analyze afternoon hours, following Foster papers

# meteorology driving STILT
# mettype <- "NAM12"
# mettype <- "WRF27"
mettype <- "HRRR"

obsdir<-"/uufs/chpc.utah.edu/common/home/lin-group4/jcl/SimCity/"  #where CH4 observations reside
winddir<-"/uufs/chpc.utah.edu/common/home/lin-group7/jcl/Transporterr_stiltread/Output"  #where wind obs and sim values are found
finalresultname <- paste0(SITE,'_CH4_UVwind_foot.wellinfo_ALL_',mettype,'.rds')
################

if(SITE=="HPL"){LAT <- 40.1434; LON <- -109.4680; zagl <- 4.06; WINDsite <- "UBHSP"; WINDsite.2015 <- "A1633"}  # UBHSP site not available in 2015, so use A1633 (RedWash) site instead in 2015
if(SITE=="ROO"){LAT <- 40.2941; LON <- -110.0090; zagl <- 4.06; WINDsite <- "QRS"}
if(SITE=="CSP"){LAT <- 40.0509; LON <- -110.0194; zagl <- 4.06; WINDsite <- "UBCSP"}

obs.all <- NULL
wind.all <- NULL
dat.all <- NULL
for(yy in 1:length(YEARs)){
  # read in observations
  objname <- paste0("SimCity_CH4_allsites_hrly_",YEARs[yy])
  print(paste("Reading in.....",objname))
  tmp <- getr(objname,path=obsdir)[,c("Time",paste0("CH4_",c(SITE,"FRU")))]
  obs.all <- rbind(obs.all,tmp)

  # read in sum(foot*wellinfo)--i.e., output from 'convolve_foot_wellinfo.r'
  #datname <- paste0(SITE,'_foot.wellinfo_',YEARs[yy],'.rds')
  datname <- paste0(SITE,'_foot.wellinfo_',YEARs[yy],'_',mettype,'.rds')
  print(paste("Reading in.....",datname))
  tmp.dat <- readRDS(datname)
  dat.all <- rbind(dat.all,tmp.dat)

  # read in HRRR winds + observed values
  windname <- paste0("HRRR_obs_",YEARs[yy],"01010000to",YEARs[yy],"12312300.RDS")
  if(YEARs[yy]=="2018")windname <- paste0("HRRR_obs_",YEARs[yy],"01010000to",YEARs[yy],"12310000.RDS")  # typo during 2018
  colnms <- c("Time",paste0(c("Usim.","Vsim.","Tsim.","Uobs.","Vobs.","Tobs."),WINDsite))
  # UBHSP site not available in 2015, so use A1633 (RedWash) site instead in 2015
  if(SITE=="HPL")colnms <- c("Time",paste0(c("Usim.","Vsim.","Tsim.","Uobs.","Vobs.","Tobs."),WINDsite),
                                    paste0(c("Usim.","Vsim.","Tsim.","Uobs.","Vobs.","Tobs."),WINDsite.2015))
  if(file.exists(paste0(winddir,"/",windname))) {
    print(paste("Reading in.....",windname))
    tmp0 <- readRDS(paste0(winddir,"/",windname))$dat
    tmp <- matrix(NA,nrow=nrow(tmp0),ncol=length(colnms))
    colnames(tmp) <- colnms
    tmp <- data.frame(tmp)
    colnms.sel <- colnames(tmp0)%in%colnms
    tmp[,colnames(tmp0)[colnms.sel]] <- tmp0[,colnms.sel]
  } else {
    print(paste("Doesn't exists:",windname,"; replaces with NAs"))
    tmp <- matrix(NA,nrow=nrow(tmp.dat),ncol=length(colnms))
    colnames(tmp) <- colnms
    tmp <- data.frame(tmp)
    tmp$Time <- tmp.dat$Time
  } # if(file.exists(paste0(winddir,"/",windname)){
  wind.all <- rbind(wind.all,tmp)

  gc()
} # for(yy in 1:length(YEARs)){

# merge observations with sum(foot*wellinfo)
result <- merge(x=obs.all,y=dat.all,by='Time',all.x=FALSE,all.y=TRUE)
result <- merge(x=wind.all,y=result,by='Time',all.x=FALSE,all.y=TRUE)
saveRDS(result,file=finalresultname)
print(paste(finalresultname,'written out'))

