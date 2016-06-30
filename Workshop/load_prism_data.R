inst# loading prism data
## assumes latest data downloaded manually from http://www.prism.oregonstate.edu/recent/
## and unzipped
# data stored in subdirectory by type
source("extract_sites.r")

# need to do this to get the projection right
# FIXME: this should be done inside extract.sites()
B = raster("./weather_data/PRISM_ppt_stable_4kmM3_198101_201506_bil/PRISM_ppt_stable_4kmM3_198101_bil.bil")
# plot(B)

# fake sites
# sites = data.frame(x=c(-100,-100.1,-99,-99.3),
#                    y=c(40,40.1,40.2,40.3))
# Pretty damn close sites. I am good. 
# sites = data.frame(lat=c("45d8'0.24\"N","45d8'0.24\"N","44d55'33.54\"N","44d58'39.24\"N"),
#                    long=c("-113d0'47.34\"E","-113d0'47.34\"E","-113d0'47.34\"E","-113d45'3.54\"E"),
#                    row.names=c("sheep.corral","reservoir.creek","McDevitt","Haynes"),
#                    stringsAsFactors=FALSE)
# sites$x = as.numeric(char2dms(sites$long))
# sites$y = as.numeric(char2dms(sites$lat))
sites = read.table("./data/AstragalusSites.txt",header=TRUE,sep="\t")
# only use the ones that represent different sites -- others have identical 
# PRISM results
# FIXME: this projection should be done inside extract.sites()
sites.sp = SpatialPoints(sites[c(1,2,5,7),3:2],CRS(projection(B)))

tmean = extract.sites(sites.sp,variable="tmean",end="201506")
tmax = extract.sites(sites.sp,variable="tmax",end="201506")
tmin = extract.sites(sites.sp,variable="tmin",end="201506")
ppt = extract.sites(sites.sp,variable="ppt",freq="M3",end="201506")

# to access the climate normals, the defaults for dataset, start and end 
# must be changed. 
tmean.norm = extract.sites(sites.sp,variable="tmean",dataset="30yr_normal",
                           start="01",end="12")
tmax.norm = extract.sites(sites.sp,variable="tmax",dataset="30yr_normal",
                           start="01",end="12")
tmin.norm = extract.sites(sites.sp,variable="tmin",dataset="30yr_normal",
                           start="01",end="12")
ppt.norm = extract.sites(sites.sp,variable="ppt",dataset="30yr_normal",
                           start="01",end="12")

# Now subtract the normal from each observation
N = nrow(tmax)
months = 0:(N-1)  %% 12 + 1
tmax.anom <- tmax[,c(1,2,5,7)] - tmax.norm[months,]
# tmax.anom = ts(tmax.anom,start=1981,frequency=12)

growing.season = months %in% 5:9
sitenames = c("Haynes Creek, ID","McDevitt Creek, ID","Reservoir Creek, MT","Sheep Corral Gulch, MT")

old.par <- par(mar=c(5,5,3,1),mfrow=c(2,2))
for (i in 1:4){
  color <- growing.season+1
  color[growing.season] <- color[growing.season] + (tmax.anom[growing.season,i]>0)
  barplot(as.vector(tmax.anom[,i]),ylab=expression(paste("Monthly ",T[max]," Anomaly ",group("[",degree*C,"]"))),
          col=c("white","blue","red")[color],border=NA,
          main=sitenames[i])
  abline(h=0,lty=2)
  axis(1,seq(1980,2015,5),at=seq(-6,length.out=8,by=72)) # center labels on growing season
}
# abline(v=seq(-12,length.out=8,by=72)) # just checking the alignment of the labels


