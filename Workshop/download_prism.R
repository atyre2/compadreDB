# downloading prism datasets directly
library(RCurl) 
# url<-c("ftp://prism.nacse.org/daily/ppt/1981/")
# x<-getURL(url, userpwd="anonymous:atyre2@unl.edu",
#           ftp.use.epsv=FALSE,dirlistonly=TRUE)
# filenames = paste(url, strsplit(x, "\r*\n")[[1]], sep = "")
# 
# url<-c("ftp://prism.nacse.org/monthly/tmean/")
# x<-getURL(url, userpwd="anonymous:atyre2@unl.edu",
#           ftp.use.epsv=FALSE,dirlistonly=TRUE)
# filenames = paste(url, strsplit(x, "\r*\n")[[1]], sep = "")
# 
# 
# data <- getURL(filenames[1],userpwd="anonymous:atyre2@unl.edu",
#                ftp.use.epsv=FALSE)
# 
# data <- getBinaryURL(filenames[1],userpwd="anonymous:atyre2@unl.edu")
# writeBin(data,"weather_data/PRISM_ppt_stable_4kmD1_19810107_bil.zip")

# http://services.nacse.org/prism/data/public/4km/tmin/1944

data <- getBinaryURL("http://services.nacse.org/prism/data/public/4km/tmin/1944")
writeBin(data,"workshop/weather/PRISM_tmin_stable_4kmM2_1944_all_bil.zip")


# OK that works, so know build up a loop of things
url_stem <- "http://services.nacse.org/prism/data/public/4km/"
years <- 1931:1972
variable <- c("tmean","ppt")
year_var <- as.matrix(expand.grid(as.character(years), variable, stringsAsFactors = FALSE))
url_tail <- paste(year_var[,2], year_var[,1], sep = "/")
url_full <- paste0(url_stem, url_tail)
filenames <- paste("workshop/weather/", year_var[,2], year_var[,1], ".bil.zip", sep="_")
for (i in 1:length(url_full)){
  data <- getBinaryURL(url_full[i])
  writeBin(data, filenames[i])
}

