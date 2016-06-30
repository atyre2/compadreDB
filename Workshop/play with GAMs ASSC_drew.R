library(mgcv)

setwd("C:/Users/ecrone01/Documents/Cornwall workshop")


# read in temp & precip data
TPdat = read.csv("Workshop/weather_astragalus.csv", header = T)
names(TPdat)[2] = "year"
names(TPdat)[3] = "month"

# Calculate 1 month SPEI
library(SPEI)

# first get potential evapotranspiration
sites <- read.table("Workshop/AstragalusSites.txt", sep = "\t", header=TRUE, stringsAsFactors = FALSE)
meanlat <- mean(sites$lat)

TPdat$PET <- thornthwaite(TPdat$tmean, meanlat)
# probably this nxt step should be done by site, currently scaling using all sites together
SPEI_1 <- with(TPdat, spei(ppt - PET, 1)) # warning! this won't work with scale > 1 -- have to do each site individually
plot(SPEI_1)
TPdat$SPEI_1 <- as.vector(SPEI_1$fitted)

# read in Baltimore checkerspot data
ASSCdat = read.csv("Workshop/Lambdas_Astragalus2.csv")
ASSCdat$year = strtoi(substr(ASSCdat$Transition, start = 6, stop = 10))

# function from Teller 2016 to aggregate data by month
aggLags<-function(datC, fun, meas, seg){
  pars <- as.list(match.call()[-1])
  measure<-datC[,as.character(pars$meas)]
  segment<-datC[,as.character(pars$seg)]
  agg<-with(datC,ftable(tapply(measure, list(datC$year, segment), fun)))
  agg<-as.data.frame(agg)
  names(agg)<-c("year",pars$seg,paste(pars$meas, pars$fun,"Now", sep=""))
  
  return(list(agg,paste(pars$fun)))
}

# function from Teller 2016 to create lagged variables
appLags<-function(agDatC,nUnits,name, fun){
  storeNames<-names(agDatC)
  for (i in 1:nUnits){
    agDatC<-as.data.frame(agDatC)
    waved<-agDatC[[3]]
    newWave<-c(rep(NA, times=i),waved) 
    newThing<-newWave[1:nrow(agDatC)]
    agDatC<-cbind(agDatC, newThing)
  }
  names(agDatC)<-c(storeNames,paste(name,fun,1:nUnits, sep=""))
  return(agDatC)
}


makeDat = function(useGDD, useBCB, monthStart){
  GDDAggWT<-aggLags(datC=useGDD, fun=sum, meas=tmean, seg=month)
  subAgg<-GDDAggWT[[1]]
  subAgg$month<-as.numeric(as.character(subAgg$month))
  subAgg<-subAgg[with(subAgg,order(year,month)),]

  GDDWaveWT<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","GDD", sep=""), fun=GDDAggWT[[2]])
  GDDWaveWT<-na.omit(GDDWaveWT)

  monthStart=monthStart
  subsWT<-subset(GDDWaveWT,GDDWaveWT$month==monthStart) #choose start month
  subsWT<-subsWT[,-which(names(subsWT)=="month")] #Get rid of "month" column
  subsWT$year<-as.numeric(as.character(subsWT$year)) #make "year" readable to merge

  #reshape to show temp over the last nUnits
  labs<-c("year","t.00", paste("t.0",1:9, sep=""),
        paste("t.",10:nUnits, sep=""))
  names(subsWT) <- labs

  # merge climate and butterfly data, following code from Teller 2016
  datag<-merge(useBCB, subsWT, by="year")
  datag$year<-as.factor(datag$year);

  ## define and enter covariates the way the fitting function wants
  tvars <- which(substr(names(datag),1,2)=="t."); 
  datag$tcovar <- as.matrix(datag[,tvars]) 

  lags <- matrix(0,nrow(datag),length(tvars)); 
  for(i in 1:ncol(lags)) lags[,i]=i; 
  datag$lags=as.matrix(lags); 
  datag$yr = strtoi(datag$year)-1990

  return(datag)
}

nUnits = 30 # number of months in the past to include in GAM
strt = 6 # starting month
dat1 = makeDat(useGDD = TPdat[TPdat$Site == "McDevitt",],useBCB = ASSCdat[ASSCdat$Site == "Devitt", ], monthStart = strt)
dat2 = makeDat(useGDD = TPdat[TPdat$Site == "Haynes",], useBCB = ASSCdat[ASSCdat$Site == "Haynes", ], monthStart = strt)
dat3 = makeDat(useGDD = TPdat[TPdat$Site == "Reservoir",],useBCB = ASSCdat[ASSCdat$Site == "Reservoir", ], monthStart = strt)
dat4 = makeDat(useGDD = TPdat[TPdat$Site == "Sheep",],useBCB = ASSCdat[ASSCdat$Site == "Sheep", ], monthStart = strt)

dat.all = rbind(dat1, dat2, dat3, dat4)
str(dat.all)
head(dat.all)

# and now, fit a GAM
gam0 = gam(log(Det.lambda) ~ -1 + Site+ Site:I(strtoi(year)-1990) + s(lags, by=tcovar, bs="cs"), data=dat.all, method="GCV.Cp",gamma=1.2, family = gaussian) 
summary(gam0)
plot(gam0, select = 1, main = "mean temperature", ylab = expression(paste("effect on annual growth rate,   ", lambda)), xlab = "time lag (months before end of demographic year)")
abline(h=0, lty=2)
rect(c(2,14,26),-1,c(8,20,32),1,col=grey(0.3,alpha =0.3), border = NA)
text(c(5,17,29),0.03,labels="Winter")
text(c(11,23),0.03,labels="Growing Season")
arrows(11,0.025,0,0.025)
text(5.5,0.02,labels="Predicted Transition")

# change makeDat to precip
makePrecipDat = function(useGDD, useBCB, monthStart){
  GDDAggWT<-aggLags(datC=useGDD, fun=sum, meas=ppt, seg=month)
  subAgg<-GDDAggWT[[1]]
  subAgg$month<-as.numeric(as.character(subAgg$month))
  subAgg<-subAgg[with(subAgg,order(year,month)),]
  
  GDDWaveWT<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","GDD", sep=""), fun=GDDAggWT[[2]])
  GDDWaveWT<-na.omit(GDDWaveWT)
  
  monthStart=monthStart
  subsWT<-subset(GDDWaveWT,GDDWaveWT$month==monthStart) #choose start month
  subsWT<-subsWT[,-which(names(subsWT)=="month")] #Get rid of "month" column
  subsWT$year<-as.numeric(as.character(subsWT$year)) #make "year" readable to merge
  
  #reshape to show temp over the last nUnits
  labs<-c("year","t.00", paste("t.0",1:9, sep=""),
          paste("t.",10:nUnits, sep=""))
  names(subsWT) <- labs
  
  # merge climate and butterfly data, following code from Teller 2016
  datag<-merge(useBCB, subsWT, by="year")
  datag$year<-as.factor(datag$year);
  
  ## define and enter covariates the way the fitting function wants
  tvars <- which(substr(names(datag),1,2)=="t."); 
  datag$tcovar <- as.matrix(datag[,tvars]) 
  
  lags <- matrix(0,nrow(datag),length(tvars)); 
  for(i in 1:ncol(lags)) lags[,i]=i; 
  datag$lags=as.matrix(lags); 
  datag$yr = strtoi(datag$year)-1990
  
  return(datag)
}

nUnits = 30 # number of months in the past to include in GAM
strt = 6 # starting month
dat1 = makePrecipDat(useGDD = TPdat[TPdat$Site == "McDevitt",],useBCB = ASSCdat[ASSCdat$Site == "Devitt", ], monthStart = strt)
dat2 = makePrecipDat(useGDD = TPdat[TPdat$Site == "Haynes",], useBCB = ASSCdat[ASSCdat$Site == "Haynes", ], monthStart = strt)
dat3 = makePrecipDat(useGDD = TPdat[TPdat$Site == "Reservoir",],useBCB = ASSCdat[ASSCdat$Site == "Reservoir", ], monthStart = strt)
dat4 = makePrecipDat(useGDD = TPdat[TPdat$Site == "Sheep",],useBCB = ASSCdat[ASSCdat$Site == "Sheep", ], monthStart = strt)

dat.ppt = rbind(dat1, dat2, dat3, dat4)
str(dat.all)
head(dat.all)

gam.ppt = gam(log(Det.lambda) ~ -1 + Site+ Site:I(strtoi(year)-1990) + s(lags, by=tcovar, bs="cs"), data=dat.ppt, method="GCV.Cp",gamma=1.2, family = gaussian) 
summary(gam.ppt)
plot(gam.ppt, select = 1, main = "precipitation", ylab = expression(paste("effect on annual growth rate,   ", lambda)), xlab = "time lag (months before end of demographic year)")
abline(h=0, lty=2)
rect(c(2,14,26),-1,c(8,20,32),1,col=grey(0.3,alpha =0.3), border = NA)
text(c(5,17,29),0.0035,labels="Winter")
text(c(11,23),0.0035,labels="Growing Season")
arrows(11,0.00325,0,0.00325)
text(5.5,0.003,labels="Predicted Transition")

gam.check(gam.ppt)
plot(fitted(gam.ppt), log(dat.ppt$Det.lambda))
abline(a=0,b=1, lty=2)
abline(lm(log(dat.ppt$Det.lambda)~fitted(gam.ppt)))

# change makeDat to SPEI
makeSPEIDat = function(useGDD, useBCB, monthStart){
  GDDAggWT<-aggLags(datC=useGDD, fun=sum, meas=SPEI_1, seg=month)
  subAgg<-GDDAggWT[[1]]
  subAgg$month<-as.numeric(as.character(subAgg$month))
  subAgg<-subAgg[with(subAgg,order(year,month)),]
  
  GDDWaveWT<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","GDD", sep=""), fun=GDDAggWT[[2]])
  GDDWaveWT<-na.omit(GDDWaveWT)
  
  monthStart=monthStart
  subsWT<-subset(GDDWaveWT,GDDWaveWT$month==monthStart) #choose start month
  subsWT<-subsWT[,-which(names(subsWT)=="month")] #Get rid of "month" column
  subsWT$year<-as.numeric(as.character(subsWT$year)) #make "year" readable to merge
  
  #reshape to show temp over the last nUnits
  labs<-c("year","t.00", paste("t.0",1:9, sep=""),
          paste("t.",10:nUnits, sep=""))
  names(subsWT) <- labs
  
  # merge climate and butterfly data, following code from Teller 2016
  datag<-merge(useBCB, subsWT, by="year")
  datag$year<-as.factor(datag$year);
  
  ## define and enter covariates the way the fitting function wants
  tvars <- which(substr(names(datag),1,2)=="t."); 
  datag$tcovar <- as.matrix(datag[,tvars]) 
  
  lags <- matrix(0,nrow(datag),length(tvars)); 
  for(i in 1:ncol(lags)) lags[,i]=i; 
  datag$lags=as.matrix(lags); 
  datag$yr = strtoi(datag$year)-1990
  
  return(datag)
}

nUnits = 30 # number of months in the past to include in GAM
strt = 6 # starting month
dat1 = makeSPEIDat(useGDD = TPdat[TPdat$Site == "McDevitt",],useBCB = ASSCdat[ASSCdat$Site == "Devitt", ], monthStart = strt)
dat2 = makeSPEIDat(useGDD = TPdat[TPdat$Site == "Haynes",], useBCB = ASSCdat[ASSCdat$Site == "Haynes", ], monthStart = strt)
dat3 = makeSPEIDat(useGDD = TPdat[TPdat$Site == "Reservoir",],useBCB = ASSCdat[ASSCdat$Site == "Reservoir", ], monthStart = strt)
dat4 = makeSPEIDat(useGDD = TPdat[TPdat$Site == "Sheep",],useBCB = ASSCdat[ASSCdat$Site == "Sheep", ], monthStart = strt)

dat.SPEI = rbind(dat1, dat2, dat3, dat4)
str(dat.all)
head(dat.all)

gam.SPEI_1 = gam(log(Det.lambda) ~ -1 + Site+ Site:I(strtoi(year)-1990) + s(lags, by=tcovar, bs="cs"), data=dat.SPEI, method="GCV.Cp",gamma=1.2, family = gaussian) 
summary(gam.SPEI_1)
plot(gam.SPEI_1, select = 1, main = "SPEI, 1 month", ylab = expression(paste("effect on annual growth rate,   log(", lambda,")")), xlab = "time lag (months before end of demographic year)")
abline(h=0, lty=2)
rect(c(2,14,26),-1,c(8,20,32),1,col=grey(0.3,alpha =0.3), border = NA)
text(c(5,17,29),0.11,labels="Winter")
text(c(11,23),0.11,labels="Growing Season")
arrows(11,0.1,0,0.1)
text(5.5,0.09,labels="Predicted Transition")
gam.check(gam.SPEI_1)

plot(fitted(gam.SPEI_1), log(dat.SPEI$Det.lambda))
abline(a=0,b=1, lty=2)
abline(lm(log(dat.SPEI$Det.lambda)~fitted(gam.SPEI_1)))

# load in vital rates of Astragalus
load("workshop/VR4.reservoir.txt")
load("workshop/VR4.haynes.txt")
load("workshop/VR4.devitt.txt")
load("workshop/VR4.sheep.txt")

# make them into data.frames
reservoir.df <- as.data.frame(VR4.reservoir,row.names = 1:nrow(VR4.reservoir))
haynes.df <- as.data.frame(VR4.haynes,row.names = 1:nrow(VR4.haynes))
devitt.df <- as.data.frame(VR4.devitt,row.names = 1:nrow(VR4.devitt))
sheep.df <- as.data.frame(VR4.sheep,row.names = 1:nrow(VR4.sheep))

# add site to each one
site_labels <- levels(ASSCdat$Site)
devitt.df$Site <- site_labels[1]
haynes.df$Site <- site_labels[2]
reservoir.df$Site <- site_labels[3]
sheep.df$Site <- site_labels[4]

# now stitch each site together with the weather data to get the lag variables
nUnits = 30 # number of months in the past to include in GAM
strt = 6 # starting month
# Change Year to year ...
names(devitt.df)[1] <- "year"
names(haynes.df)[1] <- "year"
names(reservoir.df)[1] <- "year"
names(sheep.df)[1] <- "year"

# Add one to year to account for different labelling
devitt.df$year <- strtoi(devitt.df$year) + 1
haynes.df$year <- strtoi(haynes.df$year) + 1
reservoir.df$year <- strtoi(reservoir.df$year) + 1
sheep.df$year <- strtoi(sheep.df$year) + 1

dat1 = makeSPEIDat(useGDD = TPdat[TPdat$Site == "McDevitt",],useBCB = devitt.df, monthStart = strt)
dat2 = makeSPEIDat(useGDD = TPdat[TPdat$Site == "Haynes",], useBCB = haynes.df, monthStart = strt)
dat3 = makeSPEIDat(useGDD = TPdat[TPdat$Site == "Reservoir",],useBCB = reservoir.df, monthStart = strt)
dat4 = makeSPEIDat(useGDD = TPdat[TPdat$Site == "Sheep",],useBCB = sheep.df, monthStart = strt)

dat.SPEI = rbind(dat1, dat2, dat3, dat4)


# logit transform all the probability columns
logit <- binomial()$linkfun
temp <- lapply(dat.SPEI[3:17],logit)
temp2 <- lapply(dat.SPEI[18:20],log)
dat.SPEI <- cbind(dat.SPEI[1:2],temp,temp2, dat.SPEI[21:ncol(dat.SPEI)])

# remove Dsurv and FtoD from analysis -- no variation at one or more sites
vars <- names(dat.SPEI)[c(3:15,17:20)]


lapply(vars,function(varname)hist(dat.SPEI[,varname], main = varname))
# yr is centered around 1990
core.formula <- " ~ -1 + Site + Site:yr + s(lags, by=tcovar, bs=\"cs\")"

formulas <- lapply(paste(vars,core.formula),as.formula)
list_of_gams <- lapply(formulas, function(ff) gam(ff, data=dat.SPEI, method="GCV.Cp",gamma=1.2, family = gaussian)) 

VRplot <- function(ff, driver = "SPEI, 1 month", ylim = NULL){
  variable = as.character(formula(ff))[2]
  s.table <- summary(ff)$s.table
  plot(ff, select = 1, main = variable, ylim = ylim, 
       ylab = "effect on transition rate", 
       xlab = "time lag (months before end of demographic year)")
  abline(h=0, lty=2)
  figcoord <- par("usr")
  rect(c(2,14,26),-1,c(8,20,32),1,col=grey(0.3,alpha =0.3), border = NA)
  text(c(5,17,29),0,labels="Winter")
  text(c(11,23),0,labels="Growing Season")
  arrows(11,0.95*figcoord[4],0,0.95*figcoord[4])
  text(5.5,0.9*figcoord[4],labels="Predicted Transition")
  text(0.9*figcoord[2], 0.9*figcoord[4], labels = paste("p = ",format(s.table[1,4], scientific = FALSE)))
}

pdf(file = "vital_rates_vs_SPEI.pdf")
lapply(list_of_gams,VRplot)
dev.off()

## now try precipitation
dat1 = makePrecipDat(useGDD = TPdat[TPdat$Site == "McDevitt",],useBCB = devitt.df, monthStart = strt)
dat2 = makePrecipDat(useGDD = TPdat[TPdat$Site == "Haynes",], useBCB = haynes.df, monthStart = strt)
dat3 = makePrecipDat(useGDD = TPdat[TPdat$Site == "Reservoir",],useBCB = reservoir.df, monthStart = strt)
dat4 = makePrecipDat(useGDD = TPdat[TPdat$Site == "Sheep",],useBCB = sheep.df, monthStart = strt)

dat.ppt = rbind(dat1, dat2, dat3, dat4)

# logit transform all the probability columns
logit <- binomial()$linkfun
temp <- lapply(dat.ppt[3:17],logit)
temp2 <- lapply(dat.ppt[18:20],log)
dat.ppt <- cbind(dat.ppt[1:2],temp,temp2, dat.ppt[21:ncol(dat.ppt)])

# remove Dsurv and FtoD from analysis -- no variation at one or more sites
vars <- names(dat.ppt)[c(3:15,17:20)]

ppt_gams <- lapply(formulas, function(ff) gam(ff, data=dat.ppt, method="GCV.Cp",gamma=1.2, family = gaussian)) 

pdf(file = "vital_rates_vs_ppt.pdf")
lapply(ppt_gams,VRplot)
dev.off()

