# identify individual matrices with time series > 20 years
ssdb <- subsetDB2(
              SurvivalIssue < 1, compadre)

ssmd <- ssdb$metadata
ssmd <- compadre$metadata

ssmd %>% filter(StudyDuration > 20) %>% group_by(SpeciesAccepted, MatrixPopulation) %>%
  summarize(n_ = n(), MC = first(MatrixComposite)) %>%
  select(SpeciesAccepted, n_, MC)

ssmd %>%
  filter(SpeciesAccepted == "Acer saccharum") %>%
  select(SpeciesAccepted, MatrixComposite, MatrixStartYear, MatrixPopulation)

ssmd %>% group_by(SpeciesAccepted, MatrixPopulation) %>%
  filter(MatrixComposite == "Individual", 
         MatrixTreatment == "Unmanipulated") %>%
  summarize(n_ = n(), start = min(MatrixStartYear), end = max(MatrixStartYear), duration = end - start, md = first(Authors)) %>%
  select(SpeciesAccepted, n_, start, end, duration, md) %>%
  filter(n_ > 15) 

ssmd %>% group_by(MatrixPopulation, MatrixStartYear) %>%
  filter(MatrixComposite == "Individual", 
         MatrixTreatment == "Unmanipulated",
         SpeciesAccepted == "Astragalus scaphoides")  %>%
  summarize(n_ = n()) %>%
  #select(MatrixPopulation, MatrixStartYear) %>%
  spread(MatrixPopulation, n_, fill = 0)

ssmd %>% filter(SpeciesAccepted == "Cirsium undulatum", 
                MatrixComposite == "Individual", 
                  MatrixTreatment == "Unmanipulated") %>%
  select(Authors)


ssdb <- subsetDB(compadre, Authors == "Dalgleish; Koons; Adler" & 
                   MatrixComposite == "Individual" & 
                   MatrixTreatment == "Unmanipulated") 
                  compadre)
ssdb$mat[[1]]
weather <- read.csv("Workshop/hays_weather_1931-2972.csv")

# EMXP - Extreme maximum daily precipitation
# MMNT - Monthly Mean minimum temperature
# EMXT - Extreme maximum daily temperature
# DT90 - Number days with maximum temperature greater than or equal 90.0 F
# DP10 - Number of days with greater than or equal to 1.0 inch of precipitation
# HTDD - Heating degree days
# DP01 - Number of days with greater than or equal to 0.1 inch of precipitation
# TPCP - Total precipitation
# EMNT - Extreme minimum daily temperature
# DT32 - Number days with minimum temperature less than or equal to 32.0 F
# TSNW - Total snow fall
# DT00 - Number days with minimum temperature less than or equal to 0.0 F
# DP05 - Number of days with greater than or equal to 0.5 inch of precipitation
# DX32 - Number days with maximum temperature less than or equal to 32.0 F
# CLDD - Cooling degree days
# MXSD - Maximum snow depth
# MMXT - Monthly Mean maximum temperature
# MNTM - Monthly mean temperature

library(lubridate)
weather$DATE <- ymd(weather$DATE)

## convert the missing data
pick <- weather$MNTM < -1000
weather[pick, "MNTM"] <- NA
ggplot(weather, aes(DATE, MNTM)) + geom_line()

# too short
#weather <- read.csv("workshop/daily_climate_data.csv")

#weather$date <- with(weather,ymd(paste(year, month, day, sep="-")))
                
#ggplot(weather, aes(date, max)) + geom_step()     


x <- subsetDB2(MatrixTreatment == "Unmanipulated" & MatrixComposite != "Mean" & Country == "USA", compadre)
sp <- names(table(x$metadata$SpeciesAccepted))
n <- as.vector(table(x$metadata$SpeciesAccepted))
temp <- data.frame(sp,n)
subset(temp,n > 19)

x <- subsetDB2(MatrixTreatment == "Unmanipulated" & MatrixComposite != "Mean" & SpeciesAccepted == "Trillium grandiflorum", compadre)
