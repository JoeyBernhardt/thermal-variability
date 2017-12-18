library(tidyverse)
library(lubridate)

source("~/Documents/Misc_R_scripts/NOAA_OISST_ncdf4.R")


ssts <- extractOISSTdaily("~/Desktop/sst.day.mean.2017.v2.nc","~/Desktop/lsmask.oisst.v2.nc",lonW=33.875,lonE=33.875,latS=131.125,latN=131.125,date1='2017-01-1',date2='2017-06-24')


s2 <- as.data.frame.table(ssts) %>% 
	distinct(Lat, Long, Date, Freq) %>% 
	rename(sst = Freq)

data_2016 <- extractOISSTdaily("~/Documents/NOAA/sst.day.mean.2016.v2.nc","~/Documents/NOAA/lsmask.oisst.v2.nc",lonW=131.125,lonE=131.125,latS=33.875,latN=33.875,date1='2016-01-1',date2='2016-12-31')
s2016 <- as.data.frame.table(data_2016) %>% 
	distinct(Lat, Long, Date, Freq) %>% 
	rename(sst = Freq)

data_2015 <- extractOISSTdaily("~/Documents/NOAA/sst.day.mean.2015.v2.nc","~/Documents/NOAA/lsmask.oisst.v2.nc",lonW=131.125,lonE=131.125,latS=33.875,latN=33.875,date1='2015-01-1',date2='2015-12-31')
s2015 <- as.data.frame.table(data_2015) %>% 
	distinct(Lat, Long, Date, Freq) %>% 
	rename(sst = Freq)

data_2014 <- extractOISSTdaily("~/Documents/NOAA/sst.day.mean.2014.v2.nc","~/Documents/NOAA/lsmask.oisst.v2.nc",lonW=131.125,lonE=131.125,latS=33.875,latN=33.875,date1='2014-01-1',date2='2014-12-31')
s2014 <- as.data.frame.table(data_2014) %>% 
	distinct(Lat, Long, Date, Freq) %>% 
	rename(sst = Freq)

data_2013 <- extractOISSTdaily("~/Documents/NOAA/sst.day.mean.2013.v2.nc","~/Documents/NOAA/lsmask.oisst.v2.nc",lonW=131.125,lonE=131.125,latS=33.875,latN=33.875,date1='2013-01-1',date2='2013-12-31')
s2013 <- as.data.frame.table(data_2013) %>% 
	distinct(Lat, Long, Date, Freq) %>% 
	rename(sst = Freq)

data_2012 <- extractOISSTdaily("~/Documents/NOAA/sst.day.mean.2012.v2.nc","~/Documents/NOAA/lsmask.oisst.v2.nc",lonW=131.125,lonE=131.125,latS=33.875,latN=33.875,date1='2012-01-1',date2='2012-12-31')
s2012 <- as.data.frame.table(data_2012) %>% 
	distinct(Lat, Long, Date, Freq) %>% 
	rename(sst = Freq)


all_data <- bind_rows(s2016, s2015, s2014, s2013, s2012)


all_data2 <- all_data %>% 
	mutate(Date = ymd(Date),
				 Lat = as.numeric(as.character(Lat)),
				 Long = as.numeric(as.character(Long)))
unique(all_data2$Lat)

all_data2 %>% 
	ggplot(aes(x = Date, y = sst)) + geom_line() + theme_classic()


## ok it looks like this works!
## let's try to make make it into a function
year <- 2012
extract_function <- function(year){
	year <- year[[1]]
data <- extractOISSTdaily(paste0("~/Documents/NOAA/sst.day.mean.", year, ".v2.nc"),"~/Documents/NOAA/lsmask.oisst.v2.nc",lonW=131.125,lonE=131.125,latS=33.875,latN=33.875,date1=paste0(year, '-01-1'),date2= paste0(year, '-12-31'))
data2 <- as.data.frame.table(data) %>% 
	distinct(Lat, Long, Date, Freq) %>% 
	rename(sst = Freq) %>% 
	mutate(Date = ymd(Date),
				 Lat = as.numeric(as.character(Lat)),
				 Long = as.numeric(as.character(Long)))

}

data2013f <- extract_function(2013)


years <- data.frame(year = c(2000:2017))

years_split <- years %>% 
	split(.$year)

temps_all <- years_split %>% 
	map_df(extract_function, .id = "year")

write_csv(temps_all, "Tetraselmis_experiment/data-processed/temps_all_33.875N.csv")


temps_all %>% 
ggplot(aes(x = Date, y = sst)) + geom_line() + theme_classic()

