### trying to extract the temperature data from NOAA



### now trying with Ropensci's rerddapp
library(rerddap)
library(tidyverse)
library(purrr)
library(ncdf4)

cache_delete("/Users/joeybernhardt/Library/Caches/R/rerddap/e6e4113fc99e88a79b7294717a602145.nc")
cache_delete_all(force = TRUE)
cache_details()

#### not sure what this is for ####
out <- ed_search(query = 'temperature')
SST <- str_subset(string = out$info$title, pattern = "SST, Daily Optimum Interpolation")
out$info[[out$info$title == "SST, Daily Optimum Interpolation (OI), AMSR+AVHRR, Version 2, 2002-2011, Lon+/-180"]]
info('SST, Daily Optimum Interpolation (OI), AMSR+AVHRR, Version 2, 2002-2011, Lon+/-180')

info_df <- out$info


info_df %>% 
	filter(title == "SST, Daily Optimum Interpolation (OI), AMSR+AVHRR, Version 2, 2002-2011, Lon+/-180")

griddap('ncdcOisst2AmsrAgg_LonPM180')


out <- griddap('ncdcOisst2AmsrAgg_LonPM180',
								time = c('2006-01-01', '1995-01-01'),
								latitude = c(33.875, 33.875),
								longitude = c(131.125, 131.125),
								fields = "sst")



reprex(venue = "gh", si = TRUE)


### ok so it looks like rerddap doesn't like to have to get more than 7 years of data at a time. What if we try and get only 5 years of data at a time, but then iterate over several chunks of time?


data_source <- 'ncdcOisst2AmsrAgg_LonPM180'

extract_temps <- function(df) {
	out <- griddap('ncdcOisst2AmsrAgg_LonPM180',
					time = c(time_end, time_start),
					latitude = c(33.875, 33.875),
					longitude = c(131.125, 131.125),
					fields = "sst",
					fmt = "nc")
	results <- out$data
	results
} 


time_start <- c("2010-01-01", "2008-01-01", "2006-01-01", "2004-01-01", "2002-01-01")
time_end <- c("2012-01-01", "2010-01-01", "2008-01-01", "2006-01-01", "2004-01-01")
final_year <- c("2012", "2010", "2008", "2006", "2004")
df <- data.frame(time_start, time_end, final_year) %>% 
	mutate(time_start = as.character(time_start),
				 time_end = as.character(time_end))

df_split <- df %>% 
	split(.$final_year)

df <- df[1,]

collected_data <- df_split %>% 
	map_df(extract_temps, .id = "last_time")

ndf_split

four_years_data <- four_years$data
write_csv(four_years_data, "Tetraselmis_experiment/data-processed/four_years_data2.csv")
four_years_data <- read_csv("Tetraselmis_experiment/data-processed/four_years_data2.csv")


thomas <- read_csv("data/thermal_trait_data/Thomas_2014_traits_derived_20140606.csv")

thomas2 <- thomas %>% 
	rename(topt = mu.g.opt.list) %>% 
	rename(w = mu.wlist,
				 a = mu.alist,
				 z = mu.c.opt.list,
				 b = mu.blist) %>% 
	select(topt, w, a, b, z, everything()) 
thomas3 <- thomas2 %>% 
	rename(latitude = isolation.latitude,
				 longitude = isolation.longitude) %>% 
	filter(!is.na(latitude), !is.na(longitude)) %>% 
	filter(habitat %in% c("marine", "estuarine", "saltmarsh")) %>% 
	filter(curvequal == "good") %>% 
  filter(mu.rsqrlist > 0.85)

thomas_locations <- thomas3 %>% 
	select(isolate.code, latitude, longitude)

write_csv(thomas_locations, "Tetraselmis_experiment/data-processed/thomas_locations.csv")

	

# bring in location data --------------------------------------------------


thomas_locations <- read_csv("Tetraselmis_experiment/data-processed/thomas_locations.csv")

### problematic isolates 34, 85, 86, 146, 173, 320, 344, 345, 346, 462, 463, 464, 465
## 466, 467, 470, 472, 570, 571, 603

thomas_split <- thomas_locations %>% 
	filter(isolate.code == thomas_locations$isolate.code[[91]]) %>%
	split(.$isolate.code)
time_start <- c("2004-01-01")
time_end <- c("2004-01-5")

extract_function <- function(df) {
	results <- griddap('ncdcOisst2AmsrAgg_LonPM180',
					time = c(time_end, time_start),
					latitude = c(df$latitude[[1]], df$latitude[[1]]),
					longitude = c(df$longitude[[1]], df$longitude[[1]]),
					fields = "sst")
	output <- results$data
}

temperatures <- thomas_split %>% 
	map_df(extract_function, .id = "isolate.code")
(temperatures)



write_csv(temperatures, "Tetraselmis_experiment/data-processed/daily_temperatures_1-97.csv")

write_csv(temperatures_2005, "Tetraselmis_experiment/data-processed/daily_temps_2005.csv")
write_csv(temperatures_2005_89, "Tetraselmis_experiment/data-processed/daily_temps_isolate89_2005.csv")

temperatures_above99 <- thomas_split %>% 
	map_df(extract_function, .id = "isolate.code")


write_csv(temperatures_above99, "Tetraselmis_experiment/data-processed/daily_temperatures_above99.csv")


# missing isolates --------------------------------------------------------


## now have another go at it for the 23 missing lat longs
library(rerddap)
library(purrr)
library(tidyverse)
library(ncdf4)

missing_isolates <- read_csv("Tetraselmis_experiment/data-processed/missing_isolates_NOAA.csv")

missing_split <- missing_isolates %>% 
	split(.$isolate.code)

extract_function_wider <- function(df) {
	results <- griddap('ncdcOisst2AmsrAgg_LonPM180',
										 time = c('2011-01-01', '2010-01-01'),
										 latitude = c(df$latitude[[1]], df$latitude[[1]] + 0.25),
										 longitude = c(df$longitude[[1]], df$longitude[[1]] + 0.25),
										 fields = "sst")
	output <- results$data
}



temperatures_missing <- missing_split %>% 
	map_df(extract_function_wider, .id = "isolate.code")


write_csv(temperatures_missing, "Tetraselmis_experiment/data-processed/temperatures_missing_round0.25.csv")

found_temperatures_0.25 <- temperatures_missing %>% 
	filter(!is.na(sst)) %>% 
	distinct(isolate.code) %>% 
	mutate(isolate.code = as.numeric(isolate.code))

still_missing <- setdiff(missing_isolates$isolate.code, found_temperatures_0.25$isolate.code)


still_missing_lats <- missing_isolates %>% 
	filter(isolate.code %in% still_missing)

still_missing_split <- still_missing_lats %>% 
	split(.$isolate.code)

## now extract the things missing with a wider net
extract_function_even_wider <- function(df) {
	results <- griddap('ncdcOisst2AmsrAgg_LonPM180',
										 time = c('2011-01-01', '2010-01-01'),
										 latitude = c(df$latitude[[1]], df$latitude[[1]] + 0.55),
										 longitude = c(df$longitude[[1]], df$longitude[[1]] + 0.55),
										 fields = "sst")
	output <- results$data
}



temperatures_still_missing <- still_missing_split %>% 
	map_df(extract_function_even_wider, .id = "isolate.code")

write_csv(temperatures_still_missing, "Tetraselmis_experiment/data-processed/temperatures_still_missing.csv")


found_isolates <- temperatures_still_missing %>% 
	filter(!is.na(sst)) %>% 
	distinct(isolate.code)

still_missing_again <- setdiff(still_missing_lats$isolate.code, found_isolates$isolate.code)



temperatures_still_missing_again <- missing_isolates %>% 
	filter(isolate.code %in% still_missing_again) %>% 
	split(.$isolate.code)


extract_function_even_wider_again <- function(df) {
	results <- griddap('ncdcOisst2AmsrAgg_LonPM180',
										 time = c('2011-01-01', '2010-01-01'),
										 latitude = c(df$latitude[[1]], df$latitude[[1]] + 3.5),
										 longitude = c(df$longitude[[1]], df$longitude[[1]] + 3.5),
										 fields = "sst")
	output <- results$data
}



temperatures_last_time <- temperatures_still_missing_again %>% 
	map_df(extract_function_even_wider_again, .id = "isolate.code")

temperatures_last_time <- temperatures_last_time %>% 
	filter(!is.na(sst)) 

write_csv(temperatures_last_time, "Tetraselmis_experiment/data-processed/temepratures_last_time.csv")
