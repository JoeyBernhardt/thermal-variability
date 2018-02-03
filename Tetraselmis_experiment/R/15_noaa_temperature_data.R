### trying to extract the temperature data from NOAA



### now trying with Ropensci's rerddapp
library(rerddap)
library(tidyverse)
library(purrr)
library(ncdf4)
library(cowplot)

cache_details()
cache_delete_all(force = TRUE)


#### search for data ####
?ed_search

out <- ed_search(query = 'CPC')
SST <- str_subset(string = out$info$title, pattern = "SST, Daily Optimum Interpolation")
out$info[[out$info$title == "SST, Daily Optimum Interpolation (OI), AMSR+AVHRR, Version 2, 2002-2011, Lon+/-180"]]
info('SST, Daily Optimum Interpolation (OI), AMSR+AVHRR, Version 2, 2002-2011, Lon+/-180')

info_df <- out$info

info_df %>% 
  filter(grepl("tmax", title)) %>% View


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


thomas_locations <- read_csv("Tetraselmis_experiment/data-processed/thomas_locations_nearest.csv")

thomas_locations2 <- thomas_locations %>% 
  mutate(latitude = ifelse(!is.na(nearest_lat), nearest_lat, latitude)) %>% 
  mutate(longitude = ifelse(!is.na(nearest_long), nearest_long, longitude))

### problematic isolates 34, 85, 86, 146, 173, 320, 344, 345, 346, 462, 463, 464, 465
## 466, 467, 470, 472, 570, 571, 603

thomas_split <- thomas_locations2 %>% 
	# filter(isolate.code == 34) %>% 
	split(.$isolate.code)
time_start <- c("2002-07-01")
time_end <- c("2002-12-31")

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

write_csv(temperatures, "Tetraselmis_experiment/data-processed/daily_temps_2002.csv")


### now to get temps from before 2002, need just AVHRR
## ncdcOisst2Agg_LonPM180

time_start <- c("1982-01-01")
time_end <- c("1982-12-31")
extract_function <- function(df) {
  results <- griddap('ncdcOisst2Agg_LonPM180',
                     time = c(time_end, time_start),
                     latitude = c(df$latitude[[1]], df$latitude[[1]]),
                     longitude = c(df$longitude[[1]], df$longitude[[1]]),
                     fields = "sst")
  output <- results$data
}

temperatures <- thomas_split %>% 
  map_df(extract_function, .id = "isolate.code")
write_csv(temperatures, "Tetraselmis_experiment/data-processed/daily_temps_1982.csv")

temperatures %>% 
  filter(is.na(sst)) %>% View

# plot all the freq histograms --------------------------------------------
t1982 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1982.csv")
t1983 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1983.csv")
t1984 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1984.csv")
t1985 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1985.csv")
t1986 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1986.csv")
t1987 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1987.csv")
t1988 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1988.csv")
t1989 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1989.csv")
t1990 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1990.csv")
t1991 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1991.csv")
t1992 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1992.csv")
t1993 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1993.csv")
t1994 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1994.csv")
t1995 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1995.csv")
t1996 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1996.csv")
t1997 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1997.csv")
t1998 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1998.csv")
t1999 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_1999.csv")
t2000 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_2000.csv")
t2001 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_2001.csv")
t2002 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_2002.csv")
t2003 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_2003.csv")
t2004 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_2004.csv")
t2005 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_2005.csv")
t2006 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_2006.csv")
t2007 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_2007.csv")
t2008 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_2008.csv")
t2009 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_2009.csv")
t2010 <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_2010.csv")


tsx <- bind_rows(t1982, t1983, t1984, t1985, t1986, t1987, t1988, t1989, t1990, t1991, t1992, t1993, t1994, t1995, t1996, t1997, t1998, t2003, t2004, t2005, t2006, t2007, t2008, t2009, t1999, t2000, t2001, t2002, t2010)

tsx2 <- tsx %>% 
  distinct(isolate.code, time, lat, long, .keep_all = TRUE)

write_csv(tsx2, "Tetraselmis_experiment/data-processed/OISST_data.csv")

### put the temperature and the thermal trait data together
tc <- left_join(tsx2, thomas3, by = "isolate.code")

growth_rates <- tc %>% 
  rename(temp = sst) %>% 
  group_by(isolate.code) %>% 
  mutate(growth_rate = a*exp(b*temp)*(1-((temp-z)/(w/2))^2)) %>% 
  select(isolate.code, temp, a, b, z, w, growth_rate)

growth_summary <- growth_rates %>% 
  group_by(isolate.code) %>% 
  summarise_each(funs(mean), growth_rate)


all_thermal_data <- read_csv("Tetraselmis_experiment/data-processed/all_thermal_data.csv")
mean <- all_thermal_data %>% 
  group_by(isolate.code) %>% 
  mutate(growth_rate_mean = a*exp(b*Mean)*(1-((Mean-z)/(w/2))^2)) %>% 
  select(isolate.code, growth_rate_mean, Mean)


all <- left_join(growth_summary, mean, by = "isolate.code")


all %>% 
  ggplot(aes(x = growth_rate_mean, y = growth_rate)) + geom_point() +
  geom_abline(slope = 1, intercept = 0)


results5 <- read_csv("Tetraselmis_experiment/data-processed/results5.csv")
predicted_growth_temperature2 <- read_csv("Tetraselmis_experiment/data-processed/predicted_growth_temperature2.csv")

pred_var <- predicted_growth_temperature2 %>% 
  select(isolate.code, predicted_growth_variable)

### here we have the time-averaged growth rate and the STT predicted growth rate.
all2 <- left_join(all, predicted_growth_temperature2, by = "isolate.code")
write_csv(all2, "Tetraselmis_experiment/data-processed/all2.csv")


all2 %>% 
  mutate(difference = growth_rate - predicted_growth_variable) %>% 
  ggplot(aes(x = difference)) + geom_histogram() 

tsx2 %>% 
  filter(isolate.code == 606) %>% 
  ggplot(aes(x = time, y = sst)) + geom_line()


tsx2 %>% 
  mutate(lat = round(lat, digits = 1)) %>% 
  mutate(lon = round(lon, digits = 1)) %>% 
  unite(lat_long_isolate, c("lat", "lon", "isolate.code"), sep = ".", remove = FALSE) %>% 
  mutate(region = ifelse(abs(lat) > 30, "temperate", "tropical")) %>% 
  mutate(region = ifelse(abs(lat) > 60, "polar", region)) %>% 
  ggplot(aes(x = sst, fill = region)) + geom_density() +
  facet_wrap(~ lat_long_isolate, scales = "free_y")
ggsave("Tetraselmis_experiment/figures/temp_histograms.pdf", width = 18, height = 14)
ggsave("Tetraselmis_experiment/figures/temp_density.pdf", width = 18, height = 14)


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
