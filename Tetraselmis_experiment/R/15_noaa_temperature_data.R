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
tsx2 <- read_csv("Tetraselmis_experiment/data-processed/OISST_data.csv")

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


### get SD of daily temps

sst_summary <- tsx2 %>% 
  group_by(isolate.code, lat, lon) %>% 
  summarise_each(funs(mean, sd), sst) 


# mean_temps <- 
 
temps <- seq(-2, 40, by = 0.01)

limits_fun <- function(x){
res <- tc %>% 
  group_by(isolate.code) %>% 
  mutate(cent_temp = sst - mean(sst)) %>%
  mutate(cent1 = cent_temp + x) %>% 
  mutate(growth_rate_mean = a*exp(b*cent1)*(1-((cent1-z)/(w/2))^2)) %>% 
  summarise(mean_growth = mean(growth_rate_mean)) %>% 
  mutate(temperature = x)
  return(res)
}
  

### to find limits
resp <- temps %>% 
  map_df(limits_fun, .id = "temp") %>% 
  group_by(isolate.code) %>% 
  filter(mean_growth > 0) %>% 
  group_by(isolate.code) %>% 
  top_n(n = -4, wt = mean_growth) 

### to find topt

resp_topt <- temps %>% 
  map_df(limits_fun, .id = "temp") %>% 
  group_by(isolate.code) %>% 
  # filter(mean_growth > 0) %>% 
  group_by(isolate.code) %>% 
  top_n(n = 1, wt = mean_growth) 


write_csv(resp, "Tetraselmis_experiment/data-processed/resp.csv")
write_csv(resp_topt, "Tetraselmis_experiment/data-processed/resp_topt.csv")
resp <- read_csv("Tetraselmis_experiment/data-processed/resp.csv")
low <- resp %>% 
  group_by(isolate.code) %>% 
  top_n(n = -1, wt = temperature) 

## 77 might be wrong
high <- resp %>% 
  group_by(isolate.code) %>% 
  top_n(n = 1, wt = temperature) 


high_low <- left_join(low, high, by = "isolate.code") %>% 
  rename(low = temperature.x, high = temperature.y) %>% 
  gather(key = type, value = temperature, low, high)


high_low %>% 
  ggplot(aes(x = temperature, y = 0)) + geom_point() + 
  # geom_line(data = temps_growth, aes(x = temp, y = growth_rate)) +
  facet_wrap( ~ isolate.code)

mean_temps %>% 
  filter(isolate.code == 128) %>% 
  ggplot(aes(x = mean_cent)) + geom_histogram() +
  facet_wrap( ~ isolate.code, scales = "free_y")


### let's compare topt_variable from STT with Topt_variable from averaging
resp_topt <- resp_topt %>% 
  rename(topt_averaging = temperature)

topts <- left_join(results5, resp_topt, by = "isolate.code") %>% 
  filter(mu.rsqrlist > 0.85)

topts %>% 
  ggplot(aes(x = topt_variable, y = topt_averaging)) + geom_point() +
  geom_abline(intercept = 0, slope = 1)

topts %>% 
  mutate(topt_diff = topt_averaging - topt_variable) %>% 
  mutate(topt_change = topt_averaging - topt) %>% 
  ggplot(aes(x = topt_change)) + geom_histogram()

topts %>% 
  mutate(topt_diff = topt_averaging - topt_variable) %>% 
  mutate(topt_change = topt_averaging - topt) %>% 
  ggplot(aes(x = rel.curveskew, y = topt_change, color = SD, label = isolate.code)) + geom_point(size = 3) +
  scale_color_viridis(option = "inferno") + theme_bw() + geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
  geom_hline(yintercept = 0) + ylab("Topt(var) - Topt(cons) (°C)") +
  xlab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black")) +
  theme(text = element_text(size=14, family = "Helvetica")) +
  geom_vline(xintercept = 0) + geom_point(size = 3, shape = 1, color = "black") + xlim(-0.02, 0.02)

topts %>% 
  mutate(topt_diff = topt_averaging - topt_variable) %>% 
  mutate(topt_change = topt_averaging - topt) %>% 
  ggplot(aes(x = longitude, y = latitude, color = topt_change)) +
  mapWorld +
  geom_point(size = 4) +
  geom_point(size = 4, shape = 1, color = "black") +
  # geom_point(aes(x = lon, y = lat, color = sst_sd), data = sst_summary, size = 4) +
  # scale_color_viridis(discrete = FALSE, begin = 0, end = 1, option = "inferno") +
  scale_color_gradient2(low = "purple", high = "orange")


ggplot(aes(x = lon, y = lat, color = sst_sd), data = sst_summary) +
  mapWorld +
  geom_point(size = 4) +
  geom_point(size = 4, shape = 1, color = "black") +
  scale_color_viridis(discrete = FALSE, option = "inferno")

library(sf)
library(maps)
world1 <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE, col = 1:10, wrap=c(-180,180)))
ggplot() + geom_sf(data = world1)

sst_data <- st_as_sf(sst_summary, coords = c("lon", "lat"), crs = 4326)
PROJ <- "+proj=wintri +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
sst_data <-st_transform(x = sst_data, crs = "+proj=robin")
world2 <- st_transform(world1,crs = 4326)


# map of daily sst sd -----------------------------------------------------


ggplot(sst_data)+
  scale_color_viridis(discrete = FALSE, option = "inferno", name = "SD of daily SST")+
  geom_sf(data = world1, color = "transparent", fill = "darkgrey") +
  geom_sf(geom = "point", aes(color = sst_sd), size = 4)+
  geom_sf(geom = "point", shape = 1, color = "black", size = 4)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        line = element_blank(),
        rect = element_blank())
ggsave("Tetraselmis_experiment/figures/daily_sst_points_map.pdf")


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
all2 <- read_csv("Tetraselmis_experiment/data-processed/all2.csv")

all2 %>% 
  mutate(difference = growth_rate - predicted_growth_variable) %>% 
  ggplot(aes(x = difference)) + geom_histogram() 




tc2 <- tc %>% 
  mutate(lat = round(lat, digits = 1)) %>% 
  mutate(lon = round(lon, digits = 1)) %>% 
  filter(isolate.code == 49) %>% 
  group_by(isolate.code) %>% 
  unite(lat_long_isolate, c("lat", "lon", "isolate.code"), sep = ".", remove = FALSE) %>% 
  mutate(region = ifelse(abs(lat) > 30, "temperate", "tropical")) %>% 
  mutate(region = ifelse(abs(lat) > 60, "polar", region))

tc2 %>% 
  ggplot(aes(x = sst)) + geom_density(fill = "blue") +
  stat_function(color = "blue", fun = function(sst,z,w,a,b) tc2$a[[1]]*exp(tc2$b[[1]]*sst)*(1-((sst-tc2$z[[1]])/(tc2$w[[1]]/2))^2)) +
  # facet_wrap(~ lat_long_isolate, scales = "free_y") +
  xlim(0, 40) +
  ylab("Frequency") + xlab("Daily SST") +
  scale_y_continuous(sec.axis = dup_axis(name = "Growth rate"), limits =c(0, 2)) +
  ggtitle("Skeletonema tropicum")


tc3 <- tc %>% 
  mutate(lat = round(lat, digits = 1)) %>% 
  mutate(lon = round(lon, digits = 1)) %>% 
  filter(isolate.code == 146)

### ok min SST is -1.8, max sst is 36.12, so maybe go from -2 to 37.
temp.cat <- function(x, lower = -2, upper, by = 0.25, sep = "--", above.char = "+") {
  labs <- c(paste(seq(lower, upper - by, by = by),
                  seq(lower + by, upper, by = by), sep = sep),
            paste(upper, above.char, sep = ""))
  
  cut(floor(x), breaks = c(seq(lower, upper, by = by), Inf),
      right = FALSE, labels = labs)
}


tc_bins <- tc %>% 
  group_by(isolate.code) %>% 
  do(length.con <- as.data.frame(table(temp.cat(.$sst, upper = 42)))) %>% ## use size.cat function to create 1um size class bins
  filter(Freq >0) %>% ## remove bins where there are no cells
  separate(Var1, c("lower", "upper"), sep = "--") %>% ## separate the size bins into upper and lower bounds
  mutate(temperature = as.numeric(upper)) %>% 
  mutate(frequency = Freq / 10411)

tc_all <- left_join(tc_bins, thomas3, by = "isolate.code")

tc_146 <- tc_all %>% 
  filter(isolate.code == 146)

tc_146 %>% 
  ggplot(aes(x = temperature)) + geom_bar(aes(x = temperature, y = frequency*10), stat = "identity", fill = "cadetblue") +
  stat_function(color = "cadetblue", fun = function(sst,z,w,a,b) tc3$a[[1]]*exp(tc3$b[[1]]*sst)*(1-((sst-tc3$z[[1]])/(tc3$w[[1]]/2))^2)) +
  xlim(0, 40) + ylim(0, 1.5) +
  ylab("Frequency (*10)") + xlab("Daily SST at isolation location (°C)") +
  scale_y_continuous(sec.axis = dup_axis(name = "Growth rate"), limits =c(0, 2))


### ok now let's get the TPCs onto thomas3

temp_seq<- rep(seq(-2, 40, by = 0.1), 89)
isolate.code <- sort(rep(thomas3$isolate.code, times = length(seq(-2, 40, by = 0.1))))
isolates_temps <- data_frame(isolate.code, temp_seq)


thomas_growth <- left_join(isolates_temps, thomas3) %>% 
  rename(temp = temp_seq) %>% 
  group_by(isolate.code) %>% 
  mutate(growth_rate = a*exp(b*temp)*(1-((temp-z)/(w/2))^2))

temps_growth <- left_join(thomas_growth, tc_all, by = "isolate.code") %>% 
  filter(isolate.code != 1) 

isol1 <- temps_growth %>% 
  filter(isolate.code ==1) %>% 
  select(frequency, growth_rate, everything())

bins <- temps_growth %>% 
  filter(isolate.code != 1) %>% 
  distinct(isolate.code, temperature, frequency, .keep_all = TRUE) 


ggplot() + geom_bar(aes(x = temperature, y = frequency*4),
                    stat = "identity", fill = "cadetblue", data = bins) +
  geom_hline(yintercept = 0, color = "grey")+
  geom_line(data = temps_growth, aes(x = temp, y = growth_rate)) +
  xlim(-3, 40) +
  facet_wrap( ~ isolate.code) +
  ylab("Frequency") + xlab("Daily SST at isolation location (°C)") +
  scale_y_continuous(sec.axis = dup_axis(name = "Growth rate"), limits = c(-0.5, 2)) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  geom_point(data = high_low, aes(x = temperature, y = 0), color = "red", size =1) 
  
  ggsave("Tetraselmis_experiment/figures/temps_curves_lims.pdf", height = 10, width = 13)


# see how different the tmaxes are w/two approaches -----------------------


  tmax <- high %>% 
    rename(tmax_hist = temperature)
  
 tmaxes <- left_join(results5, tmax, by = "isolate.code") %>% 
   filter(mu.rsqrlist > 0.85)
 
 tmaxes %>% 
   mutate(diff_max = temperature_max - tmax_hist) %>% 
   filter(diff_max < 5) %>% 
   ggplot(aes(x = tmax_hist, y = temperature_max)) + geom_point() +
   ylab("Tmax from STT") + xlab("Tmax from time averaging") +
   geom_abline(slope = 1, intercept = 0)
 
 tmaxes %>% 
   mutate(diff_max = temperature_max - tmax_hist) %>% 
   filter(diff_max < 5) %>% 
   mutate(rev_max_diff = tmax_hist - tmax) %>% 
   ggplot(aes(x = rel.curveskew, y = rev_max_diff, color = SD)) + geom_point(size = 3) +
   scale_color_viridis(option = "inferno") + theme_bw() +
   geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
   ylab("Tmax(var) - Tmax(cons) (°C)") +
   xlab("Curve skewness") + 
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(color="black")) +
   theme(text = element_text(size=14, family = "Helvetica")) + geom_point(size = 3, shape = 1, color = "black") + geom_hline(yintercept = 0) +
   xlim(-0.02, 0.02) + geom_vline(xintercept = 0) 
 
  
tc3 %>% 
  ggplot(aes(x = sst)) + geom_density(fill = "blue") +
  stat_function(color = "blue", fun = function(sst,z,w,a,b) tc3$a[[1]]*exp(tc3$b[[1]]*sst)*(1-((sst-tc3$z[[1]])/(tc3$w[[1]]/2))^2)) +
  # facet_wrap(~ lat_long_isolate, scales = "free_y") +
  xlim(0, 40) +
  ylab("Frequency") + xlab("Daily SST") +
  scale_y_continuous(sec.axis = dup_axis(name = "Growth rate"), limits =c(0, 1.4)) +
  ggtitle("Skeletonema costatum")
  
ggsave("Tetraselmis_experiment/figures/temp_histograms.pdf", width = 18, height = 14)
ggsave("Tetraselmis_experiment/figures/temp_density.pdf", width = 18, height = 14)


function(sst,z,w,a,b) a*exp(b*sst)*(1-((sst-z)/(w/2))^2)


write_csv(temperatures, "Tetraselmis_experiment/data-processed/daily_temperatures_1-97.csv")

write_csv(temperatures_2005, "Tetraselmis_experiment/data-processed/daily_temps_2005.csv")
write_csv(temperatures_2005_89, "Tetraselmis_experiment/data-processed/daily_temps_isolate89_2005.csv")

temperatures_above99 <- thomas_split %>% 
	map_df(extract_function, .id = "isolate.code")


write_csv(temperatures_above99, "Tetraselmis_experiment/data-processed/daily_temperatures_above99.csv")



# plotting skeletonema isolates -------------------------------------------

thomas3 %>% 
  filter(grepl("Skeletonema", speciesname)) %>% View






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
