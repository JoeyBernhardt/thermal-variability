
### attempt to detrend the seasonal data

library(tidyverse)
library(lubridate)
library(timeSeries)
library(timetk)
library(viridis)
library(cowplot)

temp1 <- read_csv("Tetraselmis_experiment/data-processed/daily_temperatures_1-97.csv")
four_years_data <- read_csv("Tetraselmis_experiment/data-processed/four_years_data.csv")
temps_all <- read_csv("Tetraselmis_experiment/data-processed/temps_all_33.875N.csv")


four_years_data %>% 
  ggplot(aes(x = sst)) + geom_histogram() + 
  facet_wrap( ~ lat, scales = "free_y")


str(four_years_data)

four_years_data %>% 
	# mutate(time = ymd_hms(time)) %>% 
	ggplot(aes(x = time, y = sst)) + geom_point() + theme_classic()


four_years <- four_years_data %>% 
	select(time, sst)


four <- timetk::tk_ts(four_years, select = sst, freq = 365)
temps <- timetk::tk_ts(temps_all, select = sst, freq = 365)

out <- decompose(temps, type = "additive")
str(out)
head(four)

five_year <- temps_all %>% 
	select(Date, sst) %>% 
	rename(day = Date) %>% 
	rename(value = sst) %>% 
	mutate(date = ymd(day)) %>% 
	select(-day) %>% 
	mutate(type = "raw data")


seasonal <- tk_tbl(out$seasonal, rename_index = "day") %>% 
	mutate(type = "seasonal") %>% 
	mutate(date = ymd(five_year$date)) %>% 
	select(-day)

trend <- tk_tbl(out$trend, rename_index = "day") %>% 
	mutate(type = "trend") %>% 
	mutate(date = ymd(five_year$date)) %>% 
	select(-day)

random <- tk_tbl(out$random, rename_index = "day") %>% 
	mutate(type = "random") %>% 
	rename(value = `x - seasonal`) %>% 
	mutate(date = ymd(five_year$date)) %>% 
	select(-day)

str(random)
str(trend)
str(seasonal)
str(five_year)

all_output <- bind_rows(seasonal, trend, random, five_year)

	
all_output %>% 
 dplyr::filter(type == "random") %>%
ggplot(aes(x = date, y = value, color = type)) + geom_line(size = 1) +
	theme_classic() + scale_color_brewer(palette = "RdBu") + ylab("Temperature (°C)") + xlab("Year")

all_output %>% 
	dplyr::filter(!is.na(value)) %>% 
	group_by(type) %>% 
	summarise_each(funs(mean, sd), value) %>% View
	ggplot(aes(x = type, y = value_sd)) + geom_point()
	
	
	all_output %>% 
		filter(type == "random") %>% 
		filter(!is.na(value)) %>% 
		mutate(abs_dev = abs(value)) %>% 
		summarise(mean_dev = mean(abs_dev))


?decompose
	