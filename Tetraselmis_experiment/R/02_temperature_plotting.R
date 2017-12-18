### ibutton temperature data plotting

# load libraries ----------------------------------------------------------

library(tidyverse)
library(purrr)
library(lubridate)
library(stringr)
library(plotrix)
library(janitor)



# load data ---------------------------------------------------------------

var24_raw <- read_csv("Tetraselmis_experiment/data-raw/ibuttons/24-variable-feb24.csv", skip = 14)
var10_raw <- read_csv("Tetraselmis_experiment/data-raw/ibuttons/10-variable-feb24.csv", skip = 14)
var30_raw <- read_csv("Tetraselmis_experiment/data-raw/ibuttons/30-variable-feb24.csv", skip = 14)


# plot it! ----------------------------------------------------------------

var24 <- var24_raw %>% 
	clean_names() %>% 
	mutate(date_time = mdy_hms(date_time)) %>% 
	mutate(treatment = "variable 24C")
	

var30 <- var30_raw %>% 
	clean_names() %>% 
	mutate(date_time = mdy_hms(date_time)) %>% 
	mutate(treatment = "variable 30C")

var10 <- var10_raw %>% 
	clean_names() %>%
	mutate(date_time = mdy_hms(date_time)) %>%
	mutate(treatment = "variable 10C")

all_temps <- bind_rows(var24, var30, var10)
	
	
all_temps %>% 
	filter(treatment == "variable 30C") %>% View
ggplot(aes(x = date_time, y = value)) + geom_point() + theme_bw() + ylab("temperature (C)") 

+
	facet_wrap( ~ treatment)



