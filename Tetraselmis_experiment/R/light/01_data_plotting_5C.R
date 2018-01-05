## data exploration and visualization


library(tidyverse)
library(lubridate)
library(broom)
library(stringr)
library(simecol)
library(purrr)

data5 <- read_csv("data-processed/5C_TT_cells.csv")


ldata2 <- data5 %>% 
	mutate(light_level = NA) %>% 
	mutate(light_level = ifelse(replicate %in% c("8", "4", "10", "1", "9"), "high_light", light_level)) %>% 
	mutate(light_level = ifelse(replicate %in% c("13", "5", "14", "11"), "low_light", light_level)) %>% 
	mutate(light_level = ifelse(replicate %in% c("7", "15", "12", "6", "2", "3"), "med_low_light", light_level)) %>%
	mutate(light_level = ifelse(replicate %in% c("16", "17", "18"), "high_N", light_level)) %>% 
	mutate(start_time = ymd_hms(start_time)) %>% 
	separate(start_time, into = c("date", "time"), sep = " ", remove = FALSE) %>% 
	mutate(date = ymd(date))

## get the days in

ldata2$start.time <- ymd("2017-05-04")
ldata2$time_since_innoc <- interval(ldata2$start.time, ldata2$date)


ldata3 <- ldata2 %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))

ldata3 %>% 
	group_by(replicate) %>% 
	ggplot(aes(x = start_time, y = cell_density, color = light_level, group = replicate)) + geom_point() + geom_line() +
	theme_bw()
