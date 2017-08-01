### attempt at figuring out a new set of predictions for the variable treatments!!


library(tidyverse)



sched <- read_csv("Tetraselmis_experiment/data-raw/light_dark_schedule.csv")



sched_sep <- sched %>% 
	gather(key = "day", value = "temperature", 4:24) %>%  
	separate(day, sep = "_", into = c("day1", "day_number"), remove = FALSE) %>% 
	select(-day1)


sched_sep %>% 
	filter(lights == "light") %>% 
	filter(day_number < 8) %>% 
	mutate(temp_time = ifelse(grepl("hot", temperature), 1, 0)) %>% 
	group_by(day) %>% 
	summarise(number_of_hot_hours = sum(temp_time)) %>% 
	summarise_each(funs(mean, min, max), number_of_hot_hours) %>% 
	mutate(fraction_hours_hot = (mean/2)/16) %>% View
	

sched %>% 
	gather(key = "day", value = "temperature", 4:24) %>% 
	filter(lights == "light") %>% 
	mutate(temp_time = ifelse(grepl("hot", temperature), 1, 0)) %>% 
	summarise(number_of_hot_hours = sum(temp_time)) %>% View
