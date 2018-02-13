
library(ggthemes)
library(tidyverse)
library(purrr)
library(lubridate)
library(stringr)
library(plotrix)
library(broom)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)

# load data ---------------------------------------------------------------

TT_raw <- read_csv("Tetraselmis_experiment/data-processed/TT_cells_round2.csv")
TT_raw_round2 <- read_csv("Tetraselmis_experiment/data-processed/TT.csv")
TT_raw_round2_extremes <- read_csv("Tetraselmis_experiment/data-processed/TT_cells_round3_extremes.csv")
cells <- read_csv("Tetraselmis_experiment/data-processed/all_cell_densities.csv")
cells_extremes <- read_csv("Tetraselmis_experiment/data-processed/cells_extremes.csv")

## ok let's try with just one temperature first

cells_extremes_sel <- cells_extremes %>% 
	mutate(temp = ifelse(temp == "o", 0, temp)) %>% 
	mutate(temp = as.numeric(temp)) %>% 
	select(temp, variability, sample_day, replicate, cell_density, time_since_innoc_hours)


cells_sel <- cells %>% 
	select(temp, variability, sample_day, replicate, cell_density, time_since_innoc_hours)
	

cells_all <- bind_rows(cells_sel, cells_extremes_sel) %>% 
	filter(variability == "c")

cells_v <- bind_rows(cells_sel, cells_extremes_sel) %>% 
	filter(variability == "v")

cells_all %>% 
	ggplot(aes(x = time_since_innoc_hours, y = cell_density)) + geom_point() +
	facet_wrap( ~ temp)

cells_v %>% 
	ggplot(aes(x = time_since_innoc_hours, y = cell_density)) + geom_point() +
	facet_wrap( ~ temp)

cells_all %>% 
	# filter(temp == 27) %>% 
	# filter(time_since_innoc_hours < 175) %>% 
	mutate(keep = "yes") %>% 
	mutate(keep = ifelse(temp == 20 & time_since_innoc_hours > 52, "no", keep)) %>% 
	mutate(keep = ifelse(temp == 24 & time_since_innoc_hours > 50, "no", keep)) %>% 
	mutate(keep = ifelse(temp == 29 & time_since_innoc_hours > 60, "no", keep)) %>% 
	mutate(keep = ifelse(temp == 5 & time_since_innoc_hours > 175, "no", keep)) %>% 
	filter(keep != "no") %>% 
	ggplot(aes(x = time_since_innoc_hours, y = cell_density)) + geom_point() +
	facet_wrap( ~ temp)


cells_exp <- cells_all %>% 
	# filter(temp == 27) %>% 
	# filter(time_since_innoc_hours < 175) %>% 
	mutate(keep = "yes") %>% 
	mutate(keep = ifelse(temp == 20 & time_since_innoc_hours > 52, "no", keep)) %>% 
	mutate(keep = ifelse(temp == 24 & time_since_innoc_hours > 50, "no", keep)) %>% 
	mutate(keep = ifelse(temp == 29 & time_since_innoc_hours > 60, "no", keep)) %>% 
	mutate(keep = ifelse(temp == 5 & time_since_innoc_hours > 175, "no", keep)) %>% 
	filter(keep != "no")

write_csv(cells_exp, "Tetraselmis_experiment/data-processed/cells_exp.csv")

### now let's pick out the exponential phase from the variable

cells_v %>% 
	ggplot(aes(x = time_since_innoc_hours, y = cell_density)) + geom_point() +
	facet_wrap( ~ temp)


cells_v_exp <- cells_v %>% 
	# filter(temp == 15) %>% 
	# filter(time_since_innoc_hours < 175) %>% 
	mutate(keep = "yes") %>% 
	mutate(keep = ifelse(temp == 20 & time_since_innoc_hours > 51, "no", keep)) %>% 
	mutate(keep = ifelse(temp == 24 & time_since_innoc_hours > 51, "no", keep)) %>% 
	# mutate(keep = ifelse(temp == 27 & time_since_innoc_hours > 50, "no", keep)) %>% 
	mutate(keep = ifelse(temp == 5 & time_since_innoc_hours > 175, "no", keep)) %>% 
	filter(keep != "no") 

write_csv(cells_v_exp, "Tetraselmis_experiment/data-processed/cells_v_exp.csv")

estimate_growth <- function(x, temperature) {
	cells_exp %>% 
	filter(temp == temperature) %>% 
		mutate(time_point = trunc(time_since_innoc_hours)) %>% 
		group_by(time_point) %>% 
	sample_n(size = x, replace = FALSE) %>% 
	group_by(temp) %>% 
	do(tidy(nls(cell_density ~ 800 * (1+a)^(time_since_innoc_hours),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>%
	mutate(temp = as.numeric(temp)) %>%
	mutate(growth_per_day = estimate*24) %>% 
	mutate(error = std.error*24)
}


samples_rep <- rep(1, 1000)
growth_0 <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 0, .id = "run") 
growth_5 <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 5, .id = "run") 
growth_10 <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 10, .id = "run") 
growth_16 <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 16, .id = "run") 
growth_20 <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 20, .id = "run") 
growth_24 <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 24, .id = "run") 
growth_27 <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 27, .id = "run") 
growth_29 <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 29, .id = "run") 
growth_32 <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 32, .id = "run") 


### bring in growth_0 from round 2,
# growth_0 <- read_csv("Tetraselmis_experiment/data-processed/growth_0_round2.csv")
growth_all <- bind_rows(growth_16, growth_10, growth_5, growth_20, growth_24, growth_27, growth_29, growth_32, growth_0)

growth_all %>% 
	group_by(temp) %>% 
	summarise_each(funs(mean, std.error), growth_per_day) %>% 
	ggplot(aes(x = temp, y = growth_per_day_mean)) + geom_point() +
	geom_errorbar(aes(ymin = growth_per_day_mean - growth_per_day_std.error, ymax = growth_per_day_mean + growth_per_day_std.error), width = 0.1)


growth_all %>% 
	group_by(temp) %>% 
	summarise(upper = quantile(growth_per_day, probs = 0.025),
						lower = quantile(growth_per_day, probs = 0.975),
						mean = mean(growth_per_day)) %>% 
	ggplot(aes(x = temp, y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)

write_csv(growth_all, "Tetraselmis_experiment/data-processed/growth_resampling.csv")

growth_all <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling.csv")
### now do the resampling for the variable 


estimate_growth <- function(x, temperature) {
	cells_v_exp %>% 
		filter(temp == temperature) %>% 
		mutate(time_point = trunc(time_since_innoc_hours)) %>% 
		group_by(time_point) %>% 
		sample_n(size = x, replace = FALSE) %>% 
		group_by(temp) %>% 
		do(tidy(nls(cell_density ~ 800 * (1+a)^(time_since_innoc_hours),
								data= .,  start=list(a=0.01),
								control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
		ungroup() %>%
		mutate(temp = as.numeric(temp)) %>%
		mutate(growth_per_day = estimate*24) %>% 
		mutate(error = std.error*24)
}


samples_rep <- rep(1, 1000)
growth_5v <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 5, .id = "run") 
growth_10v <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 10, .id = "run") 
growth_15v <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 15, .id = "run") 
growth_20v <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 20, .id = "run") 
growth_24v <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 24, .id = "run") 
growth_27v <- samples_rep %>% 
	map_df(.f = estimate_growth, temperature = 27, .id = "run") 

growth_all_v <- bind_rows(growth_15v, growth_10v, growth_5v, growth_20v, growth_24v, growth_27v)	

write_csv(growth_all_v, "Tetraselmis_experiment/data-processed/growth_resampling_v.csv")

growth_all_v %>% 
  group_by(temp) %>% 
  summarise_each(funs(mean), growth_per_day) %>% 
  ggplot(aes(x = temp, y = growth_per_day)) + geom_point()
