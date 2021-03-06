
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

cells_days_c <- read_csv("Tetraselmis_experiment/data-processed/cells_exp_mod.csv") %>% 
  mutate(days = time_since_innoc_hours/24) %>% 
  mutate(days = ifelse(days < 0, 0, days))
cells_days_v <- read_csv("Tetraselmis_experiment/data-processed/cells_days_v_mod.csv")


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

cells_c %>% 
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
cells_exp <- read_csv("Tetraselmis_experiment/data-processed/cells_exp.csv")
### now let's pick out the exponential phase from the variable

cells_days_c %>% 
  filter(temp == 5) %>% 
	ggplot(aes(x = days, y = cell_density)) + geom_point() +
	facet_wrap( ~ temp) + stat_function(fun = function(x) 800*exp(0.28956 *x))
1.130026

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
	cells_days_c %>% 
	filter(temp == temperature) %>% 
		# mutate(time_point = trunc(time_since_innoc_hours)) %>% 
		group_by(sample_group) %>% 
	sample_n(size = x, replace = FALSE) %>% 
	group_by(temp) %>% 
	do(tidy(nls(cell_density ~ 800 * exp(r*days),
							data= .,  start=list(r=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>%
	mutate(temp = as.numeric(temp))
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
	summarise_each(funs(mean, std.error), estimate) %>%
	ggplot(aes(x = temp, y = estimate_mean)) + geom_point() +
	geom_errorbar(aes(ymin = estimate_mean - estimate_std.error, ymax = estimate_mean + estimate_std.error), width = 0.1)


growth_all %>% 
	group_by(temp) %>% 
	summarise(upper = quantile(growth_per_day, probs = 0.025),
						lower = quantile(growth_per_day, probs = 0.975),
						mean = mean(growth_per_day)) %>% 
	ggplot(aes(x = temp, y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1)

write_csv(growth_all, "Tetraselmis_experiment/data-processed/growth_resampling_exp.csv")

growth_all <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling_exp.csv")
### now do the resampling for the variable 


estimate_growth <- function(x, temperature) {
  cells_days_v %>% 
    filter(temp == temperature) %>% 
    # mutate(time_point = trunc(time_since_innoc_hours)) %>% 
    group_by(sample_group) %>% 
    sample_n(size = x, replace = FALSE) %>% 
    group_by(temp) %>% 
    do(tidy(nls(cell_density ~ 800 * exp(r*days),
                data= .,  start=list(r=0.01),
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

write_csv(growth_all_v, "Tetraselmis_experiment/data-processed/growth_resampling_v_exp.csv")

growth_all_v %>% 
  group_by(temp) %>% 
  summarise_each(funs(mean), growth_per_day) %>% 
  ggplot(aes(x = temp, y = growth_per_day)) + geom_point()


### now estimating growth rates and bootstrapping with nlsBoot
library(nlstools)
?nlsBoot

## temp = 0
temperature <- 0
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_c, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_0 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)


## temp = 5

temperature <- 5
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_c, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_5 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)


## temp = 10
temperature <- 10
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_c, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_10 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)

## temp = 16
temperature <- 16
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_c, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_16 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)



## temp = 20
temperature <- 20
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_c, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_20 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)


## temp = 24
temperature <- 24
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_c, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_24 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)


## temp = 27
temperature <- 27
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_c, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_27 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)


## temp = 29
temperature <- 29
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_c, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_29 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)


## temp = 32
temperature <- 32
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_c, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_32 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)

all_estimates <- bind_rows(output_0, output_5, output_10, output_16, output_20, output_24, output_27, output_29, output_32)
write_csv(all_estimates, "Tetraselmis_experiment/data-processed/growth_estimates_boot_car.csv")

### now onto variable
## temp = 5
temperature <- 5
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_v, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_v5 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)

## temp = 10
temperature <- 10
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_v, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_v10 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)

## temp = 15
temperature <- 15
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_v, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_v15 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)


## temp = 20
temperature <- 20
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_v, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_v20 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)


## temp = 24
temperature <- 24
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_v, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_v24 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)


## temp = 27
temperature <- 27
growth <- nls(cell_density ~ 800 * exp(r*days),
              data= filter(cells_days_v, temp == temperature),  start=list(r=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))
boot_growth <- as_data_frame((bootCase(growth, B = 999)))
lower <-  quantile(boot_growth$r, probs=0.025)
upper = quantile(boot_growth$r, probs=0.975)
cis <- as_data_frame(cbind(lower, upper))
tidy_growth <- tidy(growth, conf.int = TRUE)
output_v27 <- bind_cols(cis, tidy_growth) %>% 
  mutate(temperature = temperature)

all_estimates_v <- bind_rows(output_v5, output_v10, output_v15, output_v20, output_v24, output_v27)
write_csv(all_estimates_v, "Tetraselmis_experiment/data-processed/growth_estimates_boot_car_v.csv")


### make this into a function!

estimate_growth_function <- function(df){
  growth <- nls(cell_density ~ 800 * exp(r*days),
                  data= df,  start=list(r=0.01),
                  control = nls.control(maxiter=100, minFactor=1/204800000))
  boot_growth <- as_data_frame((bootCase(growth, B = 999)))
  
  lower <-  quantile(boot_growth$r, probs=0.025)
  upper = quantile(boot_growth$r, probs=0.975)
  cis <- as_data_frame(cbind(lower, upper))
  tidy_growth <- tidy(growth, conf.int = TRUE)
  output <- bind_cols(cis, tidy_growth) %>% 
    mutate(temperature = df$temp[[1]])
  return(output)
}


cells_c_split <- cells_days_c %>% 
  split(.$temp)

cells_v_split <- cells_days_v %>% 
  split(.$temp)

all_growth_estimated <- cells_c_split %>% 
  map_df(estimate_growth_function, .id = "temperature_rep") %>% 
  mutate(temp_num = as.numeric(temperature))

all_growth_estimated_v <- cells_v_split %>% 
  map_df(estimate_growth_function, .id = "temperature_rep") %>% 
  mutate(temp_num = as.numeric(temperature))


## looking at cell sizes

cells %>% 
  filter(variability == "c") %>% 
  ggplot(aes(x = time_since_innoc_days, y = cell_density)) + geom_point() +
  facet_wrap(~ temp, scales = "free")

cells %>% 
  filter(variability == "c") %>% 
  ggplot(aes(x = time_since_innoc_days, y = cell_volume)) + geom_point() +
  facet_wrap(~ temp, scales = "free")

cells %>% 
  group_by(temp) %>% 
  filter(cell_volume < 3000, temp < 27) %>% 
  ungroup() %>% 
ggplot(aes(x = temp, y = cell_volume)) + geom_point() +
  geom_smooth(method = "lm")

cells %>% 
  group_by(temp, replicate) %>% 
  filter(cell_volume < 3000, temp < 24) %>% 
  summarise(mean_size = mean(cell_volume)) %>% 
  lm(mean_size ~ temp, data = .) %>% 
  tidy(., conf.int =  TRUE) 


cells %>% 
  group_by(temp, replicate) %>% 
  filter(cell_volume < 3000, temp < 24) %>% 
  summarise(mean_size = mean(cell_volume)) %>%
  ungroup() %>% 
  ggplot(aes(x = temp, y = mean_size)) + geom_point() +
  geom_smooth(method = "lm")

cells %>% 
  filter(variability == "c") %>% 
  ggplot(aes(x = time_since_innoc_days, y = cell_density)) + geom_point() + 
  facet_wrap( ~ temp)
