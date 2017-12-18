### Initial plotting

# load libraries ----------------------------------------------------------

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
TT_0 <- read_csv("Tetraselmis_experiment/data-processed/TT_cells_round2_0.csv")

TT <- TT_raw %>% 
	mutate(start_time = ymd_hms(start_time)) %>% 
	mutate(total_biovolume = cell_density * cell_volume) %>% 
	unite(unique, temp, variability, replicate, sample_day, sep = "-", remove = FALSE)


TT_round2 <- TT_raw_round2 %>% 
	separate(replicate, into = c("temp", "variability"), sep = -2, remove = FALSE) %>%
	mutate(cell_density = as.numeric(cell_density)) %>% 
	mutate(start_time = ymd_hms(start_time)) %>% 
	mutate(total_biovolume = cell_density * cell_volume) %>% 
	unite(unique, temp, variability, replicate, sep = "-", remove = FALSE)



### now get growth rates

TT$start.time <- ymd_hms("2017-06-07 14:46:11 UTC")

TT$time_since_innoc <- interval(TT$start.time, TT$start_time)


TT2 <- TT %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))

TT_round2$start.time <- ymd_hms("2017-06-16 11:06:11 UTC")

TT_round2$time_since_innoc <- interval(TT_round2$start.time, TT_round2$start_time)
TT2_round2 <- TT_round2 %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))






TT3 <- TT2 %>% 
	mutate(keep = NA) %>% 
	mutate(keep = ifelse(temp %in% c("24", "27") & start_time > "2017-06-10", "drop", "keep")) %>% 
	mutate(keep = ifelse(temp == "10" & variability == "v" & start_time > "2017-06-12", "drop", keep)) %>%
	mutate(keep = ifelse(temp == "10" & variability == "c" & start_time > "2017-06-13", "drop", keep)) %>%
	mutate(keep = ifelse(temp == "20" & start_time > "2017-06-10", "drop", keep)) %>%
	mutate(keep = ifelse(temp == "16" & start_time > "2017-06-11", "drop", keep)) %>%
	select(keep, everything()) %>% 
	select(-time_since_innoc) %>% 
	filter(keep == "keep")

TT3 %>% 
	group_by(temp) %>% 
	# filter(temp == 20) %>% 
	filter(unique != "20-c-4-2") %>% 
	# filter(variability == "v") %>% 
	ggplot(aes(x = start_time, y = cell_density, group = unique, color = factor(variability))) + geom_point(size = 3) +
	facet_wrap( ~ temp)

TT3$rows <- as.numeric(rownames(TT3))

TT4 <- TT3 %>% 
		mutate(keep = ifelse(rows > 328 & rows < 333, "drop", keep)) %>% ## get rid of second sampling point on Jun 11
	filter(keep == "keep")

TT4 %>% 
	group_by(temp) %>% 
	filter(unique != "20-c-4-2") %>% 
	ggplot(aes(x = start_time, y = cell_density, group = unique, color = factor(variability))) + geom_point(size = 3) +
	facet_wrap( ~ temp) + theme_bw() + ylab("cell density") + xlab("date")

## rough cut of the growth rate estimates
round1_growth <- TT4 %>% 
	group_by(temp, variability) %>% 
	do(tidy(nls(cell_density ~ 800 * (1+a)^(time_since_innoc_hours),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>%
	mutate(temp = as.numeric(temp)) %>%
	mutate(growth_per_day = estimate*24) %>% 
	mutate(error = std.error*24) %>% 
	# filter(temp > 6) %>% 
	ggplot(aes(x = temp, y = growth_per_day, color = variability)) + geom_point(size = 2, alpha = 0.5) +
	geom_errorbar(aes(ymin = growth_per_day - error, ymax = growth_per_day + error), width = 0.1) + 
	theme_bw() + geom_smooth() + ylab("intrinsic rate of increase (r) / day") + xlab("temperature (C)")


round2_growth <- TT2_round2 %>% 
	group_by(temp, variability) %>% 
	do(tidy(nls(cell_density ~ 800 * (1+a)^(time_since_innoc_hours),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>%
	mutate(temp = as.numeric(temp)) %>%
	mutate(growth_per_day = estimate*24) %>% 
	mutate(error = std.error*24) %>% 
	filter(temp != "10") %>% 
	filter(temp != "27") %>% 
	filter(temp != "5")


round1_growth_1 <- round1_growth %>% 
	filter(temp != "20") %>% 
	filter(temp != "24")


TT5 <- TT4 %>% 
	select(time_since_innoc_hours, temp, replicate, variability, cell_density)
TT6 <- TT2_round2 %>% 
	filter(temp != "10") %>% 
	filter(temp != "27") %>% 
	filter(temp != "5") %>% 
	mutate(temp = as.numeric(temp)) %>% 
	select(time_since_innoc_hours, temp, variability, cell_density)

exp_growth <- TT6 %>% 
	filter(temp == 15) %>% 
	ggplot() + geom_point(size = 4, color = "blue", alpha = 0.5, aes(x = time_since_innoc_hours, y = cell_density,  frame = time_since_innoc_hours, cumulative = TRUE)) +
	theme_classic() + ylab("Abundance (cells/mL)") + xlab("Time (hours)") +
	theme(text = element_text(size=20, family = "Helvetica")) 
library(gganimate)
library(animation)
ani.options(interval = 0.2, loop = 1)
gganimate(exp_growth, interval = 0.2, filename = "Tetraselmis_experiment/figures/exponential_growth.mp4", title_frame = FALSE, ani.width = 600, ani.height = 400, loop = FALSE)


all_densities <- bind_rows(TT5, TT6)
all_densities2 <- bind_rows(TT4, TT6)


write_csv(all_densities2, "Tetraselmis_experiment/data-processed/all_cell_densities.csv")


all_densities2 %>% 
	# filter(temp == 24) %>% 
	ggplot(aes(x = time_since_innoc_hours, y = cell_density, color = variability)) + geom_point(alpha = 0.7) + 
	facet_wrap( ~ temp) + theme_bw() + ylab("abundbance") + xlab("time (hours)")


all_growth <- all_densities2 %>% 
	mutate(keep = NA) %>% 
	mutate(keep = ifelse(temp > 19 & time_since_innoc_hours > 60, "drop", "keep")) %>% 
	# filter(keep == "keep") %>%
	# ggplot(aes(x = time_since_innoc_hours, y = cell_density, color = variability)) + geom_point() +
	# facet_wrap( ~ temp)
	group_by(temp, variability) %>% 
	do(tidy(nls(cell_density ~ 800 * (1+a)^(time_since_innoc_hours),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>%
	mutate(temp = as.numeric(temp)) %>%
	mutate(growth_per_day = estimate*24) %>% 
	mutate(error = std.error*24) 
	# filter(temp > 6) %>% 
	ggplot(aes(x = temp, y = growth_per_day, color = variability)) + geom_point(size = 2, alpha = 0.5) +
	geom_line() +
	geom_errorbar(aes(ymin = growth_per_day - error, ymax = growth_per_day + error), width = 0.1) + 
	theme_bw() +
	# geom_smooth() +
	ylab("intrinsic rate of increase (r) / day") + xlab("temperature (C)")




extremes <- read_csv("Tetraselmis_experiment/data-processed/all_r_with0.csv")
extremes2 <- extremes %>% 
	filter(temp == 0 | temp == 35) %>% 
	filter(var == "constant temperature") %>% 
	rename(variability = var) %>% 
	mutate(variability = str_replace(variability, "constant temperature", "c"))
	
growth_all <- bind_rows(all_growth, extremes2) %>% 
	mutate(growth_per_day = ifelse(is.na(growth_per_day), estimate, growth_per_day))


growth_all %>% 
	ggplot(aes(x = temp, y = growth_per_day, color = variability)) + geom_point(size = 2, alpha = 0.5) +
	geom_errorbar(aes(ymin = growth_per_day - error, ymax = growth_per_day + error), width = 0.1) + 
	geom_line() +
	theme_bw() +
	ylab("intrinsic rate of increase (r) / day") + xlab("temperature (C)")



## add predictions from the constant temperatures

View(growth_all)

prediction5V <- 0.5*((growth_all$growth_per_day[growth_all$temp == 0 & growth_all$variability == "c"][[1]]) + (growth_all$growth_per_day[growth_all$temp == 10 & growth_all$variability == "c"][[1]]))
prediction15V <- 0.5*((growth_all$growth_per_day[growth_all$temp == 10 & growth_all$variability == "c"][[1]]) + (growth_all$growth_per_day[growth_all$temp == 20 & growth_all$variability == "c"][[1]]))
prediction15V <- 0.5*((growth_all$growth_per_day[growth_all$temp == 10 & growth_all$variability == "c"][[1]]) + (growth_all$growth_per_day[growth_all$temp == 20 & growth_all$variability == "c"][[1]]))
prediction10V <- 0.5*((growth_all$growth_per_day[growth_all$temp == 5 & growth_all$variability == "c"][[1]]) + (growth_all$growth_per_day[growth_all$temp == 16 & growth_all$variability == "c"][[1]]))


temp <- 5

prediction5Vdf <- data.frame(prediction5V, temp) %>% 
	mutate(treatment = "prediction") %>% 
	rename(growth_per_day = prediction5V)
temp15 <- 15
prediction15Vdf <- data.frame(prediction15V, temp15) %>% 
	mutate(treatment = "prediction") %>% 
	rename(growth_per_day = prediction15V) %>% 
	rename(temp = temp15)


temp10 <- 10
prediction10Vdf <- data.frame(prediction10V, temp10) %>% 
	mutate(treatment = "prediction") %>% 
	rename(growth_per_day = prediction10V) %>% 
	rename(temp = temp10)

write_csv(growth_all, "/Users/student/Downloads/Thermal_variability-master/Tetraselmis_experiment/data-processed/growth_estimates_round3.csv")
growth_all <- read_csv("Tetraselmis_experiment/data-processed/growth_estimates_round3.csv")

growth_all2 <- growth_all %>% 
	mutate(treatment = "experiment")

growth_all3 <- bind_rows(growth_all2, prediction5Vdf, prediction15Vdf, prediction10Vdf)

growth_all3 %>% 
	ggplot(aes(x = temp, y = growth_per_day, color = variability)) + geom_point(size = 3, alpha = 0.5) +
	geom_errorbar(aes(ymin = growth_per_day - error, ymax = growth_per_day + error), width = 0.1) + 
	geom_line() +
	theme_bw() +
	ylab("intrinsic rate of increase (r) / day") + xlab("temperature (C)") + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(text=element_text(family="Helvetica", size=16)) +
	geom_rect(aes(xmin=0, xmax=14.72, ymin=-Inf, ymax=Inf), alpha = 0.02, colour = NA) +
	annotate("text", x = 7, y = 1.7, label = "predicted positive \n effect of variability") +
	annotate("text", x = 21, y = 1.7, label = "predicted negative \n effect of variability")
ggsave("Tetraselmis_experiment/figures/growth_round3.pdf")


### what is the activation energy of growth rate?

growth_all3 %>% 
	filter(variability == "c") %>% 
	filter(temp < 28, temp > 0) %>% 
	mutate(inverse_temp = (1/(.00008617*(temp+273.15)))) %>%
	# ggplot(aes(x = inverse_temp, y = log(estimate))) + geom_point() + geom_smooth(method = "lm")
	# ungroup() %>% 
	do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View




### Next step is to fit the TPC from Mridul's code

### now let's make that 27C plot


T27 <- TT2 %>% 
	filter(temp == 27) %>% 
	select(-time_since_innoc) 

### hacky way of adding background shading
T27 %>% 
filter(variability == "v") %>% 
	ggplot(aes(x = time_since_innoc_hours, y = cell_density, group = unique)) + geom_point(size = 3) + 
theme_bw() +
	annotate("rect", xmin = 7, xmax = (7+12), ymin = 0, ymax = 29000,
											alpha = .2, fill = "blue") +
	annotate("rect", xmin = 19, xmax = (19+12), ymin = 0, ymax = 29000,
					 alpha = .2, fill = "red") +
	annotate("rect", xmin = 31, xmax = (31+12), ymin = 0, ymax = 29000,
					 alpha = .2, fill = "blue") +
	annotate("rect", xmin = 43, xmax = (43+12), ymin = 0, ymax = 29000,
					 alpha = .2, fill = "red")+
	annotate("rect", xmin = 55, xmax = (55+12), ymin = 0, ymax = 29000,
					 alpha = .2, fill = "blue") +
	annotate("rect", xmin = 67, xmax = (67+12), ymin = 0, ymax = 29000,
					 alpha = .2, fill = "red") +
	annotate("rect", xmin = 79, xmax = (79+12), ymin = 0, ymax = 29000,
					 alpha = .2, fill = "blue") +
	annotate("rect", xmin = 91, xmax = (91+12), ymin = 0, ymax = 29000,
					 alpha = .2, fill = "red") +
	ylab("cell density, cells/ml") + xlab("time, hours")



T27 %>% 
	mutate(var_temp_night = NA) %>% 
	mutate(time_of_day_minute = minute(start_time)) %>% 
	mutate(time_of_day_hour = hour(start_time)) %>%
	select(time_of_day_minute, time_of_day_hour, everything()) %>% 
	mutate(var_temp_night = ifelse(time_of_day_hour >= 21 | time_of_day_hour < 10, "cold", "hot")) %>%
	select(var_temp_night, everything()) %>% View
	
	
#### now let's get it in the right shape to fit the logistic model
	
	
TTg <- TT2 %>% 
	mutate(hours = floor(time_since_innoc_hours)) 

write_csv(TTg, "Tetraselmis_experiment/data-processed/TTg.csv")


TT_raw_round2_extremes

TT_raw_round2_extremes$start.time <- ymd_hms("2017-06-23 15:07:07 UTC")

TT_raw_round2_extremes$time_since_innoc <- interval(TT_raw_round2_extremes$start.time, TT_raw_round2_extremes$start_time)


growth_32 <- TT_raw_round2_extremes %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1)) %>% 
	filter(temp == 32) %>% 
	group_by(temp) %>% 
	do(tidy(nls(cell_density ~ 800 * (1+a)^(time_since_innoc_hours),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temp = as.numeric(temp)) %>%
	mutate(growth_per_day = estimate*24) %>% 
	mutate(error = std.error*24) 

write_csv(growth_32, "Tetraselmis_experiment/data-processed/growth_32.csv")


cells_extremes <- TT_raw_round2_extremes %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))
	
write_csv(cells_extremes, "Tetraselmis_experiment/data-processed/cells_extremes.csv")


### now get growth rates for 0

TT_0$start.time <- ymd_hms("2017-03-18 16:17:46 UTC")

TT_0$time_since_innoc <- interval(TT_0$start.time, TT_0$start_time)



cells_0 <- TT_0 %>%
	select(-variability) %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))


cells_0 %>%
	group_by(temperature) %>% 
		do(tidy(nls(cell_density ~ 800 * (1+a)^(time_since_innoc_hours),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>%
	mutate(temp = as.numeric(temperature)) %>%
	mutate(growth_per_day = estimate*24) %>% 
	mutate(error = std.error*24)

cells_0 %>% 
	filter(time_since_innoc_hours < 5) %>% 
	summarise(mean_density = mean(cell_density))


estimate_growth <- function(x, temperature) {
	cells_0 %>% 
		rename(temp = temperature) %>% 
		filter(temp == temperature) %>% 
		mutate(time_point = trunc(time_since_innoc_hours)) %>% 
		group_by(time_point) %>% 
		sample_n(size = x, replace = FALSE) %>% 
		group_by(temp) %>% 
		do(tidy(nls(cell_density ~ 1000 * (1+a)^(time_since_innoc_hours),
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

growth_0b <- growth_0 %>% 
	distinct(growth_per_day, .keep_all = TRUE)

write_csv(growth_0, "Tetraselmis_experiment/data-processed/growth_0_round2.csv")
