### Now looking at the daily temperature data!


library(tidyverse)
library(lubridate)
library(viridis)


## bring in data

temp1 <- read_csv("Tetraselmis_experiment/data-processed/daily_temperatures_1-97.csv")
temp2 <- read_csv("Tetraselmis_experiment/data-processed/daily_temperatures_above99.csv")
temp3 <- read_csv("Tetraselmis_experiment/data-processed/temperatures_missing_round0.25.csv")
temp4 <- read_csv("Tetraselmis_experiment/data-processed/temperatures_still_missing.csv")
temp5 <- read_csv("Tetraselmis_experiment/data-processed/temepratures_last_time.csv")
missing <- read_csv("Tetraselmis_experiment/data-processed/missing_isolates_NOAA.csv")
results5 <- read_csv("Tetraselmis_experiment/data-processed/results5.csv")



ts <- bind_rows(temp1, temp2)
identical(unique(ts$isolate.code), results5$isolate.code)




results5 %>% 
  filter(mu.rsqrlist > 0.85) %>% View

length(unique(temp1$isolate.code))
length(unique(temp2$isolate.code))
intersect(unique(temp1$isolate.code), unique(temp2$isolate.code))


close_temps <- bind_rows(temp1, temp2) %>% 
	filter(!is.na(sst)) %>% 
	distinct(isolate.code)

length(setdiff(results5$isolate.code, close_temps$isolate.code)) ### ok looks like we are missing 23 isolates


temp4_isolates <- temp4 %>% 
	filter(!is.na(sst)) %>% 
	distinct(isolate.code)

temp5_isolates <- temp5 %>% 
	filter(!is.na(sst)) %>% 
	distinct(isolate.code)

temp3_isolates <- temp3 %>% 
	filter(!is.na(sst)) %>% 
	distinct(isolate.code)

temp2_isolates <- temp2 %>% 
	filter(!is.na(sst)) %>% 
	distinct(isolate.code)

temp1_isolates <- temp1 %>% 
	filter(!is.na(sst)) %>% 
	distinct(isolate.code)

all_isolates<- bind_rows(temp3_isolates, temp4_isolates, temp5_isolates)

setdiff(missing$isolate.code, all_isolates$isolate.code) ## ok it looks like isolates 344, 345, 346 are missing. but that's ok for now.


all_temps <- bind_rows(temp1, temp2, temp3, temp4, temp5)


otemps <- bind_rows(temp1, temp2)

otemps_sd <- otemps %>% 
	mutate(month = month(time)) %>% 
	filter(!is.na(sst)) %>% 
	group_by(isolate.code, month, time) %>% 
	summarise(temp = mean(sst)) %>%
	group_by(isolate.code, month) %>% 
	summarise_each(funs(mean, sd), temp) 

otemps_id <- otemps %>% 
	distinct(isolate.code, .keep_all = TRUE) %>% 
	select(isolate.code, lat, lon)

otemps_annual <- otemps %>% 
	group_by(isolate.code) %>% 
	summarise(annual_mean = mean(sst)) %>% 
	filter(!is.na(annual_mean))


otemps_lat <- left_join(otemps_sd, otemps_id)

otemps_all <- left_join(otemps_lat, otemps_annual)

otemps_all %>% 
	mutate(hemisphere = ifelse(lat > 0, "northern", "southern")) %>% 
	mutate(absolute_latitude = abs(lat)) %>% 
	ggplot(aes(x = absolute_latitude, y = temp_sd, color = annual_mean)) + geom_point(size = 4) + theme_classic() + 
	ylab("SST SD within each month") +
	scale_color_viridis(option = "inferno") +
	geom_point(size = 4, shape = 1, color = "black") +
	facet_wrap( ~ hemisphere) + xlab("absolute latitude")
	
	


str(otemps)


temp_summary <- all_temps %>% 
	filter(!is.na(sst)) %>% 
	group_by(isolate.code, time) %>% 
	summarise(mean_temp = mean(sst)) %>% 
	group_by(isolate.code) %>% 
	summarise_each(funs(mean, sd), mean_temp) 

results6 <- left_join(results5, temp_summary, by = "isolate.code")


results6 %>% 
	filter(minqual == "good", maxqual == "good") %>% 
	filter(mu.n > 4) %>% 
	# filter(isolate.code %in% greater_90) %>% 
	# filter(rel.curveskew < 0) %>% 
	ggplot(aes(x = rel.curveskew, y = diff, color = mean_temp_sd, label = isolate.code)) + geom_point(size = 2) +
	# geom_label()
	scale_color_viridis(option = "inferno") + theme_bw() + geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
	geom_hline(yintercept = 0) + ylab("Topt(cons) - Topt(var) (°C)") +
	xlab("Curve skewness") + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	geom_vline(xintercept = 0)

# make new predictions with daily SD --------------------------------------

derivative <- function(f, x, ..., order = i, delta = 0.1, sig = 6) {
	# Numerically computes the specified order derivative of f at x
	vals <- matrix(NA, nrow = order + 1, ncol = order + 1)
	grid <- seq(x - delta/2, x + delta/2, length.out = order + 1)
	vals[1, ] <- sapply(grid, f, ...) - f(x, ...)
	for (i in 2:(order + 1)) {
		for (j in 1:(order - i + 2)) {
			stepsize <- grid[i + j - 1] - grid[i + j - 2]
			vals[i, j] <- (vals[i - 1, j + 1] - vals[i - 1, j])/stepsize
		}
	}
	return(signif(vals[order + 1, 1], sig))
}

tpc<-function(x){
	res<-all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)
	res
}

## step 2 
variable_predictions <- function(x) {
	y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(all$daily_SD^2)
}


predict_function <- function(data) {
	all <- data
	x <- seq(0, 45, by = 0.1)
	tpc<-function(x){
		res<-all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)
		res
	}
	## step 2 
	variable_predictions <- function(x) {
		y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(all$daily_SD^2)
	}
	
	predicted_growth_variable <- sapply(x, variable_predictions)
	predicted_growth_variable2 <- data.frame(x, predicted_growth_variable) %>% 
		rename(temperature = x, 
					 growth.rate = predicted_growth_variable) %>% 
		top_n(n = 1, wt = growth.rate)
	# data.frame(all$isolate.code, predicted_growth_variable2)
}


breadth_function <- function(data) {
	all <- data
	x <- seq(0, 45, by = 0.1)
	tpc<-function(x){
		res<-all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)
		res
	}
	## step 2 
	variable_predictions <- function(x) {
		y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(all$daily_SD^2)
	}
	
	predicted_growth_variable <- sapply(x, variable_predictions)
	predicted_growth_variable2 <- data.frame(x, predicted_growth_variable) %>% 
		rename(temperature = x, 
					 growth.rate = predicted_growth_variable) %>% 
		filter(growth.rate > 0) %>% 
		top_n(n = -20, wt = growth.rate) %>% 
		summarise_each(funs(min, max), temperature)
	# data.frame(all$isolate.code, predicted_growth_variable2)
}


all_split <- results6 %>% 
	filter(!is.na(topt)) %>% 
	rename(daily_SD = mean_temp_sd) %>% 
	select(1:44, daily_SD) %>% 
	split(.$isolate.code)

results_dailySD <- all_split %>% 
	map_df(predict_function, .id = "isolate_code")


results_dailySD2 <- results_dailySD %>% 
	rename(topt_variable_daily = temperature,
				 rmax_daily = growth.rate,
				 isolate.code = isolate_code) %>% 
	mutate(isolate.code = as.integer(isolate.code))


results_all <- left_join(results6, results_dailySD2) %>% 
	rename(daily_SD = mean_temp_sd)

results_all %>% 
	filter(minqual == "good", maxqual == "good") %>% 
	filter(mu.n > 4) %>% 
	mutate(diff_daily = topt - topt_variable_daily) %>% 
	# filter(isolate.code %in% greater_90) %>% 
	# filter(rel.curveskew < 0) %>% 
	ggplot(aes(x = rel.curveskew, y = diff_daily, color = daily_SD, label = isolate.code)) + geom_point(size = 2) +
	# geom_label()
	scale_color_viridis(option = "inferno") + theme_bw() + geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
	geom_hline(yintercept = 0) + ylab("Topt(cons) - Topt(var) (°C)") +
	xlab("Curve skewness") + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	geom_vline(xintercept = 0)


results_all %>% 
	ggplot(aes(x = latitude, y = daily_SD)) + geom_point(color = "blue") +
	# geom_abline(slope = 1, intercept = 0) + 
	# geom_smooth(method = "lm") + theme_classic() +
	ylab("SST SD calculated with daily data") + xlab("latitude") +
	theme_classic()

hist(results_all$daily_SD)
hist(results_all$SD)

 all_temps %>% 
	filter(!is.na(sst)) %>% 
	group_by(isolate.code, time, lat) %>% 
	summarise(mean_temp = mean(sst)) %>% 
 	filter(isolate.code %in% close_temps$isolate.code) %>% 
 	mutate(hemisphere = ifelse(lat > 0, "northern hemisphere", "southern hemisphere")) %>% 
 	mutate(absolute_latitude = abs(lat)) %>% 
	ggplot(aes(x = time, y = mean_temp, color = absolute_latitude, group = isolate.code)) + geom_line() +
	facet_wrap( ~ hemisphere) + 
 	theme_classic() +
 	theme(strip.background = element_rect(colour="white", fill="white")) +
 	scale_color_viridis(option = "inferno") + ylab("SST") + xlab("Date")
