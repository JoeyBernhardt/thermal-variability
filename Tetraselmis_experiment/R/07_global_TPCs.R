## gaol for this script -- recreate Fig 2A in Thomas et al. 2012 Science,
## ask, how does predicted Topt including variability alter the match between
## Topt and mean annual temperature?


library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(fuzzyjoin)
library(purrr)
library(viridis)

thomas <- read_csv("data/thermal_trait_data/Thomas_2014_traits_derived_20140606.csv")
therm_raw <- read_csv("Tetraselmis_experiment/data-processed/global_therm.csv")

therm <- therm_raw %>% 
	filter(Biome == "marine")


### ok so what we need to do here is start by taking one of these curves, get the mean
### and sd of the temperature where it is, and calculate the new new Topt!


thomas %>% 
	filter(habitat != "freshwater") %>% View

hist(thomas$rel.curveskew)

## lat 28, long -79.270

thomas2 <- thomas %>% 
	rename(topt = mu.g.opt.list) %>% 
	rename(w = mu.wlist,
				 a = mu.alist,
				 z = mu.c.opt.list,
				 b = mu.blist) %>% 
	select(topt, w, a, b, z, everything()) 


therm2 <- therm %>% 
	rename(latitude = Lat, 
				 longitude = Lon)

thomas3 <- thomas2 %>% 
	rename(latitude = isolation.latitude,
				 longitude = isolation.longitude) %>% 
	filter(!is.na(latitude), !is.na(longitude)) %>% 
	filter(habitat %in% c("marine", "estuarine", "saltmarsh"))


therm_round <- therm %>% 
	mutate(latround = round(Lat, digits = 0)) %>%
	mutate(longround = round(Lon, digits = 0)) 

thomas_round <- thomas2 %>% 
	mutate(latround = round(isolation.latitude, digits = 0)) %>% 
	mutate(longround = round(isolation.longitude, digits = 0))

g <- geo_left_join(thomas3, therm2, max_dist = 200, method = "geo")

 thermal_sd <- g %>% 
 	select(starts_with("lat"), starts_with("long"), Mean, isolate.code, everything()) %>% 
 	mutate(lat_diff = abs(latitude.x - latitude.y),
 				 long_diff = abs(longitude.x - longitude.y)) %>% 
 	select(lat_diff, long_diff, everything()) %>%
 	mutate(total_diff = lat_diff + long_diff) %>% 
 	select(total_diff, everything()) %>% 
 	group_by(isolate.code) %>% 
 	top_n(wt = total_diff, n = -1) %>% 
 	select(SD, everything()) %>% 
 	group_by(isolate.code) %>% 
 	summarise_each(funs(mean), SD, Mean)

write_csv(thermal_sd, "Tetraselmis_experiment/data-processed/thermal_sd.csv")
write_csv(thermal_sd, "Tetraselmis_experiment/data-processed/thermal_sd_aug.csv")

thermal_sd <- read_csv("Tetraselmis_experiment/data-processed/thermal_sd_aug.csv")

all_thermal_data <- left_join(thermal_sd, thomas3, by = c("isolate.code"))




## let's plot one curve
## 592, 324, 614, 613, 317, 209, 538, 607...77 could work for the figures
curve <- all_thermal_data %>% 
	filter(isolate.code == 194)

tpc1<-function(x){
	res<-curve$a[1]*exp(curve$b[1]*x)*(1-((x-curve$z[1])/(curve$w[1]/2))^2)
	res
}

##c
# 24.5505
# 0.0653433
# 15.1353
# 33.9425
# 0.0092910
# 0.0944360
# 0.981972
# 0.00343878
# 10
tpc1<-function(x){
	res<-1.6*exp(-0.07*x)*(1-((x-23)/(36/2))^2)
	res
}

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
right_skew <- p + stat_function(fun = tpc1, color = "black", size = 2) +xlim(0, 40) + ylim(0, 0.5) + 
	theme_bw() + ylab("Growth rate") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line.y = element_line(color="black"),
				axis.line.x = element_line(color="black")) +
	theme(text = element_text(size=14, family = "Helvetica")) +
	xlab(expression("Temperature (" *degree * "C)"))
ggsave("Tetraselmis_experiment/figures/example_right_skew_TPC.pdf")


# derivative --------------------------------------------------------------

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



### ok let's generalize this!

## step 1

tpc<-function(x){
	res<-all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)
	res
}

## step 2 
variable_predictions <- function(x) {
	y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(all$SD^2)
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
		y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(all$SD^2)
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
		y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(all$SD^2)
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



all_split <- all_thermal_data %>% 
	filter(!is.na(topt)) %>% 
	split(.$isolate.code)

results_breadth <- all_split %>% 
	map_df(breadth_function, .id = "isolate_code")


write_csv(results_breadth, "Tetraselmis_experiment/data-processed/Topt_variable_global_breadth.csv")
write_csv(results, "Tetraselmis_experiment/data-processed/Topt_variable_global.csv")
results <- read_csv("Tetraselmis_experiment/data-processed/Topt_variable_global.csv")
results_breadth <- read_csv("Tetraselmis_experiment/data-processed/Topt_variable_global_breadth.csv")


all_variable <- left_join(results, results_breadth)

results2 <- all_variable %>% 
	mutate(isolate.code = as.integer(isolate_code))

results3 <- left_join(thomas3, results2, by = c("isolate.code"))

results4 <- left_join(results3, thermal_sd, by = c("isolate.code"))

results5 <- results4 %>% 
	rename(topt_variable = temperature) %>% 
	filter(curvequal == "good") %>% 
	mutate(diff = topt - topt_variable) %>% 
	filter(habitat %in% c("marine", "estuarine", "saltmarsh"))


results6 <- results4 %>% 
	rename(topt_variable = temperature) 

results5 %>%
	# filter(rel.curveskew<0) %>%
	mutate(skew_sign = ifelse(rel.curveskew > 0, 1, 2)) %>% 
	ggplot(aes(x = topt, y = Mean, color = rel.curveskew)) + geom_point(size = 2) +
	geom_abline(slope = 1, intercept = 0) + theme_bw()  + ylim(0, 40) + scale_color_viridis(option = "inferno") +
	xlim(0, 40) +ylab("Topt under variable conditions") + xlab("Topt under constant conditions")

results5 %>%
	# filter(rel.curveskew<0) %>% 
	rename(`curve skew` = rel.curveskew) %>% 
	ggplot(aes(x = topt, y = Mean, color = `curve skew`, size = SD)) + geom_point() +
	geom_abline(slope = 1, intercept = 0) + theme_bw()  + ylim(0, 30) + scale_color_viridis(option = "inferno") +
	xlim(0, 30) +ylab("Mean environmental temperature") + xlab("Topt under constant conditions") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=14, family = "Helvetica")) 
ggsave("Tetraselmis_experiment/figures/env_temp_v_Topt.pdf")

## skew and SD best explain variation in the effect of variability
m1 <- lm(diff ~ rel.curveskew + SD, data = results5)
m2 <- lm(diff ~ rel.curveskew, data = results5)
m3 <- lm(diff ~ SD, data = results5)
AIC(m1)
AIC(m2)
AIC(m3)
summary(m1)

mean(abs(results5$diff), na.rm = TRUE)
hist(results5$diff)

results5 %>% 
	ggplot(aes(x = SD, y = rel.curveskew)) + geom_point() + geom_smooth()

cor(use = "na.or.complete", x = abs(results5$latitude), y = results5$rel.curveskew)

topt_skew <- results5 %>% 
	# mutate(diff = topt_variable - topt) %>% 
	ggplot(aes(x = rel.curveskew, y = diff, color = SD, label = isolate.code)) + geom_point(size = 2) +
	# geom_label()
	scale_color_viridis(option = "inferno") + theme_bw() + geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
	geom_hline(yintercept = 0) + ylab("Topt(cons) - Topt(var) (°C)") +
	xlab("Curve skewness") + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	geom_vline(xintercept = 0)
ggsave("Tetraselmis_experiment/figures/variability_effect_skew.pdf")
ggsave("Tetraselmis_experiment/figures/variability_effect_skew.png", width = 5, height = 4)

### now let's look at the mins and maxes

results5 %>% 
	select(isolate.code, tmin, tmax, temperature_min, temperature_max, SD) %>% 
	mutate(max_diff = tmax - temperature_max) %>% 
	mutate(min_diff = tmin - temperature_min) %>% 
	# gather(key = type, value = temperature, 2:5) %>% 
	# ggplot(aes(x = type, y = temperature)) + geom_boxplot()
	select(isolate.code, SD, max_diff, min_diff) %>% 
	gather(key = type, value = difference, 3:4) %>% 
	ggplot(aes(x = type, y = difference, color = SD)) + geom_point() +
	scale_color_viridis(option = "inferno")

ct_max <- results5 %>% 
	select(isolate.code, tmin, tmax, temperature_min, temperature_max, SD, w, rel.curveskew) %>% 
	mutate(breadth_constant = tmax-tmin) %>%
	mutate(max_diff = tmax - temperature_max) %>% 
	mutate(breadth_variable = temperature_max - temperature_min) %>% 
	mutate(breadth_diff = breadth_constant - breadth_variable) %>% 
	ggplot(aes(x = rel.curveskew, y = max_diff, color = SD)) + geom_point(size = 2) +
	scale_color_viridis(option = "inferno") + theme_bw() +
	# geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
 ylab("CTmax(cons) - CTmax(var) (°C)") +
	xlab("Curve skewness") + 
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) + ylim(0, 4)
ggsave("Tetraselmis_experiment/figures/CT_max_variability_effect_skew.png", width = 5, height = 4)




# giant plot --------------------------------------------------------------

library(cowplot)

## c plot comes from line 294 in 09_bootstrap_plots.R
plot3 <- plot_grid(c,topt_skew, map_plot,ct_max, labels = c("A", "B", "C", "D"), align = "v", ncol = 2, nrow = 2)
save_plot("Tetraselmis_experiment/figures/figure3.png", plot3, base_width = 12, base_height = 8.1)

results5 %>% 
	ggplot(aes(x = Mean, y = topt)) + geom_point() +
	geom_abline(slope = 1, intercept = 0) + theme_bw()  + ylim(0, 40)

results5 %>% 
	mutate(diff = topt_variable - topt) %>% 
	ggplot(aes(x = topt, y = diff, color = rel.curveskew, label = isolate.code)) + geom_point(size = 3) + 
	geom_hline(yintercept = 0) + theme_bw() +  scale_color_viridis(option = "inferno") +ylim(-4, 4)

results6 %>% 
	mutate(diff = topt_variable - topt) %>% 
	ggplot(aes(x = SD, y = diff, color = rel.curveskew, label = isolate.code)) + geom_point(size = 3) + 
	geom_hline(yintercept = 0) + theme_bw() +  scale_color_viridis(option = "inferno") +ylim(-4, 4)

results5 %>% 
	mutate(diff = topt_variable - topt) %>% 
	ggplot(aes(x = SD, y = diff, color = rel.curveskew, label = isolate.code)) + geom_point(size = 3) + 
	geom_hline(yintercept = 0) + theme_bw() +  scale_color_viridis() +ylim(-4, 4)

results6 %>% 
	mutate(diff = topt_variable - topt) %>% 
	ggplot(aes(x = Mean, y = SD, color = rel.curveskew, label = isolate.code)) + geom_point(size = 3) + 
	theme_bw() +  scale_color_viridis(option = "inferno") + geom_smooth()



results5 %>% 
	ggplot(aes(x = Mean, y = SD, label = isolate.code)) + geom_point()


results5 %>% 
	gather(key = type, value = temperature, topt, topt_variable) %>% 
	ggplot(aes(x = temperature, y = Mean, color = SD)) + geom_point() +
	geom_abline(slope = 1, intercept = 0) + theme_bw() + scale_color_viridis() +
	facet_wrap( ~ type)


# new question: what is the predicted difference in r, given the l --------

r_predicted <- function(data) {
	all <- data
	x <- seq(0, 45, by = 0.1)
	tpc<-function(x){
		res<-all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)
		res
	}
	## step 2 
	variable_predictions <- function(x) {
		y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(all$SD^2)
	}
	
	predicted_growth_variable <- sapply(x, variable_predictions)
	predicted_growth_variable2 <- data.frame(x, predicted_growth_variable) %>% 
		rename(temperature = x, 
					 growth.rate = predicted_growth_variable) 
	# data.frame(all$isolate.code, predicted_growth_variable2)
}

all_split <- all_thermal_data %>% 
	filter(!is.na(topt)) %>% 
	split(.$isolate.code)

predicted_r <- all_split %>% 
	map_df(r_predicted, .id = "isolate_code")

predicted_r_variable <- predicted_r %>% 
	rename(growth_variable = growth.rate)

write_csv(predicted_r_variable, "Tetraselmis_experiment/data-processed/predicted_r_variable.csv")
## now get constant prediction

r_predicted_constant <- function(data) {
	all <- data
	x <- seq(0, 45, by = 0.1)
	tpc<-function(x){
		res<-all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)
		res
	}
	## step 2 
	constant_predictions <- function(x) {
		y <- tpc(x)
	}
	
	predicted_growth_constant <- sapply(x, constant_predictions)
	predicted_growth_constant2 <- data.frame(x, predicted_growth_constant) %>% 
		rename(temperature = x, 
					 growth.rate = predicted_growth_constant) 
	# data.frame(all$isolate.code, predicted_growth_variable2)
}

all_split <- all_thermal_data %>% 
	filter(!is.na(topt)) %>% 
	split(.$isolate.code)

predicted_r_constant <- all_split %>% 
	map_df(r_predicted_constant, .id = "isolate_code") 

predicted_r_constant <- predicted_r_constant %>% 
	rename(growth_constant = growth.rate)

write_csv(predicted_r_constant, "Tetraselmis_experiment/data-processed/predicted_r_constant.csv")

all_predicted_growth <- left_join(predicted_r_constant, predicted_r_variable, by = c("isolate_code", "temperature"))

str(all_predicted_growth2)
all_predicted_growth2 <- all_predicted_growth %>% 
	rename(Mean = temperature) %>% 
	rename(isolate.code = isolate_code) %>% 
	mutate(isolate.code = as.integer(isolate.code))

thermal_data_round <- all_thermal_data %>% 
	mutate(Mean = round(Mean, digits = 1)) 

predicted_growth_temperature <- left_join(thermal_data_round, all_predicted_growth2, by = c("Mean", "isolate.code"))

predicted_growth_temperature %>% 
	mutate(growth_diff = growth_constant - growth_variable) %>% 
	ggplot(aes(x = rel.curveskew, y = growth_diff, color = SD)) + geom_point(size = 2) +
	scale_color_viridis(option = "inferno") + theme_bw() +
	# geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
	ylab("predicted growth under constant - predicted growth under variable") +
	xlab("curve skew") + 
	theme_bw() +
	geom_hline(yintercept = 0) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica"))

# ok let’s try this again -------------------------------------------------

sub77 <- all_thermal_data %>% 
	filter(isolate.code == 77)

r_predicted_cons_var <- function(data) {
	all <- data
	x <- seq(0, 45, by = 0.1)
	tpc<-function(x){
		res<-all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)
		res
	}
	## step 2 
	variable_predictions <- function(x) {
		y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(all$SD^2)
	}
	
	t <- all$Mean[[1]]
	predicted_growth_variable <- sapply(t, variable_predictions)
	predicted_growth_variable2 <- data.frame(t, predicted_growth_variable) %>% 
		rename(temperature = t, 
					 predicted_growth_variable = predicted_growth_variable) 
	
	constant_predictions <- function(x) {
		y <- tpc(x)
	}
	
	predicted_growth_constant <- sapply(t, constant_predictions)
	
	predicted_growth_constant2 <- data.frame(predicted_growth_variable2, predicted_growth_constant) 
	
	# data.frame(all$isolate.code, predicted_growth_variable2)
}

all_split <- all_thermal_data %>% 
	filter(!is.na(topt)) %>% 
	filter(!is.na(SD)) %>% 
	filter(!is.na(Mean)) %>% 
	split(.$isolate.code)

predicted_r_var_cons <- all_split %>% 
	map_df(r_predicted_cons_var, .id = "isolate_code")

pred_r_vc <- predicted_r_var_cons %>% 
	rename(isolate.code = isolate_code) %>% 
	mutate(isolate.code = as.integer(isolate.code))
predicted_growth_temperature2 <- left_join(all_thermal_data, pred_r_vc, by = "isolate.code")

write_csv(predicted_growth_temperature2, "Tetraselmis_experiment/data-processed/predicted_growth_temperature2.csv")

predicted_growth_temperature2 <- read_csv("Tetraselmis_experiment/data-processed/predicted_growth_temperature2.csv")


predicted_growth_temperature2 %>% 
	mutate(growth_diff = predicted_growth_constant - predicted_growth_variable) %>% 
	filter(growth_diff < 0.5 ) %>% 
	mutate(temp_diff = topt- Mean) %>% 
	ggplot(aes(x = temp_diff, y = growth_diff, color = rel.curveskew)) + geom_point(size = 2) +
	scale_color_viridis(option = "inferno") + theme_bw() +
	# geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
	ylab("predicted growth under constant - predicted growth under variable") +
	xlab("Temperature difference") + 
	theme_bw() +
	geom_hline(yintercept = 0) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica"))
