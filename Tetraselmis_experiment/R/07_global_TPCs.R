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

# g <- geo_left_join(thomas3, therm2, max_dist = 200, method = "geo")

# thermal_sd <- g %>% 
# 	select(starts_with("lat"), starts_with("long"), Mean, isolate.code, everything()) %>% 
# 	mutate(lat_diff = abs(latitude.x - latitude.y),
# 				 long_diff = abs(longitude.x - longitude.y)) %>% 
# 	select(lat_diff, long_diff, everything()) %>%
# 	mutate(total_diff = lat_diff + long_diff) %>% 
# 	select(total_diff, everything()) %>% 
# 	group_by(isolate.code) %>% 
# 	top_n(wt = total_diff, n = -1) %>% 
# 	select(SD, everything()) %>% 
# 	group_by(isolate.code) %>% 
# 	summarise_each(funs(mean), SD, Mean)

write_csv(thermal_sd, "Tetraselmis_experiment/data-processed/thermal_sd.csv")

thermal_sd <- read_csv("Tetraselmis_experiment/data-processed/thermal_sd.csv")

thermal_sd %>% 
	ggplot(aes(x = Mean, y = SD)) + geom_point()

all_thermal_data <- left_join(thermal_sd, thomas3, by = c("isolate.code"))




## let's plot one curve
## 592, 324, 614, 613, 317, 209, 538, 607
curve <- thomas2 %>% 
	filter(isolate.code == 40)
tpc<-function(x){
	res<-curve$a[1]*exp(curve$b[1]*x)*(1-((x-curve$z[1])/(curve$w[1]/2))^2)
	res
}

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + stat_function(fun = tpc, color = "black", size = 2) +xlim(10, 40) + ylim(0, 0.5) + 
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

x <- seq(0, 45, by = 0.1)
variable_predictions <- function(x) {
	y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(2.176435^2)
}


predicted_growth_variable <- sapply(x, variable_predictions)
predicted_growth_variable2 <- data.frame(x, predicted_growth_variable) %>% 
	rename(temperature = x, 
				 growth.rate = predicted_growth_variable) %>% 
	top_n(n = 1, wt = growth.rate)



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
data.frame(all$isolate.code, predicted_growth_variable2)
}


all_split <- all_thermal_data %>% 
	# filter(isolate.code < 15) %>% 
	split(.$isolate.code)

results <- all_split %>% 
	map_df(predict_function, .id = "isolate_code")

write_csv(results, "Tetraselmis_experiment/data-processsed/Topt_variable_global.csv")
results <- read_csv("Tetraselmis_experiment/data-processsed/Topt_variable_global.csv")


results2 <- results %>% 
	mutate(isolate.code = as.integer(isolate_code))

results3 <- left_join(thomas3, results2, by = c("isolate.code"))

results4 <- left_join(results3, thermal_sd, by = c("isolate.code"))

results5 <- results4 %>% 
	rename(topt_variable = temperature) %>% 
	filter(curvequal == "good") %>% 
	mutate(diff = topt_variable - topt) %>% 
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

results5 %>% 
	mutate(diff = topt_variable - topt) %>% 
	ggplot(aes(x = rel.curveskew, y = diff, color = SD, label = isolate.code)) + geom_point(size = 2) +
	# geom_label()
	scale_color_viridis(option = "inferno") + theme_bw() + geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
	geom_hline(yintercept = 0) + ylab("Effect of thermal variability \nTopt constant - Topt variable (C)") +
	xlab("Curve skewness") + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=14, family = "Helvetica")) 
ggsave("Tetraselmis_experiment/figures/variability_effect_skew.pdf")

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
