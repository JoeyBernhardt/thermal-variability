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
library(broom)
library(cowplot)

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


# import thermal sd -------------------------------------------------------


thermal_sd <- read_csv("Tetraselmis_experiment/data-processed/thermal_sd_aug.csv")

all_thermal_data <- left_join(thermal_sd, thomas3, by = c("isolate.code"))


write_csv(all_thermal_data, "Tetraselmis_experiment/data-processed/all_thermal_data.csv")

## let's plot one curve
## 592, 324, 614, 613, 317, 209, 538, 607...77 could work for the figures

## let's look at 550, 132, 140, 337 for good examples of the effect of skew
curve <- all_thermal_data %>% 
	filter(isolate.code == 94)

tpc1<-function(x){
	res<-curve$a[1]*exp(curve$b[1]*x)*(1-((x-curve$z[1])/(curve$w[1]/2))^2)
	res
}

### negative skew example
curve_neg <- all_thermal_data %>% 
	filter(isolate.code == 550)

tpc_neg<-function(x){
	res<-curve_neg$a[1]*exp(curve_neg$b[1]*x)*(1-((x-curve_neg$z[1])/(curve_neg$w[1]/2))^2)
	res
}


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
neg_skew_ex <- p + stat_function(fun = tpc_neg, color = "black", size = 1) +xlim(0, 33) + 
	ylim(0, curve_neg$mu.g.opt.val.list+0.05) + 
	theme_bw() + ylab("") +
	geom_vline(xintercept = curve_neg$Mean, color = "grey", linetype = "dashed", size = 1) +
	geom_vline(xintercept = curve_neg$topt) +
	geom_rect(mapping=aes(xmin=curve_neg$Mean-curve_neg$SD, xmax=curve_neg$Mean+curve_neg$SD, ymin=0, ymax=curve_neg$mu.g.opt.val.list+0.05), alpha=0.3)+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line.y = element_line(color="black"),
				axis.line.x = element_line(color="black")) +
	theme(text = element_text(size=14, family = "Helvetica")) +
	xlab(expression("Temperature (" *degree * "C)"))
ggsave("Tetraselmis_experiment/figures/pos_skew_ex.png", width = 3, height = 2)


### positive skew example, possibly 140, 214, 337, 35, 292, 77, 221, 211
curve_pos <- all_thermal_data %>% 
	filter(isolate.code == 211)

tpc_pos<-function(x){
	res<-curve_pos$a[1]*exp(curve_pos$b[1]*x)*(1-((x-curve_pos$z[1])/(curve_pos$w[1]/2))^2)
	res
}


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
pos_skew_ex <- p + stat_function(fun = tpc_pos, color = "black", size = 1) +xlim(0, 20) + 
	ylim(0, curve_pos$mu.g.opt.val.list+0.05) + 
	theme_bw() + ylab("") +
	geom_vline(xintercept = curve_pos$Mean, color = "grey", linetype = "dashed", size = 1) +
	geom_vline(xintercept = curve_pos$topt) +
	geom_rect(mapping=aes(xmin=curve_pos$Mean-curve_pos$SD, xmax=curve_pos$Mean+curve_pos$SD, ymin=0, ymax=curve_pos$mu.g.opt.val.list+0.05), alpha=0.3)+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line.y = element_line(color="black"),
				axis.line.x = element_line(color="black")) +
	theme(text = element_text(size=14, family = "Helvetica")) +
	xlab(expression("Temperature (" *degree * "C)")) 


plot3_examples <- plot_grid(neg_skew_ex, pos_skew_ex, ncol = 1, nrow = 2)
save_plot("Tetraselmis_experiment/figures/figure3_examples.png", plot3_examples, base_width = 2.5, base_height = 4)


curves_selected <- all_thermal_data %>% 
	filter(isolate.code %in% c("211", "550"))


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
right_skew <- p + stat_function(fun = tpc1, color = "black", size = 2) +xlim(0, 28) + 
	ylim(0, curve$mu.g.opt.val.list+0.05) + 
	theme_bw() + ylab("Growth rate") +
	geom_vline(xintercept = curve$Mean, color = "grey", linetype = "dashed", size = 1.5) +
	geom_vline(xintercept = curve$topt) +
	geom_rect(mapping=aes(xmin=curve$Mean-curve$SD, xmax=curve$Mean+curve$SD, ymin=0, ymax=curve$mu.g.opt.val.list+0.05), alpha=0.3)+
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



# Comparison of var and cons TPCs -----------------------------------------

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
	mutate(breadth_constant = tmax-tmin) %>%
	mutate(max_diff = tmax - temperature_max) %>% 
	mutate(breadth_variable = temperature_max - temperature_min) %>% 
	mutate(breadth_diff = breadth_constant - breadth_variable) %>% 
	filter(habitat %in% c("marine", "estuarine", "saltmarsh"))

write_csv(results5, "Tetraselmis_experiment/data-processed/results5.csv")

# read in results5 --------------------------------------------------------


results5 <- read_csv("Tetraselmis_experiment/data-processed/results5.csv")
results6 <- results5 %>% 
	# rename(topt_variable = temperature) %>% 
	mutate(rev_diff = topt_variable - topt)

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


# model selection ---------------------------------------------------------


## skew and SD best explain variation in the effect of variability
m1 <- lm(diff ~ rel.curveskew + SD, data = results5)
m1 <- lm(max_diff ~ rel.curveskew + SD, data = results5)
m1 <- lm(breadth_diff ~ rel.curveskew + SD, data = results5)
m2 <- lm(diff ~ rel.curveskew, data = results5)
m3 <- lm(diff ~ SD, data = results5)
AIC(m1)
AIC(m2)
AIC(m3)
summary(m1)
tidy(m1, conf.int = TRUE)

mean(abs(results5$diff), na.rm = TRUE)
hist(results5$diff)

results5 %>% 
	ggplot(aes(x = SD, y = rel.curveskew)) + geom_point() + geom_smooth()

cor(use = "na.or.complete", x = abs(results5$latitude), y = results5$rel.curveskew)


# plots -------------------------------------------------------------------


topt_skew <- results6 %>% 
	filter(minqual == "good", maxqual == "good") %>% 
 filter(mu.n > 4) %>% 
	# filter(isolate.code %in% greater_90) %>% 
	# filter(rel.curveskew < 0) %>% 
	ggplot(aes(x = rel.curveskew, y = rev_diff, color = SD, label = isolate.code)) + geom_point(size = 3) +
	# geom_label()
	scale_color_viridis(option = "inferno") + theme_bw() + geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
	geom_hline(yintercept = 0) + ylab("Topt(var) - Topt(cons) (°C)") +
	xlab("") + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=14, family = "Helvetica")) +
	geom_vline(xintercept = 0) + geom_point(size = 3, shape = 1, color = "black") + xlim(-0.02, 0.02)
ggsave("Tetraselmis_experiment/figures/variability_effect_skew.pdf")
ggsave("Tetraselmis_experiment/figures/variability_effect_skew.png", width = 5, height = 3)

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
	mutate(rev_max_diff = temperature_max - tmax) %>% 
	mutate(breadth_variable = temperature_max - temperature_min) %>% 
	mutate(breadth_diff = breadth_constant - breadth_variable) %>% 
	ggplot(aes(x = rel.curveskew, y = rev_max_diff, color = SD)) + geom_point(size = 3) +
	scale_color_viridis(option = "inferno") + theme_bw() +
	geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
	ylab("Tmax(var) - Tmax(cons) (°C)") +
	xlab("Curve skewness") + 
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=14, family = "Helvetica")) + geom_point(size = 3, shape = 1, color = "black") + geom_hline(yintercept = 0) +
	xlim(-0.02, 0.02) + geom_vline(xintercept = 0) 
ggsave("Tetraselmis_experiment/figures/CT_max_variability_effect_skew.png", width = 5, height = 3)


skewness_plots <- plot_grid(topt_skew, ct_max, align = "v", ncol = 1, nrow = 2)
save_plot("Tetraselmis_experiment/figures/skewness_plots.png", skewness_plots, base_width = 4.5, base_height = 5)

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

predicted_r_constant <- read_csv("Tetraselmis_experiment/data-processed/predicted_r_constant.csv")
predicted_r_variable <- read_csv("Tetraselmis_experiment/data-processed/predicted_r_variable.csv")

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
	der <- derivative(f = tpc, x = t, order = 2)
	
	predicted_growth_constant2 <- data.frame(predicted_growth_variable2, predicted_growth_constant, der) 
	
	# data.frame(all$isolate.code, predicted_growth_variable2)
}

all_split <- all_thermal_data %>% 
	filter(!is.na(topt)) %>% 
	filter(!is.na(SD)) %>% 
	filter(!is.na(Mean)) %>% 
	split(.$isolate.code)

predicted_r_var_cons_der <- all_split %>% 
	map_df(r_predicted_cons_var, .id = "isolate_code")

pred_r_vc <- predicted_r_var_cons_der %>% 
	rename(isolate.code = isolate_code) %>% 
	mutate(isolate.code = as.integer(isolate.code))
predicted_growth_temperature2 <- left_join(all_thermal_data, pred_r_vc, by = "isolate.code")

write_csv(predicted_growth_temperature2, "Tetraselmis_experiment/data-processed/predicted_growth_temperature2.csv")

predicted_growth_temperature2 <- read_csv("Tetraselmis_experiment/data-processed/predicted_growth_temperature2.csv")

names(predicted_growth_temperature2)

r_diff <- 
	predicted_growth_temperature2 %>% 
	filter(curvequal == "good", minqual == "good", maxqual == "good") %>% 
	mutate(growth_diff_rel = ((predicted_growth_constant - predicted_growth_variable)/mu.g.opt.val.list)*100) %>% 
	mutate(growth_diff = predicted_growth_constant - predicted_growth_variable) %>% 
	mutate(rev_growth_diff = predicted_growth_variable - predicted_growth_constant) %>% 
	mutate(temp_diff = topt- Mean) %>% 
	mutate(skew_dir = ifelse(rel.curveskew<0, "negative skew", "positive skew")) %>% 
	filter(growth_diff_rel < 60) %>%
	ggplot(aes(x = temp_diff, y = rev_growth_diff, color = SD, label = isolate.code)) + geom_point(size = 3) +
	scale_color_viridis(option = "inferno") + theme_bw() +
	facet_wrap( ~ skew_dir) +
	# geom_smooth(color = "grey", alpha = 0.3, size = 0.4) +
	ylab("r(var) - r(cons)") +
	xlab("TOpt - Mean temperature at isolation location (°C)") + 
	theme_bw() +
	geom_hline(yintercept = 0) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +geom_vline(xintercept = 0) + 
	ylim(-0.4, 0.02) +
	# theme(legend.position = "none") +
	theme(strip.background = element_rect(colour="white", fill="white")) +
	geom_point(size = 3, shape = 1, color = "black")
ggsave("Tetraselmis_experiment/figures/r_diff_plot.png", width = 6, height = 3)
	
## let's look at 550, 132, 140, 337

plot3b <- plot_grid(c,topt_skew, r_diff,ct_max, map_plot, pos_skew_ex, labels = c("A", "B", "C", "D", "E"), align = "v", ncol = 2, nrow = 3)
save_plot("Tetraselmis_experiment/figures/figure3b.png", plot3b, base_width = 12, base_height = 12)


# What is the value of the deriv at MAT? ----------------------------------

variability_results <- predicted_growth_temperature2 %>% 
	# filter(curvequal == "good") %>% 
	mutate(growth_diff_rel = ((predicted_growth_constant - predicted_growth_variable)/mu.g.opt.val.list)*100) %>% 
	mutate(growth_diff = predicted_growth_constant - predicted_growth_variable) %>% 
	mutate(temp_diff = topt- Mean) %>% 
	mutate(skew_dir = ifelse(rel.curveskew<0, "negative skew", "positive skew")) 

write_csv(variability_results, "Tetraselmis_experiment/data-processed/variability_results.csv")

variability_results <- read_csv("Tetraselmis_experiment/data-processed/variability_results.csv")
derivative_sign <- read_csv("Tetraselmis_experiment/data-processed/derivative_sign.csv")
results5 <- read_csv("Tetraselmis_experiment/data-processed/results5.csv")

variability_results %>% 
	filter(growth_diff_rel < 60) %>%
	# filter(curvequal == "good", minqual == "good", maxqual == "good") %>% 
	# filter(der > -0.10) %>% 
	# filter(skew_dir == "negative skew" ) %>% 
	# filter(growth_diff < 0 & temp_diff < 0) %>% 
	ggplot(aes(x = growth_diff_rel, y = der, color = temp_diff, shape = skew_dir, size = SD, label = isolate.code)) + geom_point() +
	scale_color_viridis(option = "inferno") +
	theme_bw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	ylab("P''") + xlab("r(cons)- r(var)")
ggsave("Tetraselmis_experiment/figures/derivative_r_diff.png")

## the answer seems to be that it's both skewness and derivative value that are driving the pattern, 
## and there is some overlap between these variables. I.e. there is are no curves with a positive 
## 2nd derivative that are negatively skewed, and within a value of second deriv. there are negative and positive
## skews that both occur

## For the topt results -- well none of the curves will be accelerating at Topt, BUT, I can look to see 
## at how many of the curves have positive vs. negative second derivatives at any point along the curve


hist(predicted_growth_temperature2$der)

der_phyto <- variability_results %>% 
	select(isolate.code, der)

all_phyto <- left_join(der_phyto, results5)
all_phyto2 <- left_join(derivative_sign, results5)

all_phyto2 %>%
	mutate(derivative_sign = ifelse(derivative_sign == "negative", "decelerating", "accelerating")) %>% 
	# mutate(diff = topt_variable - topt) %>% 
	ggplot(aes(x = rel.curveskew, y = diff, color = derivative_sign, label = isolate.code)) + geom_point(size = 2) +
	# geom_label()
	# scale_color_viridis(option = "inferno") +
	theme_bw() + geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
	geom_hline(yintercept = 0) + ylab("Topt(cons) - Topt(var) (°C)") +
	xlab("Curve skewness") + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	geom_vline(xintercept = 0)
ggsave("Tetraselmis_experiment/figures/topt_diff_derivative_sign.png")


all_phyto2 %>%
	mutate(skew_dir = ifelse(rel.curveskew<0, "negative skew", "positive skew")) %>% 
	filter(!is.na(skew_dir)) %>% 
	mutate(derivative_sign = ifelse(derivative_sign == "negative", "decelerating", "accelerating")) %>% 
	ggplot(aes(x = rel.curveskew, y = diff, shape = derivative_sign, color = SD)) + geom_jitter(size = 3) + theme_bw() +
	geom_hline(yintercept = 0) +
	scale_color_viridis(option = "inferno") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	ylab("Topt(cons) - Topt(var)") + xlab("Curve skewness") + geom_vline(xintercept = 0)
ggsave("Tetraselmis_experiment/figures/topt_diff_derivative_sign_sd.png")


all_phyto2 %>%
	mutate(skew_dir = ifelse(rel.curveskew<0, "negative skew", "positive skew")) %>% 
	filter(mu.rsqrlist > 0.85) %>% 
	filter(!is.na(skew_dir)) %>% 
	filter(curvequal == "good", maxqual == "good", minqual == "good") %>% 
	mutate(derivative_sign = ifelse(derivative_sign == "negative", "decelerating", "accelerating")) %>% 
	ggplot(aes(x = rel.curveskew, y = diff, color = SD)) + geom_point(size = 3) + theme_bw() +
	geom_smooth(method = "lm", color = "black") +
	geom_hline(yintercept = 0) +
	scale_color_viridis(option = "inferno") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	ylab("Topt(cons) - Topt(var)") + xlab("Curve skewness") + geom_vline(xintercept = 0)
ggsave("Tetraselmis_experiment/figures/topt_diff_rsq85.png", width = 6, height = 4)


all_phyto2 %>%
	mutate(derivative_sign = ifelse(derivative_sign == "negative", "fully decelerating", "sometimes accelerating")) %>% 
	mutate(skew_dir = ifelse(rel.curveskew<0, "negative skew", "positive skew")) %>% 
	filter(!is.na(skew_dir)) %>% 
	ggplot(aes(x = rel.curveskew, y = derivative_sign)) + geom_point(size = 5, alpha = 0.3) + 
	geom_vline(xintercept = 0) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	ylab("Curve shape") + xlab("Curve skewness")
ggsave("Tetraselmis_experiment/figures/skew_shape_dotchart.png")


all_phyto %>% 
	ggplot(aes(x = der, y = diff, color = SD)) + geom_point() +
	geom_point(size = 2) +
	# geom_label()
	scale_color_viridis(option = "inferno") + theme_bw() + geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
	geom_hline(yintercept = 0) + ylab("Topt(cons) - Topt(var) (°C)") +
	xlab("Derivative") + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	geom_vline(xintercept = 0)


variability_results %>% 
	ggplot(aes(x = latitude, y = mu.g.opt.val.list)) + geom_point() +
	theme_bw() +ylab("r(max)") + xlab("Latitude")

library(visreg)

m1 <- lm(growth_diff_rel ~ rel.curveskew + SD + der, data = variability_results)
m2 <- lm(growth_diff_rel ~ der, data = variability_results)
m3 <- lm(growth_diff_rel ~ der + rel.curveskew, data = variability_results)
summary(m1)
summary(m2)
AIC(m1)
AIC(m2)
AIC(m3)

visreg(m1)
