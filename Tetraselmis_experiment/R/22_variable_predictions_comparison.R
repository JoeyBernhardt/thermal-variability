

## goal: estimate the growth rate over time...with the daily temperature data
### so get time averaged growth rate

library(tidyverse)
library(janitor)
library(lubridate)
library(cowplot)


curve_data <- read_csv("Tetraselmis_experiment/data-processed/curve_data_20140606.csv") %>% 
	clean_names()

growth_raw <- read_csv("Tetraselmis_experiment/data-processed/growth_data_20140606.csv")
all_thermal_data <- read_csv("Tetraselmis_experiment/data-processed/all_thermal_data.csv") %>% 
	clean_names

temperature <- read_csv("Tetraselmis_experiment/data-processed/temps_all_33.875N.csv")


all_thermal_data %>% 
	filter(source == "Dokai Bay, Japan" , curvequal == "good") %>% View

curve <- all_thermal_data %>% 
	# filter(grepl("Detonula", speciesname)) %>% 
	filter(isolate_code == 462) 


tpc1<-function(x){
	res<-curve$a[1]*exp(curve$b[1]*x)*(1-((x-curve$z[1])/(curve$w[1]/2))^2)
	res
}


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + stat_function(fun = tpc1, color = "black", size = 2) +xlim(35, 35.5) + 
	# ylim(0, 1.5) + 
	theme_bw() + ylab("Growth rate") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line.y = element_line(color="black"),
				axis.line.x = element_line(color="black")) +
	theme(text = element_text(size=14, family = "Helvetica")) +
	xlab(expression("Temperature (" *degree * "C)")) + geom_hline(yintercept = 0)



# now take time averaged growth rate --------------------------------------

tpc <- function(df){
	x <- df$sst
	res<-curve$a[1]*exp(curve$b[1]*x)*(1-((x-curve$z[1])/(curve$w[1]/2))^2)
	data.frame(temperature = x, growth_rate = res)
}



temp_split <- temperature %>% 
	split(.$Date)
	
growth_time  <- temp_split %>% 
	map_df(tpc, .id = "date") 


growth_time2 <- growth_time %>% 
	mutate(date = ymd(date)) 

growth_time2 %>% 
	ggplot(aes(x = date, y = growth_rate)) + geom_line() +
	theme_classic()

str(temperature)
temperature2 <- temperature %>% 
	rename(date = Date)

temps_growth <- left_join(temperature2, growth_time2, by = "date") %>% 
	gather(key = type, value = value, temperature:growth_rate)

temps_growth_wide <- left_join(temperature2, growth_time2, by = "date")

mean(temperature$sst)

temp_plot <- temps_growth_wide %>% 
	ggplot(aes(x = date, y = temperature)) + geom_line() +
	theme_classic() + geom_hline(yintercept = 24.3521) +
	geom_hline(yintercept = 9.7, color = "blue") + geom_hline(yintercept = 35.7, color = "red") +
	geom_hline(yintercept = 18.727, color = "grey", linetype = "dashed") 


growth_plot <- temps_growth_wide %>% 
	ggplot(aes(x = date, y = growth_rate)) + geom_line() +
	theme_classic() 
	
plot_grid(temp_plot, growth_plot, nrow=2)


temps_growth %>% 
	ggplot(aes(x = date, y = value, color = type)) + geom_line() +
	theme_classic() + facet_wrap( ~ type, scales = "free")

time_integrated_prediction <- mean(growth_time$growth_rate)


## now the approximation
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


tpc1 <- function(x){
	x <- curve$mean
	res<-curve$a[1]*exp(curve$b[1]*x)*(1-((x-curve$z[1])/(curve$w[1]/2))^2)
	res
}


x <- seq(0, 45, by = 0.1)
x <- 15
## step 2 
variable_predictions <- function(x) {
	y <- tpc1(x) + derivative(f = tpc1, x = x, order = 2)*0.5*(curve$sd^2)
}

predicted_growth_variable <- sapply(x, variable_predictions)


variable_predictions <- function(data) {
	data <- curve
	x <- data$mean
	SD <- data$sd
	y <- tpc1(x) + derivative(f = tpc1, x = x, order = 2)*0.5*(SD^2)
	y
}



variable_predictions(curve)

tpc1(curve$mean)

all <- curve

x <- curve$mean
tpc<-function(x){
	res<-all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)
	res
}

## step 2 
variable_predictions <- function(x) {
	x <- data$mean
	y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(all$sd^2)
	y
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


approx_prediction <- variable_predictions(curve)
no_var_prediction <- tpc(curve$mean)
time_integrated_prediction <- mean(growth_time$growth_rate)

compare <- data.frame(approx = approx_prediction, no_var = no_var_prediction, time_int = time_integrated_prediction)

compare %>% 
	rename(time_integration = time_int,
				 approximation = approx,
				 no_variability = no_var) %>% 
	gather(key = "prediction type", value = "estimated growth rate") %>% 
	ggplot(aes(x = reorder(`prediction type`, `estimated growth rate`), y = `estimated growth rate`)) + geom_point() +
	xlab("prediction type")
