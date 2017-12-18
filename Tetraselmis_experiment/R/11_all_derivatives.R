library(tidyverse)


all_thermal_data <- read_csv("Tetraselmis_experiment/data-processed/all_thermal_data.csv")


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



get_derivatives <- function(data) {
	all <- data
	x <- seq(all$tmin,all$tmax, by = 0.1)
	tpc<-function(x){
		res<-all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)
		res
	}
	
	der <- sapply(x, derivative, f = tpc, order = 2)
	derivatives <- data.frame(x, der) 
}

all_split <- all_thermal_data %>% 
	filter(!is.na(topt)) %>% 
	filter(!is.na(SD)) %>% 
	filter(!is.na(Mean)) %>% 
	split(.$isolate.code)


isolate1 <- all_thermal_data %>% 
	filter(isolate.code == 1)

all_derivatives <- all_split %>% 
	map_df(get_derivatives, .id = "isolate_code")

write_csv(all_derivatives, "Tetraselmis_experiment/data-processed/all_derivatives.csv")


derivative_sign <- all_derivatives %>% 
	group_by(isolate_code) %>% 
	summarise_each(funs(max, min), der) %>% 
	mutate(derivative_sign = ifelse(der_max > 0, "positive", "negative"))

write_csv(derivative_sign, "Tetraselmis_experiment/data-processed/derivative_sign.csv")
