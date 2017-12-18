library(tidyverse)

### now try and make fake right and left skewed TPCs!

fits_constant_variable <- read_csv("Tetraselmis_experiment/data-processed/fits_constant_variable.csv")

curve_variable<-function(x){
	res<-fits_constant_variable$a.list[2]*exp(fits_constant_variable$b.list[2]*x)*(1-((x-fits_constant_variable$z.list[2])/(fits_constant_variable$w.list[2]/2))^2)
	res
}



briere <- function(x, Tmin, Tmax, m, c){
	res <- c*(x*(x-Tmin)*(((Tmax-x))^(1/m)))
	res
}

mod_briere <- function(x, Tmin, Tmax, m, c){
	res <- c*(-x + Tmin + Tmax)*(Tmax - x)*((x - Tmin)^(1/m))
	res
}

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 
p + 
	# stat_function(fun = briere, args = list(Tmin = 0, Tmax = 30, m = 2, c = 0.0005)) +
	stat_function(fun = mod_briere, args = list(Tmin = 0, Tmax = 32, m = 3, c = 0.0004), color = "yellow", size = 2) +
	stat_function(fun = mod_briere, args = list(Tmin = 0, Tmax = 32, m = 2, c = 0.0004), color = "black") +
	stat_function(fun = mod_briere, args = list(Tmin = 0, Tmax = 32, m = 1.5, c = 0.0004), color = "red") +
	stat_function(fun = mod_briere, args = list(Tmin = 0, Tmax = 32, m = 1.3, c = 0.0004), color = "blue") +
	stat_function(fun = mod_briere, args = list(Tmin = 0, Tmax = 32, m = 1.35, c = 0.0004), color = "grey") +
	stat_function(fun = mod_briere, args = list(Tmin = 0, Tmax = 32, m = 1.1, c = 0.0004), color = "green") +
	stat_function(fun = mod_briere, args = list(Tmin = 0, Tmax = 32, m = 1.2, c = 0.0004), color = "pink") +
	stat_function(fun = mod_briere, args = list(Tmin = 0, Tmax = 32, m = 1, c = 0.0004), color = "purple") +
	 xlim(0, 31) + ylim(0, 3) + theme_bw() + xlab("Temperature (°C)") + ylab("Performance")

