library(tidyverse)

fits_variable <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v.csv")

### goal: find the tmin and tmaxes of the bootstrapped curves

breadth_function <- function(data) {
  all <- data
  x <- seq(-5, 55, by = 0.1)
  tpc<-function(x){
    res<-all$a.list[1]*exp(all$b.list[1]*x)*(1-((x-all$z.list[1])/(all$w.list[1]/2))^2)
    res
  }
  ## step 2 
  variable_predictions <- function(x) {
    y <- tpc(x)
  }
  
  predicted_growth_variable <- sapply(x, variable_predictions)
  curve.id.list <- rep(x = all$curve.id.list, times = length(x))
  predicted_growth_variable2 <- data.frame(x, predicted_growth_variable, curve.id.list) %>% 
    rename(temperature = x, 
           growth.rate = predicted_growth_variable) %>% 
    filter(growth.rate > 0) %>% 
    group_by(curve.id.list) %>% 
    top_n(n = -10, wt = growth.rate) %>% 
    summarise_each(funs(min, max), temperature)
  # data.frame(all$isolate.code, predicted_growth_variable2)
}


all_split <- fits_variable %>% 
  split(.$curve.id.list)



results_breadth <- all_split %>% 
  map_df(breadth_function, .id = "curve_id")

### next step, get rid of all curves that estimate a tmin lower than -5

results_above_minus5 <- results_breadth %>% 
  filter(temperature_min > -5)

curves_above_minus5 <- results_above_minus5$curve.id.list


### now re-do above process, but with smaller increments

