### limiting resampling exploring


library(tidyverse)



d <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample.csv")
d10 <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_10000.csv")

fits_split <- d %>% 
  # sample_n(size = 500, replace = FALSE) %>% 
  filter(rsqr.list > 0.98) %>% 
  split(.$curve.id.list)

prediction_function <- function(curve1) {
  x <- seq(-3, 38, 0.1)
  predictions <- curve1$a.list[[1]]*exp(curve1$b.list[[1]]*x)*(1-((x-curve1$z.list[[1]])/(curve1$w.list[[1]]/2))^2)
  data.frame(x, predictions)
}

## make predictions for each set of parameter estimates. 
all_predictions <- fits_split %>% 
  map_df(prediction_function) 

## get the upper and lower limits of the predicted growth rates at each value of x (temperature)
boot_limits <- all_predictions %>% 
  group_by(x) %>% 
  summarise(q2.5=quantile(predictions, probs=0.025),
            q97.5=quantile(predictions, probs=0.975),
            mean = mean(predictions)) 


fits_split10 <- d10 %>% 
  # sample_n(size = 500, replace = FALSE) %>% 
  filter(rsqr.list > 0.98) %>% 
  split(.$curve.id.list)

prediction_function <- function(curve1) {
  x <- seq(-3, 38, 0.1)
  predictions <- curve1$a.list[[1]]*exp(curve1$b.list[[1]]*x)*(1-((x-curve1$z.list[[1]])/(curve1$w.list[[1]]/2))^2)
  data.frame(x, predictions)
}

## make predictions for each set of parameter estimates. 
all_predictions10 <- fits_split10 %>% 
  map_df(prediction_function) 




boot_limits_10 <- all_predictions10 %>% 
  group_by(x) %>% 
  summarise(q2.5=quantile(predictions, probs=0.025),
            q97.5=quantile(predictions, probs=0.975),
            mean = mean(predictions)) 


