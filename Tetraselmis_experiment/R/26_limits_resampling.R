library(tidyverse)
library(cowplot)

### get the upper and lower 95% on Tmax

fits_constant <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample.csv")

## ok now take the fits, and make prediction curves and then take the 97.5 and 2.5% CIs

## split up the fits df by curve id
fits_split <- fits %>% 
	# filter(a.list > 0, z.list > 0) %>%
	split(.$curve.id.list)

prediction_function <- function(curve1) {
	x <- seq(-15, 38, 0.01)
	predictions <- curve1$a.list[[1]]*exp(curve1$b.list[[1]]*x)*(1-((x-curve1$z.list[[1]])/(curve1$w.list[[1]]/2))^2)
	data.frame(x, predictions)
}

## make predictions for each set of parameter estimates. 
all_predictions <- fits_split %>% 
	map_df(prediction_function) 

all_predictions %>% 
	group_by(x) %>% 
	summarise(q2.5=quantile(predictions, probs=0.025),
						q97.5=quantile(predictions, probs=0.975),
						mean = mean(predictions)) %>%
	# filter(x < 0) %>% 
	ggplot(aes(x = x, y = q97.5)) + geom_point() + geom_hline(yintercept = 0) +
	xlim(-10, 0) + ylim(-0.01, 0.01)
## tmax = 31.25, 32.92
##tmin = -1.25, 

### find the upper and lower boundary curves for the constant conditions

fits_constant %>% 
  mutate(q2.5=quantile(topt.list, probs=0.025),
         q97.5=quantile(topt.list, probs=0.975),
         q50=quantile(topt.list, probs=0.50),
         mean = mean(topt.list),
         median = median(topt.list)) %>% 
  filter(topt.list > 24.57) %>% 
  filter(topt.list <24.575) %>% View

## ok curve 695 is the lower 2.5, curve 82 is the upper 97.5, 468 is the middle curve. 


fits_constant %>% 
  filter(curve.id.list %in% c(695, 82, 468)) %>% View


### now for variable

fits <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v.csv")
fits_c <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params.csv")
fits_v <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_v.csv")
## ok now take the fits, and make prediction curves and then take the 97.5 and 2.5% CIs

## split up the fits df by curve id
fits_split <- fits %>% 
	filter(a.list >0) %>% 
	split(.$curve.id.list)

prediction_function <- function(curve1) {
	x <- seq(-10, 38, 0.01)
	predictions <- curve1$a.list[[1]]*exp(curve1$b.list[[1]]*x)*(1-((x-curve1$z.list[[1]])/(curve1$w.list[[1]]/2))^2)
	data.frame(x, predictions)
}

## make predictions for each set of parameter estimates. 
all_predictions <- fits_split %>% 
	map_df(prediction_function) 

all_predictions %>% 
	filter(x < 0) %>% 
  ggplot(aes(x = x, y = predictions) + geom_point()

## 31.90, -4.47

preds <- all_predictions %>% 
	group_by(x) %>% 
	summarise(q2.5=quantile(predictions, probs=0.025),
						q97.5=quantile(predictions, probs=0.975),
						mean = mean(predictions))



### tmax = 28.66, 32.77


ggplot(aes(x = x, y = q2.5)) + geom_point() + geom_hline(yintercept = 0) +
	xlim(30.5, 34)



### let's just get the actual fitted 0

fits_v <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_v.csv")
fits_c <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params.csv")


### or do this by taking the tmax associated with the highest and lowest 95% CI on Topt

fits <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v.csv")

fits %>% 
  mutate(q2.5=quantile(topt.list, probs=0.025),
            q97.5=quantile(topt.list, probs=0.975),
         q50=quantile(topt.list, probs=0.50),
            mean = mean(topt.list),
         median = median(topt.list)) %>% 
  filter(topt.list > 21.468) %>% 
  filter(topt.list <21.4684) %>% View
  
### 521 is the lower limit curve, curve.id 81 is the upper limit curve, curve.id 773 is the median

## ok so now let's get the tmaxes on these 

### here are the 2.5, 50 and 97.5% quantile curves, used to get the 95% CI on the parameter estimates. 
fits_constant %>% 
  filter(curve.id.list %in% c(81, 521, 773)) %>% View

fits_constant %>% 
  filter(curve.id.list == 521) %>% View
fits_constant %>% 
  filter(curve.id.list == 773) %>% View


fits_variable <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v.csv")



  x <- seq(-10, 38, 0.01)
  predictions <- fits$a.list[fits$curve.id.list==521]*exp(fits$b.list[fits$curve.id.list==521]*x)*(1-((x-fits$z.list[fits$curve.id.list==521])/(fits$w.list[fits$curve.id.list==521]/2))^2)
  upper_curve <- data.frame(x, predictions)
  
  
  x <- seq(-10, 38, 0.01)
  predictions <- fits$a.list[fits$curve.id.list==773]*exp(fits$b.list[fits$curve.id.list==773]*x)*(1-((x-fits$z.list[fits$curve.id.list==773])/(fits$w.list[fits$curve.id.list==773]/2))^2)
  mid_curve <- data.frame(x, predictions)
  
  upper_curve %>% 
    filter(x < 5) %>% View
  
  ##tmax is 30.36, tmin is 0.44
  
  
  x <- seq(-20, 38, 0.01)
  lower_predictions <- fits$a.list[fits$curve.id.list==81]*exp(fits$b.list[fits$curve.id.list==81]*x)*(1-((x-fits$z.list[fits$curve.id.list==81])/(fits$w.list[fits$curve.id.list==81]/2))^2)
  lower_curve <- data.frame(x, lower_predictions)
  
  lower_curve %>% 
    filter(x < 5) %>% View
## tmax is 31.08

  
  growth_all_v <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling_v.csv")
  fits_v <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v.csv")
  
  curve_variable_resamp<-function(x){
    res<-fits_v$a.list[1]*exp(fits_v$b.list[1]*x)*(1-((x-fits_v$z.list[1])/(fits_v$w.list[1]/2))^2)
    res
  }
  
  
  lower_curve %>% 
    ggplot(aes(x = x, y = lower_predictions)) + geom_line() + xlim(-10, 35) + ylim(-0.05, 2) +
    geom_line(aes(x = x, y = predictions), data = upper_curve) +
    geom_line(aes(x = x, y = predictions), data = mid_curve, color = "red") +
    geom_vline(xintercept = 20.27451) +
    geom_vline(xintercept = 21.46802) +
    geom_vline(xintercept = 22.52605) +
    geom_point(aes(x = temp, y = growth_per_day), data = growth_all_v) +
    stat_function(fun = curve_variable_resamp, color = "blue", size = 2)
  
  fits$w.list[fits$curve.id.list == 81]
  fits$w.list[fits$curve.id.list == 773]
  fits$w.list[fits$curve.id.list == 521]
  
  
  fits %>% 
    filter(topt.list < 30) %>% 
    ggplot(aes(x = topt.list)) + geom_histogram() +
    geom_vline(xintercept = 20.27451) +
    geom_vline(xintercept = 21.46802) +
    geom_vline(xintercept = 22.52605)
  
  all_predictions %>% 
    ggplot(aes(x = x, y = predictions)) + geom_point()
