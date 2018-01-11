### limiting resampling exploring


library(tidyverse)
library(colormap)
library(cowplot)
library(rootSolve)


d <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample.csv")
d10 <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_10000_alltemps.csv") %>% 
  filter(rsqr.list > 0.98)
v10 <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_10000_v.csv") %>% 
  filter(rsqr.list > 0.98)
table(v10$n.list)
table(d$n.list)


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

## now for 10,000
fits_split10 <- d10 %>% 
  filter(rsqr.list > 0.98) %>% 
  # sample_n(size = 1790, replace = FALSE) %>% 
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

fits10 <- d10 %>% 
  rename(z = z.list,
         a = a.list, 
         b = b.list, 
         w = w.list)


cf10 <- fits10 %>%
  filter(rsqr.list > 0.98) %>% 
  group_by(curve.id.list) %>% str()
  mutate(tmax = ifelse(length(uniroot.all(function(x) nbcurve(x, z, w, a, b),c(topt.list,150)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w,  a, b),c(topt.list,150)))) %>% 
  mutate(tmin = ifelse(length(uniroot.all(function(x) nbcurve(x, z,w, a, b),c(-1.8,topt.list)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w, a, b),c(-1.8,topt.list))))


fits_10 <- cf10 %>% 
  filter(!is.na(tmax)) %>% 
  mutate(tmin = ifelse(is.na(tmin), -1.8, tmin)) %>% 
  ungroup()


fits_above_freezing_constant_10 <- fits_10 %>% 
  mutate(breadth = tmax - tmin) %>%
  summarise(w_low=quantile(w, probs=0.025),
            w_high=quantile(w, probs=0.975),
            w_mean=mean(w),
            tmin_low=quantile(tmin, probs=0.025),
            tmin_high=quantile(tmin, probs=0.975),
            tmin_mean=mean(tmin),
            tmax_low=quantile(tmax, probs=0.025),
            tmax_high=quantile(tmax, probs=0.975),
            tmax_mean=mean(tmax),
            breadth_low=quantile(breadth, probs=0.025),
            breadth_high=quantile(breadth, probs=0.975),
            breadth_mean=mean(breadth),
            topt_low=quantile(topt.list, probs=0.025),
            topt_high=quantile(topt.list, probs=0.975),
            topt_mean=mean(topt.list)) %>% 
  gather(key = "metric", value = "value_incomplete") 
write_csv(fits_above_freezing_constant_10, "Tetraselmis_experiment/data-processed/fits_above_freezing_constant_10.csv")



# now for variable --------------------------------------------------------

fits_splitv <- v10 %>% 
  # sample_n(size = 500, replace = FALSE) %>% 
  filter(rsqr.list > 0.98) %>% 
  split(.$unique_id)

prediction_function <- function(curve1) {
  x <- seq(-3, 38, 0.1)
  predictions <- curve1$a.list[[1]]*exp(curve1$b.list[[1]]*x)*(1-((x-curve1$z.list[[1]])/(curve1$w.list[[1]]/2))^2)
  data.frame(x, predictions)
}

## make predictions for each set of parameter estimates. 
all_predictionsv <- fits_splitv%>% 
  map_df(prediction_function) 

## get the upper and lower limits of the predicted growth rates at each value of x (temperature)
boot_limitsv <- all_predictionsv %>% 
  group_by(x) %>% 
  summarise(q2.5=quantile(predictions, probs=0.025),
            q97.5=quantile(predictions, probs=0.975),
            mean = mean(predictions)) 

## now for 10,000
fits_split10v <- v10 %>% 
  filter(rsqr.list > 0.98) %>% 
  split(.$unique_id)

prediction_function <- function(curve1) {
  x <- seq(-3, 38, 0.1)
  predictions <- curve1$a.list[[1]]*exp(curve1$b.list[[1]]*x)*(1-((x-curve1$z.list[[1]])/(curve1$w.list[[1]]/2))^2)
  data.frame(x, predictions)
}

## make predictions for each set of parameter estimates. 
all_predictions10v <- fits_split10v %>% 
  map_df(prediction_function) 




boot_limits_10v <- all_predictions10v %>% 
  group_by(x) %>% 
  summarise(q2.5=quantile(predictions, probs=0.025),
            q97.5=quantile(predictions, probs=0.975),
            mean = mean(predictions)) 

fits10v <- v10 %>% 
  rename(z = z.list,
         a = a.list, 
         b = b.list, 
         w = w.list)


cf10v <- fits10v %>%
  filter(rsqr.list > 0.98) %>% 
  group_by(unique_id) %>% 
  mutate(tmax = ifelse(length(uniroot.all(function(x) nbcurve(x, z, w, a, b),c(topt.list,150)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w,  a, b),c(topt.list,150)))) %>% 
  mutate(tmin = ifelse(length(uniroot.all(function(x) nbcurve(x, z,w, a, b),c(-1.8,topt.list)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w, a, b),c(-1.8,topt.list))))




fits_10v <- cf10v %>% 
  filter(!is.na(tmax)) %>% 
  mutate(tmin = ifelse(is.na(tmin), -1.8, tmin)) %>% 
  ungroup()


fits_above_freezing_constant_10v <- fits_10v %>% 
  mutate(breadth = tmax - tmin) %>%
  summarise(w_low=quantile(w, probs=0.025),
            w_high=quantile(w, probs=0.975),
            w_mean=mean(w),
            tmin_low=quantile(tmin, probs=0.025),
            tmin_high=quantile(tmin, probs=0.975),
            tmin_mean=mean(tmin),
            tmax_low=quantile(tmax, probs=0.025),
            tmax_high=quantile(tmax, probs=0.975),
            tmax_mean=mean(tmax),
            breadth_low=quantile(breadth, probs=0.025),
            breadth_high=quantile(breadth, probs=0.975),
            breadth_mean=mean(breadth),
            topt_low=quantile(topt.list, probs=0.025),
            topt_high=quantile(topt.list, probs=0.975),
            topt_mean=mean(topt.list)) %>% 
  gather(key = "metric", value = "value_incomplete") 
write_csv(fits_above_freezing_constant_10, "Tetraselmis_experiment/data-processed/fits_above_freezing_variable_10.csv")

### now just with full 9 temps
fits_split9 <- d10 %>% 
  filter(rsqr.list > 0.98, n.list == 9) %>% 
  split(.$curve.id.list)

prediction_function <- function(curve1) {
  x <- seq(-3, 38, 0.1)
  predictions <- curve1$a.list[[1]]*exp(curve1$b.list[[1]]*x)*(1-((x-curve1$z.list[[1]])/(curve1$w.list[[1]]/2))^2)
  data.frame(x, predictions)
}

## make predictions for each set of parameter estimates. 
all_predictions9 <- fits_split9 %>% 
  map_df(prediction_function) 




boot_limits_9 <- all_predictions9 %>% 
  group_by(x) %>% 
  summarise(q2.5=quantile(predictions, probs=0.025),
            q97.5=quantile(predictions, probs=0.975),
            mean = mean(predictions)) 

fits9 <-  d10 %>% 
  filter(rsqr.list > 0.98, n.list == 9) %>% 
  rename(z = z.list,
         a = a.list, 
         b = b.list, 
         w = w.list)


cf9 <- fits9 %>%
  group_by(curve.id.list) %>% 
  mutate(tmax = ifelse(length(uniroot.all(function(x) nbcurve(x, z, w, a, b),c(topt.list,150)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w,  a, b),c(topt.list,150)))) %>% 
  mutate(tmin = ifelse(length(uniroot.all(function(x) nbcurve(x, z,w, a, b),c(-1.8,topt.list)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w, a, b),c(-1.8,topt.list))))


fits_9 <- cf9 %>% 
  filter(!is.na(tmax)) %>% 
  mutate(tmin = ifelse(is.na(tmin), -1.8, tmin)) %>% 
  ungroup()


fits_above_freezing_constant_9 <- fits_9 %>% 
  mutate(breadth = tmax - tmin) %>%
  summarise(w_low=quantile(w, probs=0.025),
            w_high=quantile(w, probs=0.975),
            w_mean=mean(w),
            tmin_low=quantile(tmin, probs=0.025),
            tmin_high=quantile(tmin, probs=0.975),
            tmin_mean=mean(tmin),
            tmax_low=quantile(tmax, probs=0.025),
            tmax_high=quantile(tmax, probs=0.975),
            tmax_mean=mean(tmax),
            breadth_low=quantile(breadth, probs=0.025),
            breadth_high=quantile(breadth, probs=0.975),
            breadth_mean=mean(breadth),
            topt_low=quantile(topt.list, probs=0.025),
            topt_high=quantile(topt.list, probs=0.975),
            topt_mean=mean(topt.list)) %>% 
  gather(key = "metric", value = "value_complete") 

left_join(fits_above_freezing_constant_10, fits_above_freezing_constant_9) %>% View


#### now with variable
v <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v_10000_2547.csv")

fits_splitv <- v %>% 
  filter(rsqr.list > 0.98, z.list >0) %>% 
  split(.$curve.id.list)

prediction_function <- function(curve1) {
  x <- seq(-3, 38, 0.1)
  predictions <- curve1$a.list[[1]]*exp(curve1$b.list[[1]]*x)*(1-((x-curve1$z.list[[1]])/(curve1$w.list[[1]]/2))^2)
  data.frame(x, predictions)
}

## make predictions for each set of parameter estimates. 
all_predictionsv <- fits_splitv %>% 
  map_df(prediction_function) 




boot_limits_v <- all_predictionsv %>% 
  group_by(x) %>% 
  summarise(q2.5=quantile(predictions, probs=0.025),
            q97.5=quantile(predictions, probs=0.975),
            mean = mean(predictions)) 

fitsv <-  v %>% 
  filter(rsqr.list > 0.98, z.list >0) %>% 
  rename(z = z.list,
         a = a.list, 
         b = b.list, 
         w = w.list)


cfv <- fitsv %>%
  group_by(curve.id.list) %>% 
  mutate(tmax = ifelse(length(uniroot.all(function(x) nbcurve(x, z, w, a, b),c(topt.list,150)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w,  a, b),c(topt.list,150)))) %>% 
  mutate(tmin = ifelse(length(uniroot.all(function(x) nbcurve(x, z,w, a, b),c(-1.8,topt.list)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w, a, b),c(-1.8,topt.list))))


fits_v <- cfv %>% 
  filter(!is.na(tmax)) %>% 
  mutate(tmin = ifelse(is.na(tmin), -1.8, tmin)) %>% 
  ungroup()

cfv %>% 
  filter(!is.na(tmin)) %>% View

fits_above_freezing_v <- fits_v %>% 
  mutate(breadth = tmax - tmin) %>%
  summarise(w_low=quantile(w, probs=0.025),
            w_high=quantile(w, probs=0.975),
            w_mean=mean(w),
            tmin_low=quantile(tmin, probs=0.025),
            tmin_high=quantile(tmin, probs=0.975),
            tmin_mean=mean(tmin),
            tmax_low=quantile(tmax, probs=0.025),
            tmax_high=quantile(tmax, probs=0.975),
            tmax_mean=mean(tmax),
            breadth_low=quantile(breadth, probs=0.025),
            breadth_high=quantile(breadth, probs=0.975),
            breadth_mean=mean(breadth),
            topt_low=quantile(topt.list, probs=0.025),
            topt_high=quantile(topt.list, probs=0.975),
            topt_mean=mean(topt.list)) %>% 
  gather(key = "metric", value = "value_complete") 


thomas %>% 
  filter(used.for.tmin.analysis == 1) %>% View
