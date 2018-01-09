### limiting resampling exploring


library(tidyverse)
library(colormap)
library(cowplot)
library(rootSolve)


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

ic <- colormap(colormap = colormaps$viridis, nshades = 8, format = "hex",
               alpha = 1, reverse = FALSE)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 
p + 
  geom_vline(xintercept = 24.5681, color = ic[3], alpha = 0.5) +
  geom_hline(yintercept = 0, color = "grey") +
  # geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits, fill = ic[3], alpha = 1) +
  geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_10, fill = ic[4], alpha = 0.5) +
  coord_cartesian(ylim = c(-0.2, 1.7), xlim = c(0, 32))
  


fits10 <- d10 %>% 
  rename(z = z.list,
         a = a.list, 
         b = b.list, 
         w = w.list)


cf10 <- fits10 %>%
  group_by(curve.id.list) %>% 
  mutate(tmax = ifelse(length(uniroot.all(function(x) nbcurve(x, z, w, a, b),c(topt.list,150)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w,  a, b),c(topt.list,150)))) %>% 
  mutate(tmin = ifelse(length(uniroot.all(function(x) nbcurve(x, z,w, a, b),c(-1.8,topt.list)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w, a, b),c(-1.8,topt.list))))


fits_10 <- cf10 %>% 
  filter(!is.na(tmax)) %>% 
  mutate(tmin = ifelse(is.na(tmin), -1.8, tmin)) %>% 
  ungroup()


fits_above_freezing_constant_10 <- fits_real %>% 
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
            topt_mean=mean(topt.list)) 
write_csv(fits_above_freezing_constant_10, "Tetraselmis_experiment/data-processed/fits_above_freezing_constant_10.csv")
