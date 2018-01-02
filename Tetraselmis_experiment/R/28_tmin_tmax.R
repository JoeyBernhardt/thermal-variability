
### Find tmin and tmax
library(rootSolve)
library(tidyverse)
library(cowplot)

fits_constant <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample.csv")

fits2 <- fits_constant %>% 
  rename(z = z.list,
         a = a.list, 
         b = b.list, 
         w = w.list)

nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

cf4 <- fits2 %>%
  group_by(curve.id.list) %>% 
  mutate(tmax = ifelse(length(uniroot.all(function(x) nbcurve(x, z, w, a, b),c(topt.list,150)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w,  a, b),c(topt.list,150)))) %>% 
  mutate(tmin = ifelse(length(uniroot.all(function(x) nbcurve(x, z,w, a, b),c(-1.8,topt.list)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w, a, b),c(-1.8,topt.list))))



fits_real_constant <- cf4 %>% 
  mutate(tmin = ifelse(is.na(tmin), -1.8, tmin)) %>% 
  ungroup()

cf4 %>% 
  filter(!is.na(tmin)) %>% View

fits_real_constant %>% 
  ggplot(aes(x = tmin)) + geom_histogram(bins = 10)


write_csv(fits_real_constant, "Tetraselmis_experiment/data-processed/fits_real_constant.csv")


fits_above_freezing_constant <- fits_real_constant %>% 
  mutate(breadth = tmax - tmin) %>%
  summarise(w_low=quantile(w, probs=0.025),
            w_high=quantile(w, probs=0.975),
            tmin_low=quantile(tmin, probs=0.025),
            tmin_high=quantile(tmin, probs=0.975),
            tmax_low=quantile(tmax, probs=0.025),
            tmax_high=quantile(tmax, probs=0.975),
            breadth_low=quantile(breadth, probs=0.025),
            breadth_high=quantile(breadth, probs=0.975)) 

View(fits_above_freezing_constant)
### Now for the variable curves


fits_variable <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v.csv")

fits3 <- fits_variable %>% 
  rename(z = z.list,
         a = a.list, 
         b = b.list, 
         w = w.list)


cf5 <- fits3 %>%
  group_by(curve.id.list) %>% 
  mutate(tmax = ifelse(length(uniroot.all(function(x) nbcurve(x, z, w, a, b),c(topt.list,150)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w,  a, b),c(topt.list,150)))) %>% 
  mutate(tmin = ifelse(length(uniroot.all(function(x) nbcurve(x, z,w, a, b),c(-2,topt.list)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w, a, b),c(-2,topt.list))))


fits_real_variable <- cf5 %>% 
  filter(!is.na(tmin)) %>%
  filter(!is.na(tmax)) %>% 
  ungroup()

write_csv(fits_real_variable, "Tetraselmis_experiment/data-processed/fits_real_variable.csv")

fits_real_variable %>% 
  summarise(q2.5=quantile(w, probs=0.025),
            q97.5=quantile(w, probs=0.975),
            mean = mean(w)) %>% View

fits_above_freezing_variable <- fits_real_variable %>% 
  mutate(breadth = tmax - tmin) %>%
  summarise(w_low=quantile(w, probs=0.025),
            w_high=quantile(w, probs=0.975),
            tmin_low=quantile(tmin, probs=0.025),
            tmin_high=quantile(tmin, probs=0.975),
            tmax_low=quantile(tmax, probs=0.025),
            tmax_high=quantile(tmax, probs=0.975)) 
  
View(fits_above_freezing_variable)
