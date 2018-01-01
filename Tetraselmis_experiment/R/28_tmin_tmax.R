
### Find tmin and tmax
library(rootSolve)

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
  mutate(tmin = ifelse(length(uniroot.all(function(x) nbcurve(x, z,w, a, b),c(-2,topt.list)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w, a, b),c(-2,topt.list))))



fits_real_constant <- cf4 %>% 
  filter(!is.na(tmin)) %>% 
  ungroup()


fits_real_constant %>% 
  summarise(q2.5=quantile(w, probs=0.025),
            q97.5=quantile(w, probs=0.975),
            mean = mean(w)) %>% View


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
  ungroup()

fits_real_variable %>% 
  summarise(q2.5=quantile(w, probs=0.025),
            q97.5=quantile(w, probs=0.975),
            mean = mean(w)) %>% View

fits_real_variable %>% 
  mutate(breadth = tmax - tmin) %>% View
  

