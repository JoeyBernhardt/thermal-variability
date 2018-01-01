
### Find tmin and tmax
library(rootSolve)

fits <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample.csv")
for (i in 1:length(curve.id.list)){
  print(i)
  nbcurve.tmax<-function(x){
    nb<-nbcurve(x,fits$z.list[i],fits$w.list[i],fits$a.list[i],fits$b.list[i])
    nb
  }
  fits$tmax[i]<-uniroot.all(nbcurve.tmax,c(fits$topt.list[i],150))[1]
  fits$tmin[i]<-uniroot.all(nbcurve.tmax,c(-2,fits$topt.list[i]))[1] 
}



fits_real_constant <- fits %>% 
  filter(!is.na(tmin)) 


fits_real_constant %>% 
  summarise(q2.5=quantile(w.list, probs=0.025),
            q97.5=quantile(w.list, probs=0.975),
            mean = mean(w.list)) %>% View


### Now for the variable curves


fits <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v.csv")
for (i in 1:length(curve.id.list)){
  print(i)
  nbcurve.tmax<-function(x){
    nb<-nbcurve(x,fits$z.list[i],fits$w.list[i],fits$a.list[i],fits$b.list[i])
    nb
  }
  fits$tmax[i]<-uniroot.all(nbcurve.tmax,c(fits$topt.list[i],150))[1]
  fits$tmin[i]<-uniroot.all(nbcurve.tmax,c(-2,fits$topt.list[i]))[1] 
}



fits_real_variable <- fits %>% 
  filter(!is.na(tmin)) 


fits_real_variable %>% 
  summarise(q2.5=quantile(w.list, probs=0.025),
            q97.5=quantile(w.list, probs=0.975),
            mean = mean(w.list)) %>% View

fits_real_variable %>% 
  mutate(breadth = tmax - tmin) %>%
  

