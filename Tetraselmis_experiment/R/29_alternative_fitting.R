
library(bbmle)
library(rootSolve)
library(purrr)
library(tidyverse)
library(broom)


growth_all <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling.csv")
growth_all_v <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling_v.csv")

dat.full <- growth_all %>%
  group_by(temp) %>% 
  summarise(growth.rate = mean(growth_per_day)) %>% 
  rename(temperature = temp) %>% 
  mutate(curve.id = "1") %>% 
  select(curve.id, temperature, growth.rate) 





nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

# Create new vector of unique curve.id values
curve.id.list<-unique(dat.full$curve.id)	

# Create empty vectors to populate with parameter values and trait estimates
z.list<-rep(NA, length(curve.id.list))				#Parameter 'z'
w.list<-rep(NA, length(curve.id.list))				#Parameter 'w', which is the temperature niche width
a.list<-rep(NA, length(curve.id.list))				#Parameter 'a'
b.list<-rep(NA, length(curve.id.list))				#Parameter 'b'
topt.list<-rep(NA, length(curve.id.list))			#Topt, Optimum temperature for growth
maxgrowth.list<-rep(NA, length(curve.id.list))		#Maximum growth rate (i.e. growth rate at Topt)
rsqr.list<-rep(NA, length(curve.id.list))			#R^2 values for all fits
s.list<-rep(NA, length(curve.id.list))				#Error
n.list<-rep(NA, length(curve.id.list))				#Number of growth rate measurements used in the curve


  
  # Take a subset of the data corressponding to the ith curve.id.list value
  dat<-subset(dat.full,dat.full$curve.id==curve.id.list[1])
  
  # guess starting values for parameters 'z' and 'w'
  z.guess<-mean(dat$temperature[dat$growth.rate==max(dat$growth.rate)])		#starting estimates for 'z'
  w.guess<-diff(range(dat$temperature))										#starting estimates for niche width
  
  ## This loop fits the model using a range of different starting guesses. We choose the best one using AIC. This helps find good solutions even if there are
  # convergence problems.
  # Starting estimates for parameters 'a' and 'b' use a plausible range but with broadly spaced estimates to speed up fitting. 
  avals<-seq(-0.2,1.2,0.05)		
  bvals<-seq(-0.2,0.3,0.05)
  mod.list<-list()
  AIC.list<-c()
  
  for(ia in 1:length(avals)){
    for(ib in 1:length(bvals)){
      a.guess<-avals[ia]
      b.guess<-bvals[ib]
      res2<-try(fit<-mle2(dat$growth.rate~dnorm(mean=nbcurve(dat$temperature,z=z,w=w,a=a,b=b),sd=s),start=list(z=z.guess,w=w.guess,a=a.guess,b=b.guess,s=0.3),
                          skip.hessian=TRUE,data=dat))
      if(class(res2)!="try-error"){
        mod.list<-append(mod.list,fit)
        AIC.list<-append(AIC.list,AIC(fit))
      }
    }
  }

  coeffs <- mod.list %>% 
  map_df(.f = tidy, .id = "id") %>% 
  group_by(id) %>% 
  select(term, estimate) %>% 
  spread(key = term, value = estimate)

## now get the tmin and tmaxes
cf3 <- coeffs %>% 
   mutate(tmax = ifelse(length(uniroot.all(function(x) nbcurve(x, z, w, a, b),c(15,150)))==0, NA,
                        uniroot.all(function(x) nbcurve(x, z, w,  a, b),c(15,150)))) %>% 
   mutate(tmin = ifelse(length(uniroot.all(function(x) nbcurve(x, z,w, a, b),c(-1.8,10)))==0, NA,
                        uniroot.all(function(x) nbcurve(x, z, w, a, b),c(-1.8,10))))


aics <- as.data.frame(AIC.list) %>% 
  mutate(id = rownames(.))

### pick only fits that have tmins above -2C
best_mod <- left_join(cf3, aics, by = "id") %>% 
  filter(!is.na(tmin)) %>% 
  ungroup() %>% 
  dplyr::top_n(., n = -1, wt = AIC.list)
  
  
### plot it!
curveid <- 1224  
  
  p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
  p + stat_function(fun = function(x) nbcurve(x, cf3$z[cf3$id == curveid], cf3$w[cf3$id == curveid], cf3$a[cf3$id == curveid], cf3$b[cf3$id == curveid]), size = 1) +
    ylim(-15, 3) + xlim(-15, 55)
  
expected <-nbcurve(dat$temperature,best_mod$z,best_mod$w,best_mod$a,best_mod$b)
rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)  

  
  # Use the curve fit to find Topt and the estimated maximum growth rate (i.e. growth rate at Topt)
  grfunc<-function(x){
    -nbcurve(x,best_mod$z,best_mod$w,best_mod$a,best_mod$b)
  }
  optinfo<-optim(c(x=best_mod$z[[1]]),grfunc)
  opt<-optinfo$par[[1]]
  maxgrowth<- -optinfo$value
  
 params_constant <- data.frame(opt, maxgrowth, best_mod) %>% 
   mutate(treatment = "constant")
 params_variable <- data.frame(opt, maxgrowth, best_mod) 
 
 params_variable <- params_variable %>% 
   mutate(treatment = "variable")
 
 all_params_above_freezing <- bind_rows(params_constant, params_variable) 
   write_csv(all_params_above_freezing, "Tetraselmis_experiment/data-processed/all_params_above_freezing.csv")
   
   
all_params_above_freezing <- read_csv("Tetraselmis_experiment/data-processed/all_params_above_freezing.csv")   
View(all_params_above_freezing)   
   
   
 