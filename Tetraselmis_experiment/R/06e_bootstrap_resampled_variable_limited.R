### Parametric bootstrapping, now with resampled growth rates


# load packages -----------------------------------------------------------


library(purrr)
library(bbmle)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggthemes)
library(broom)
library(tidyverse)
library(rootSolve)


# load data ---------------------------------------------------------------


params_raw <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_v.csv") ## estimated TPC parameters for constant conditions
growth_sum <- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary_v.csv") ## empirically observed growth rates

ps <- read_csv("Tetraselmis_experiment/data-processed/all_params_above_freezing.csv")
ps_v <- ps %>% 
  filter(treatment == "variable")

params <- ps_v %>% 
  select(z, w, a, b) %>%
  gather(key = term, value = estimate) 

# params <- params_raw %>% 
#   select(z.list, w.list, a.list, b.list) %>%
#   rename(z = z.list, 
#          w = w.list, 
#          a = a.list, 
#          b = b.list) %>% 
#   gather(key = term, value = estimate) 



growth2 <- growth_sum
# get the fitted values along the TPC (i.e. Eqn S.1 in Thomas et a --------

## this is the TPC function
nbcurve1<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

predicted_data <- function(temp) {
  predicted_growth <- nbcurve1(temp, params$estimate[params$term== "z"], 
                               params$estimate[params$term== "w"],
                               params$estimate[params$term== "a"],
                               params$estimate[params$term== "b"])
  return(predicted_growth)
}
temps <- growth2 %>% 
  select(temp) %>% 
  arrange(temp) %>% 
  split(.$temp)


predictions <- temps %>% 
  map_df(predicted_data, .id = "temperature") %>% 
  rename(predicted_growth = temp)


# Now pick out new values for growth rate (S.3) ---------------------------------
## at first I had sd 

sample_size <- 1

EqnS.3 <- function(sample_size){
  x_5 <- rnorm(n = sample_size, sd = ((growth2$sd[growth2$temp == 5])^2*(((9-1)/rchisq(n = 1, df = 8))^0.5))^0.5,
               mean = predictions$predicted_growth[predictions$temperature ==5])
  x_10 <- rnorm(n = sample_size, sd = ((growth2$sd[growth2$temp == 10])^2*(((9-1)/rchisq(n = 1, df = 8))^0.5))^0.5,
                mean = predictions$predicted_growth[predictions$temperature ==10])
  x_15 <- rnorm(n = sample_size, sd = ((growth2$sd[growth2$temp == 15])^2*(((9-1)/rchisq(n = 1, df = 8))^0.5))^0.5,
                mean = predictions$predicted_growth[predictions$temperature ==15])
  x_20 <- rnorm(n = sample_size, sd = ((growth2$sd[growth2$temp == 20])^2*(((9-1)/rchisq(n = 1, df = 8))^0.5))^0.5,
                mean = predictions$predicted_growth[predictions$temperature ==20])
  x_24 <- rnorm(n = sample_size, sd = ((growth2$sd[growth2$temp == 24])^2*(((9-1)/rchisq(n = 1, df = 8))^0.5))^0.5,
                mean = predictions$predicted_growth[predictions$temperature ==24])
  x_27 <- rnorm(n = sample_size, sd = ((growth2$sd[growth2$temp == 27])^2*(((9-1)/rchisq(n = 1, df = 8))^0.5))^0.5,
                mean = predictions$predicted_growth[predictions$temperature ==27])
  data.frame(x_5, x_10, x_15, x_20, x_24, x_27)
}

samples <- rep(1, 1000)

## generate all our new synthetic datasets to which we will fit our TPCs
dat.full <- samples %>% 
  map_df(EqnS.3, .id = "curve.id") %>% 
  gather(key = "temperature", value = "growth.rate", starts_with("x")) %>% 
  separate(temperature, into = c("x", "temperature")) %>% 
  select(-x) %>% 
  filter(growth.rate >=0) %>% 
  mutate(temperature = as.numeric(temperature))

## store a mini dataframe for plotting later
data_full <- dat.full %>% 
  filter(curve.id == "1") %>% 
  mutate(temperature = as.numeric(temperature))

# fit the data! -----------------------------------------------------------


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

### Loop through all curve.id.list values to estimate parameters for all curves
### on my computer this takes QUITE a while! This code is taken directly from
### Mridul Thomas, http://mridulkthomas.weebly.com/data--code.html

## to skip this step, just read in the outputs of this step:
## ("Tetraselmis_experiment/data-processed/boot_fits.csv")


for(i in 1:length(curve.id.list)){
  print(i)
  
  # Take a subset of the data corressponding to the ith curve.id.list value
  dat<-subset(dat.full,dat.full$curve.id==curve.id.list[i])
  
  # guess starting values for parameters 'z' and 'w'
  z.guess<-mean(dat$temperature[dat$growth.rate==max(dat$growth.rate)])		#starting estimates for 'z'
  w.guess<-diff(range(dat$temperature))										#starting estimates for niche width
  
  ## This loop fits the model using a range of different starting guesses. We choose the best one using AIC. This helps find good solutions even if there are
  # convergence problems.
  # Starting estimates for parameters 'a' and 'b' use a plausible range but with broadly spaced estimates to speed up fitting. 
  avals<-seq(-0.2,1.2,0.1)		
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
  bestmod <- left_join(cf3, aics, by = "id") %>% 
    filter(!is.na(tmin)) %>% 
    ungroup() %>% 
    dplyr::top_n(., n = -1, wt = AIC.list)
  
  cfs<-c(bestmod$z[[1]], bestmod$w[[1]], bestmod$a[[1]], bestmod$b[[1]])
  expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
  rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)
  
  # # Identify the best model from the list and save coefficients and R^2 values
  # if(!is.null(AIC.list)){
  # 	bestmodind<-which(AIC.list==min(AIC.list))
  # 	if(length(bestmodind)>1){
  # 		bestmodind<-sample(bestmodind,1)
  # 	}
  # 	bestmod<-mod.list[[bestmodind]]
  # 	cfs<-coef(bestmod)
  # 	expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
  # 	rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)
  # }
  # 
  # # If the quick fit yielded poor results (low R^2), try a more thorough search through parameter space
  # if(rsqr<0.95){
  # 	avals<-seq(-0.2,1.2,0.02)
  # 	bvals<-seq(-0.2,0.3,0.02)
  # 	mod.list<-list()
  # 	AIC.list<-c()
  # 	for(ia in 1:length(avals)){
  # 		for(ib in 1:length(bvals)){
  # 			a.guess<-avals[ia]
  # 			b.guess<-bvals[ib]
  # 			res2<-try(fit<-mle2(dat$growth.rate~dnorm(mean=nbcurve(dat$temperature,z=z,w=w,a=a,b=b),sd=s),start=list(z=z.guess,w=w.guess,a=a.guess,b=b.guess,s=0.3),
  # 													skip.hessian=TRUE,data=dat))
  # 			if(class(res2)!="try-error"){
  # 				mod.list<-append(mod.list,fit)
  # 				AIC.list<-append(AIC.list,AIC(fit))
  # 			}
  # 		}
  # 	}
  # 	# Identify the best model from the list and save coefficients and R^2 values
  # 	bestmodind<-which(AIC.list==min(AIC.list))
  # 	if(length(bestmodind)>1){
  # 		bestmodind<-sample(bestmodind,1)
  # 	}
  # 	
  # 	bestmod<-mod.list[[bestmodind]]
  # 	cfs<-coef(bestmod)
  # 	expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
  # 	rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)
  # }
  
  
  # Use the curve fit to find Topt and the estimated maximum growth rate (i.e. growth rate at Topt)
  grfunc<-function(x){
    -nbcurve(x,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
  }
  optinfo<-optim(c(x=cfs[[1]]),grfunc)
  opt<-optinfo$par[[1]]
  maxgrowth<- -optinfo$value
  
  #stash results		
  rsqr.list[i]<-rsqr
  z.list[i]<-cfs[[1]]
  w.list[i]<-cfs[[2]]
  a.list[i]<-cfs[[3]]
  b.list[i]<-cfs[[4]]
  # s.list[i]<-cfs[[5]]
  topt.list[i]<-opt
  maxgrowth.list[i]<-maxgrowth
  n.list[i]<-length(dat$temperature)
}
fits_limited_variable <-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,n.list)
write_csv(fits_limited_variable, "Tetraselmis_experiment/data-processed/boot_fits_resample_limited_variable.csv")

fits_limited_variable <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_limited_variable.csv")
fits_variable <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v.csv")


fits_limited_variable %>% 
  filter(rsqr.list > 0.98) %>% View


nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

fits_limited_variable %>%
  filter(rsqr.list > 0.98) %>% 
  rename(z = z.list, 
         w = w.list, 
         a = a.list, 
         b = b.list) %>% 
 group_by(curve.id.list) %>%
  mutate(tmax = uniroot.all(function(x) nbcurve(x, z, w, a, b), c(topt.list,150))) %>% 
  mutate(tmin = uniroot.all(function(x) nbcurve(x, z,w, a, b),c(-1.8,10))) %>% View


## now get the upper and lower limits on the curve

fits_split <- fits_variable %>% 
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
boot_limits_variable_limited <- all_predictions %>% 
  group_by(x) %>% 
  summarise(q2.5=quantile(predictions, probs=0.025),
            q97.5=quantile(predictions, probs=0.975),
            mean = mean(predictions)) 
write_csv(boot_limits_variable_limited, "Tetraselmis_experiment/data-processed/boot_limits_variable_limited.csv")

boot_limits_variable <- all_predictions %>% 
  group_by(x) %>% 
  summarise(q2.5=quantile(predictions, probs=0.025),
            q97.5=quantile(predictions, probs=0.975),
            mean = mean(predictions)) 
write_csv(boot_limits_variable, "Tetraselmis_experiment/data-processed/boot_limits_variable.csv")



### let's plot it with the fitted curve
library(colormap)
ic <- colormap(colormap = colormaps$inferno, nshades = 72, format = "hex",
               alpha = 1, reverse = FALSE)

curve_variable_resamp<-function(x){
  res<-ps_v$a[1]*exp(ps_v$b[1]*x)*(1-((x-ps_v$z[1])/(ps_v$w[1]/2))^2)
  res
}



variable_predictions_points <- read_csv("Tetraselmis_experiment/data-processed/variable_predictions_points.csv")

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 
b <- p + 
  geom_ribbon(aes(x = temperature, ymin = growth.rate.lower, ymax = growth.rate.upper, linetype = NA), fill = ic[60], alpha = 0.5, data = variable_predictions_points) +
  geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_variable_limited, fill = ic[10], alpha = 0.5) +
  # 	theme_bw() +
  #   labs(y = expression ("Population growth rate"~day^-1))+
  stat_function(fun = curve_variable_resamp, color = ic[10], size = 1) +
  geom_point(aes(x = temp, y = mean), data = growth_sum_v, color = ic[10], size = 2) +
  geom_errorbar(aes(x = temp, ymin = lower, ymax = upper), data = growth_sum_v, width = 0.1, color = ic[10]) +
  geom_point(aes(x = temp, y = mean), data = growth_sum_v, shape = 1, color = "black", size = 2) +
  xlab("") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_cartesian(ylim = c(-0.5, 1.8), xlim= c(0, 32)) +
  # ylim(-0.5, 1.8) +
  theme_classic() +
  # xlim(0, 32) +
  theme(text = element_text(size=14)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_vline(xintercept = 21.41687, color = ic[10], linetype = "dashed", alpha = 0.7)


