
### Bootstrap the variable curve now!

# load packages -----------------------------------------------------------


library(purrr)
library(bbmle)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)


# load data ---------------------------------------------------------------


params <- read_csv("Tetraselmis_experiment/data-processed/variable_TPC_params.csv") ## estimated TPC parameters for constant conditions
growth <- read_csv("Tetraselmis_experiment/data-processed/growth_estimates_round3.csv") ## empirically observed growth rates



growth2 <- growth %>% 
	filter(variability == "v") %>% 
	# rename(growth.rate = estimate) %>% 
	mutate(sd = std.error*sqrt(4)) %>% 
	select(temp, estimate, sd)
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

samples <- rep(1, 2)

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
	
	# Identify the best model from the list and save coefficients and R^2 values
	if(!is.null(AIC.list)){
		bestmodind<-which(AIC.list==min(AIC.list))
		if(length(bestmodind)>1){
			bestmodind<-sample(bestmodind,1)
		}
		bestmod<-mod.list[[bestmodind]]
		cfs<-coef(bestmod)
		expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
		rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)
	}
	
	# If the quick fit yielded poor results (low R^2), try a more thorough search through parameter space
	if(rsqr<0.95){
		avals<-seq(-0.2,1.2,0.02)
		bvals<-seq(-0.2,0.3,0.02)
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
		# Identify the best model from the list and save coefficients and R^2 values
		bestmodind<-which(AIC.list==min(AIC.list))
		if(length(bestmodind)>1){
			bestmodind<-sample(bestmodind,1)
		}
		
		bestmod<-mod.list[[bestmodind]]
		cfs<-coef(bestmod)
		expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
		rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)
	}
	
	
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
	s.list[i]<-cfs[[5]]
	topt.list[i]<-opt
	maxgrowth.list[i]<-maxgrowth
	n.list[i]<-length(dat$temperature)
}

fits<-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list)

write_csv(fits, "Tetraselmis_experiment/data-processed/boot_fits_variable.csv")

fits <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_variable.csv")

## ok now take the fits, and make prediction curves and then take the 97.5 and 2.5% CIs

## split up the fits df by curve id
fits_split <- fits %>% 
	split(.$curve.id.list)

prediction_function <- function(curve1) {
	x <- seq(0, 32, 1)
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
write_csv(boot_limits, "Tetraselmis_experiment/data-processed/boot_limits_variable.csv")	

## plot it! (bands look super skinny now??). I think this is what we are after.
boot_limits <- read_csv("Tetraselmis_experiment/data-processed/boot_limits_variable.csv")
	
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5), fill = "grey", alpha = 0.7, data = boot_limits) +
	geom_line(aes(x = x, y = mean), data = boot_limits) +
	theme_bw() + ylab("growth rate") + xlab("Temperature")


### now try with the upper and lower limits of the paramaters themselves (I don't think this is the right thing to do here)

boot_params <- fits %>% 
	select(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list) %>% 
	gather(key = parameter, value = estimate, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list) %>%
	group_by(parameter) %>% 
	summarise(q2.5=quantile(estimate, probs=0.025),
						q97.5=quantile(estimate, probs=0.975),
						mean = mean(estimate))

bootcurve_upper<-function(x){
	res<-boot_params$q97.5[boot_params$parameter == "a.list"]*exp(boot_params$q97.5[boot_params$parameter == "b.list"]*x)*(1-((x-boot_params$q97.5[boot_params$parameter == "z.list"])/(boot_params$q97.5[boot_params$parameter == "w.list"]/2))^2)
	res
}

bootcurve_lower<-function(x){
	res<-boot_params$q2.5[boot_params$parameter == "a.list"]*exp(boot_params$q2.5[boot_params$parameter == "b.list"]*x)*(1-((x-boot_params$q2.5[boot_params$parameter == "z.list"])/(boot_params$q2.5[boot_params$parameter == "w.list"]/2))^2)
	
	res
}

bootcurve_mean<-function(x){
	res<-boot_params$mean[boot_params$parameter == "a.list"]*exp(boot_params$mean[boot_params$parameter == "b.list"]*x)*(1-((x-boot_params$mean[boot_params$parameter == "z.list"])/(boot_params$mean[boot_params$parameter == "w.list"]/2))^2)
	
	res
}

## now bands look super wide? 
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p + geom_point(aes(x = temperature, y = growth.rate), data = data_full, size = 0.5, alpha = 0) +
	stat_function(fun = bootcurve_lower, color = "red", size = 1) + 
	stat_function(fun = bootcurve_upper, color = "blue", size = 1) + 
	stat_function(fun = bootcurve_mean, color = "black", size = 1) + ylim(0, 0.1)

