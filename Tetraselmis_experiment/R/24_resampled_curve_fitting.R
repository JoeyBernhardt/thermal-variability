
library(purrr)
library(bbmle)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)


growth_all <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling_exp.csv")
growth_all_v <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling_v_exp.csv") %>% 
  select(-growth_per_day) %>% 
  select(-error)



growth_all_v %>% 
  ggplot(aes(x = temp, y = estimate)) + geom_point()


dat.full <- growth_all_v %>%
	group_by(temp) %>% 
	summarise(growth.rate = mean(estimate)) %>% 
	rename(temperature = temp) %>% 
	mutate(curve.id = "1") %>% 
	select(curve.id, temperature, growth.rate) 
str(dat.full)

# freeze <- data.frame(curve.id = "1", temperature = -1.8, growth.rate = 0)

dat.full <- bind_rows(dat.full, freeze)


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

write_csv(fits, "Tetraselmis_experiment/data-processed/resampling_TPC_params_exp.csv")
write_csv(fits, "Tetraselmis_experiment/data-processed/resampling_TPC_params_v_exp.csv") ## this is for the fitting for the variable dataset


pc <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_exp.csv")
dc <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params.csv")
tc <- read_csv("Tetraselmis_experiment/data-processed/time_resampling_fits.csv")

bind_rows(pc, dc) %>% View

nbcurvec<-function(temp,z,w,a,b){
  res<-pc$a.list[[1]]*exp(pc$b.list[[1]]*temp)*(1-((temp-pc$z.list[[1]])/(pc$w.list[[1]]/2))^2)
  res
}
etpc<-function(temp,z,w,a,b){
  res<-tc$a[[1]]*exp(tc$b[[1]]*temp)*(1-((temp-tc$z[[1]])/(tc$w[[1]]/2))^2)
  res
}
dtpc<-function(temp,z,w,a,b){
  res<-dc$a.list[[1]]*exp(dc$b.list[[1]]*temp)*(1-((temp-dc$z.list[[1]])/(dc$w.list[[1]]/2))^2)
  res
}

nbcurve2<-function(x){
	res<-cfs[3]*exp(cfs[4]*x)*(1-((x-cfs[1])/(cfs[2]/2))^2)
	res
}
growth_all <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling_exp.csv")

growth_sum <- growth_all %>% 
	group_by(temp) %>% 
	summarise(lower = quantile(estimate, probs = 0.025),
						upper = quantile(estimate, probs = 0.975),
						mean = mean(estimate),
						sd = sd(estimate)) 

write_csv(growth_sum, "Tetraselmis_experiment/data-processed/resampled_growth_rates_summary.csv")

growth_sum_v <- growth_all_v %>% 
	group_by(temp) %>% 
	summarise(lower = quantile(estimate, probs = 0.025),
						upper = quantile(estimate, probs = 0.975),
						mean = mean(estimate),
						sd = sd(estimate)) 

write_csv(growth_sum_v, "Tetraselmis_experiment/data-processed/resampled_growth_rates_summary_v.csv")

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p + geom_point(aes(x = temp, y = estimate), data = growth_all, size = 0.05) +
	stat_function(fun = nbcurvec, size = 1) +
  stat_function(fun = etpc, size = 1, color = "purple") +
  stat_function(fun = dtpc, size = 1, color = "blue") +
	geom_point(aes(x = temp, y = mean), data = growth_sum) +
	geom_errorbar(aes(ymin = lower, ymax = upper, x = temp), width = 0.1, data = growth_sum) +
	# geom_errorbar(aes(ymin = mean-sd , ymax = mean+sd, x = temp), width = 0.1, data = growth_sum) +
	theme_classic(base_family = 'GillSans') + xlim(-8, 32) +
	ylab("Exponential growth rate") + xlab("Temperature (C)") +
	theme(text=element_text(size=14)) + geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = -1.8)


fits %>%
  rename(z = z.list,
         a = a.list, 
         b = b.list, 
         w = w.list) %>% 
  group_by(curve.id.list) %>% 
  mutate(tmax = ifelse(length(uniroot.all(function(x) nbcurve(x, z, w, a, b),c(topt.list,150)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w,  a, b),c(topt.list,150)))) %>% 
  mutate(tmin = ifelse(length(uniroot.all(function(x) nbcurve(x, z,w, a, b),c(-5,topt.list)))==0, NA,
                       uniroot.all(function(x) nbcurve(x, z, w, a, b),c(-5,topt.list)))) %>% View

### now let's fit the TPC to the resampled TPC for the variable treatment







