 ### try fitting a gaussian x gompertz function to our data



# ggcurve<-function(temp,z,w,a,b){
# 	res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
# 	res
# }
# 
# ggcurve<-function(rmax, B, Tb, Trmax, alpha){
# 	res<- rmax*(exp(-exp(B(Tb - Trmax)-8) - alpha(Tb - Trmax)^2))
# 	res
# }

## w is rmax, z is B, a is Trmax, b is alpha

nbcurve<-function(temp,z,w,a,b){
	res<-a*(exp(-exp(b(temp - z)-8) - w(temp - z)^2))
	res
}


library(tidyverse)
library(stringr)
library(rootSolve)
library(bbmle)
all_r <- read_csv("Tetraselmis_experiment/data-processed/growth_estimates_round3.csv")

all_r %>% 
	filter(temp < 31) %>% 
	ggplot(aes(x = temp, y = estimate, color = variability)) + geom_point() + geom_line()

temperature <- 32
growth.rate <- 0.001
curve.id <- "c"
growth_32 <- data.frame(temperature, growth.rate, curve.id)

data_full <- all_r %>% 
	# filter(var == "variable temperature") %>% 
	select(-temperature) %>% 
	rename(temperature = temp) %>% 
	rename(growth.rate = estimate) %>% 
	# mutate(var = str_replace(var, "constant temperature" ,"1")) %>% 
	# mutate(var = str_replace(var, "variable temperature" ,"2")) %>% 
	rename(`thermal environment` = variability) %>% 
	select(`thermal environment`, temperature, growth.rate) %>% 
	filter(temperature < 33)

dat.full_raw <- all_r %>% 
	# filter(variability == "c") %>% 
	select(-temperature) %>% 
	rename(temperature = temp) %>% 
	rename(growth.rate = estimate) %>% 
	mutate(var = str_replace(variability, "c" ,"1")) %>% 
	mutate(var = str_replace(variability, "v" ,"2")) %>% 
	rename(curve.id = variability) %>% 
	select(curve.id, temperature, growth.rate) %>% 
	filter(temperature < 33)

dat.full <- bind_rows(dat.full_raw, growth_32) %>% 
	mutate(curve.id = as.factor(curve.id)) %>% 
	filter(curve.id == "c")

str(dat.full)


#### from Mridul's code, get the best fits for both of the TPCs

nbcurve<-function(temp,z,w,a,b){
	res<-a*(exp(-exp((b*(temp - z)-8)) - (w*(temp - z)^2)))
	res
}


# Create new vector of unique curve.id values
curve.id.list<-unique(dat.full$curve.id)	

# Create empty vectors to populate with parameter values and trait estimates
z.list<-rep(NA, length(curve.id.list))				#Parameter 'z'
w.list<-rep(NA, length(curve.id.list))				#Parameter 'w', which is the temperature niche width
a.list<-rep(NA, length(curve.id.list))				#Parameter 'a'
b.list<-rep(NA, length(curve.id.list))				#Parameter 'b'
# topt.list<-rep(NA, length(curve.id.list))			#Topt, Optimum temperature for growth
maxgrowth.list<-rep(NA, length(curve.id.list))		#Maximum growth rate (i.e. growth rate at Topt)
rsqr.list<-rep(NA, length(curve.id.list))			#R^2 values for all fits
s.list<-rep(NA, length(curve.id.list))				#Error
n.list<-rep(NA, length(curve.id.list))				#Number of growth rate measurements used in the curve

# Loop through all curve.id.list values to estimate parameters for all curves

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

			
			res2<-try(fit<-mle2(dat$growth.rate~dnorm(mean=nbcurve(dat$temperature,z=z,w=w,a=a,b=b),sd=s),start=list(z=27,w=1,a=0.075,b=1,s=0.3), 
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
	
	

	
	
	#stash results		
	rsqr.list[i]<-rsqr
	z.list[i]<-cfs[[1]]
	w.list[i]<-cfs[[2]]
	a.list[i]<-cfs[[3]]
	b.list[i]<-cfs[[4]]
	# s.list[i]<-cfs[[5]]
	# topt.list[i]<-opt
	# maxgrowth.list[i]<-maxgrowth
	n.list[i]<-length(dat$temperature)
}

fits_constant_variable2 <-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) 
write_csv(fits_constant_variable, "Tetraselmis_experiment/data-processed/gaussian_gompertz_fit.csv")

fits_constant_variable <- read_csv("Tetraselmis_experiment/data-processed/gaussian_gompertz_fit.csv")
res2 <- fits_constant_variable
nbcurve2<-function(temp,z,w,a,b){
	res<-res2$a.list*(exp(-exp((res2$b.list*(temp - res2$z.list)-8)) - (res2$w.list*(temp - res2$z.list)^2)))
	res
}


dat %>%
	ggplot(aes(x = temperature, y = growth.rate)) + geom_point() +
	stat_function(fun = nbcurve2) + theme_bw() + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) + 
	ylab("Growth rate") + xlab("Temperature (°C)")
ggsave("Tetraselmis_experiment/figures/gaussian_gompertz_fit.png")


cfs<-coef(res2)
expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
(rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2))




# now let’s run the GG curve through to predictions -----------------------

derivative <- function(f, x, ..., order = i, delta = 0.1, sig = 6) {
	# Numerically computes the specified order derivative of f at x
	vals <- matrix(NA, nrow = order + 1, ncol = order + 1)
	grid <- seq(x - delta/2, x + delta/2, length.out = order + 1)
	vals[1, ] <- sapply(grid, f, ...) - f(x, ...)
	for (i in 2:(order + 1)) {
		for (j in 1:(order - i + 2)) {
			stepsize <- grid[i + j - 1] - grid[i + j - 2]
			vals[i, j] <- (vals[i - 1, j + 1] - vals[i - 1, j])/stepsize
		}
	}
	return(signif(vals[order + 1, 1], sig))
}

nbcurve2<-function(temp,z,w,a,b){
	res<-res2$a.list*(exp(-exp((res2$b.list*(temp - res2$z.list)-8)) - (res2$w.list*(temp - res2$z.list)^2)))
	res
}


x <- seq(0, 35, by = 0.1)

second_derivative <- sapply(x, derivative, f = nbcurve2, order = 2)

sder <- data.frame(x, second_derivative)


x <- seq(0, 32, by = 0.01)

variable_predictions <- function(x) {
	y <- nbcurve2(x) + derivative(f = nbcurve2, x = x, order = 2)*0.5*25
}

predicted_growth_variable <- sapply(x, variable_predictions)
predicted_growth_variable2 <- data.frame(x, predicted_growth_variable) %>% 
	rename(temperature = x, 
				 growth.rate = predicted_growth_variable)


predicted_growth_variable2 %>% 
	ggplot(aes(x = temperature, y = growth.rate)) + geom_line(color = "grey", size = 1) +xlim(0, 30) + ylim(0, 0.075) +
	stat_function(fun = nbcurve2, color = "black", size = 1) + theme_bw() +
	geom_point(aes(x = temperature, y = growth.rate, group = `thermal environment`, color = `thermal environment`), data = data_full, size = 1) +
	scale_color_manual(values = c("black", "grey")) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) + 
	ylab("Growth rate") + xlab("Temperature (°C)")
ggsave("Tetraselmis_experiment/figures/GG_fits_w_data.png", width = 8, height = 4)
