
library(tidyverse)
library(bbmle)

fits_c <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params.csv")
fits_c <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_5000_exp.csv")
fits_c <- read_csv("Tetraselmis_experiment/data-processed/ctpc.csv") %>% 
  rename(a.list = a,
         z.list = z,
         b.list = b,
         w.list = w)

### Goal: find inflection point of the constant temperature TPC

derivative <- function(f, x, ..., order = 1, delta = 0.1, sig = 6) {
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


curve_constant_resamp<-function(x){
	res<-fits_c$a.list[1]*exp(fits_c$b.list[1]*x)*(1-((x-fits_c$z.list[1])/(fits_c$w.list[1]/2))^2)
	res
}


par(mfrow = c(1, 2))
x <- seq(14.716, 14.717, by = 0.001)
for (i in 2:2) {
	plot(x, sapply(x, derivative, f = curve_constant_resamp(), order = i), type = "l", xlab = "x", 
			 ylab = "y", main = paste("Order", i, "Derivative"))
	abline(h = 0, col = "gray60")
}





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
###

x <- seq(16, 18, by = 0.001)


deriv_2 <- function(x) {
	y <- derivative(f = curve_constant_resamp, x = x, order = 2)
}

deriv_value <- sapply(x, deriv_2)
variable_upper2 <- data.frame(x, deriv_value) %>% 
	rename(temperature = x)


variable_upper2 %>% 
	filter(deriv_value < 0.001) %>% View ### inflection point of constant curve is 16.763
	ggplot(aes(x = temperature, y = deriv_value)) + geom_line() +
	geom_hline(yintercept = 0) + xlim(16.5, 17) + ylim(-0.01, 0.01)
	
	
	### find the inflection points of the upper and lower curves
	
	boot_limits_constant <- read_csv("Tetraselmis_experiment/data-processed/boot_limits_constant_resample_10k.csv")
	limits_c <- read_csv("Tetraselmis_experiment/data-processed/limits_c_direct.csv")

	
	dat.full <- limits_c %>% 
	  gather(key = type, value = y, 2:4) %>% 
	  rename(curve.id = type, 
	         growth.rate = y) %>% 
	  mutate(curve.id = ifelse(curve.id == "mean", 1, curve.id),
	         curve.id = ifelse(curve.id == "q2.5", 2, curve.id),
	         curve.id = ifelse(curve.id == "q97.5", 3, curve.id))
	
	
	
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
	
	write_csv(fits, "Tetraselmis_experiment/data-processed/upper_lower_curve_constant_fits_direct.csv")
	write_csv(fits, "Tetraselmis_experiment/data-processed/upper_lower_curve_constant_fits.csv")
	
	curve_lower<-function(x){
	  res<-fits$a.list[curve.id.list == 2]*exp(fits$b.list[curve.id.list == 2]*x)*(1-((x-fits$z.list[curve.id.list == 2])/(fits$w.list[curve.id.list == 2]/2))^2)
	  res
	}
	
	curve_upper<-function(x){
	  res<-fits$a.list[curve.id.list == 3]*exp(fits$b.list[curve.id.list == 3]*x)*(1-((x-fits$z.list[curve.id.list == 3])/(fits$w.list[curve.id.list == 3]/2))^2)
	  res
	}
	
	curve_mean<-function(x){
	  res<-fits$a.list[curve.id.list == 1]*exp(fits$b.list[curve.id.list == 1]*x)*(1-((x-fits$z.list[curve.id.list == 1])/(fits$w.list[curve.id.list == 1]/2))^2)
	  res
	}
	
	deriv_lower <- function(x) {
	  y <- derivative(f = curve_lower, x = x, order = 2)
	}
	
	deriv_lower <- sapply(x, deriv_lower)
	constant_lower <- data.frame(x, deriv_lower) %>% 
	  rename(temperature = x) %>% 
	  filter(deriv_lower < 0.0001) 
	
	## inflection of the lower curve is at 15.669, 16.71 (indirect) and 16.813, 16.764, 16.761 (direct) (update may12 2018)
	
	deriv_upper <- function(x) {
	  y <- derivative(f = curve_upper, x = x, order = 2)
	}
	
	deriv_upper <- sapply(x, deriv_upper)
	constant_upper <- data.frame(x, deriv_upper) %>% 
	  rename(temperature = x) %>% 
	  filter(deriv_upper < 0.0001)
	
	
	deriv_mean <- function(x) {
	  y <- derivative(f = curve_mean, x = x, order = 2)
	}
	
	deriv_mean <- sapply(x, deriv_mean)
	constant_mean <- data.frame(x, deriv_mean) %>% 
	  rename(temperature = x) %>% 
	  filter(deriv_upper < 0.0001)
	
	
	## inflection of the lower curve is at 16.933, 17.124
	
	
	
	
	p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 
TPC <- p +
  # geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 16.96, color = "grey")+
  stat_function(fun = curve_constant_resamp, size = 1.5) +
  ylim(0, 1.8) + xlim(-2, 32) +
  ylab("P(T)") + xlab(expression("")) 

deriv <- p +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 16.96, color = "grey") +
  ylab("P''(T)") + xlab(expression("Temperature (" *degree * "C)")) +
  xlim(-2, 32) + ylim(-0.05, 0.05) +
  geom_line(aes(x = temperature, y = deriv_value), data = variable_upper2, size = 1.5)
  

plot1 <- plot_grid(TPC, deriv, labels = c("A", "B"), align = "v", nrow = 2)
ggsave(plot1, file = "Tetraselmis_experiment/figures/figure1.png", width = 6, height = 7)

