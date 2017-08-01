library(tidyverse)
library(readr)
library(dplyr)
library(purrr)
library(RCurl)
library(tidyr)
library(gridExtra)

all_r <- read_csv("Tetraselmis_experiment/data-processed/growth_estimates_round3.csv")
params <- read_csv("Tetraselmis_experiment/data-processed/constant_TPC_params.csv")

constant_growth <- all_r %>% 
	filter(variability == "c") %>% 
	select(temp, estimate, std.error) %>% 
	mutate(std.error = ifelse(is.na(std.error), 0, std.error)) %>% 
	rename(growth.rate = estimate)

temp <- 32
growth.rate <- 0.001
std.error <- 0.00001
growth.32 <- data.frame(temp, growth.rate, std.error)

all_cons <- bind_rows(constant_growth, growth.32)

all_cons2 <- all_cons %>% 
	mutate(sd = std.error*sqrt(4)) %>% 
	filter(temp != 35) %>% 
	mutate(row = rownames(.)) %>% 
	filter(row != 8)
write_csv(all_cons2, "Tetraselmis_experiment/data-processed/constant_growth_rates.csv")

dat.full_raw <- all_r %>% 
	filter(variability == "c") %>% 
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
	select(-curve.id)


### ok now let's try to make a function out of this!!

nbcurve_pred<-function(x){
	res<-sample(rnorm(1000, mean=params$estimate[[3]], sd=params$std.error[[3]]), size = 1)*exp(sample(rnorm(1000, mean=params$estimate[[4]], sd=params$std.error[[4]]), size = 1)*x)*(1-((x-sample(rnorm(1000, mean=params$estimate[[1]], sd=params$std.error[[1]]), size = 1))/(sample(rnorm(1000, mean=params$estimate[[2]], sd=params$std.error[[2]]), size = 1)/2))^2)
	res
}
##abzw
nbcurve_constant<-function(x){
	res<-params$estimate[[3]]*exp(params$estimate[[4]]*x)*(1-((x-params$estimate[[1]])/(params$estimate[[2]]/2))^2)
	res
}


resampling <- function(a){
	x <- seq(0, 32, 1)
	predictions <- sample(rnorm(1000, mean=params$estimate[[3]], sd=params$std.error[[3]]), size = a)*exp(sample(rnorm(1000, mean=params$estimate[[4]], sd=params$std.error[[4]]), size = a)*x)*(1-((x-sample(rnorm(1000, mean=params$estimate[[1]], sd=params$std.error[[1]]), size = a))/(sample(rnorm(1000, mean=params$estimate[[2]], sd=params$std.error[[2]]), size = a)/2))^2)
	dataframe <- data.frame(x, predictions)
}

as <- rep(1, 50000)
boostrapped <- as %>% 
	map_df(resampling, .id = "run_id") 

write_csv(boostrapped, "Tetraselmis_experiment/data-processed/constant_TPC_bootstrapped.csv")	

bootsummary <- boostrapped %>% 
	group_by(x) %>% 
	summarise(q2.5=quantile(predictions, probs=0.025),
						q97.5=quantile(predictions, probs=0.975),
						mean = mean(predictions))

write_csv(bootsummary, "Tetraselmis_experiment/data-processed/bootsummary.csv")

# plot a ------------------------------------------------------------------

plota <- all_cons %>% 
ggplot(aes(x = temp, y = growth.rate)) + geom_point() +
	geom_errorbar(aes(ymin = growth.rate - std.error*1.96, ymax = growth.rate + std.error*1.96), width = 0.1)+
	geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5), fill = "grey", alpha = 0.3, data = bootsummary) +
	stat_function(fun = nbcurve_constant, size = 1, color = "black") + theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line.y = element_line(color="black"),
				axis.line.x = element_line(color="black")) +
	ylab("exponential growth rate (r)") +
	theme(text = element_text(size=14, family = "Helvetica")) +
	xlab(expression("Temperature (" *degree * "C)")) + ylim(0, 0.1)
ggsave("Tetraselmis_experiment/figures/TPC_constant_werror.png", width = 4, height = 3)


# plot c ------------------------------------------------------------------


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
plotc <- p + geom_point(aes(x = temperature, y = growth.rate), data = constant, size = 0.05, alpha = 0) +
	stat_function(fun = nbcurve_constant, color = "black", size = 2) +
	stat_function(fun = nbcurve_variable, color = "grey", size = 2) +
	ylim(0, 0.1) + theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line.y = element_line(color="black"),
				axis.line.x = element_line(color="black"))+
	ylab("exponential growth rate (r)") +
	theme(text = element_text(size=14, family = "Helvetica")) +
	xlab(expression("Temperature (" *degree * "C)")) 
	


boostrapped %>% 
	# filter(run_id < 1000) %>% 
	ggplot(aes(x = x, y = predictions, group = run_id)) + geom_line(size = 0.5)


### now onto the derivatives

dat.full <- bootsummary %>% 
	gather(key = curve.id, value = growth.rate, 2, 3, 4) %>% 
	rename(temperature = x)
	


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
			res2<-try(fit<-mle2(dat$growth.rate~dnorm(mean=nbcurve(dat$temperature,z=z,w=w,a=a,b=b),sd=s),start=list(z=z.guess,w=w.guess,a=a.guess,b=b.guess,s=0.3),
													skip.hessian=FALSE,data=dat))
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
	
	
	#Save .png plot with fitted curve. File is saved with the curve.id.list value as the name
	# png(paste(curve.id.list[i],'.png',sep=''))
	plot(dat$growth.rate~dat$temperature,ylim=c(pmin(0,min(dat$growth.rate)),max(dat$growth.rate)+(0.2)*max(dat$growth.rate)),main=curve.id.list[1],
			 xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
	plot(dat$growth.rate~dat$temperature,ylim=c(pmin(0,min(dat$growth.rate)),max(dat$growth.rate)+(0.2)*max(dat$growth.rate)),main=curve.id.list[2],
			 xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
	
	plot(dat$growth.rate~dat$temperature,ylim=c(pmin(0,min(dat$growth.rate)),max(dat$growth.rate)+(0.2)*max(dat$growth.rate)),main=curve.id.list[i],
			 xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
	curve(nbcurve(x,cfs[1],cfs[2],cfs[3],cfs[4]),col='red', lwd=2,add=TRUE)
	dev.off()
	
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

boot_fits <-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) 

write_csv(boot_fits, "Tetraselmis_experiment/data-processed/boot_fits.csv")

bootcurve_upper<-function(x){
	res<-fits$a.list[2]*exp(fits$b.list[2]*x)*(1-((x-fits$z.list[2])/(fits$w.list[2]/2))^2)
	res
}

bootcurve_lower<-function(x){
	res<-fits$a.list[1]*exp(fits$b.list[1]*x)*(1-((x-fits$z.list[1])/(fits$w.list[1]/2))^2)
	res
}


bootcurve_mean<-function(x){
	res<-fits$a.list[3]*exp(fits$b.list[3]*x)*(1-((x-fits$z.list[3])/(fits$w.list[3]/2))^2)
	res
}

### function for getting the derivatives
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

## now get the predictions

prediction0_upper <- bootcurve_upper(0) + derivative(f = bootcurve_upper, x = 0, order = 2)*0.5*25
prediction0_mean <- bootcurve_mean(0) + derivative(f = bootcurve_mean, x = 0, order = 2)*0.5*25
prediction0_lower <- bootcurve_lower(0) + derivative(f = bootcurve_lower, x = 0, order = 2)*0.5*25
prediction5_upper <- bootcurve_upper(5) + derivative(f = bootcurve_upper, x = 5, order = 2)*0.5*25
prediction5_mean <- bootcurve_mean(5) + derivative(f = bootcurve_mean, x = 5, order = 2)*0.5*25
prediction5_lower <- bootcurve_lower(5) + derivative(f = bootcurve_lower, x = 5, order = 2)*0.5*25
prediction10_upper <- bootcurve_upper(10) + derivative(f = bootcurve_upper, x = 10, order = 2)*0.5*25
prediction10_mean <- bootcurve_mean(10) + derivative(f = bootcurve_mean, x = 10, order = 2)*0.5*25
prediction10_lower <- bootcurve_lower(10) + derivative(f = bootcurve_lower, x = 10, order = 2)*0.5*25
prediction15_upper <- bootcurve_upper(15) + derivative(f = bootcurve_upper, x = 15, order = 2)*0.5*25
prediction15_mean <- bootcurve_mean(15) + derivative(f = bootcurve_mean, x = 15, order = 2)*0.5*25
prediction15_lower <- bootcurve_lower(15) + derivative(f = bootcurve_lower, x = 15, order = 2)*0.5*25
prediction20_upper <- bootcurve_upper(20) + derivative(f = bootcurve_upper, x = 20, order = 2)*0.5*25
prediction20_mean <- bootcurve_mean(20) + derivative(f = bootcurve_mean, x = 20, order = 2)*0.5*25
prediction20_lower <- bootcurve_lower(20) + derivative(f = bootcurve_lower, x = 20, order = 2)*0.5*25
prediction24_upper <- bootcurve_upper(24) + derivative(f = bootcurve_upper, x = 24, order = 2)*0.5*25
prediction24_mean <- bootcurve_mean(24) + derivative(f = bootcurve_mean, x = 24, order = 2)*0.5*25
prediction24_lower <- bootcurve_lower(24) + derivative(f = bootcurve_lower, x = 24, order = 2)*0.5*25
prediction27_upper <- bootcurve_upper(27) + derivative(f = bootcurve_upper, x = 27, order = 2)*0.5*25
prediction27_mean <- bootcurve_mean(27) + derivative(f = bootcurve_mean, x = 27, order = 2)*0.5*25
prediction27_lower <- bootcurve_lower(27) + derivative(f = bootcurve_lower, x = 27, order = 2)*0.5*25
prediction29_upper <- bootcurve_upper(29) + derivative(f = bootcurve_upper, x = 29, order = 2)*0.5*25
prediction29_mean <- bootcurve_mean(29) + derivative(f = bootcurve_mean, x = 29, order = 2)*0.5*25
prediction29_lower <- bootcurve_lower(29) + derivative(f = bootcurve_lower, x = 29, order = 2)*0.5*25


temperature <- c(0, 0, 0, 5, 5, 5, 10, 10, 10, 15, 15, 15, 20, 20, 20, 24, 24, 24, 27, 27, 27, 29, 29, 29)
type <- c("upper","mean", "lower", "upper","mean", "lower", "upper","mean", "lower", "upper", "mean", "lower", "upper", "mean", "lower", "upper","mean", "lower", "upper", "mean", "lower", "upper","mean", "lower")
prediction <- c(prediction0_upper, prediction0_mean, prediction0_lower, prediction5_upper, prediction5_mean, prediction5_lower, prediction10_upper, prediction10_mean, prediction10_lower,
								prediction15_upper, prediction15_mean, prediction15_lower, prediction20_upper, prediction20_mean, prediction20_lower,
								prediction24_upper, prediction24_mean, prediction24_lower, 
								prediction27_upper, prediction27_mean, 	prediction27_lower,
								prediction30_upper, prediction30_mean, 	prediction30_lower)
predictions <- data.frame(type, temperature, prediction)

predictions_wide <- predictions %>% 
	spread(key = type, value = prediction)

r_estimates <- all_r %>% 
	filter(variability == "v") %>%
	select(1:6) %>% 
	rename(temperature = temp) 

obs_pred <- left_join(predictions_wide, r_estimates)

obs_pred %>% 
	ggplot(aes(x = temperature, y = estimate)) + geom_point(size = 2) +
	geom_errorbar(aes(ymin = estimate - std.error*1.96, ymax = estimate + std.error*1.96), width = 0.1) +
	# geom_point(aes(x = temperature, y = prediction, color = type), data = predictions, size = 2) +
	geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.3) + theme_bw() +
	stat_function(fun = nbcurve_variable, color = "black") +
	stat_function(fun = nbcurve_constant, color = "black", linetype = "dotted") +
	xlim(0, 32) + ylim(0, 0.078) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("exponential growth rate (r)") +
	theme(text = element_text(size=14, family = "Helvetica")) +
	xlab(expression("Temperature (" *degree * "C)")) +
	geom_line(aes(x = temperature, y = growth.rate), data = predicted_growth_variable2, color = "grey38")



# now can we make a line that represents predictions? ---------------------



x <- seq(0, 32, by = 0.001)

variable_predictions <- function(x) {
	y <- nbcurve_constant(x) + derivative(f = nbcurve_constant, x = x, order = 2)*0.5*25
}

predicted_growth_variable <- sapply(x, variable_predictions)
predicted_growth_variable2 <- data.frame(x, predicted_growth_variable) %>% 
	rename(temperature = x, 
				 growth.rate = predicted_growth_variable)

variable_predictions_upper <- function(x) {
	y <- bootcurve_upper(x) + derivative(f = bootcurve_upper, x = x, order = 2)*0.5*25
}

variable_upper <- sapply(x, variable_predictions_upper)
variable_upper2 <- data.frame(x, variable_upper) %>% 
	rename(temperature = x, 
				 growth.rate.upper = variable_upper)

variable_predictions_lower <- function(x) {
	y <- bootcurve_lower(x) + derivative(f = bootcurve_lower, x = x, order = 2)*0.5*25
}

variable_lower <- sapply(x, variable_predictions_lower)
variable_lower2 <- data.frame(x, variable_lower) %>% 
	rename(temperature = x, 
				 growth.rate.lower = variable_lower)

variable_predictions_points <- left_join(variable_lower2, variable_upper2) %>% 
	filter(growth.rate.lower >= 0, growth.rate.upper >=0)

# plot b ------------------------------------------------------------------

plotb <- obs_pred %>% 
	ggplot(aes(x = temperature, y = estimate)) + geom_point(size = 0.5) +
	geom_ribbon(aes(x = temperature, ymin = growth.rate.lower, ymax = growth.rate.upper),
							inherit.aes = FALSE, data = variable_predictions_points, fill = "grey", alpha = 0.5) +
	geom_errorbar(aes(ymin = estimate - std.error*1.96, ymax = estimate + std.error*1.96), width = 0.1) +
	ylab("exponential growth rate (r)") +
	theme(text = element_text(size=14, family = "Helvetica")) +
	xlab(expression("Temperature (" *degree * "C)")) +
	theme_bw() +
	geom_line(aes(x = temperature, y = growth.rate), data = predicted_growth_variable2, color = "grey", size = 1.2) +
	# stat_function(fun = nbcurve_variable, color = "black", size = 1.2) +
	geom_point(aes(x = temperature, y = estimate), data = obs_pred, size = 1.5) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line.y = element_line(color="black"),
				axis.line.x = element_line(color="black")) +
	xlim(0, 32) + ylim(0, 0.1)
ggsave("Tetraselmis_experiment/figures/variable_growth_boot.pdf")
ggsave("Tetraselmis_experiment/figures/variable_growth_boot.png", width = 4, height =3)	

plots <- grid.arrange(plota, plotb, plotc, nrow = 1, ncol =3)
ggsave(plots, file = "Tetraselmis_experiment/figures/figure2.png", width = 8, height = 3)
