## plotting the confidence intervals, are the curves different??



# load packages -----------------------------------------------------------

library(tidyverse)
library(bbmle)
library(cowplot)
library(colormap)


# load data ---------------------------------------------------------------
all_r <- read_csv("Tetraselmis_experiment/data-processed/growth_estimates_round3.csv")
growth_constant <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling.csv")
# boot_limits_variable <- read_csv("Tetraselmis_experiment/data-processed/boot_limits_variable.csv")
# boot_limits_constant <- read_csv("Tetraselmis_experiment/data-processed/boot_limits_constant.csv")
boot_limits_constant <- read_csv("Tetraselmis_experiment/data-processed/boot_limits_constant_resample.csv") ### these are for the new resampled curve
boot_limits_variable <- read_csv("Tetraselmis_experiment/data-processed/boot_limits_constant_resample_v.csv")
fits_constant_variable <- read_csv("Tetraselmis_experiment/data-processed/fits_constant_variable.csv")

variable_growth <- all_r %>% 
	filter(variability == "v")

constant_growth <- all_r %>% 
	filter(variability == "c") %>% 
	filter(temp< 33)



variable2 <- boot_limits_variable %>% 
	mutate(q2.5 = ifelse(q2.5 < -0.01, -0.01, q2.5))

## plot!

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5), fill = "grey", alpha = 0.7, data = variable2) +
	geom_line(aes(x = x, y = mean), data = boot_limits_variable) +
	geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5), fill = "grey", alpha = 0.7, data = boot_limits_constant) +
	geom_line(aes(x = x, y = mean), data = boot_limits_constant) +
	theme_bw() + ylab("growth rate/hr") + ylim(-0.01, 0.075) +
	xlab("Temperature (°C)") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=14, family = "Helvetica")) 
ggsave("Tetraselmis_experiment/figures/TPC_w_CIs.pdf")


## ok here are the upper and lower confidence lines around the constant curve.
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p +geom_line(aes(x = x, y = q2.5), data = boot_limits_constant) +
	geom_line(aes(x = x, y = q97.5), data = boot_limits_constant) +
	theme_bw() + ylab("growth rate/hr") + ylim(-0.01, 0.075) +
	xlab("Temperature (°C)") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=14, family = "Helvetica")) 

## next step is to fit the TPC to them, so that we can get their second derivatives.
### here we fit the TPC to the upper and lower edges of the bootstrapped TPCs
dat.full <- boot_limits_constant %>% 
	gather(key = curve.id, value = growth.rate, 2:4) %>% 
	rename(temperature = x) 

dat.full %>% 
	ggplot(aes(x = temperature, y = growth.rate, group = curve.id)) + geom_line()


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

fits <-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) 
write_csv(fits, "Tetraselmis_experiment/data-processed/boot_upper_lower_fits.csv")
fits <- read_csv("Tetraselmis_experiment/data-processed/boot_upper_lower_fits.csv")

## now make the plots for predictions

# 
# fits <- read_csv("Tetraselmis_experiment/data-processed/fits_constant_variable.csv")
# fits <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params.csv") ## new fitted params from resampling approach

bootcurve_upper<-function(x){
	res<-fits$a.list[2]*exp(fits$b.list[2]*x)*(1-((x-fits$z.list[2])/(fits$w.list[2]/2))^2)
	res
}

bootcurve_lower<-function(x){
	res<-fits$a.list[1]*exp(fits$b.list[1]*x)*(1-((x-fits$z.list[1])/(fits$w.list[1]/2))^2)
	res
}




###
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
###

x <- seq(0, 33, by = 0.001)


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

write_csv(variable_predictions_points, "Tetraselmis_experiment/data-processed/variable_predictions_points.csv")
variable_predictions_points <- read_csv("Tetraselmis_experiment/data-processed/variable_predictions_points.csv")

curve_variable<-function(x){
	res<-fits_constant_variable$a.list[2]*exp(fits_constant_variable$b.list[2]*x)*(1-((x-fits_constant_variable$z.list[2])/(fits_constant_variable$w.list[2]/2))^2)
	res
}

curve_constant<-function(x){
	res<-fits_constant_variable$a.list[1]*exp(fits_constant_variable$b.list[1]*x)*(1-((x-fits_constant_variable$z.list[1])/(fits_constant_variable$w.list[1]/2))^2)
	res
}


variable_predictions_points <- read_csv("Tetraselmis_experiment/data-processed/variable_predictions_points.csv")


growth_sum_v <- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary_v.csv")
growth_sum<- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary.csv")

fits_v <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_v.csv")
fits_c <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params.csv")

curve_variable_resamp<-function(x){
	res<-fits_v$a.list[1]*exp(fits_v$b.list[1]*x)*(1-((x-fits_v$z.list[1])/(fits_v$w.list[1]/2))^2)
	res
}

curve_constant_resamp<-function(x){
	res<-fits_c$a.list[1]*exp(fits_c$b.list[1]*x)*(1-((x-fits_c$z.list[1])/(fits_c$w.list[1]/2))^2)
	res
}
## plot the predictions

ic <- colormap(colormap = colormaps$inferno, nshades = 72, format = "hex",
							 alpha = 1, reverse = FALSE)


## plot A
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 
a <- p + geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_constant, fill = ic[20], alpha = 0.5) + theme_bw()+
	stat_function(fun = curve_constant_resamp, color = ic[20], size = 1) +
	geom_point(aes(x = temp, y = mean), data = growth_sum, color = ic[20], size = 2) + geom_errorbar(aes(x = temp, ymin = lower, ymax = upper), data = growth_sum, width = 0.1, color =ic[20]) +
	geom_point(aes(x = temp, y = mean), data = growth_sum, shape = 1, size = 2) +
	xlab("") + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black"),
				panel.border = element_rect(colour = "black", fill=NA, size=1)) +
	ylim(-0.5, 1.8) +
	theme_classic() + xlim(0, 32) +
	theme(text = element_text(size=14)) +
	geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
	labs(y = expression ("Population growth rate"~day^-1)) +
	geom_vline(xintercept = 24.5681, color = ic[20], linetype = "dashed", alpha = 0.7)
ggsave("Tetraselmis_experiment/figures/constant_tpc_figure_with_inflection.png", width = 5, height = 3)
ggsave("Tetraselmis_experiment/figures/constant_tpc.png", width = 5, height = 3)




## plot B
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 
b <- p + geom_ribbon(aes(x = temperature, ymin = growth.rate.lower, ymax = growth.rate.upper, linetype = NA), fill = ic[60], alpha = 0.5, data = variable_predictions_points) +
# 	geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_variable, fill = ic[10], alpha = 0.5) +
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
	ylim(-0.5, 1.8) +
	theme_classic() + xlim(0, 32) +
	theme(text = element_text(size=14)) +
	geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
	geom_vline(xintercept = 21.41687, color = ic[10], linetype = "dashed", alpha = 0.7)
ggsave("Tetraselmis_experiment/figures/variable_predictions_data.pdf")
ggsave("Tetraselmis_experiment/figures/variable_predictions_data_prediction_band.png", width = 5, height = 3)
ggsave("Tetraselmis_experiment/figures/variable_predictions_data_prediction_band_points.png", width = 5, height = 3)

## plot C




p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 
c <- p + 
	theme_bw() +
	geom_vline(xintercept = 24.5681, color = ic[20], linetype = "dashed", alpha = 0.5) +
	geom_vline(xintercept = 21.41687, color = ic[10], linetype = "dashed", alpha = 0.5)+
	stat_function(fun = curve_constant_resamp, color = ic[20], size = 1) +
	stat_function(fun = curve_variable_resamp, color = ic[10], size = 1)  + xlim(0, 33) + 
	xlab("Temperature (°C)") + ylab("") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	 			panel.background = element_blank(),
	 			axis.line = element_line(color="black"),
				panel.border = element_rect(colour = "black", fill=NA, size=1)) +
	ylim(-0.5, 1.8) +
	theme_classic() + xlim(0, 32) +
	theme(text = element_text(size=14)) +
	labs(y = expression ("Population growth rate"~day^-1))+
	geom_hline(yintercept = 0, color= "grey", linetype = "dotted") +
	geom_vline(xintercept = 16.96, color = "grey", linetype = "dotted") +
	geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_variable, fill = ic[10], alpha = 0.5) +
	geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_constant, fill = ic[20], alpha = 0.5)
ggsave("Tetraselmis_experiment/figures/variable_predictions_data_prediction_together.png", width = 5, height = 3)

	
	

plots <- plot_grid(a, b, c, labels = c("A", "B", "C"), align = "v", nrow = 3)
ggsave(plots, file = "Tetraselmis_experiment/figures/figure2_resampling_color_tall.png", width = 5, height = 10.15)
ggsave(plots, file = "Tetraselmis_experiment/figures/figure2_resampling_color_tall.pdf", width = 5, height = 10.15)


library(RColorBrewer)
View(palette(viridis(6)))
display.brewer.pal(n = 8, name = 'Dark2')


### Now let's try with all the curves and points overlaid


p + 
	theme_bw() +
	stat_function(fun = curve_constant_resamp, color = "#43a2ca", size = 1) +
	stat_function(fun = curve_variable_resamp, color = "#2ca25f", size = 1)  + xlim(0, 33) + ylim(-0.5, 1.7) +
	xlab("Temperature (°C)") + ylab("") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black"),
				panel.border = element_rect(colour = "black", fill=NA, size=1)) +
	ylim(-0.5, 1.8) +
	theme_classic(base_family = 'GillSans') + xlim(0, 32) +
	theme(text = element_text(size=14)) +
	geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
	geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_variable, fill = "#2ca25f", alpha = 0.5) +
	geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_constant, fill = "#43a2ca", alpha = 0.5) +
	geom_ribbon(aes(x = temperature, ymin = growth.rate.lower, ymax = growth.rate.upper, linetype = NA), fill = "orange", alpha = 0.5, data = variable_predictions_points) +
	geom_point(aes(x = temp, y = mean), data = growth_sum_v, color = "#2ca25f", size = 1) +
	geom_errorbar(aes(x = temp, ymin = lower, ymax = upper), data = growth_sum_v, width = 0.1, color = "#2ca25f") +
	geom_point(aes(x = temp, y = mean), data = growth_sum_v, shape = 1, color = "black", size = 2) +
	geom_point(aes(x = temp, y = mean), data = growth_sum, color = "#43a2ca", size = 2) + geom_errorbar(aes(x = temp, ymin = lower, ymax = upper), data = growth_sum, width = 0.1, color ="#43a2ca") +
	geom_point(aes(x = temp, y = mean), data = growth_sum, shape = 1, size = 2) +
	geom_vline(xintercept = 16.96, color = "grey", linetype = "dotted")

library(modelr)


temperature <- seq(-2, 33, 0.1)


curve_variable_resamp_points <- function(x){
	res<-fits_v$a.list[1]*exp(fits_v$b.list[1]*x)*(1-((x-fits_v$z.list[1])/(fits_v$w.list[1]/2))^2)
	res
}

curve_constant_resamp<-function(x){
	res<-fits_c$a.list[1]*exp(fits_c$b.list[1]*x)*(1-((x-fits_c$z.list[1])/(fits_c$w.list[1]/2))^2)
	res
}

x <- seq(-5, 33, 0.01)
variable_prediction <- function(x) {
	predictions <- fits_v$a.list[1]*exp(fits_v$b.list[1]*x)*(1-((x-fits_v$z.list[1])/(fits_v$w.list[1]/2))^2)
	data.frame(x, predictions)
}

constant_prediction <- function(x) {
	predictions <- fits_c$a.list[1]*exp(fits_c$b.list[1]*x)*(1-((x-fits_c$z.list[1])/(fits_c$w.list[1]/2))^2)
	data.frame(x, predictions)
}

variation_points <- x %>% 
	map_df(variable_prediction)

constant_points <- x %>% 
	map_df(constant_prediction)

constant_points %>% 
	filter(x < 10) %>% View

variation_points %>% 
	filter(x < 10) %>% View


### let's find the 95% CI on the Topt and rmax


cboot <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample.csv")
vboot <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v.csv")


cboot %>% 
	filter(a.list > 0) %>% 
	summarise(q2.5=quantile(maxgrowth.list, probs=0.025),
						q97.5=quantile(maxgrowth.list, probs=0.975),
						mean = mean(maxgrowth.list)) %>% View

cboot %>% 
	filter(a.list > 0, z.list > 0) %>% 
	summarise(q2.5=quantile(w.list, probs=0.025),
						q97.5=quantile(w.list, probs=0.975),
						mean = mean(w.list)) %>% View

# 33.59421
# 41.85766
# 36.77638


vboot %>%
	filter(a.list > 0, z.list > 0) %>% 
	summarise(q2.5=quantile(w.list, probs=0.025),
						q97.5=quantile(w.list, probs=0.975),
						mean = mean(w.list)) %>% View

# 29.37522
# 56.5662
# 37.0481

vboot %>% 
	filter(a.list > 0) %>% 
	summarise(q2.5=quantile(topt.list, probs=0.025),
						q97.5=quantile(topt.list, probs=0.975),
						mean = mean(topt.list)) %>% View

cboot %>% 
	filter(a.list > 0) %>% 
	summarise(q2.5=quantile(topt.list, probs=0.025),
						q97.5=quantile(topt.list, probs=0.975),
						mean = mean(topt.list)) %>% View

curve229 <- cboot %>% 
	filter(curve.id.list == 458)

fits_c <- curve229

curve_constant_1<-function(x){
	res<-fits_c$a.list[1]*exp(fits_c$b.list[1]*x)*(1-((x-fits_c$z.list[1])/(fits_c$w.list[1]/2))^2)
	res
}
 
p + 
	theme_bw() +
	stat_function(fun = curve_constant_1, color = "#43a2ca", size = 1) + geom_vline(xintercept = fits_c$topt.list[1]) +
	geom_hline(yintercept = fits_c$maxgrowth.list[1]) + geom_hline(yintercept = 0) + xlim(-82, -30)  + ylim(-0.01, 0.01)

## tmax = 32.0

p + 
	theme_bw() +
	stat_function(fun = curve_constant_resamp, color = "#43a2ca", size = 1) + xlim(31.9, 31.925) + geom_vline(xintercept = fits_c$topt.list[1]) +
	geom_hline(yintercept = fits_c$maxgrowth.list[1]) + geom_hline(yintercept = 0) +ylim(-0.01, 0.01)

##-4.5, 31.90
	