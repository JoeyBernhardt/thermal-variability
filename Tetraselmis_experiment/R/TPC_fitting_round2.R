


# TPC fitting with round 2 data -------------------------------------------
library(tidyverse)
library(stringr)
library(rootSolve)
library(bbmle)
library(cowplot)
library(broom)

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
	mutate(curve.id = as.factor(curve.id))

str(dat.full)


#### from Mridul's code, get the best fits for both of the TPCs
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

fits_constant_variable <-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) 
summary(bestmod)
constant_TPC_params <- tidy(bestmod)
write_csv(constant_TPC_params, "Tetraselmis_experiment/data-processed/constant_TPC_params.csv")
write_csv(fits_constant_variable, "Tetraselmis_experiment/data-processed/fits_constant_variable.csv")

fits$tmax<-fits$z.list+(fits$w.list/2)
fits$tmin<-fits$z.list-(fits$w.list/2)

### alternate way of getting CTmin and CTmax
for (i in 1:length(curve.id.list)){
	 print(i)
	nbcurve.tmax<-function(x){
	nb<-nbcurve(x,fits$z.list[i],fits$w.list[i],fits$a.list[i],fits$b.list[i])
	nb
	}
	fits$tmax[i]<-uniroot.all(nbcurve.tmax,c(fits$topt.list[i],150))[1]
	fits$tmin[i]<-uniroot.all(nbcurve.tmax,c(-2,fits$topt.list[i]))[1] 
	 }


### now make the plots!
nbcurve_constant<-function(x){
	res<-fits_constant_variable$a.list[1]*exp(fits_constant_variable$b.list[1]*x)*(1-((x-fits_constant_variable$z.list[1])/(fits_constant_variable$w.list[1]/2))^2)
	res
}



nbcurve_variable<-function(x){
	res<-fits_constant_variable$a.list[2]*exp(fits_constant_variable$b.list[2]*x)*(1-((x-fits_constant_variable$z.list[2])/(fits_constant_variable$w.list[2]/2))^2)
	res
}

data_full2 <- dat.full %>% 
	mutate(curve.id = ifelse(curve.id == "c", "constant", "variable")) %>% 
	rename(`thermal environment` = curve.id) 
	

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p + geom_point(aes(x = temperature, y = growth.rate, color = `thermal environment`), data = data_full2, size = 0.5, alpha = 1) +
	stat_function(fun = nbcurve_constant, color = "black", size = 2) +
	stat_function(fun = nbcurve_variable, color = "grey", size = 2) +
 geom_point(aes(x = temperature, y = growth.rate, color = `thermal environment`), data = data_full2, size = 2, alpha = 0.9) + scale_color_manual(values = c("black", "grey"))+
	# geom_rect(aes(xmin=0, xmax=16.59, ymin=-Inf, ymax=Inf), alpha = 0.1) +
	ylim(0, 0.07) + theme_bw() +
	geom_vline(xintercept = 16.59, linetype = "dashed", color = "grey") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Intrinsic growth rate (r)") +
	theme(text = element_text(size=14, family = "Helvetica")) +
	xlab(expression("Temperature (" *degree * "C)")) + theme(legend.position = c(0.2, 0.85))
ggsave("Tetraselmis_experiment/figures/figure2_TPC_with_data.pdf")

# get the derivatives of the constant TPC ---------------------------------

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

derivative2 <- function(x, f = nbcurve_constant,  order = 2, delta = 0.1, sig = 6) {
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


x <- seq(0, 35, by = 0.001)

second_derivative <- sapply(x, derivative, f = nbcurve_constant, order = 2)

sder <- data.frame(x, second_derivative)


x <- seq(0, 35, by = 0.001)

first_derivative <- sapply(x, derivative, f = nbcurve_constant, order = 1)

fder <- data.frame(x, first_derivative)


par(mfrow = c(1, 2))
grid <- seq(24.549, 24.55, by = 0.001)
for (i in 2:2) {
	plot(grid, sapply(grid, derivative, f = nbcurve_constant, order = 1), type = "l", xlab = "x", 
			 ylab = "y", main = paste("Order", i, "Derivative"))
	abline(h = 0, col = "gray60")
}

## 16.591 is where the second deriv crosses the 0 line
## 24.55 is where the first deriv crosses the 0 line
### let's plot our predictions

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

TPC_plot <- p + 
	# geom_point(aes(x = temperature, y = growth.rate, color = curve.id), data = dat.full, size = 0.5, alpha = 1) +
	stat_function(fun = nbcurve_constant, color = "black", size = 1.3) + xlim(0, 32) + theme_bw() + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank()) + 
	theme(axis.title.y=element_text(face="italic"))+
	ylab("P") + xlab(expression("Temperature (" *degree * "C)")) +
	theme(text = element_text(size=14, family = "Helvetica")) + geom_vline(xintercept = 16.59, color = "grey")+
	geom_vline(xintercept = 24.55, color = "grey", linetype = "dashed")  +
	annotate("text", x = 11.7, y = 0.055, label = "Inflection point", size = 3) +
	annotate("text", x = 24, y = 0.070, label = "Pmax", size = 3) +
	annotate("segment", x = 14, xend = 16, y = 0.05, yend = 0.045, colour="black", size=0.5) +ylim(0, 0.075)

arrow=arrow(type = "closed")

sd_plot <- sder %>% 
	ggplot(aes(x = x, y = second_derivative)) + geom_line(size = 1.3) +
	xlim(0, 32) + theme_bw() + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank()) + 
	ylab("P''") + xlab(expression("Temperature (" *degree * "C)")) + geom_hline(yintercept = 0) +
	theme(text = element_text(size=14, family = "Helvetica")) +
	theme(axis.title.y=element_text(face="italic")) + geom_vline(xintercept = 16.59, color = "grey")+
	geom_vline(xintercept = 24.55, color = "grey", linetype = "dashed")


fd_plot <- fder %>% 
	ggplot(aes(x = x, y = first_derivative)) + geom_line(size = 1.3) +
	xlim(0, 32) + theme_bw() + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank()) + 
	ylab("P'") + xlab(expression("Temperature (" *degree * "C)")) + geom_hline(yintercept = 0) +
	theme(text = element_text(size=14, family = "Helvetica")) +
	theme(axis.title.y=element_text(face="italic")) + geom_vline(xintercept = 16.59, color = "grey")+
	geom_vline(xintercept = 24.55, color = "grey", linetype = "dashed")


## now plot them all together

predictions_plot <- plot_grid(TPC_plot, fd_plot, sd_plot, labels = c("A", "B", "C"), align = "v", ncol = 1)
save_plot("Tetraselmis_experiment/figures/figure1_predictions.pdf", predictions_plot, ncol = 1, base_height = 7, base_width = 4)


# let’s get the second derivative at a specific point ---------------------
## variable temps are 5, 10, 15, 20, 24, 29

prediction5 <- nbcurve_constant(5) + derivative(f = nbcurve_constant, x = 5, order = 2)*0.5*25
prediction10 <- nbcurve_constant(10) + derivative(f = nbcurve_constant, x = 10, order = 2)*0.5*25
prediction15 <- nbcurve_constant(15) + derivative(f = nbcurve_constant, x = 15, order = 2)*0.5*25
prediction20 <- nbcurve_constant(20) + derivative(f = nbcurve_constant, x = 20, order = 2)*0.5*25
prediction24 <- nbcurve_constant(24) + derivative(f = nbcurve_constant, x = 24, order = 2)*0.5*25
prediction27 <- nbcurve_constant(27) + derivative(f = nbcurve_constant, x = 27, order = 2)*0.5*25


temperature <- c(5, 10, 15, 20, 24, 27)
prediction <- c(prediction5, prediction10, prediction15, prediction20, prediction24, prediction27)
predictions <- data.frame(temperature, prediction)


# plot the predictions and the observations -------------------------------


all_r %>% 
	filter(variability == "v") %>% 
	rename(temperature = temp) %>% 
	ggplot(aes(x = temperature, y = estimate)) + geom_point(size = 2) +
	geom_errorbar(aes(ymin = estimate - std.error*1.96, ymax = estimate + std.error*1.96), width = 0.1) +
	geom_point(aes(x = temperature, y = prediction), color = "blue", data = predictions, size = 2)


### ok so it looks like the 2nd derivative switches from positive to negative at about 14.72


# extra code, I forget what this does -------------------------------------



### ok now let's try to estimate what we would predict for the variable temps just based on the data we collected

expected_v10 <- data.frame(prediction = (0.4033370 + 0.0968503)/2) %>% 
	mutate(temp = 10) %>% 
	mutate(var = "expectation")
expected_v5 <- data.frame(prediction = (0.3092085 + 0.0000000)/2) %>% 
	mutate(temp = 5) %>% 
	mutate(var = "expectation")
expected_v24 <- data.frame(prediction =(0.4649875 + 0.6368179)/2) %>% 
	mutate(temp = 24) %>% 
	mutate(var = "expectation")

expected_all <- bind_rows(expected_v10, expected_v24, expected_v5)

all <- bind_rows(all_r, expected_all) %>% 
	mutate(estimate = ifelse(is.na(estimate), prediction, estimate))


all %>% 
	# filter(var %in% c("variable temperature", "expectation")) %>% 
	ggplot(aes(x = temp, y = estimate, color = var)) + geom_point(size = 2) +
	geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.1)

### update on the comparison between variable predicted and expected
##5v
expected0 <- nbcurve2(0)[["a"]]*(1-0.58)
expected10 <- nbcurve2(10)[["a"]]*(0.58)

##10v
expected5 <- nbcurve2(5)[["a"]]*(1-0.43)
expected15 <- nbcurve2(15)[["a"]]*(0.43)

## 20v
expected15 <- nbcurve2(15)[["a"]]*(1-0.54)
expected25 <- nbcurve2(25)[["a"]]*(0.54)

##24v
expected19 <- nbcurve2(19)[["a"]]*(1-0.51)
expected29 <- nbcurve2(29)[["a"]]*(0.51)

## 27v
expected22 <- nbcurve2(22)[["a"]]*(1-0.51)
expected32 <- nbcurve2(32)[["a"]]*(0.51)

## 30v
expected25 <- nbcurve2(25)[["a"]]*(1-0.51)
expected35 <- nbcurve2(35)[["a"]]*(0.51)

expected_5v_2 <- data.frame(prediction = (expected0+expected10)) %>% 
	mutate(temp = 5) %>% 
	mutate(var = "expectation")

expected_10v_2 <- data.frame(prediction = (expected5+expected15)) %>% 
	mutate(temp = 10) %>% 
	mutate(var = "expectation")

expected_20v_2 <- data.frame(prediction = (expected15+expected25)) %>% 
	mutate(temp = 20) %>% 
	mutate(var = "expectation")

expected_24v_2 <- data.frame(prediction = (expected19+expected29)) %>% 
	mutate(temp = 24) %>% 
	mutate(var = "expectation")

expected_27v_2 <- data.frame(prediction = (expected22+expected32)) %>% 
	mutate(temp = 27) %>% 
	mutate(var = "expectation")


expected_30v_2 <- data.frame(prediction = (expected25+expected35)) %>% 
	mutate(temp = 30) %>% 
	mutate(var = "expectation")

all_expected <- bind_rows(expected_5v_2, expected_10v_2, expected_20v_2, expected_24v_2, expected_27v_2, expected_30v_2)

all_expected %>% 
	ggplot(aes(x = temp, y = prediction)) + geom_point() + geom_smooth()


all_2 <- bind_rows(all_r, all_expected) %>% 
	mutate(estimate = ifelse(is.na(estimate), prediction, estimate))


all_2 %>%
	filter(temp != 10) %>% 
	filter(var != "constant temperature") %>% 
	# filter(var %in% c("variable temperature", "expectation")) %>% 
	ggplot(aes(x = temp, y = estimate, color = var)) + geom_point(size = 2) +
	geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.1) + geom_smooth()