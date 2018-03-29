### attempt at plotting and fitting curves to the thermal variability growth rate data
### april 27 2017
### next steps: need to get a better fit for the constant temperature curve. Right now it's underestimating the peak of the curve




library(tidyverse)
library(stringr)
library(rootSolve)
library(bbmle)
library(extrafont)
loadfonts()

all_r1 <- read_csv("Tetraselmis_experiment/data-processed/all_r_with0.csv")
all_r <- read_csv("Tetraselmis_experiment/data-processed/growth_estimates_round3.csv")

all_r %>% 
	filter(temp < 31) %>% 
	ggplot(aes(x = temp, y = estimate, color = variability)) + geom_point() + geom_line()


data_full <- all_r %>% 
	# filter(var == "variable temperature") %>% 
	select(-temperature) %>% 
	rename(temperature = temp) %>% 
	rename(growth.rate = estimate) %>%
  mutate(growth.rate = growth.rate*24) %>% 
	# mutate(var = str_replace(var, "constant temperature" ,"1")) %>% 
	# mutate(var = str_replace(var, "variable temperature" ,"2")) %>% 
	rename(`thermal environment` = variability) %>% 
	select(`thermal environment`, temperature, growth.rate) %>% 
	filter(temperature < 33)

dat.full <- all_r %>% 
	filter(variability == "c") %>% 
	select(-temperature) %>% 
	rename(temperature = temp) %>% 
	rename(growth.rate = estimate) %>% 
  mutate(growth.rate = growth.rate*24) %>% 
	mutate(var = str_replace(variability, "constant temperature" ,"1")) %>% 
	mutate(var = str_replace(variability, "variable temperature" ,"2")) %>% 
	rename(curve.id = variability) %>% 
	select(curve.id, temperature, growth.rate) %>% 
	filter(temperature < 33)


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
	
	
	#Save .png plot with fitted curve. File is saved with the curve.id.list value as the name
	png(paste(curve.id.list[i],'.png',sep=''))
	plot(dat$growth.rate~dat$temperature,ylim=c(pmin(0,min(dat$growth.rate)),max(dat$growth.rate)+(0.2)*max(dat$growth.rate)),main=curve.id.list[1],
			 xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
	plot(dat$growth.rate~dat$temperature,ylim=c(pmin(0,min(dat$growth.rate)),max(dat$growth.rate)+(0.2)*max(dat$growth.rate)),main=curve.id.list[2],
			 xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)

	plot(dat$growth.rate~dat$temperature,ylim=c(pmin(0,min(dat$growth.rate)),max(dat$growth.rate)+(0.2)*max(dat$growth.rate)),main=curve.id.list[i],
			 xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
	curve(nbcurve(x,cfs[1],cfs[2],cfs[3],cfs[4]),col='red', lwd=2,add=TRUE)
	
	
	
	curve(nbcurve(x,17.03711,30.94049,0.2599170,0.04544617),col='red', lwd=4,add=TRUE)
	curve(nbcurve(x,16.57264,35.61432,0.1203635,0.07903706),col='black', lwd=4,add=TRUE)
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

fits2<-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) ## constant
fits3<-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list) ## variable

### now make the plots!
nbcurve1<-function(x){
	res<-fits2$z.list[1]*exp(fits2$w.list[1]*x)*(1-((x-fits2$a.list[1])/(fits2$b.list[1]/2))^2)
	res
}


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p + geom_point(aes(x = temperature, y = growth.rate), data = dat.full, size = 0.5, alpha = 0) +
	# stat_function(fun = nbcurve1, color = "red", size = 2) +
	stat_function(fun = nbcurve2, color = "black", size = 1) +
	stat_function(fun = nbcurve3, color = "blue", size = 2) +
	geom_point(aes(x = temperature, y = growth.rate, color = `thermal environment`), data = data_full, size = 2, alpha = 0.9) + scale_color_manual(values = c("black", "red")) 


trans = yearmon_trans()
k <- 8.62 * 10^(-5)

temp_arr_trans <- function(x) {(1/(.00008617*(x+273.15)))}
temp_arr_trans(15)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = 0))

dat.full2 <- dat.full %>% 
  mutate(inverse_temp = (1/(.00008617*(temperature+273.15))))

p + geom_point(aes(x = inverse_temp, y = growth.rate), data = dat.full2, size = 0.5, alpha = 0) +
  stat_function(fun = nbcurve2, color = "black", size = 1) +
  xlab("Temperature (°C)") +
  ylab(bquote('Exponential growth rate'*~day^-1*'')) +
  theme_bw(base_family = "Arial", base_size = 12) +
  scale_x_continuous(sec.axis = sec_axis(~(temp_arr_trans(.))), limits = c(0, 32)) + 
   ylim(-0.2, 1.7) +
  geom_hline(yintercept = 0, color = "grey") +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  theme_bw() +
  theme(text = element_text(size=12, family = "Arial"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Temperature (1/kT)") 
ggsave("Tetraselmis_experiment/figures/TPC_old_k-temp.pdf", width = 5, height = 3.5)


nbcurve2<-function(x){
	res<-cfs[3]*exp(cfs[4]*x)*(1-((x-cfs[1])/(cfs[2]/2))^2)
	res
}

uniroot.all(function(x) nbcurve2(x),c(20,150))
                    

nbcurve3<-function(x){
	res<-0.1203635*exp(0.0899376*x)*(1-((x-13.57264)/(41.61432/2))^2)
	res
}
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p + geom_point(aes(x = temperature, y = growth.rate), data = dat.full, size = 2, alpha = 1) +
	stat_function(fun = nbcurve2, color = "red", size = 2) +
	# stat_function(fun = nbcurve2, color = "black", size = 2) +
	# stat_function(fun = nbcurve3, color = "blue", size = 2) +
	geom_point(aes(x = temperature, y = growth.rate, color = `thermal environment`), data = data_full, size = 2, alpha = 0.9) + scale_color_manual(values = c("black", "red")) +
	geom_rect(aes(xmin=0, xmax=14.72, ymin=-Inf, ymax=Inf), alpha = 0.4) + ylim(0, 0.1)



data_full <- all_r %>% 
	# filter(var == "variable temperature") %>% 
	select(-temperature) %>% 
	rename(temperature = temp) %>% 
	rename(growth.rate = estimate) %>% 
	# mutate(var = str_replace(var, "constant temperature" ,"1")) %>% 
	# mutate(var = str_replace(var, "variable temperature" ,"2")) %>% 
	rename(`thermal environment` = var) %>% 
	select(`thermal environment`, temperature, growth.rate) 

xvals<-seq(0,35,0.001)

df <- as.data.frame(xvals)

df_preds_5 <- df %>% 
	mutate(prediction = (nbcurve3(xvals-5) + nbcurve3(xvals+5))/2) %>% 
	mutate(prediction_10 = (nbcurve3(xvals-10) + nbcurve3(xvals+10))/2) %>% 
	mutate(tpc = (nbcurve3(xvals) + nbcurve3(xvals))/2) %>% 
	mutate(diff = prediction - tpc) %>% 
	mutate(effect_of_variability = ifelse(diff > 0, "beneficial", "detrimental")) %>% 
	rename(temperature = xvals)

df_preds_2 <- df %>% 
	mutate(prediction = (nbcurve2(xvals-5) + nbcurve2(xvals+5))/2) %>% 
	mutate(prediction_10 = (nbcurve2(xvals-10) + nbcurve2(xvals+10))/2) %>% 
	mutate(tpc = (nbcurve2(xvals) + nbcurve2(xvals))/2) %>% 
	mutate(diff = prediction - tpc) %>% 
	mutate(effect_of_variability = ifelse(diff > 0, "beneficial", "detrimental")) %>% 
	rename(temperature = xvals)




df_preds_10 <- df %>% 
	mutate(prediction = (nbcurve2(xvals-10) + nbcurve2(xvals+10))/2) %>% 
	rename(temperature = xvals)

### predictions
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p + geom_point(aes(x = temperature, y = growth.rate), data = data_full, size = 0.5, alpha = 0) +
	stat_function(fun = nbcurve1, color = "red", size = 2) +
	stat_function(fun = nbcurve2, color = "black", size = 2) +
	# stat_function(fun = nbcurve3, color = "green", size = 2) +
	geom_point(aes(x = temperature, y = growth.rate, color = `thermal environment`), data = data_full, size = 2, alpha = 0.9) + scale_color_manual(values = c("black", "red")) +
	geom_line(aes(x = temperature, y = prediction), data = df_preds_2, size = 2, color = "blue") +
	# geom_line(aes(x = temperature, y = prediction), data = df_preds_2, size = 2, color = "orange") +
	# geom_line(aes(x = temperature, y = prediction_10), data = df_preds_5, size = 2, color = "orange") +
	theme_bw() + xlab("temperature (C)") + ylab("intrinsic growth rate (r)") +
	geom_rect(aes(xmin=0, xmax=17.42, ymin=-Inf, ymax=Inf), alpha = 0.4) + ylim(-0.2, 0.9)
ggsave("Tetraselmis_experiment/figures/TPCs-with-data.pdf")


### variable treatments (blue is predicted, red is observed)
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p + geom_point(aes(x = temperature, y = growth.rate), data = data_full, size = 0.5, alpha = 0) +
	stat_function(fun = nbcurve1, color = "red", size = 2) +
	stat_function(fun = nbcurve2, color = "black", size = 2) +
	# stat_function(fun = nbcurve3, color = "green", size = 2) +
	geom_line(aes(x = temperature, y = prediction), data = df_preds_2, size = 2, color = "blue") +
	geom_point(aes(x = temperature, y = growth.rate, color = `thermal environment`), data = data_full, size = 2, alpha = 0.9) + scale_color_manual(values = c("black", "red")) +
	theme_bw() + xlab("temperature (C)") + ylab("intrinsic growth rate (r)") +
	geom_rect(aes(xmin=0, xmax=17.42, ymin=-Inf, ymax=Inf), alpha = 0.4) + ylim(-0.2, 0.9)
ggsave("Tetraselmis_experiment/figures/fitted-TPCs-with-data.pdf")





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



par(mfrow = c(1, 2))
grid <- seq(14.716, 14.717, by = 0.001)
for (i in 2:2) {
	plot(grid, sapply(grid, derivative, f = nbcurve2, order = i), type = "l", xlab = "x", 
			 ylab = "y", main = paste("Order", i, "Derivative"))
	abline(h = 0, col = "gray60")
}


### ok so it looks like the 2nd derivative switches from positive to negative at about 14.72

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

