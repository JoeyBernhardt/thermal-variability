## now let's try to fit the curve with a modified beta curve


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



## a, b, z, w, 

# nbcurve<-function(temp,z,w,a,b){
# 	res<-a*(exp(-exp((b*(temp - z)-8)) - (w*(temp - z)^2)))
# 	res
# }


nbcurve <- function(temp, a, b, c, d, e) {
	res <- a*((((temp - b) + ((c*(d-1))/(d+e-2)))/c)^(d-1))(((1-((temp - b) + ((c*(d-1))/(d+e-2)))/c)^(e-1)))/ ((((d-1)/(d+e-2))^(d-1))(((e-1)/(d+e-2))^(e-1)))
res	
}


curve.id.list<-unique(dat.full$curve.id)	

# Create empty vectors to populate with parameter values and trait estimates
a.list<-rep(NA, length(curve.id.list))				#Parameter 'z'
b.list<-rep(NA, length(curve.id.list))				#Parameter 'w', which is the temperature niche width
c.list<-rep(NA, length(curve.id.list))				#Parameter 'a'
d.list<-rep(NA, length(curve.id.list))				#Parameter 'b'
e.list<-rep(NA, length(curve.id.list))				#Parameter 'b'
# topt.list<-rep(NA, length(curve.id.list))			#Topt, Optimum temperature for growth
maxgrowth.list<-rep(NA, length(curve.id.list))		#Maximum growth rate (i.e. growth rate at Topt)
rsqr.list<-rep(NA, length(curve.id.list))			#R^2 values for all fits
s.list<-rep(NA, length(curve.id.list))				#Error
n.list<-rep(NA, length(curve.id.list))				#Number of growth rate measurements used in the curve



for(i in 1:length(curve.id.list)){
	print(i)
	
	# Take a subset of the data corressponding to the ith curve.id.list value
	dat<-subset(dat.full,dat.full$curve.id==curve.id.list[i])
	
	# guess starting values for parameters 'z' and 'w'
	# z.guess<-mean(dat$temperature[dat$growth.rate==max(dat$growth.rate)])		#starting estimates for 'z'
	# w.guess<-diff(range(dat$temperature))										#starting estimates for niche width
	
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
			
			
			res2<-try(fit<-mle2(dat$growth.rate~dnorm(mean=nbcurve(dat$temperature,a=a, b=b,c=c,d=d, e = e),sd=s),start=list(a=1,b=1,c=1,d=1,e = 1, s=0.3), 
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
		expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]], cfs[[5]])
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
		expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]], cfs[[5]])
		rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)
	}
	
	
	
	
	
	#stash results		
	rsqr.list[i]<-rsqr
	a.list[i]<-cfs[[1]]
	b.list[i]<-cfs[[2]]
	c.list[i]<-cfs[[3]]
	d.list[i]<-cfs[[4]]
	e.list[i]<-cfs[[5]]
	# s.list[i]<-cfs[[5]]
	# topt.list[i]<-opt
	# maxgrowth.list[i]<-maxgrowth
	n.list[i]<-length(dat$temperature)
}

