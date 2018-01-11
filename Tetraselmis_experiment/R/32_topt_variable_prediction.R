library(tidyverse)
library(bbmle)
library(cowplot)
library(colormap)

variable_predictions_points <- read_csv("Tetraselmis_experiment/data-processed/variable_predictions_points.csv")

dat.full <- variable_predictions_points %>% 
  gather(key = curve.id, value = growth.rate, 2:3) 

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
write_csv(fits, "Tetraselmis_experiment/data-processed/predicted_parms_variable.csv")

curve_fit_low<-function(x){
  res<-fits$a.list[1]*exp(fits$b.list[1]*x)*(1-((x-fits$z.list[1])/(fits$w.list[1]/2))^2)
  res
}

curve_fit_high<-function(x){
  res<-fits$a.list[2]*exp(fits$b.list[2]*x)*(1-((x-fits$z.list[2])/(fits$w.list[2]/2))^2)
  res
}

dat.full %>% 
  ggplot(aes(x = temperature, y = growth.rate, group = curve.id)) + geom_line() +
  stat_function(fun = curve_fit_high, color = "red") +
  stat_function(fun = curve_fit_low, color = "green") +
  ylim(0, 1.5)
  


