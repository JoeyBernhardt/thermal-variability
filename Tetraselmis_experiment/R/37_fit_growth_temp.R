

library(purrr)
library(bbmle)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)


growth_all <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling.csv")
growth_all_v <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling_v.csv")

cells_exp <- read_csv("Tetraselmis_experiment/data-processed/cells_exp.csv")

dat.full <- growth_all %>%
  group_by(temp) %>% 
  summarise(growth.rate = mean(growth_per_day)) %>% 
  rename(temperature = temp) %>% 
  mutate(curve.id = "1") %>% 
  select(curve.id, temperature, growth.rate) 
str(dat.full)

freeze <- data.frame(curve.id = "1", temperature = -1.8, growth.rate = 0)

dat.full <- bind_rows(dat.full, freeze)


dat <- cells_exp %>% 
  rename(temperature = temp) %>% 
  rename(hours = time_since_innoc_hours)

z.guess<-mean(dat$temperature[dat$growth.rate==max(dat$growth.rate)])		#starting estimates for 'z'
w.guess<-diff(range(dat$temperature))					

nbcurve<-function(temp,z,w,a,b, hours){
  res<-800*exp(((1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^hours))
  res
}



avals<-seq(-0.2,1.2,0.1)		
bvals<-seq(-0.2,0.3,0.05)
mod.list<-list()
AIC.list<-c()

for(ia in 1:length(avals)){
  for(ib in 1:length(bvals)){
    a.guess<-avals[ia]
    b.guess<-bvals[ib]
    res2<-try(fit<-mle2(dat$cell_density~dnorm(mean=nbcurve(dat$temperature,z=z,w=w,a=a,b=b,dat$hours),sd=s),start=list(z=16,w=30,a=0.2,b=0.1,s=0.3),
                        skip.hessian=TRUE,data=dat))
    
    
    if(class(res2)!="try-error"){
      mod.list<-append(mod.list,fit)
      AIC.list<-append(AIC.list,AIC(fit))
    }
  }
}

res2<-try(fit<-mle2(dat$cell_density~dnorm(mean=nbcurve(dat$temperature,z=z,w=w,a=a,b=b, hours = hours),sd=s),start=list(z=16,w=30,a=0.2,b=0.1,s=0.3),
                    skip.hessian=TRUE,data=dat))



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


fits_vals <- fits_real_constant %>% 
  summarise_each(funs(mean, sd), z, w, a, b)

a*exp(b*temp)*(1-((temp-z)/(w/2))^2)

library(minpack.lm)

cells_days <- cells_exp %>% 
  mutate(days = time_since_innoc_hours/24) %>% 
  mutate(days = ifelse(days < 0, 0, days))

cells_days %>% 
  ggplot(aes(x = days, y = cell_density, color = factor(temp))) + geom_point() 



fit_c <- nlsLM(cell_density ~ 800 * (1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^(time_since_innoc_hours),
    data= cells_days,  start=list(z=16,w=36,a=0.2,b=0.1),
    control = nls.control(maxiter=1000, minFactor=1/204800000))

avals<-seq(-0.2,1.2,0.02)
bvals<-seq(-0.2,0.3,0.02)


df <-  expand.grid(a = avals, b = bvals) %>% 
  mutate(unique_id = rownames(.))

fit_d <- nlsLM(cell_density ~ 800 * (1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^(days),
               data= cells_days,  
               start=list(z=14.4,w=35,a=-0.2, b=-0.2),
               lower = c(0, 0, -0.2, -0.2),
               upper = c(30, 40, 1.2, 0.3),
               control = nls.control(maxiter=1024, minFactor=1/204800000))


fit_growth <- function(df){
  res <- try(nlsLM(cell_density ~ 800 * (1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^(days),
      data= cells_days,  
      start=list(z=14.4,w=35,a= df$a[[1]], b=df$b[[1]]),
      lower = c(0, 0, -0.2, -0.2),
      upper = c(30, 40, 1.2, 0.3),
      control = nls.control(maxiter=1024, minFactor=1/204800000)))
  if(class(res)!="try-error"){
    out1 <- tidy(res) %>% 
      select(estimate, term) %>% 
      spread(key = term, value = estimate)
    out2 <- glance(res)
    
  }
  all <- bind_cols(out1, out2)
  all
}


df_split <- df %>% 
  split(.$unique_id)


output <- df_split %>%
  map_df(fit_growth, .id = "run") 
  
tp <- output %>% 
  top_n(n = -1, wt = AIC) 

cfd <- (coef(fit_d))

tpc1<-function(x){
  res<-(cfd[["a"]]*exp(cfd[["b"]]*x)*(1-((x-cfd[["z"]])/(cfd[["w"]]/2))^2))
  res
}

tpc2<-function(x){
  res<-(tp$a[[1]]*exp(tp$b[[1]]*x)*(1-((x-tp$z[[1]])/(tp$w[[1]]/2))^2))
  res
}

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + stat_function(fun = tpc2, color = "black", size = 1) +xlim(-4, 33) + 
  ylim(-2, 5) + geom_hline(yintercept = 0) +
  stat_function(fun = curve_constant_resamp, color = ic[3], size = 1.5, alpha = 0.7)
  
nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}
nbcurve1<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

grfunc<-function(x){
  -nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]])
}

?optim
optinfo<-optim(c(x=tp$z[[1]]),grfunc)
opt <-optinfo$par[[1]]
maxgrowth <- -optinfo$value

library(rootSolve)
uniroot.all(function(x) nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]]),c(opt,150))
uniroot.all(function(x) nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]]),c(-3,opt))


tidy(fit_d) %>% View
tidy(fit_d) %>% View

days <- 2
temp <- 27
growth_fun_27 <- function(x, temp = 27){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}

growth_fun_5 <- function(x, temp = 5){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}

growth_fun_10 <- function(x, temp = 10){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}
growth_fun_24 <- function(x, temp = 24){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}

growth_fun_16 <- function(x, temp = 16){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}

growth_fun_29 <- function(x, temp = 29){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}

growth_fun_32 <- function(x, temp = 32){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}
growth_fun_0 <- function(x, temp = 0){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}
growth_fun_20 <- function(x, temp = 20){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}


library(viridis)
library(colormap)
library(cowplot)
ic <- colormap(colormap = colormaps$viridis, nshades = 9, format = "hex",
               alpha = 1, reverse = FALSE)

cells_days %>% 
  ggplot(aes(x = days, y = cell_density, color = factor(temp))) + geom_point(size = 2) +
  geom_point(size = 2, shape= 1, color = "black") +
  stat_function(fun = growth_fun_27, color = ic[7], size = 1) +ylim(0, 30000) +
  stat_function(fun = growth_fun_5, color = ic[2], size = 1) +
  stat_function(fun = growth_fun_10, color = ic[3], size = 1) +
  stat_function(fun = growth_fun_24, color = ic[6], size = 1) +
  stat_function(fun = growth_fun_16, color = ic[4], size = 1) +
  stat_function(fun = growth_fun_29, color = ic[8], size = 1) +
  stat_function(fun = growth_fun_32, color = ic[9], size = 1) +
  stat_function(fun = growth_fun_0, color = ic[1], size = 1) +
  stat_function(fun = growth_fun_20, color = ic[5], size = 1) +
   scale_color_viridis(discrete = TRUE) +
  ylab("Cell density") + xlab("Days")



cf <- (coef(fit_c))


expected<-nbcurve(dat$temperature,cfd[[1]],cfd[[2]],cfd[[3]],cfd[[4]])
rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)


fits_c$a.list[1]*exp(fits_c$b.list[1]*x)*(1-((x-fits_c$z.list[1])/(fits_c$w.list[1]/2))^2)

tpch<-function(x){
  res<-(cf[["a"]]*exp(cf[["b"]]*x)*(1-((x-cf[["z"]])/(cf[["w"]]/2))^2))*24
  res
}



p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + stat_function(fun = tpch, color = "black", size = 1) +xlim(-2, 33) + 
  ylim(-2, 1.6)

temp <- 24
r <- cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)

