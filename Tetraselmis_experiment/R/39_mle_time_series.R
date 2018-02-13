
library(bbmle)
library(tidyverse)

cells_days_v_mod <- read_csv("Tetraselmis_experiment/data-processed/cells_days_v_mod.csv")
cells_exp_raw <- read_csv("Tetraselmis_experiment/data-processed/cells_exp_mod.csv")



# get data set up ---------------------------------------------------------

cells_exp <- cells_exp_raw %>% 
  mutate(days = time_since_innoc_hours/24) %>% 
  mutate(days = ifelse(days < 0,0, days)) %>% 
  select(temp, cell_density, days) 

cells_v <- cells_days_v_mod %>% 
  select(temp, cell_density, days)

cells <- cells_exp 

# define functions --------------------------------------------------------


growth_function <- function(a, b, z, w, temp, days){
  res <- 800 * (1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^(days)
  res
}


nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}


# Set up guesses ----------------------------------------------------------


z.guess <- 15
w.guess <- 30
a.guess <- 0.2
b.guess <- 0.1

avals<-seq(-0.2,1.2,0.05)
bvals<-seq(-0.2,0.3,0.05)
zvals<-seq(12,16,0.5)
mod.list<-list()
AIC.list<-c()


for(ia in 1:length(avals)){
  for(ib in 1:length(bvals)){
    for(iz in 1:length(zvals)){
    a.guess<-avals[ia]
    b.guess<-bvals[ib]
    z.guess<-zvals[iz]
res2<-try(fit<-mle2(cells$cell_density~dnorm(mean=growth_function(temp = cells$temp,days = cells$days, 
                                                               z=z,w=w,a=a,b=b),sd=s),
                    start=list(z=z.guess,w=w.guess,a=a.guess,b=b.guess,s=0.3),
                    skip.hessian=TRUE,data=cells))
if(class(res2)!="try-error"){
  mod.list<-append(mod.list,fit)
  AIC.list<-append(AIC.list,AIC(fit))
}
  }
}}

if(!is.null(AIC.list)){
  bestmodind<-which(AIC.list==min(AIC.list))
  if(length(bestmodind)>1){
    bestmodind<-sample(bestmodind,1)
  }
  bestmod<-mod.list[[bestmodind]]
  cfs<-coef(bestmod)
}


grfunc<-function(x){
  -nbcurve(x,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
}

optinfo<-optim(c(x=cfs[[1]]),grfunc)
opt<-optinfo$par[[1]]
maxgrowth<- -optinfo$value

z.list<-cfs[[1]]
w.list<-cfs[[2]]
a.list<-cfs[[3]]
b.list<-cfs[[4]]
s.list<-cfs[[5]]
topt.list<-opt
maxgrowth.list<-maxgrowth

fits <-data.frame(topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,s.list) 

