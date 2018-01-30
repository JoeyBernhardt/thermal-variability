

library(tidyverse)
library(cowplot)
library(DescTools)


ts <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_isolate11_2005.csv")
# ts <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_isolate89_2005.csv")
results5 <- read_csv("Tetraselmis_experiment/data-processed/results5.csv")

isolate_11 <- results5 %>% 
  filter(isolate.code == 11)
d <- isolate_11


ts %>% 
  ggplot(aes(x = time, y = sst)) + geom_line()


ts_89 %>% 
  ggplot(aes(x = sst)) + geom_density()


growth <- ts %>% 
  rename(temp = sst) %>% 
  mutate(growth_rate = d$a*exp(d$b*temp)*(1-((temp-d$z)/(d$w/2))^2)) %>% 
  mutate(day = as.numeric(rownames(.)))

write_csv(growth, "Tetraselmis_experiment/data-processed/growth_toy.csv")

AUC(x = growth$day, y = growth$growth_rate)/366

lin_avg <- ts %>% 
  rename(temp = sst) %>% 
  mutate(growth_rate = d$a*exp(d$b*temp)*(1-((temp-d$z)/(d$w/2))^2)) %>% 
  summarise_each(funs(mean, sd), growth_rate)

library(tidyverse)

## read in data

growth <- read.csv("Tetraselmis_experiment/data-processed/growth_toy.csv")

## make a function that approximates the time series
data_fun <- approxfun(growth$day, growth$growth_rate, method="constant", 0, 0) 

## this is equivalent to equation 2.3 in Vasseur et al. 2014 (I think), here we have 366 days
out <- integrate(data_fun, lower = 1, upper = 366, subdivisions = 2000) ## integrate over the time series

## the value of equation 2.3 gives the total 'area under the curve' i.e. the cumulative growth rate
out$value

### but what makes more sense to me is this (which divides by number of days)
out$value/365


AUC(x, y)

plot(x, y)

x <- seq(10, 20, by = 0.1)
y <- sapply(x, data_fun)


growth %>% 
  ggplot(aes(x = day, y = growth_rate)) + geom_line() +
  stat_function(fun = data_fun, color = "red")


nbcurve<-function(temp,z,w,a,b){
  res<-d$a*exp(d$b*temp)*(1-((temp-d$z)/(d$w/2))^2)
  res
}




temp <- mean(ts$sst)
res <- d$a*exp(d$b*temp)*(1-((temp-d$z)/(d$w/2))^2)


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 
p + 
  stat_function(fun = nbcurve) + xlim(-2, 10) + ylim(0, 1) +
  geom_vline(xintercept = mean(ts$sst)) +
  geom_vline(xintercept = mean(ts$sst) + sd(ts$sst), color = "green") +
  geom_vline(xintercept = mean(ts$sst) - sd(ts$sst), color = "green") +
  geom_vline(xintercept = max(ts$sst), color = "orange") +
  geom_vline(xintercept = min(ts$sst), color = "orange") +
  ylab("Population growth rate") + xlab("Temperature (°C)") +
  geom_hline(yintercept = res[[1]], color = "blue") + 
  geom_hline(yintercept = lin_avg$growth_rate_mean)



nbcurve1<-function(temp,z,w,a,b){
  d$a*exp(d$b*temp)*(1-((temp-d$z)/(d$w/2))^2)
}


integrate(nbcurve1, lower = 1, upper = 10)  
nbcurve1(1)
(nbcurve1(1) + nbcurve1(10))/2
  