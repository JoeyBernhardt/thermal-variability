

library(tidyverse)
library(cowplot)


ts <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_isolate11_2005.csv")
ts <- read_csv("Tetraselmis_experiment/data-processed/daily_temps_isolate89_2005.csv")


results5 <- read_csv("Tetraselmis_experiment/data-processed/results5.csv")

isolate_11 <- results5 %>% 
  filter(isolate.code == 11)
d <- isolate_11
d <- results5 %>% 
  filter(isolate.code == 89)

ts_89 %>% 
  ggplot(aes(x = time, y = sst)) + geom_line()


ts_89 %>% 
  ggplot(aes(x = sst)) + geom_density()


nbcurve<-function(temp,z,w,a,b){
  res<-d$a*exp(d$b*temp)*(1-((temp-d$z)/(d$w/2))^2)
  res
}


lin_avg <- ts %>% 
  rename(temp = sst) %>% 
  mutate(growth_rate = d$a*exp(d$b*temp)*(1-((temp-d$z)/(d$w/2))^2)) %>% 
  summarise_each(funs(mean, sd), growth_rate)

temp <- mean(ts$sst)
res <- d$a*exp(d$b*temp)*(1-((temp-d$z)/(d$w/2))^2)


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 
p + 
  stat_function(fun = nbcurve) + xlim(-2, 40) + ylim(0, 1) +
  geom_vline(xintercept = mean(ts$sst)) +
  geom_vline(xintercept = mean(ts$sst) + sd(ts$sst), color = "green") +
  geom_vline(xintercept = mean(ts$sst) - sd(ts$sst), color = "green") +
  geom_vline(xintercept = max(ts$sst), color = "orange") +
  geom_vline(xintercept = min(ts$sst), color = "orange") +
  ylab("Population growth rate") + xlab("Temperature (°C)") +
  geom_hline(yintercept = res[[1]], color = "blue") + 
  geom_hline(yintercept = lin_avg$growth_rate_mean)
  
  