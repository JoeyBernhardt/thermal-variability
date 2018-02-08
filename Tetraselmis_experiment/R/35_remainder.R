
library(tidyverse)
library(purrr)
library(viridis)
library(broom)
library(cowplot)


### goal, find the 3rd derivative of each TPC at the mean temperature

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



thomas <- read_csv("data/thermal_trait_data/Thomas_2014_traits_derived_20140606.csv")
thomas2 <- thomas %>% 
  rename(topt = mu.g.opt.list) %>% 
  rename(w = mu.wlist,
         a = mu.alist,
         z = mu.c.opt.list,
         b = mu.blist) %>% 
  select(topt, w, a, b, z, everything()) 


thomas3 <- thomas2 %>% 
  rename(latitude = isolation.latitude,
         longitude = isolation.longitude) %>% 
  filter(!is.na(latitude), !is.na(longitude)) %>% 
  filter(habitat %in% c("marine", "estuarine", "saltmarsh")) %>% 
  filter(curvequal == "good") %>% 
  filter(mu.rsqrlist > 0.85)


all <- left_join(thomas3, sst_summary)


## derivative(f = tpc, x = x, order = 2)

TPC <- function(x) all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)

all %>%
 mutate(third_deriv = derivative(f = function(x) a*exp(b*x)*(1-((x-z)/(w/2))^2)), x = sst_mean, order = 3) %>% View 

derivative(TPC, x = all$sst_mean[1], order = 3)

all_split <- all %>% 
  split(.$isolate.code)

third_deriv_function <- function(df){
  TPC <- function(x) df$a[1]*exp(df$b[1]*x)*(1-((x-df$z[1])/(df$w[1]/2))^2)
  deriv_out <- derivative(TPC, x = df$sst_mean[1], order = 3)
 output <- data.frame(deriv_out)
}

third_derivs <- all_split %>% 
  map_df(third_deriv_function, .id = "isolate.code")

third_derivs %>% 
  mutate(abs_deriv = abs(deriv_out)) %>% View

third_derivs %>% 
  ggplot(aes(x = deriv_out)) + geom_histogram()
