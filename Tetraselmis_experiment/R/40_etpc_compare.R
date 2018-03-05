library(cowplot)
library(tidyverse)

etpc_fits <- read_csv("Tetraselmis_experiment/data-processed/time_resampling_fits.csv")
ctpc_fits <- read_csv("Tetraselmis_experiment/data-processed/ctpc.csv")
vtpc_fits <- read_csv("Tetraselmis_experiment/data-processed/time_resampling_fits_v.csv")
c_fits <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_exp.csv")



c_fits <- c_fits %>% 
  rename(z = z.list,
         a = a.list, 
         b = b.list, 
         w = w.list)



etpc <-function(x){
  res<-(etpc_fits$a[[1]]*exp(etpc_fits$b[[1]]*x)*(1-((x-etpc_fits$z[[1]])/(etpc_fits$w[[1]]/2))^2))
  res
}

ctpc <-function(x){
  res<-(ctpc_fits$a[[1]]*exp(ctpc_fits$b[[1]]*x)*(1-((x-ctpc_fits$z[[1]])/(ctpc_fits$w[[1]]/2))^2))
  res
}

vtpc_fits <- read_csv("Tetraselmis_experiment/data-processed/vtpc.csv") ## new vtpc

vtpc <-function(x){
  res<-(vtpc_fits$a[[1]]*exp(vtpc_fits$b[[1]]*x)*(1-((x-vtpc_fits$z[[1]])/(vtpc_fits$w[[1]]/2))^2))
  res
}
 

tpc_c <-function(x){
  res<-(c_fits$a[[1]]*exp(c_fits$b[[1]]*x)*(1-((x-c_fits$z[[1]])/(c_fits$w[[1]]/2))^2))
  res
}


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + 
  stat_function(fun = etpc, color = "black", size = 1) +
  stat_function(fun = tpc_c, color = "green", size = 1) +
  xlim(-3, 33) + 
  ylim(-1, 2) + geom_hline(yintercept = 0) + ylab("Exponential growth rate (per day)") +
  xlab("Temperature (°C)")


library(rootSolve)
nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

uniroot.all(function(x) nbcurve(x, z = etpc_fits$z[[1]],w = etpc_fits$w[[1]],a = etpc_fits$a[[1]], b = etpc_fits$b[[1]]),c(opt,150))
uniroot.all(function(x) nbcurve(x, z = c_fits$z[[1]],w = c_fits$w[[1]],a = c_fits$a[[1]],b = c_fits$b[[1]]),c(opt,150))

uniroot.all(function(x) nbcurve(x, z = etpc_fits$z[[1]],w = etpc_fits$w[[1]],a = etpc_fits$a[[1]],b = etpc_fits$b[[1]]),c(-15,opt))
uniroot.all(function(x) nbcurve(x, z = c_fits$z[[1]],w = c_fits$w[[1]],a = c_fits$a[[1]],b = c_fits$b[[1]]),c(-15,opt))


tpc_fit <- c_fits

e_grfunc<-function(x){
  -nbcurve(x, z = etpc_fits$z[[1]],w = etpc_fits$w[[1]],a = etpc_fits$a[[1]],b = etpc_fits$b[[1]])
}
e_optinfo<-optim(c(x=etpc_fits$z[[1]]),e_grfunc)
e_opt <-e_optinfo$par[[1]]
e_maxgrowth <- -e_optinfo$value


grfunc<-function(x){
  -nbcurve(x, z = c_fits$z[[1]],w = c_fits$w[[1]],a = c_fits$a[[1]],b = c_fits$b[[1]])
}
optinfo<-optim(c(x=c_fits$z[[1]]),grfunc)
opt <-optinfo$par[[1]]
maxgrowth <- -optinfo$value



# now plot on CIs for the etpc approach -----------------------------------

bs_v <- read_csv("Tetraselmis_experiment/data-processed/bootstrap_time_series_fitsv.csv")
bs_v <- read_csv("Tetraselmis_experiment/data-processed/bootstrap_time_series_fitsv_bounds.csv")
bs_v <- read_csv("Tetraselmis_experiment/data-processed/nls_boot.csv") %>% 
  mutate(replicate = rownames(.))

bs_c <- read_csv("Tetraselmis_experiment/data-processed/nls_boot_c.csv") %>% 
  mutate(replicate = rownames(.))
bs_c <- read_csv("Tetraselmis_experiment/data-processed/bootstrap_time_series_fits.csv")

x <- seq(-3, 32, by = 0.01)


predict_tpc <- function(x) {
  y <- 0.5*(etpc(x + 5) + etpc(x - 5))
}


predictions <- sapply(x, predict_tpc)
predictions <- data.frame(x, predictions) %>% 
  rename(temperature = x, 
         growth = predictions)




all_preds_average <- predictions



prediction_function <- function(df) {
  tpc <-function(x){
    res<-(df$a[[1]]*exp(df$b[[1]]*x)*(1-((x-df$z[[1]])/(df$w[[1]]/2))^2))
    res
  }
  
  pred <- function(x) {
    y <- tpc(x)
  }
  
  x <- seq(-3, 33, by = 0.1)
  
  preds <- sapply(x, pred)
  preds <- data.frame(x, preds) %>% 
    rename(temperature = x, 
           growth = preds)
}


bs_v_split <- bs_v %>%
  # filter(z != 30) %>% 
  # filter(b > 0, a != 1, z != 20, z != 0, b != 0) %>%
  # filter(b < 0.12) %>% 
  split(.$replicate)

all_preds_v <- bs_v_split %>% 
  map_df(prediction_function, .id = "replicate")


limits_v <- all_preds_v %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            mean = mean(growth))

bs_split <- bs_c %>% 
  split(.$replicate)

all_preds_c <- bs_split %>% 
  map_df(prediction_function, .id = "replicate")

limits_c <- all_preds_c %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            mean = mean(growth))

bs_v2 <- bs_v %>% 
  mutate(replicate = as.character(replicate))

all_variable <- left_join(all_preds_v, bs_v2, by = "replicate")


hist(bs_v$a)
hist(bs_v$b)
hist(bs_v$w)
hist(bs_v$z)


# NLA predictions ---------------------------------------------------------

prediction_NLA <- function(df) {
  tpc <-function(x){
    res<-(df$a[[1]]*exp(df$b[[1]]*x)*(1-((x-df$z[[1]])/(df$w[[1]]/2))^2))
    res
  }
  
  pred <- function(x) {
    y <-  0.5*(tpc(x + 5) + tpc(x - 5))
  }
  
  x <- seq(-3, 32, by = 0.1)
  
  preds <- sapply(x, pred)
  preds <- data.frame(x, preds) %>% 
    rename(temperature = x, 
           growth = preds)
}


all_preds_NLA <- bs_split %>% 
  map_df(prediction_NLA, .id = "replicate")



limits_prediction <- all_preds_NLA %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            mean = mean(growth))


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + 
  # stat_function(fun = vtpc, color = "orange", size = 1) +
  # stat_function(fun = tpc_c, color = "green", size = 1) +
  geom_line(aes(x = temperature, y = growth, group = replicate), color = "cadetblue", data = all_preds_c, alpha = 0.2) +
  stat_function(fun = ctpc, color = "black", size = 1) +
  stat_function(fun = vtpc, color = "orange", size = 1) +
  # # geom_line(aes(x = temperature, y = growth, group = replicate), color = "orange", data = all_preds_v, alpha = 0.1) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_v, fill = "orange", alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_c, fill = "cadetblue", alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
              fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
  # geom_point(aes(x = temp, y = mean), data = gs, size = 2, color = "cadetblue") +
  geom_point(aes(x = temp, y = mean), data = gsv, size = 2, color = "orange") +
  # geom_errorbar(aes(ymin = lower, ymax = upper, x = temp), data = gs, color = "cadetblue", width = 0.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = temp), data = gsv, color = "orange", width = 0.2) +
  # geom_line(aes(x = temperature, y = growth), data = all_preds_average, color = "purple", size = 0.5) +
  # geom_line(aes(x = temperature, y = growth), data = all_preds_NLA, color = "purple", size = 0.5) +
  geom_hline(yintercept = 0) + ylab("Exponential growth rate (per day)") +
  xlab("Temperature (°C)") + coord_cartesian(xlim = c(-3,33), ylim = c(-1, 1.6))
ggsave("Tetraselmis_experiment/figures/nls_boot_figure2.png", width = 5, height = 4)


# now bring in the growth rates -------------------------------------------


growth_sum <- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary.csv") ## empirically observed growth rates
growth_sum_v <- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary_v.csv") ## empirically observed growth rates


gs <- growth_sum %>% 
  mutate(treatment = "constant")

gsv <- growth_sum_v %>% 
  mutate(treatment = "variable")

gs_all <- bind_rows(gs, gsv)





library(viridis)

p + 
  geom_line(aes(x = temperature, y = growth, group = replicate, color = b), data = all_variable, alpha = 1, size = 1.5) +
  scale_color_viridis(discrete = FALSE) +  xlim(-3, 33) +
  ylim(-1, 2) 

