library(cowplot)
library(tidyverse)
library(purrr)
library(rootSolve)

# etpc_fits <- read_csv("Tetraselmis_experiment/data-processed/time_resampling_fits.csv")
ctpc_fits <- read_csv("Tetraselmis_experiment/data-processed/ctpc.csv")
# vtpc_fits <- read_csv("Tetraselmis_experiment/data-processed/time_resampling_fits_v.csv")
vtpc_fits <- read_csv("Tetraselmis_experiment/data-processed/vtpc.csv") ## new vtpc
# c_fits <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_exp.csv")


# 
# c_fits <- c_fits %>% 
#   rename(z = z.list,
#          a = a.list, 
#          b = b.list, 
#          w = w.list)
# 
# 
# 
etpc <-function(x){
  res<-(etpc_fits$a[[1]]*exp(etpc_fits$b[[1]]*x)*(1-((x-etpc_fits$z[[1]])/(etpc_fits$w[[1]]/2))^2))
  res
}

ctpc <-function(x){
  res<-(ctpc_fits$a[[1]]*exp(ctpc_fits$b[[1]]*x)*(1-((x-ctpc_fits$z[[1]])/(ctpc_fits$w[[1]]/2))^2))
  res
}



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
  # stat_function(fun = etpc, color = "black", size = 1) +
  stat_function(fun = ctpc, color = "black", size = 1) +
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

# bs_v <- read_csv("Tetraselmis_experiment/data-processed/bootstrap_time_series_fitsv.csv")
# bs_v <- read_csv("Tetraselmis_experiment/data-processed/bootstrap_time_series_fitsv_bounds.csv")

## bring in bootstrap TPCs
bs_v <- read_csv("Tetraselmis_experiment/data-processed/nls_boot.csv") %>% 
  mutate(replicate = rownames(.))

bs_c <- read_csv("Tetraselmis_experiment/data-processed/nls_boot_c.csv") %>% 
  mutate(replicate = rownames(.))



# bs_c1 <- read_csv("Tetraselmis_experiment/data-processed/bootstrap_time_series_fits.csv")

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


write_csv(limits_c, "Tetraselmis_experiment/data-processed/limits_c_direct.csv")
bs_v2 <- bs_v %>% 
  mutate(replicate = as.character(replicate))

all_variable <- left_join(all_preds_v, bs_v2, by = "replicate")

# NLA predictions ---------------------------------------------------------

prediction_NLA <- function(df) {
  tpc <-function(x){
    res<-(df$a[[1]]*exp(df$b[[1]]*x)*(1-((x-df$z[[1]])/(df$w[[1]]/2))^2))
    res
  }
  
  pred <- function(x) {
    y <-  0.5*(tpc(x + 5) + tpc(x - 5))
  }
  
  x <- seq(-20, 32, by = 0.01)
  
  preds <- sapply(x, pred)
  preds <- data.frame(x, preds) %>% 
    rename(temperature = x, 
           growth = preds)
}


### this contains all the predictions using NLA
all_preds_NLA <- bs_split %>% 
  map_df(prediction_NLA, .id = "replicate")



limits_prediction <- all_preds_NLA %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            q50=quantile(growth, probs=0.5),
            mean = mean(growth),
            median = median(growth))

write_csv(limits_prediction, "Tetraselmis_experiment/data-processed/limits_prediction.csv")

limits_prediction <- read_csv("Tetraselmis_experiment/data-processed/limits_prediction.csv")

### update April 17 2018, ok now let's get the NLA predictions for the variable treatment
## but based on the constant curve derived from the indirect approach

fits_c <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_5000_exp.csv")
fits_c <- fits_c %>%
  rename(a = a.list,
         b = b.list,
         w = w.list,
         z = z.list) %>% 
  filter(maxgrowth.list < 2)


fits_c_split <- fits_c %>% 
  split(.$curve.id.list)


all_preds_NLA_indirect <- fits_c_split %>% 
  map_df(prediction_NLA, .id = "replicate")



limits_prediction_indirect <- all_preds_NLA_indirect %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            q50=quantile(growth, probs=0.5),
            mean = mean(growth),
            median = median(growth))
write_csv(limits_prediction_indirect, "Tetraselmis_experiment/data-processed/limits_prediction_indirect.csv")

### ok now let's get the actual predicted (mean) under variable conditions

prediction_NLA <-  prediction_NLA(ctpc_fits)


### here's the scale transition theory prediction
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


prediction_STT <- function(df) {
  tpc <-function(x){
    res<-(df$a[[1]]*exp(df$b[[1]]*x)*(1-((x-df$z[[1]])/(df$w[[1]]/2))^2))
    res
  }
  
  pred <- function(x) {
    y <-  y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*25
  }
  
  x <- seq(-2, 32, by = 0.01)
  
  preds <- sapply(x, pred)
  preds <- data.frame(x, preds) %>% 
    rename(temperature = x, 
           growth = preds)
}


all_preds_STT <- bs_split %>% 
  map_df(prediction_STT, .id = "replicate")



limits_prediction_STT <- all_preds_STT %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            q50=quantile(growth, probs=0.5),
            mean = mean(growth),
            median = median(growth))

write_csv(limits_prediction_STT, "Tetraselmis_experiment/data-processed/limits_prediction_STT.csv")
# now bring in the growth rates -------------------------------------------


growth_sum <- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary.csv") ## empirically observed growth rates
growth_sum_v <- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary_v.csv") ## empirically observed growth rates


gs <- growth_sum %>% 
  mutate(treatment = "constant")

gsv <- growth_sum_v %>% 
  mutate(treatment = "variable")

gs_all <- bind_rows(gs, gsv)



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
  geom_errorbar(aes(ymin = lower, ymax = upper, x = temp), data = gsv, color = "red", width = 0.2) +
  # geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, x = temp), data = gsv, color = "orange", width = 0.2) +
  # geom_line(aes(x = temperature, y = growth), data = all_preds_average, color = "purple", size = 0.5) +
  # geom_line(aes(x = temperature, y = growth), data = all_preds_NLA, color = "purple", size = 0.5) +
  geom_hline(yintercept = 0) + ylab("Exponential growth rate (per day)") +
  xlab("Temperature (°C)") + coord_cartesian(xlim = c(-3,33), ylim = c(-1, 1.6))
ggsave("Tetraselmis_experiment/figures/nls_boot_figure2.png", width = 5, height = 4)



# plot this again ---------------------------------------------------------

all_estimates <- read_csv("Tetraselmis_experiment/data-processed/growth_estimates_boot_car.csv")
all_estimates_v <- read_csv("Tetraselmis_experiment/data-processed/growth_estimates_boot_car_v.csv")

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
panel_a_np <- p + 
  geom_line(aes(x = temperature, y = growth, group = replicate), color = "cadetblue", data = all_preds_c, alpha = 0.2) +
  stat_function(fun = ctpc, color = "black", size = 1) +
  stat_function(fun = vtpc, color = "orange", size = 1) +
   geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_v, fill = "orange", alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_c, fill = "cadetblue", alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
              fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates_v, color = "black", width = 0.2) +
  geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = "orange") +
  geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = "black", shape = 1) +
  geom_hline(yintercept = 0) + ylab("") +
  xlab("") + coord_cartesian(xlim = c(-2,33), ylim = c(-0.1, 1.6))
ggsave("Tetraselmis_experiment/figures/nls_boot_figure2_with_points.png", width = 5, height = 3)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
panel_b_no <- p + 
  # geom_line(aes(x = temperature, y = growth, group = replicate), color = "cadetblue", data = all_preds_c, alpha = 0.2) +
  # stat_function(fun = ctpc, color = "black", size = 1) +
  # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_c, fill = "cadetblue", alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
              fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
  # geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates, color = "black", width = 0.2) +
  # geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = "cadetblue") +
  # geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = "black", shape = 1) +
  geom_hline(yintercept = 0) + ylab("") +
  xlab("") + coord_cartesian(xlim = c(-2,33), ylim = c(-0.1, 1.6))
ggsave("Tetraselmis_experiment/figures/nls_boot_figure2A_with_points.png", width = 5, height = 3)
ggsave("Tetraselmis_experiment/figures/predictions.png", width = 5, height = 3)

## plot together

figure2_vasseur <- plot_grid(panel_b_no, panel_a_np, labels = c("A", "B"), align = "v", ncol = 1)
save_plot("Tetraselmis_experiment/figures/figure2_vasseur_no_points.png", figure2_vasseur, ncol = 1, base_height = 7, base_width = 6.2)


library(viridis)

p + 
  geom_line(aes(x = temperature, y = growth, group = replicate, color = b), data = all_variable, alpha = 1, size = 1.5) +
  scale_color_viridis(discrete = FALSE) +  xlim(-3, 33) +
  ylim(-1, 2) 


# get the roots and topt --------------------------------------------------


library(rootSolve)
nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

## get Tmax
tmax_c <- uniroot.all(function(x) nbcurve(x, z = ctpc_fits$z[[1]],w = ctpc_fits$w[[1]],a = ctpc_fits$a[[1]], b = ctpc_fits$b[[1]]),c(10,150))
tmax_v <- uniroot.all(function(x) nbcurve(x, z = vtpc_fits$z[[1]],w = vtpc_fits$w[[1]],a = vtpc_fits$a[[1]], b = vtpc_fits$b[[1]]),c(10,150))


## get Tmin
tmin_c <- uniroot.all(function(x) nbcurve(x, z = ctpc_fits$z[[1]],w = ctpc_fits$w[[1]],a = ctpc_fits$a[[1]], b = ctpc_fits$b[[1]]),c(-15,10))
tmin_v <- uniroot.all(function(x) nbcurve(x, z = vtpc_fits$z[[1]],w = vtpc_fits$w[[1]],a = vtpc_fits$a[[1]], b = vtpc_fits$b[[1]]),c(-15,10))

## get the upper and lower bounds on Tmax from the nls bootstrapping 

bs_split <- bs_c %>% 
  split(.$replicate)

bs_v_split <- bs_v %>%
  split(.$replicate)

get_tmax <- function(df){
  uniroot.all(function(x) nbcurve(x, z = df$z[[1]],w = df$w[[1]],a = df$a[[1]], b = df$b[[1]]),c(10,150))
}

tmaxes_c <- bs_split %>% 
  map(get_tmax) %>% 
  unlist() %>% 
  as_data_frame() %>% 
  rename(tmax = value) %>% 
  mutate(replicate = rownames(.)) %>% 
  mutate(treatment = "constant")

tmaxes_v <- bs_v_split %>% 
  map(get_tmax) %>% 
  unlist() %>% 
  as_data_frame() %>% 
  rename(tmax = value) %>% 
  mutate(replicate = rownames(.)) %>% 
  mutate(treatment = "variable")

all_tmaxes <- bind_rows(tmaxes_c, tmaxes_v)

tmax_summ <- all_tmaxes %>% 
  group_by(treatment) %>% 
  summarise(lower=quantile(tmax, probs=0.025),
            upper=quantile(tmax, probs=0.975),
            mean = mean(tmax),
            median = median(tmax)) %>% 
  mutate(param = "tmax")

tmax_summ %>% 
  ggplot(aes(x = treatment, y = mean)) + geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)

## get tmins

get_tmin <- function(df){
  uniroot.all(function(x) nbcurve(x, z = df$z[[1]],w = df$w[[1]],a = df$a[[1]], b = df$b[[1]]),c(-40,5))
}

tmins_c <- bs_split %>% 
  map(get_tmin) %>% 
  unlist() %>% 
  as_data_frame() %>% 
  rename(tmin = value) %>% 
  mutate(replicate = rownames(.)) %>% 
  mutate(treatment = "constant")

tmins_v <- bs_v_split %>% 
  map(get_tmin) %>% 
  unlist() %>% 
  as_data_frame() %>% 
  rename(tmin = value) %>% 
  mutate(replicate = rownames(.)) %>% 
  mutate(treatment = "variable")

all_tmins <- bind_rows(tmins_c, tmins_v)

tmin_summ <- all_tmins %>% 
  group_by(treatment) %>% 
  summarise(lower=quantile(tmin, probs=0.025),
            upper=quantile(tmin, probs=0.975),
            mean = mean(tmin),
            median = median(tmin)) %>% 
  mutate(param = "tmin")

tmin_summ %>% 
  ggplot(aes(x = treatment, y = mean)) + geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)


### get topts

get_topt <- function(df){
grfunc<-function(x){
  -nbcurve(x, z = df$z[[1]],w = df$w[[1]],a = df$a[[1]],b = df$b[[1]])
}
optinfo<-optim(c(x=df$z[[1]]),grfunc)
opt <-c(optinfo$par[[1]])
maxgrowth <- c(-optinfo$value)
results <- data.frame(topt = opt, rmax = maxgrowth)
return(results)
}

topts_c <- bs_split %>% 
  map_df(get_topt, .id = "replicate") %>% 
  mutate(treatment = "constant")

topts_v <- bs_v_split %>% 
  map_df(get_topt, .id = "replicate") %>% 
  mutate(treatment = "variable")

all_topts <- bind_rows(topts_c, topts_v)

topts_summ <- all_topts %>% 
  group_by(treatment) %>% 
  summarise(lower=quantile(topt, probs=0.025),
            upper=quantile(topt, probs=0.975),
            mean = mean(topt),
            median = median(topt)) %>% 
  mutate(param = "topt")

rmax_summ <- all_topts %>% 
  group_by(treatment) %>% 
  summarise(lower=quantile(rmax, probs=0.025),
            upper=quantile(rmax, probs=0.975),
            mean = mean(rmax),
            median = median(rmax)) %>% 
  mutate(param = "rmax")


topts_summ %>% 
  ggplot(aes(x = treatment, y = mean)) + geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)





### now need to get the upper and lower limits predictions

## try with approx fun
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
 p + 
  # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
               # fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
   stat_function(fun = lower_line, color = "red") +
   # stat_function(fun = dtpc, color = "pink") +
   stat_function(fun = ftpc, color = "orange") +
   stat_function(fun = upper_line, color = "blue") +
   geom_hline(yintercept = 0) + xlim(0, 32) + ylim(0, 1.6) +
   # geom_vline(xintercept =  tmax_lower, color = "red") +
   # geom_vline(xintercept =  21.95, color = "green") +
   # geom_vline(xintercept =  21.90, color = "blue") +
   # geom_vline(xintercept =  21.88, color = "red") +
   stat_function(fun = mean_line, color = "green")
 ggsave("Tetraselmis_experiment/figures/prediction_nla_colors.png")


   
   
   ### new option: find the first derivative, and find where that function crosses 0

   dtpc <-function(x) ctpc_fits$a[[1]]*exp(ctpc_fits$b[[1]]*x)*(1-((x-ctpc_fits$z[[1]])/(ctpc_fits$w[[1]]/2))^2)
   ftpc <- function(x) 0.5*(dtpc(x + 5) + dtpc(x - 5))

   #### ok trying this again!!
   
   bs_c <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_5000_exp.csv") %>% 
     mutate(replicate = rownames(.))
   fits_c <- bs_c %>%
     rename(a = a.list,
            b = b.list,
            w = w.list,
            z = z.list)
   
   bs_split <- fits_c %>%
     filter(maxgrowth.list < 2) %>% 
     split(.$curve.id.list)
   
   get_topts <- function(df){
   
  dtpc <-function(x) df$a[[1]]*exp(df$b[[1]]*x)*(1-((x-df$z[[1]])/(df$w[[1]]/2))^2)
   ftpc <- function(x) 0.5*(dtpc(x + 5) + dtpc(x - 5))
   
   grfunc_ftpc<-function(x){
     -ftpc(x)
   }
   
   optinfo<-optim(c(df$z[[1]]),grfunc_ftpc)
   opt <-c(optinfo$par[[1]])
   maxgrowth <- c(-optinfo$value)
   results_ftpc <- data.frame(topt = opt, rmax = maxgrowth)
   return(results_ftpc)
   }
   
   
   topts_predicted <- bs_split %>% 
     map_df(get_topts, .id = "replicate") 
   
   topt_lims <- topts_predicted %>% 
     filter(topt < 30) %>% 
     summarise(lower=quantile(topt, probs=0.025),
               upper=quantile(topt, probs=0.975),
               mean = mean(topt),
               median = median(topt)) %>% 
     # gather(key = limit_type, value = value) %>% 
     mutate(param = "predicted_topt") %>% 
     mutate(treatment = "variable")
   
   rmax_lims <- topts_predicted %>% 
     filter(topt < 30) %>% 
     summarise(lower=quantile(rmax, probs=0.025),
               upper=quantile(rmax, probs=0.975),
               mean = mean(rmax),
               median = median(rmax)) %>% 
     # gather(key = limit_type, value = value) %>% 
     mutate(param = "predicted_rmax") %>% 
     mutate(treatment = "variable")
   

   
  # get predicted tmax ------------------------------------------------------

   get_tmaxes <- function(df){
     
     dtpc <-function(x) df$a[[1]]*exp(df$b[[1]]*x)*(1-((x-df$z[[1]])/(df$w[[1]]/2))^2)
     ftpc <- function(x) 0.5*(dtpc(x + 5) + dtpc(x - 5))
     
     tmax <- uniroot.all(ftpc,c(10,150))
     return(data.frame(tmax))
   }
   
   
   tmax_predicted <- bs_split %>% 
     map_df(get_tmaxes, .id = "replicate") 
   
   tmax_lims <- tmax_predicted %>% 
     summarise(lower=quantile(tmax, probs=0.025),
               upper=quantile(tmax, probs=0.975),
               mean = mean(tmax),
               median = median(tmax)) %>% 
     # gather(key = limit_type, value = value) %>% 
     mutate(param = "predicted_tmax") %>% 
     mutate(treatment = "variable")
   
   # get predicted tmin ------------------------------------------------------
   
   get_tmins <- function(df){
     
     dtpc <-function(x) df$a[[1]]*exp(df$b[[1]]*x)*(1-((x-df$z[[1]])/(df$w[[1]]/2))^2)
     ftpc <- function(x) 0.5*(dtpc(x + 5) + dtpc(x - 5))
     
     tmin <- uniroot.all(ftpc,c(-70,10))
     return(data.frame(tmin))
   }
   
   tmin_predicted <- bs_split %>% 
     map_df(get_tmins, .id = "replicate") 
   
   tmin_lims <- tmin_predicted %>% 
     # mutate(tmin = ifelse(tmin < 1.8, 1.8, tmin)) %>% 
     summarise(lower=quantile(tmin, probs=0.025),
               upper=quantile(tmin, probs=0.975),
               mean = mean(tmin),
               median = median(tmin)) %>% 
     # gather(key = limit_type, value = value) %>% 
     mutate(param = "predicted_tmin") %>% 
     mutate(treatment = "variable")
   
  
   
 
   
   crit_temps_c <- left_join(tmins_c, tmaxes_c) %>% 
     mutate(tmin = ifelse(tmin <-1.8, -1.8, tmin)) %>% 
     mutate(width = tmax - tmin)
   
   
   width_summ_c <- crit_temps_c %>% 
     summarise(lower=quantile(width, probs=0.025),
               upper=quantile(width, probs=0.975),
               mean = mean(width),
               median = median(width)) %>% 
     # gather(key = limit_type, value = value) %>% 
     mutate(param = "width") %>% 
     mutate(treatment = "constant")
   
   
   crit_temps_v <- left_join(tmins_v, tmaxes_v) %>% 
     mutate(tmin = ifelse(tmin <-1.8, -1.8, tmin)) %>% 
     mutate(width = tmax - tmin)
   
   
   width_summ_v <- crit_temps_v %>% 
     summarise(lower=quantile(width, probs=0.025),
               upper=quantile(width, probs=0.975),
               mean = mean(width),
               median = median(width)) %>% 
     # gather(key = limit_type, value = value) %>% 
     mutate(param = "width") %>% 
     mutate(treatment = "variable")
   
  
     
   
   
   all_critical_temps %>% 
     ggplot(aes(x = param, y = mean, color = treatment)) + geom_point() +
     geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
     facet_wrap( ~ param, scales = "free")
   
   
   all_critical_temps %>% 
     filter(grepl("tmin", param), treatment == "variable") %>% 
     ggplot(aes(x = param, y = mean, color = treatment)) + geom_point() +
     geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)
   
   ### find the limits on w
   
   w_summ_c <- bs_c %>% 
     summarise(lower=quantile(w, probs=0.025),
               upper=quantile(w, probs=0.975),
               mean = mean(w),
               median = median(w)) %>% 
     # gather(key = limit_type, value = value) %>% 
     mutate(param = "w") %>% 
     mutate(treatment = "constant")
   
   w_summ_v <- bs_v %>% 
     summarise(lower=quantile(w, probs=0.025),
               upper=quantile(w, probs=0.975),
               mean = mean(w),
               median = median(w)) %>% 
     # gather(key = limit_type, value = value) %>% 
     mutate(param = "w") %>% 
     mutate(treatment = "variable")
   
   
   all_critical_temps <- bind_rows(w_summ_c, w_summ_v, topts_summ, tmin_summ, tmax_summ, topt_lims, tmin_lims, tmax_lims, rmax_summ, rmax_lims, width_summ_c, width_summ_v)
   
   all_critical_temps %>% 
     mutate_at(.funs = round, digits = 2, .vars = 1:4) %>% View
   
   write_csv(all_critical_temps, "Tetraselmis_experiment/data-processed/all_critical_temps.csv")
   
   
   
   
# big plot with critical temps --------------------------------------------

   crit_temps <- all_critical_temps %>% 
     filter(param %in% c("topt", "tmax"))
   
   crit_temps_predicted <- all_critical_temps %>% 
     filter(param %in% c("predicted_topt", "predicted_tmax"))
   
   p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
   panel_a <- 
    p + 
     # stat_function(fun = ctpc, color = "black", size = 1) +
     # stat_function(fun = vtpc, color = "orange", size = 1) +
     # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_v, fill = "orange", alpha = 0.5) +
     # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_c, fill = "cadetblue", alpha = 0.7) +
     geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
                 fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
     # geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates_v, color = "black", width = 0.2) +
     # geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = "orange") +
     # geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = "black", shape = 1) +
     # # geom_point(aes(x = mean, y = 1.58), data = filter(crit_temps, treatment == "constant"), color = "cadetblue", shape = 17) +
     # geom_point(aes(x = mean, y = 1.65), data = filter(all_critical_temps, param == "predicted_topt"), color = "black", shape = 17) +
     # geom_errorbarh(aes(xmin = lower, xmax = upper, y = 1.65, x = mean), data = filter(all_critical_temps, param == "predicted_topt"), height = 0.1, color = "black", shape = 17) +
     # geom_point(aes(x = mean, y = 1.65), data = filter(all_critical_temps, param == "predicted_tmax"), color = "black", shape = 17) +
     # geom_errorbarh(aes(xmin = lower, xmax = upper, y = 1.65, x = mean), data = filter(all_critical_temps, param == "predicted_tmax"), height = 0.1, color = "black", shape = 17) +
     # geom_point(aes(x = mean, y = 1.58), data = filter(crit_temps, treatment == "variable"), color = "orange", shape = 17) +
     # geom_point(aes(x = -1.5, y = mean), data = filter(all_critical_temps, param == "predicted_rmax"), color = "black", shape = 17) +
     # geom_errorbar(aes(ymin = lower, ymax = upper, x = -1.5), data = filter(all_critical_temps, param == "predicted_rmax"), height = 0.1, color = "black", shape = 17) +
     # geom_point(aes(x = -1.5, y = mean), data = filter(all_critical_temps, param == "rmax", treatment == "variable"), color = "orange", shape = 17) +
     # geom_errorbar(aes(ymin = lower, ymax = upper, x =-1.5), data = filter(all_critical_temps, param == "rmax", treatment == "variable"), height = 0.1, color = "orange", shape = 17) +
     # geom_point(aes(x = -1.5, y = mean), data = filter(all_critical_temps, param == "rmax", treatment == "constant"), color = "cadetblue", shape = 17) +
     # geom_errorbar(aes(ymin = lower, ymax = upper, x =-1.5), data = filter(all_critical_temps, param == "rmax", treatment == "constant"), height = 0.1, color = "cadetblue", shape = 17) +
     # geom_errorbarh(aes(xmin = lower, xmax = upper, y = 1.58, x = mean), data = filter(crit_temps, treatment == "constant"), height = 0.1, color = "cadetblue", shape = 17) +
     # geom_errorbarh(aes(xmin = lower, xmax = upper, y = 1.58, x = mean), data = filter(crit_temps, treatment == "variable"), height = 0.1, color = "orange", shape = 17) +
     # geom_point(aes(x = mean, y = 1.58), data = crit_temps, color = "black", shape = 2) +
     geom_hline(yintercept = 0) + ylab("") +
     xlab("") + coord_cartesian(xlim = c(-2,33), ylim = c(-0.1, 1.6)) 
   ggsave("Tetraselmis_experiment/figures/nls_boot_figure2_with_points.png", width = 6.5, height = 4)
   ggsave("Tetraselmis_experiment/figures/nls_boot_variable.png", width = 6.5, height = 4)
   ggsave("Tetraselmis_experiment/figures/variable_prediction.png", width = 6.5, height = 4)
   p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
   panel_b_no <- p + 
     # geom_line(aes(x = temperature, y = growth, group = replicate), color = "cadetblue", data = all_preds_c, alpha = 0.2) +
     stat_function(fun = ctpc, color = "black", size = 1) +
     geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_c, fill = "cadetblue", alpha = 0.5) +
     geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
                 fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
     geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates, color = "black", width = 0.2) +
     geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = "cadetblue") +
     geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = "black", shape = 1) +
     # geom_errorbarh(aes(xmin = lower, xmax = upper, y = 1.64, x = mean), data = filter(crit_temps, treatment == "constant"), height = 0.1, color = "cadetblue") +
     # geom_point(aes(x = mean, y = 1.64), data = filter(crit_temps, treatment == "constant"), color = "cadetblue", shape = 17) +
     # geom_point(aes(x = mean, y = 1.64), data = filter(crit_temps, treatment == "constant"), color = "black", shape = 2) +
     geom_hline(yintercept = 0) + ylab("") +
     xlab("") + coord_cartesian(xlim = c(-2,33), ylim = c(-0.1, 1.6)) 
   ggsave("Tetraselmis_experiment/figures/nls_boot_constant.png", width = 6.5, height = 4)
   ggsave("Tetraselmis_experiment/figures/nls_boot_constant_w_pred.png", width = 6.5, height = 4)
   figure2_direct <- plot_grid(panel_b_no, panel_a, labels = c("A", "B"), align = "v", ncol = 1)
   save_plot("Tetraselmis_experiment/figures/figure2_direct.png", figure2_direct, ncol = 1, base_height = 7, base_width = 6.2)
   
   
   constant_lims <- p +
     geom_errorbarh(aes(xmin = lower, xmax = upper, y = 0, x = mean), height = 0.05, data = filter(crit_temps, treatment == "constant"), height = 0.05, color = "cadetblue") +
     geom_point(aes(x = mean, y = 0), data = filter(crit_temps, treatment == "constant"), color = "cadetblue", shape = 17) +
     geom_point(aes(x = mean, y = 0), data = filter(crit_temps, treatment == "constant"), color = "black", shape = 2) +
   geom_errorbarh(aes(xmin = lower, xmax = upper, y = 0, x = mean), data = filter(crit_temps, treatment == "constant"), height = 0.05, color = "cadetblue") +
   geom_point(aes(x = mean, y = 0), data = filter(crit_temps, treatment == "constant"), color = "cadetblue", shape = 17) +
   geom_point(aes(x = mean, y = 0), data = filter(crit_temps, treatment == "constant"), color = "black", shape = 2) +
   geom_point(aes(x = mean, y = 0), data = filter(crit_temps, treatment == "constant"), color = "cadetblue", shape = 17) +
   geom_point(aes(x = mean, y = 0.05), data = filter(all_critical_temps, param == "predicted_topt"), color = "black", shape = 17) +
   geom_errorbarh(aes(xmin = lower, xmax = upper, y = 0.05, x = mean), data = filter(all_critical_temps, param == "predicted_topt"), height = 0.05, color = "black", shape = 17) +
   geom_point(aes(x = mean, y = 0.05), data = filter(all_critical_temps, param == "predicted_tmax"), color = "black", shape = 17) +
   geom_errorbarh(aes(xmin = lower, xmax = upper, y = 0.05, x = mean), data = filter(all_critical_temps, param == "predicted_tmax"), height = 0.05, color = "black", shape = 17) +
   geom_point(aes(x = mean, y = 0), data = filter(crit_temps, treatment == "variable"), color = "orange", shape = 17) +
   # geom_point(aes(x = -1.5, y = mean), data = filter(all_critical_temps, param == "predicted_rmax"), color = "black", shape = 17) +
   # geom_errorbar(aes(ymin = lower, ymax = upper, x = -1.5), data = filter(all_critical_temps, param == "predicted_rmax"), height = 0.1, color = "black", shape = 17) +
   # geom_point(aes(x = -1.5, y = mean), data = filter(all_critical_temps, param == "rmax", treatment == "variable"), color = "orange", shape = 17) +
   # geom_errorbar(aes(ymin = lower, ymax = upper, x =-1.5), data = filter(all_critical_temps, param == "rmax", treatment == "variable"), height = 0.1, color = "orange", shape = 17) +
   # geom_point(aes(x = -1.5, y = mean), data = filter(all_critical_temps, param == "rmax", treatment == "constant"), color = "cadetblue", shape = 17) +
   # geom_errorbar(aes(ymin = lower, ymax = upper, x =-1.5), data = filter(all_critical_temps, param == "rmax", treatment == "constant"), height = 0.1, color = "cadetblue", shape = 17) +
   geom_errorbarh(aes(xmin = lower, xmax = upper, y = 0, x = mean), data = filter(crit_temps, treatment == "constant"), height = 0.05, color = "cadetblue", shape = 17) +
   geom_errorbarh(aes(xmin = lower, xmax = upper, y = 0, x = mean), data = filter(crit_temps, treatment == "variable"), height = 0.05, color = "orange", shape = 17) +
   geom_point(aes(x = mean, y = 0), data = crit_temps, color = "black", shape = 2) +
     coord_cartesian(xlim = c(-2,33), ylim = c(0.05, 1.66)) + ylab("") + xlab("")
   
   figure2_limits <- plot_grid(constant_lims, labels = c("A"), align = "v", ncol = 1)
     save_plot("Tetraselmis_experiment/figures/figure2_limits.png", figure2_limits, ncol = 1, base_height = 7, base_width = 6.2)
   
   