library(cowplot)
library(tidyverse)

etpc_fits <- read_csv("Tetraselmis_experiment/data-processed/time_resampling_fits.csv")
ctpc_fits <- read_csv("Tetraselmis_experiment/data-processed/ctpc.csv")
vtpc_fits <- read_csv("Tetraselmis_experiment/data-processed/time_resampling_fits_v.csv")
vtpc_fits <- read_csv("Tetraselmis_experiment/data-processed/vtpc.csv") ## new vtpc
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

### ok now let's get the actual predicted (mean) under variable conditions

prediction_NLA <-  prediction_NLA(ctpc_fits)


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
  # geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates_v, color = "black", width = 0.2) +
  # geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = "orange") +
  # geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = "black", shape = 1) +
  geom_hline(yintercept = 0) + ylab("") +
  xlab("") + coord_cartesian(xlim = c(-2,33), ylim = c(-0.1, 1.6))
ggsave("Tetraselmis_experiment/figures/nls_boot_figure2_with_points.png", width = 5, height = 3)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
panel_b_no <- p + 
  geom_line(aes(x = temperature, y = growth, group = replicate), color = "cadetblue", data = all_preds_c, alpha = 0.2) +
  stat_function(fun = ctpc, color = "black", size = 1) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_c, fill = "cadetblue", alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
              fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
  # geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates, color = "black", width = 0.2) +
  # geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = "cadetblue") +
  # geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = "black", shape = 1) +
  geom_hline(yintercept = 0) + ylab("") +
  xlab("") + coord_cartesian(xlim = c(-2,33), ylim = c(-0.1, 1.6))
ggsave("Tetraselmis_experiment/figures/nls_boot_figure2A_with_points.png", width = 5, height = 3)


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
uniroot.all(function(x) nbcurve(x, z = ctpc_fits$z[[1]],w = ctpc_fits$w[[1]],a = ctpc_fits$a[[1]], b = ctpc_fits$b[[1]]),c(10,150))
uniroot.all(function(x) nbcurve(x, z = vtpc_fits$z[[1]],w = vtpc_fits$w[[1]],a = vtpc_fits$a[[1]], b = vtpc_fits$b[[1]]),c(10,150))


## get Tmin
uniroot.all(function(x) nbcurve(x, z = ctpc_fits$z[[1]],w = ctpc_fits$w[[1]],a = ctpc_fits$a[[1]], b = ctpc_fits$b[[1]]),c(-15,10))
uniroot.all(function(x) nbcurve(x, z = vtpc_fits$z[[1]],w = vtpc_fits$w[[1]],a = vtpc_fits$a[[1]], b = vtpc_fits$b[[1]]),c(-15,10))

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
            median = median(tmax))

tmax_summ %>% 
  ggplot(aes(x = treatment, y = mean)) + geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)

## get tmins

get_tmin <- function(df){
  uniroot.all(function(x) nbcurve(x, z = df$z[[1]],w = df$w[[1]],a = df$a[[1]], b = df$b[[1]]),c(-10,5))
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
            median = median(tmin))

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
            median = median(topt))

topts_summ %>% 
  ggplot(aes(x = treatment, y = mean)) + geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)


### now need to get the upper and lower limits predictions

## try with approx fun
lower_line <- approxfun(limits_prediction$temperature, limits_prediction$q2.5)
upper_line <- approxfun(limits_prediction$temperature, limits_prediction$q97.5)
mean_line <- approxfun(limits_prediction$temperature, limits_prediction$mean)
median_line <- approxfun(limits_prediction$temperature, limits_prediction$median)
prediction_line <- approxfun(prediction_NLA$temperature, prediction_NLA$growth)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
 p + 
  # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
              # fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
   stat_function(fun = lower_line, color = "red") +
   stat_function(fun = upper_line, color = "blue") +
   geom_hline(yintercept = 0) + xlim(21.85, 22) + ylim(1.26, 1.295) +
   # geom_vline(xintercept =  tmax_lower, color = "red") +
   geom_vline(xintercept =  21.95, color = "green") +
   geom_vline(xintercept =  21.90, color = "blue") +
   geom_vline(xintercept =  21.88, color = "red") +
   stat_function(fun = mean_line, color = "green")
 

uniroot.all(function(x) nbcurve(x, z = df$z[[1]],w = df$w[[1]],a = df$a[[1]], b = df$b[[1]]),c(-10,5)) 
 (tmax_lower <- uniroot.all(lower_line,c(20,35)))
 (tmax_upper <- uniroot.all(upper_line,c(20,35)))
 (tmax_mean <- uniroot.all(mean_line,c(20,35)))
 
 tmin_upper <- uniroot.all(upper_line,c(-20,10))
 tmin_lower <- uniroot.all(lower_line,c(-20,10))
 tmin_mean <- uniroot.all(mean_line,c(-20,10))
 
 ### topt predicted
 
 lower_line
 grfunc<-function(x){
   -lower_line(x)
 }
   optinfo<-optim(c(x = ctpc_fits$z[[1]]),grfunc)
   opt <-c(optinfo$par[[1]])
   maxgrowth <- c(-optinfo$value)
   (results_lower <- data.frame(topt = opt, rmax = maxgrowth))
results_lower$limit <- "lower"
   
grfunc_upper<-function(x){
     -upper_line(x)
   }
   optinfo<-optim(c(x=ctpc_fits$z[[1]]),grfunc_upper)
   opt <-c(optinfo$par[[1]])
   maxgrowth <- c(-optinfo$value)
   (results_upper <- data.frame(topt = opt, rmax = maxgrowth))
   results_upper$limit <- "upper"
   
   
   grfunc_mean<-function(x){
     -mean_line(x)
   }
   
  optinfo<-optim(c(ctpc_fits$z[[1]]),grfunc_mean)
   opt <-c(optinfo$par[[1]])
   maxgrowth <- c(-optinfo$value)
   (results_mean <- data.frame(topt = opt, rmax = maxgrowth))
   results_mean$limit <- "mean"   
   
   
   grfunc_median<-function(x){
     -median_line(x)
   }
   
   optinfo<-optim(c(ctpc_fits$z[[1]]),grfunc_median)
   opt <-c(optinfo$par[[1]])
   maxgrowth <- c(-optinfo$value)
   (results_median <- data.frame(topt = opt, rmax = maxgrowth))
   results_median$limit <- "median"   
   
   grfunc_prediction<-function(x){
     -prediction_line(x)
   }
   
   optinfo<-optim(c(ctpc_fits$z[[1]]),grfunc_prediction)
   opt <-c(optinfo$par[[1]])
   maxgrowth <- c(-optinfo$value)
   (results_prediction <- data.frame(topt = opt, rmax = maxgrowth))
   results_prediction$limit <- "prediction"   
   
   
   ### ok let's put all the predictions and observations together
  
   topt_predictions <- bind_rows(results_lower, results_upper, results_mean, results_median, results_prediction) 
   
   
   topt_predictions %>% 
     ggplot(aes(x = limit, y = topt)) + geom_point()
   
   limits_prediction %>% 
     top_n(mean, n = 1) %>% View
   
   limits_prediction %>% 
     top_n(q2.5, n = 1) %>% View
   
   limits_prediction %>% 
     gather(key = limit, value = value, 2:6) %>% 
     group_by(limit) %>% 
     top_n(value, n = 1) %>% View
   
   
   ### new option: find the first derivative, and find where that function crosses 0

   dtpc <-function(x) ctpc_fits$a[[1]]*exp(ctpc_fits$b[[1]]*x)*(1-((x-ctpc_fits$z[[1]])/(ctpc_fits$w[[1]]/2))^2)

   deriv((y ~ sin(cos(x) * y)), c("x","y"), func = TRUE)
   
   pred <- function(x) 0.5*(ctpc(x + 5) + ctpc(x - 5))

   # fderiv <- approxfun(deriv(0.5*(ctpc_fits$a[[1]]*exp(ctpc_fits$b[[1]]*x+5)*(1-((x+5-ctpc_fits$z[[1]])/(ctpc_fits$w[[1]]/2))^2) + 
                # ctpc_fits$a[[1]]*exp(ctpc_fits$b[[1]]*x-5)*(1-((x-5-ctpc_fits$z[[1]])/(ctpc_fits$w[[1]]/2))^2)), c("x")))
   
   fderiv <- deriv(0.5*(dtpc(x + 5) + dtpc(x - 5)), c("x"), func = TRUE)
   
   uniroot(fderiv, c(1, 40))

   
   p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
   p + 
     # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
     # fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
     stat_function(fun = fderiv, color = "red") + xlim(0, 30) + ylim(0, 1)
   
   
   
   