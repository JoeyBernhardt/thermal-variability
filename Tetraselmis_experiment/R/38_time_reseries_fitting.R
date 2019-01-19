
library(cowplot)
library(broom)
library(purrr)
library(tidyverse)
library(minpack.lm)
library(nlstools)
library(viridis)

cells_exp <- read_csv("Tetraselmis_experiment/data-processed/cells_exp.csv")
cells_v_exp <- read_csv("Tetraselmis_experiment/data-processed/cells_v_exp.csv")

cells_days_v <- read_csv("Tetraselmis_experiment/data-processed/cells_days_v_mod.csv") %>% 
  mutate(environment = "variable")
cells_days_c <- read_csv("Tetraselmis_experiment/data-processed/cells_exp_mod.csv") %>% 
  mutate(days = time_since_innoc_hours/24) %>% 
  mutate(days = ifelse(days < 0, 0, days)) %>% 
  mutate(environment = "constant")

cells_days_c %>% 
  ggplot(aes(x = days, y = cell_density)) + geom_point() +
  facet_wrap( ~ temp)

### make a plot of time series

all_cells_plot <- bind_rows(cells_days_c, cells_days_v)

all_cells_plot %>%
  ggplot(aes(x = days, y = cell_density, color = temp)) + geom_point() +
  facet_wrap( ~ environment) + scale_color_viridis() +
  stat_function(fun = growth_function(temp = 5))
  

growth_function <- function(temp, days)(800 * exp((results$a*exp(results$b*temp)*(1-((temp-results$z)/(results$w/2))^2))*(days)))

avals<-seq(-0.2,1.2,0.02)
bvals<-seq(-0.2,0.3,0.02)

avals<-seq(0,0.5,0.07)
bvals<-seq(0,0.16,0.08)
zvals<-seq(10,15,1)
wvals<-seq(34,37,1)


df <-  expand.grid(a = avals, b = bvals, z = zvals, w = wvals) %>% 
  mutate(unique_id = rownames(.)) %>% 
  mutate(sample_size = 1)


##cell_density ~ 800 * exp(r*days)

fit_growth <- function(df){
  res <- try(nlsLM(cell_density ~ 800 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
                   data= cells_days_c,  
                   start=list(z=df$z[[1]],w=df$w[[1]],a= df$a[[1]], b=df$b[[1]]),
                   lower = c(z = 0, w= 0, a = -0.2, b = 0),
                   upper = c(z = 20, w= 80,a =  0.5, b = 0.15),
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

output_v <- df_split %>%
  map_df(fit_growth, .id = "run") 


results_v <- output_v %>% 
  filter(df.residual > 5) %>%
  top_n(n = -1, wt = AIC)
results <- output %>% 
  filter(df.residual > 5) %>%
  top_n(n = -1, wt = AIC)


write_csv(results, "Tetraselmis_experiment/data-processed/time_resampling_fits.csv")
write_csv(results_v, "Tetraselmis_experiment/data-processed/time_resampling_fits_v.csv")


results <- read_csv("Tetraselmis_experiment/data-processed/time_resampling_fits.csv")

fitc <- nlsLM(cell_density ~ 800 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
             data= cells_days_c,  
             start=list(z= results$z[[1]],w= results$w[[1]],a= results$a[[1]], b= results$b[[1]]),
             lower = c(z = 0, w= 0, a = -0.2, b = 0),
             upper = c(z = 20, w= 80,a =  0.5, b = 0.15),
             control = nls.control(maxiter=1024, minFactor=1/204800000))

summary(fitc)
best_fit_c <- coef(fitc)
nlsResiduals(fitc)
nls_boot_c <- nlsBoot(fitc, niter = 999)
nls_boot_coefs_c <- as_data_frame(nls_boot_c$coefboot)
ctpc <- as_data_frame(best_fit_c) %>% 
  rownames_to_column(.) %>% 
  spread(key = rowname, value = value)

write_csv(nls_boot_coefs_c, "Tetraselmis_experiment/data-processed/nls_boot_c.csv")
write_csv(ctpc, "Tetraselmis_experiment/data-processed/ctpc.csv")
nls_boot_coefs_c <- read_csv("Tetraselmis_experiment/data-processed/nls_boot_c.csv")
ctpc <- read_csv("Tetraselmis_experiment/data-processed/ctpc.csv")
# bootstrap variable ---------------------------------------------------------------
library(nlstools)
results_v <- read_csv("Tetraselmis_experiment/data-processed/time_resampling_fits_v.csv")

fit <- nlsLM(cell_density ~ 800 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
                 data= cells_days_v,  
                 start=list(z= 12.14402,w= 36.23217,a= 0.2672314, b=0.08407602),
                 lower = c(z = 0, w= 0, a = -0.2, b = 0),
                 upper = c(z = 20, w= 80,a =  0.5, b = 0.15),
                 control = nls.control(maxiter=1024, minFactor=1/204800000))

best_fit_v <- coef(fit)

nls_boot <- nlsBoot(fit, niter = 999)
nls_boot_coefs <- as_data_frame(nls_boot$coefboot)
vtpc <- as_data_frame(best_fit_v) %>% 
  rownames_to_column(.) %>% 
  spread(key = rowname, value = value)
write_csv(nls_boot_coefs, "Tetraselmis_experiment/data-processed/nls_boot.csv")
write_csv(vtpc, "Tetraselmis_experiment/data-processed/vtpc.csv")

View(nls_boot$coefboot)
nls_boot_coefs_v <- read_csv("Tetraselmis_experiment/data-processed/nls_boot.csv")
vtpc <- read_csv("Tetraselmis_experiment/data-processed/vtpc.csv")
tp <- results_v

grfunc<-function(x){
  -nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]])
}

### find limits
optinfo<-optim(c(x=tp$z[[1]]),grfunc)
opt <-optinfo$par[[1]]
maxgrowth <- -optinfo$value

library(rootSolve)
uniroot.all(function(x) nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]]),c(opt,150))
uniroot.all(function(x) nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]]),c(-15,opt))





# plot both curves --------------------------------------------------------

fc <- read_csv("Tetraselmis_experiment/data-processed/time_resampling_fits.csv")
fv <- read_csv("Tetraselmis_experiment/data-processed/time_resampling_fits_v.csv")

tpc_c <-function(x){
  res<-(fc$a[[1]]*exp(fc$b[[1]]*x)*(1-((x-fc$z[[1]])/(fc$w[[1]]/2))^2))
  res
}

tpc_v <-function(x){
  res<-(fv$a[[1]]*exp(fv$b[[1]]*x)*(1-((x-fv$z[[1]])/(fv$w[[1]]/2))^2))
  res
}

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + 
  stat_function(fun = tpc_c, color = "black", size = 1) +
  xlim(-3, 33) + 
  ylim(0, 2) + geom_hline(yintercept = 0) +
  stat_function(fun = tpc_v, color = "grey", size = 1) + ylab("Growth rate") + xlab("Temperature (°C)") +
  # geom_line(aes(x = temperature, y = growth, group = replicate), color = "cadetblue", data = all_preds, alpha = 0.2) +
  # geom_line(aes(x = temperature, y = growth, group = replicate), color = "purple", data = all_preds_v, alpha = 0.2) +
  # geom_line(aes(x = temperature, y = growth, group = replicate), color = "orange", data = all_preds_NLA, alpha = 0.2) +
  # geom_ribbon(aes(x = temperature, ymin = prediction_lower, ymax = prediction_upper, linetype = NA), fill = "transparent", alpha = 0.01,
              # data = all_preds_average, linetype = "dashed", color = "red", size = 0.5) +
  # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_v, fill = "purple", alpha = 0.5) +
  geom_vline(xintercept = 27) +
  coord_cartesian()
  



# let’s get the NLA prediction on there -----------------------------------


x <- seq(-3, 32, by = 0.01)


predict_tpc <- function(x) {
  y <- 0.5*(ctpc(x + 5) + ctpc(x - 5))
}


predictions <- sapply(x, predict_tpc)
predictions <- data.frame(x, predictions) %>% 
  rename(temperature = x, 
         growth = predictions)




all_preds_average <- predictions



bs <- read_csv("Tetraselmis_experiment/data-processed/bootstrap_time_series_fitsb.csv")
vs <- read_csv("Tetraselmis_experiment/data-processed/bootstrap_time_series_fitsv.csv")

bs2 <- bs %>% 
  select(replicate, z, w, a, b)

vs2 <- vs %>% 
  filter(a != 1.2, w != 45, a!= -0.2, z != 0) %>% 
  filter(b > 0, z > 15) %>% 
  distinct(z, .keep_all = TRUE) %>% 
  select(replicate, z, w, a, b)

# make predictions --------------------------------------------------------

df <- filter(bs2, replicate ==1)



prediction_function <- function(df) {
  tpc <-function(x){
  res<-(df$a[[1]]*exp(df$b[[1]]*x)*(1-((x-df$z[[1]])/(df$w[[1]]/2))^2))
  res
}

pred <- function(x) {
  y <- tpc(x)
}

x <- seq(-3, 32, by = 0.1)

preds <- sapply(x, pred)
preds <- data.frame(x, preds) %>% 
  rename(temperature = x, 
         growth = preds)
}


bs_split <- bs2 %>% 
  split(.$replicate)

vs_split <- vs2 %>% 
  split(.$replicate)

all_preds <- bs_split %>% 
  map_df(prediction_function, .id = "replicate")

all_preds_v <- vs_split %>% 
  map_df(prediction_function, .id = "replicate")

# Make NLA prediction -----------------------------------------------------

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


limits <- all_preds_NLA %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            mean = mean(growth))

limits_v <- all_preds_v %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            mean = mean(growth))
