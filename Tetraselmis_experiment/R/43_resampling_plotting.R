
library(tidyverse)
library(cowplot)

growth_sum <- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary.csv") ## empirically observed growth rates
growth_sum_v <- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary_v.csv") ## empirically observed growth rates


gs <- growth_sum %>% 
  mutate(treatment = "constant")

gsv <- growth_sum_v %>% 
  mutate(treatment = "variable")

gs_all <- bind_rows(gs, gsv)


gs_all %>% 
  ggplot(aes(x = temp, y = mean, color = treatment)) + geom_point()




fits_variable <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v_exp.csv")

fits_c <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_1000_exp.csv")

fits_c <- fits_c %>%
  rename(a = a.list,
         b = b.list,
         w = w.list,
         z = z.list)

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


fits_c_split <- fits_c %>%
  # filter(z != 30) %>% 
  # filter(b > 0, a != 1, z != 20, z != 0, b != 0) %>%
  # filter(b < 0.12) %>% 
  split(.$curve.id.list)

all_preds_c <- fits_c_split %>% 
  map_df(prediction_function, .id = "replicate")

v_fits <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_v_exp.csv")

tpc_v <-function(x){
  res<-v_fits$a.list[[1]]*exp(v_fits$b.list[[1]]*x)*(1-((x-v_fits$z.list[[1]])/(v_fits$w.list[[1]]/2))^2)
  res
}

tpc_c <-function(x){
  res<-fits_c$a.list[[1]]*exp(fits_c$b.list[[1]]*x)*(1-((x-fits_c$z.list[[1]])/(fits_c$w.list[[1]]/2))^2)
  res
}

limits_c <- all_preds_c %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            mean = mean(growth))


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + 
  # geom_line(aes(x = temperature, y = growth, group = replicate), data = all_preds_c, alpha = 0.2) +
  stat_function(fun = tpc_c, color = "purple") +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_c, fill = "orange", alpha = 0.5) +
  geom_point(aes(x = temp, y = mean), data = gs) +
 coord_cartesian(xlim = c(-2, 32), ylim = c(-0.5, 2))
