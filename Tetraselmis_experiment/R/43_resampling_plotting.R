
library(tidyverse)
library(cowplot)
library(extrafont)
loadfonts()
growth_sum <- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary.csv") ## empirically observed growth rates
growth_sum_v <- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary_v.csv") ## empirically observed growth rates


gs <- growth_sum %>% 
  mutate(treatment = "constant")

gsv <- growth_sum_v %>% 
  mutate(treatment = "variable")

gs_all <- bind_rows(gs, gsv)


gs_all %>% 
  ggplot(aes(x = temp, y = mean, color = treatment)) + geom_point()






fits_c <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_5000_exp.csv")
# fits_c <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_1000_exp.csv")
fits_v <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v_exp3.csv")
fits_v_1000 <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v_exp2.csv")
fits_c <- fits_c %>%
  rename(a = a.list,
         b = b.list,
         w = w.list,
         z = z.list)

fits_v <- fits_v %>%
  rename(a = a.list,
         b = b.list,
         w = w.list,
         z = z.list)

fits_v_1000 <- fits_v_1000 %>%
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

fits_v_split <- fits_v %>%
  split(.$curve.id.list)

fits_v_split_1000 <- fits_v_1000 %>%
  split(.$curve.id.list)

all_preds_v <- fits_v_split %>% 
  map_df(prediction_function, .id = "replicate")

all_preds_v_1000 <- fits_v_split_1000 %>% 
  map_df(prediction_function, .id = "replicate")

limits_c <- all_preds_c %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            mean = mean(growth))

limits_v <- all_preds_v %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            mean = mean(growth))

limits_v_1000 <- all_preds_v_1000 %>% 
  group_by(temperature) %>% 
  summarise(q2.5=quantile(growth, probs=0.025),
            q97.5=quantile(growth, probs=0.975),
            mean = mean(growth))

v_fits <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_v_exp.csv")
c_fits <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_exp.csv")

tpc_v <-function(x){
  res<-v_fits$a.list[[1]]*exp(v_fits$b.list[[1]]*x)*(1-((x-v_fits$z.list[[1]])/(v_fits$w.list[[1]]/2))^2)
  res
}

tpc_c <-function(x){
  res<-c_fits$a.list[[1]]*exp(c_fits$b.list[[1]]*x)*(1-((x-c_fits$z.list[[1]])/(c_fits$w.list[[1]]/2))^2)
  res
}


all_estimates <- read_csv("Tetraselmis_experiment/data-processed/growth_estimates_boot_car.csv")
all_estimates_v <- read_csv("Tetraselmis_experiment/data-processed/growth_estimates_boot_car_v.csv")
limits_prediction <- read_csv("Tetraselmis_experiment/data-processed/limits_prediction.csv")
limits_prediction_STT <- read_csv("Tetraselmis_experiment/data-processed/limits_prediction_STT.csv")
variable_predictions_points <- read_csv("Tetraselmis_experiment/data-processed/variable_predictions_points.csv")

library(colormap)
ic <- colormap(colormap = colormaps$viridis, nshades = 8, format = "hex",
               alpha = 1, reverse = FALSE)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + 
  # geom_line(aes(x = temperature, y = growth, group = replicate), data = all_preds_v, alpha = 0.2) +
  stat_function(fun = tpc_c, color = "cadetblue") +
  stat_function(fun = tpc_v, color = "red") +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_c, fill = "cadetblue", alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_v, fill = "orange", alpha = 0.8) +
  # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_v_1000, fill = "red", alpha = 0.5) +
  # geom_point(aes(x = temp, y = mean), data = gs, color = "cadetblue") +
  geom_point(aes(x = temp, y = mean), data = gsv, color = "red") +
  geom_point(aes(x = temp, y = mean), data = gs, color = "cadetblue") +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates_v, color = "black", width = 0.2) +
  geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = "orange") +
  geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = "black", shape = 1) +
 coord_cartesian(xlim = c(-2, 32), ylim = c(-0.5, 1.6))


panel_a <- p + 
  stat_function(fun = tpc_c, color = ic[3], size = 1.5) +
  # stat_function(fun = tpc_v, color = "orange") +
  # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_v, fill = "orange", alpha = 0.5) +
  # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_c, fill = ic[3], alpha = 0.5) +
  # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
              # fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates, color = "black", width = 0.2) +
  geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = ic[3]) +
  geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = "black", shape = 1) +
  # geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates_v, color = "black", width = 0.2) +
  # geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = "orange") +
  # geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = "black", shape = 1) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
              fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction_STT,
              fill = "transparent", alpha = 0.01, linetype = "dashed", color = "darkgrey", size = 0.5) +
  geom_hline(yintercept = 0) + ylab("") +
  xlab("") + coord_cartesian(xlim = c(-2,33), ylim = c(-0.1, 1.6)) +
  labs(y = expression ("Population growth rate"~day^-1))

panel_b <- p + 
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_v, fill = ic[5], alpha = 0.5) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_c, fill = ic[3], alpha = 0.5) +
  stat_function(fun = tpc_c, color = ic[3], size = 1.5) +
  stat_function(fun = tpc_v, color = ic[5], size = 1.5) +
  # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
  # fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
  # geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates, color = "black", width = 0.2) +
  # geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = "cadetblue") +
  # geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = "black", shape = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates_v, color = "black", width = 0.2) +
  geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = ic[5]) +
  geom_point(aes(x = temperature, y = estimate), data = all_estimates_v, size = 2, color = "black", shape = 1) +
  geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction,
              fill = "transparent", alpha = 0.01, linetype = "dashed", color = "black", size = 0.5) +
  # geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_prediction_STT,
  #             fill = "transparent", alpha = 0.01, linetype = "dashed", color = "darkgrey", size = 0.5) +
  geom_hline(yintercept = 0) + ylab("") +
  xlab("") + coord_cartesian(xlim = c(-2,33), ylim = c(-0.1, 1.6)) +
  labs(y = expression ("Population growth rate"~day^-1))

plots <- plot_grid(panel_a, panel_b, labels = c("A", "B"), align = "v", nrow = 2)
ggsave(plots, file = "Tetraselmis_experiment/figures/figure2_indirect.png", width = 6, height = 7)


# fig for k-temp-supp -----------------------------------------------------

temp_arr_trans <- function(x) {(1/(.00008617*(x+273.15)))}

p + 
  geom_hline(yintercept = 0, color = "grey") +
  stat_function(fun = tpc_c, color ="black", size = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = temperature), data = all_estimates, color = "black", width = 0.2) +
  geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = "grey") +
  geom_point(aes(x = temperature, y = estimate), data = all_estimates, size = 2, color = "black", shape = 1) +
  geom_hline(yintercept = 0) + ylab("") +
  xlab("Temperature (°C)") + coord_cartesian(xlim = c(0,33), ylim = c(-0.1, 1.6)) +
  scale_x_continuous(sec.axis = sec_axis(~(temp_arr_trans(.))), limits = c(0, 32)) + 
  labs(y = expression ("Exponential growth rate"~day^-1))+
  xlab("Temperature (°C)") +
  ylab(bquote('Exponential growth rate'*~day^-1*'')) +
  theme_bw(base_family = "Arial", base_size = 12) +
  scale_x_continuous(sec.axis = sec_axis(~(temp_arr_trans(.))), limits = c(0, 32)) + 
  ylim(-0.2, 1.7) +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  theme_bw() +
  theme(text = element_text(size=12, family = "Arial"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("Temperature (1/kT)") 
ggsave("Tetraselmis_experiment/figures/k-temp-supp-TPC.pdf", width = 6, height = 4)
