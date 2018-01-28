

library(tidyverse)
library(colormap)

boot_limits_constant <- read_csv("Tetraselmis_experiment/data-processed/boot_limits_constant_resample_10k.csv") ### these are for the new resampled curve
boot_limits_variable <- read_csv("Tetraselmis_experiment/data-processed/boot_limits_constant_resample_v10k.csv")

fits_v <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_v.csv")
fits_c <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params.csv")
growth_sum_v <- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary_v.csv")
growth_sum<- read_csv("Tetraselmis_experiment/data-processed/resampled_growth_rates_summary.csv")
variable_predictions_points <- read_csv("Tetraselmis_experiment/data-processed/variable_predictions_points.csv")


fits1 <- read_csv("Tetraselmis_experiment/data-processed/boot_upper_lower_fits2.csv")
# fits2 <- read_csv("Tetraselmis_experiment/data-processed/boot_upper_lower_fits.csv")

fits <- fits1
bootcurve_upper<-function(x){
  res<-fits$a.list[2]*exp(fits$b.list[2]*x)*(1-((x-fits$z.list[2])/(fits$w.list[2]/2))^2)
  res
}

bootcurve_lower<-function(x){
  res<-fits$a.list[1]*exp(fits$b.list[1]*x)*(1-((x-fits$z.list[1])/(fits$w.list[1]/2))^2)
  res
}


curve_constant_resamp<-function(x){
  res<-fits_c$a.list[1]*exp(fits_c$b.list[1]*x)*(1-((x-fits_c$z.list[1])/(fits_c$w.list[1]/2))^2)
  res
}

curve_variable_resamp<-function(x){
  res<-fits_v$a.list[1]*exp(fits_v$b.list[1]*x)*(1-((x-fits_v$z.list[1])/(fits_v$w.list[1]/2))^2)
  res
}


preds <- growth_sum_v %>% 
  mutate(prediction = 0.5*(curve_constant_resamp(temp -5) + curve_constant_resamp(temp + 5))) %>% 
  mutate(prediction_upper = 0.5*(bootcurve_upper(temp -5) + bootcurve_upper(temp + 5))) %>% 
  mutate(prediction_lower = 0.5*(bootcurve_lower(temp -5) + bootcurve_lower(temp + 5))) 


ic <- colormap(colormap = colormaps$viridis, nshades = 8, format = "hex",
               alpha = 1, reverse = FALSE)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 
p + 
  geom_hline(yintercept = 0, color = "darkgrey") +
  # geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = inflection, color = "grey", size = 0.5) +
  # stat_function(fun = curve_constant_resamp, color = ic[3], size = 1.5, alpha = 0.7) +
  # geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_constant, fill = ic[3], alpha = 0.5) +
  # geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_variable, fill = ic[5], alpha = 0.3) +
  labs(y = expression ("Population growth rate"~day^-1))+
  # stat_function(fun = curve_variable_resamp, color = ic[5], size = 1.5, alpha = 0.7) +
  geom_point(aes(x = temp, y = mean), data = growth_sum_v, color = ic[5], size = 2.5) +
  geom_errorbar(aes(x = temp, ymin = lower, ymax = upper), data = growth_sum_v, width = 0.1, color = ic[5]) +
   xlab("") + 
  geom_errorbar(aes(x = temp, ymin = lower, ymax = upper), data = growth_sum_v, width = 0.2, color = ic[5]) +
  geom_point(aes(x = temp, y = prediction), data = preds, color = "orange", size = 2.5) +
  geom_errorbar(aes(x = temp, ymin = prediction_lower, ymax = prediction_upper), data = preds, width = 0.1, color = "orange") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_cartesian(ylim = c(-0.5, 1.8), xlim= c(0, 32)) +
  theme_classic() +
  theme(text = element_text(size=14)) +
  geom_ribbon(aes(x = temperature, ymin = growth.rate.lower, ymax = growth.rate.upper, linetype = NA), fill = "transparent", alpha = 0.01, data = variable_predictions_points, linetype = "dashed", color = "black", size = 0.5) +
  coord_cartesian(ylim = c(-0.2, 1.7), xlim = c(0, 32))
ggsave("Tetraselmis_experiment/figures/averaging_preds2.pdf")



x <- seq(-2, 32, by = 0.01)


variable_predictions_upper <- function(x) {
  y <- 0.5*(bootcurve_upper(x + 5) + bootcurve_upper(x - 5))
}
 

variable_upper <- sapply(x, variable_predictions_upper)
variable_upper2 <- data.frame(x, variable_upper) %>% 
  rename(temperature = x, 
         prediction_upper = variable_upper)


variable_predictions_lower <- function(x) {
  y <- 0.5*(bootcurve_lower(x + 5) + bootcurve_lower(x - 5))
}


variable_lower <- sapply(x, variable_predictions_lower)
variable_lower2 <- data.frame(x, variable_lower) %>% 
  rename(temperature = x, 
         prediction_lower = variable_lower)


all_preds_average <- left_join(variable_lower2, variable_upper2)

library(colormap)
variable_predictions_points <- read_csv("Tetraselmis_experiment/data-processed/variable_predictions_points.csv")

ic <- colormap(colormap = colormaps$viridis, nshades = 8, format = "hex",
               alpha = 1, reverse = FALSE)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 
p + 
  geom_hline(yintercept = 0, color = "darkgrey") +
   labs(y = expression ("Population growth rate"~day^-1))+
  stat_function(fun = curve_variable_resamp, color = ic[5], size = 1.5, alpha = 0.7) +
  stat_function(fun = curve_constant_resamp, color = ic[3], size = 1.5, alpha = 0.7) +
  geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_constant, fill = ic[3], alpha = 0.5) +
  geom_ribbon(aes(x = x, ymin = q2.5, ymax = q97.5, linetype=NA), data = boot_limits_variable, fill = ic[5], alpha = 0.3) +
  geom_point(aes(x = temp, y = mean), data = growth_sum_v, color = ic[5], size = 2.5) +
  geom_errorbar(aes(x = temp, ymin = lower, ymax = upper), data = growth_sum_v, width = 0.1, color = ic[5]) +
  xlab("") + 
  geom_errorbar(aes(x = temp, ymin = lower, ymax = upper), data = growth_sum_v, width = 0.2, color = ic[5]) +
  geom_point(aes(x = temp, y = prediction), data = preds, color = "orange", size = 2.5) +
  geom_errorbar(aes(x = temp, ymin = prediction_lower, ymax = prediction_upper), data = preds, width = 0.1, color = "orange") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_cartesian(ylim = c(-0.5, 1.8), xlim= c(0, 32)) +
  theme_classic() +
  theme(text = element_text(size=14)) +
  geom_ribbon(aes(x = temperature, ymin = prediction_lower, ymax = prediction_upper, linetype = NA), fill = "transparent", alpha = 0.01, data = all_preds_average, linetype = "dashed", color = "black", size = 0.5) +
  coord_cartesian(ylim = c(-0.2, 1.7), xlim = c(0, 32))
  # geom_ribbon(aes(x = temperature, ymin = growth.rate.lower, ymax = growth.rate.upper, linetype = NA), fill = "transparent", alpha = 0.01, data = variable_predictions_points, linetype = "dashed", color = "red", size = 0.5)
  
ggsave("Tetraselmis_experiment/figures/average_prediction.png")
