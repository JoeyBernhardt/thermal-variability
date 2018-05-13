

params_raw <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_exp.csv")

nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

## get Tmax
tmax_c <- uniroot.all(function(x) nbcurve(x, z = c_fits$z.list[[1]],w = c_fits$w.list[[1]],a = c_fits$a.list[[1]], b = c_fits$b.list[[1]]),c(10,150))
tmax_v <- uniroot.all(function(x) nbcurve(x, z = v_fits$z.list[[1]],w = v_fits$w.list[[1]],a = v_fits$a.list[[1]], b = v_fits$b.list[[1]]),c(10,150))



fits_c <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_5000_exp.csv")
# fits_c <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_1000_exp.csv")
fits_v <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v_exp3.csv")
# fits_v_1000 <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v_exp2.csv")
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

bs_split <- fits_c %>% 
  filter(maxgrowth.list < 2) %>% 
  split(.$curve.id.list)

bs_v_split <- fits_v %>%
  filter(maxgrowth.list < 2) %>% 
  split(.$curve.id.list)


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


rmax_summ <- fits_v %>% 
  filter(maxgrowth.list < 2) %>% 
  summarise(lower=quantile(maxgrowth.list, probs=0.025),
            upper=quantile(maxgrowth.list, probs=0.975),
            mean = mean(maxgrowth.list),
            median = median(maxgrowth.list)) %>% 
  mutate(param = "rmax")

rmax_summ_c <- fits_c %>% 
  filter(maxgrowth.list < 2) %>% 
  summarise(lower=quantile(maxgrowth.list, probs=0.025),
            upper=quantile(maxgrowth.list, probs=0.975),
            mean = mean(maxgrowth.list),
            median = median(maxgrowth.list)) %>% 
  mutate(param = "rmax")


library(colormap)
ic <- colormap(colormap = colormaps$viridis, nshades = 8, format = "hex",
               alpha = 1, reverse = FALSE)

points <- read_csv("Tetraselmis_experiment/data-processed/pred_obs_points.csv") %>% 
  mutate(obs_type = as.factor(obs_type))

?levels

points$obs_type <- factor(points$obs_type, c("rmax", "Topt", "Tmax"))
points %>% 
  unite(col = env_pred, sep = " ", environment, predicted, remove = FALSE) %>% 
  ggplot(aes(x = env_pred, y = obs, color = environment)) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + facet_wrap( ~ obs_type, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_manual(values = c(ic[3], ic[5])) +
  theme(strip.background = element_rect(colour="white", fill="white")) + xlab("") + ylab("")
ggsave("Tetraselmis_experiment/figures/pred_obs_dot_plot.pdf", width = 10, height = 4)



all_critical_temps <- read_csv("Tetraselmis_experiment/data-processed/all_critical_temps.csv")
library(stringr)

all_ct <- all_critical_temps %>% 
  mutate(predicted = ifelse(grepl("predicted", param), "predicted", "observed")) %>% 
  mutate(param = str_replace(param, "predicted_", "")) %>% 
  mutate(param = str_replace(param, "topt", "Topt")) %>% 
  mutate(param = str_replace(param, "tmax", "Tmax")) %>% 
  filter(param %in% c("rmax", "Topt", "Tmax"))

all_ct$param <- factor(all_ct$param, c("rmax", "Topt", "Tmax"))

all_ct2 <- all_ct %>% 
  unite(col = env_pred, sep = " ", treatment, predicted, remove = FALSE) %>% 
  rename(environment = treatment) %>% 
  mutate(approach = "Direct")

points$approach <- "Indirect"

points2 <- points %>% 
  unite(col = env_pred, sep = " ", environment, predicted, remove = FALSE) %>% 
  rename(mean = obs) %>% 
  rename(param = obs_type)

all_points <- bind_rows(points2, all_ct2)

all_points %>%   
ggplot(aes(x = env_pred, y = mean, color = environment, shape = predicted)) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + facet_wrap( ~ param + approach, scales = "free", ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_manual(values = c(ic[3], ic[5])) +
  theme(strip.background = element_rect(colour="white", fill="white")) + xlab("") + ylab("") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.text=element_text(size=14))
ggsave("Tetraselmis_experiment/figures/pred_obs_dot_plot_all.pdf", width = 8, height = 8)
