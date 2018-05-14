

cells_days_v <- read_csv("Tetraselmis_experiment/data-processed/cells_days_v_mod.csv") %>% 
  mutate(environment = "variable")
cells_days_c <- read_csv("Tetraselmis_experiment/data-processed/cells_exp_mod.csv") %>% 
  mutate(days = time_since_innoc_hours/24) %>% 
  mutate(days = ifelse(days < 0, 0, days)) %>% 
  mutate(environment = "constant")

all_cells_plot <- bind_rows(cells_days_c, cells_days_v)

all_cells_plot %>% 
  distinct(temp, variability) %>% View

cfc <- read_csv("Tetraselmis_experiment/data-processed/ctpc.csv")
cfv <- read_csv("Tetraselmis_experiment/data-processed/vtpc.csv")

800 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days))

growth_fun_27c <- function(x, temp = 27){
  res <- 800 * exp((cfc[["a"]]*exp(cfc[["b"]]*temp)*(1-((temp-cfc[["z"]])/(cfc[["w"]]/2))^2))*(x))
  res
}
growth_fun_27v <- function(x, temp = 27){
  res <- 800 * exp((cfv[["a"]]*exp(cfv[["b"]]*temp)*(1-((temp-cfv[["z"]])/(cfv[["w"]]/2))^2))*(x))
  res
}

growth_fun_5c <- function(x, temp = 5){
  res <- 800 * exp((cfc[["a"]]*exp(cfc[["b"]]*temp)*(1-((temp-cfc[["z"]])/(cfc[["w"]]/2))^2))*(x))
  res
}

growth_fun_5v <- function(x, temp = 5){
  res <- 800 * exp((cfv[["a"]]*exp(cfv[["b"]]*temp)*(1-((temp-cfv[["z"]])/(cfv[["w"]]/2))^2))*(x))
  res
}

growth_fun_10c <- function(x, temp = 10){
  res <- 800 * exp((cfc[["a"]]*exp(cfc[["b"]]*temp)*(1-((temp-cfc[["z"]])/(cfc[["w"]]/2))^2))*(x))
  res
}
growth_fun_10v <- function(x, temp = 10){
  res <- 800 * exp((cfv[["a"]]*exp(cfv[["b"]]*temp)*(1-((temp-cfv[["z"]])/(cfv[["w"]]/2))^2))*(x))
  res
}
growth_fun_24c <- function(x, temp = 24){
  res <- 800 * exp((cfc[["a"]]*exp(cfc[["b"]]*temp)*(1-((temp-cfc[["z"]])/(cfc[["w"]]/2))^2))*(x))
  res
}

growth_fun_24v <- function(x, temp = 24){
  res <- 800 *  (1+(cfv[["a"]]*exp(cfv[["b"]]*temp)*(1-((temp-cfv[["z"]])/(cfv[["w"]]/2))^2)))^(x)
  res
}

growth_fun_16c <- function(x, temp = 16){
  res <- 800 * exp((cfc[["a"]]*exp(cfc[["b"]]*temp)*(1-((temp-cfc[["z"]])/(cfc[["w"]]/2))^2))*(x))
  res
}

growth_fun_15v <- function(x, temp = 15){
  res <- 800 * exp((cfv[["a"]]*exp(cfv[["b"]]*temp)*(1-((temp-cfv[["z"]])/(cfv[["w"]]/2))^2))*(x))
  res
}

growth_fun_29c <- function(x, temp = 29){
  res <- 800 * exp((cfc[["a"]]*exp(cfc[["b"]]*temp)*(1-((temp-cfc[["z"]])/(cfc[["w"]]/2))^2))*(x))
  res
}

growth_fun_32c <- function(x, temp = 32){
  res <- 800 * exp((cfc[["a"]]*exp(cfc[["b"]]*temp)*(1-((temp-cfc[["z"]])/(cfc[["w"]]/2))^2))*(x))
  res
}
growth_fun_0c <- function(x, temp = 0){
  res <- 800 * exp((cfc[["a"]]*exp(cfc[["b"]]*temp)*(1-((temp-cfc[["z"]])/(cfc[["w"]]/2))^2))*(x))
  res
}
growth_fun_20c <- function(x, temp = 20){
  res <- 800 * exp((cfc[["a"]]*exp(cfc[["b"]]*temp)*(1-((temp-cfc[["z"]])/(cfc[["w"]]/2))^2))*(x))
  res
}
growth_fun_20v <- function(x, temp = 20){
  res <- 800 * exp((cfv[["a"]]*exp(cfv[["b"]]*temp)*(1-((temp-cfv[["z"]])/(cfv[["w"]]/2))^2))*(x))
  res
}


library(viridis)
library(colormap)
library(cowplot)
ic <- colormap(colormap = colormaps$viridis, nshades = 9, format = "hex",
               alpha = 1, reverse = FALSE)

all_cells_plot %>% 
  rename(Environment = variability) %>% 
  mutate(Environment = ifelse(Environment == "c", "constant", "variable")) %>% 
  rename(Temperature = temp) %>% 
  ggplot(aes(x = days, y = cell_density, color = factor(Temperature), shape = Environment)) +
  geom_point(size = 4) +
  # # geom_point(size = 4, shape= variability, color = "black") +
  stat_function(fun = growth_fun_27v, color = ic[7], size = 1) +
  stat_function(fun = growth_fun_27c, color = ic[7], size = 1, linetype = "dashed") +
  stat_function(fun = growth_fun_29c, color = ic[8], size = 1, linetype = "dashed") +
  stat_function(fun = growth_fun_16c, color = ic[4], size = 1, linetype = "dashed") +
  stat_function(fun = growth_fun_15v, color = ic[4], size = 1) +
  stat_function(fun = growth_fun_5c, color = ic[2], size = 1, linetype = "dashed") +
  stat_function(fun = growth_fun_5v, color = ic[2], size = 1) +
  stat_function(fun = growth_fun_20c, color = ic[5], size = 1, linetype = "dashed") +
  stat_function(fun = growth_fun_20v, color = ic[5], size = 1) +
   ylim(0, 30000) +
  stat_function(fun = growth_fun_10c, color = ic[3], size = 1, linetype = "dashed") +
  stat_function(fun = growth_fun_10v, color = ic[3], size = 1) +
  stat_function(fun = growth_fun_24c, color = ic[6], size = 1, linetype = "dashed") +
  stat_function(fun = growth_fun_24v, color = ic[6], size = 1) +
  stat_function(fun = growth_fun_32c, color = ic[9], size = 1, linetype = "dashed") +
  stat_function(fun = growth_fun_0c, color = ic[1], size = 1, linetype = "dashed") +
  scale_color_viridis(discrete = TRUE, name = "Temperature") +
  ylab("Population abundance (cells/mL)") + xlab("Days")
ggsave("Tetraselmis_experiment/figures/cell_abundance_days.pdf", width = 10, height = 6)

