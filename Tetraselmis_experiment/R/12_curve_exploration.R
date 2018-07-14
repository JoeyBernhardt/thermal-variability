library(tidyverse)
library(janitor)
library(purrr)
library(cowplot)

growth_raw <- read_csv("Tetraselmis_experiment/data-processed/growth_data_20140606.csv")
all_thermal_data <- read_csv("Tetraselmis_experiment/data-processed/all_thermal_data.csv") %>% 
	clean_names
curve_data <- read_csv("Tetraselmis_experiment/data-processed/curve_data_20140606.csv") %>% 
	clean_names()

growth <- growth_raw %>% 
	clean_names()


# all_growth <- left_join(growth, all_thermal_data, by = c("curve_code" = "isolate_code"))



all1 <- left_join(curve_data, all_thermal_data) %>% 
	filter(!is.na(topt))


all2 <- left_join(growth, all1)



growth %>% 
	ggplot(aes(x = temperature, y = growth_rate_mu, group = curve_code, color = curve_code)) + geom_line()


pos_skew <- all2 %>% 
	filter(rel_curveskew > 0)

pos_curves <- unique(pos_skew$curve_code)

all2 %>% 
	filter(!is.na(topt)) %>% 
	# filter(curve_code %in% pos_curves) %>% 
	filter(curvequal == "good", maxqual == "good", minqual == "good") %>% 
	filter(mu_n > 4) %>% 
	ggplot(aes(x = temperature, y = growth_rate_mu, group = curve_code, color = curve_code)) + geom_point() +
	geom_line() + facet_wrap( ~ isolate_code, scales = "free") +theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=10, family = "Helvetica")) +
	theme(strip.background = element_rect(colour="white", fill="white"))
ggsave("Tetraselmis_experiment/figures/all_curves.png", width = 12, height = 12)


curve <- all2 %>% 
	# filter(grepl("Detonula", speciesname)) %>% 
	filter(isolate_code == 49) 


tpc1<-function(x){
	res<-curve$a[1]*exp(curve$b[1]*x)*(1-((x-curve$z[1])/(curve$w[1]/2))^2)
	res
}

library(plotrix)
curve %>% 
	group_by(temperature) %>% 
	summarise_each(funs(mean, std.error), growth_rate_mu) %>% 
	ggplot(aes(x = temperature, y = growth_rate_mu_mean)) + geom_point(size = 2) +
	geom_errorbar(aes(ymin = growth_rate_mu_mean - growth_rate_mu_std.error, ymax = growth_rate_mu_mean + growth_rate_mu_std.error))+
	stat_function(fun = tpc1) + theme_bw() + xlim(0, 50) 
	
	
#### figure out how many actual temperatures each data set has ####


temp_numbers <- all2 %>% 
	group_by(isolate_code) %>% 
summarise(number_of_temps = length(unique(temperature)))


library(modelr)
??modelr

all2 %>% 
	filter(isolate_code == 52) %>% View



tpc1<-function(x){
	res<-curve$a[1]*exp(curve$b[1]*x)*(1-((x-curve$z[1])/(curve$w[1]/2))^2)
	res
}

predict_function <- function(data){
data_sub <- data
x <- seq(min(data_sub$temperature), max(data_sub$temperature), by = 0.5)
y <- data$a[1]*exp(data$b[1]*x)*(1-((x-data$z[1])/(data$w[1]/2))^2)
df <- data.frame(x, y)
df
}

sub <- all2 %>% 
	filter(isolate_code == 52)

all_split <- all2 %>% 
	split(.$isolate_code)


predictions_df <- all_split %>% 
	map_df(predict_function, .id = "isolate_code")

predictions_df <- predictions_df %>% 
	mutate(isolate_code = as.integer(isolate_code))

write_csv(predictions_df, "Tetraselmis_experiment/data-processed/TPC_predictions.csv")

all3 <- left_join(predictions_df, all2, by = "isolate_code")

write_csv(all3, "Tetraselmis_experiment/data-processed/global_TPCs.csv")

all3 <- read_csv("Tetraselmis_experiment/data-processed/global_TPCs.csv")
all3 %>% 
	filter(!is.na(topt)) %>% 
  filter(mu_n > 4) %>% 
	filter(curve_code %in% pos_curves) %>% 
	mutate(skew_dir = ifelse(rel_curveskew<0, "negative skew", "positive skew")) %>% 
	filter(curvequal == "good", maxqual == "good", minqual == "good") %>% 
	filter(mu_rsqrlist > 0.85) %>% 
	filter(mu_n > 4) %>% 
	distinct(isolate_code, .keep_all = TRUE) %>% View
	ggplot(aes(x = temperature, y = growth_rate_mu, color = skew_dir)) + geom_point() +
	geom_line(aes(x = x, y = y)) + facet_wrap( ~ isolate_code, scales = "free_y") +theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=10, family = "Helvetica")) +
	theme(strip.background = element_rect(colour="white", fill="white"))
ggsave("Tetraselmis_experiment/figures/all_curves_fits.png", width = 12, height = 12)
ggsave("Tetraselmis_experiment/figures/all_curves_fits_skew.png", width = 12, height = 12)
ggsave("Tetraselmis_experiment/figures/all_curves_fits_skew_rsq95.png", width = 12, height = 12)
ggsave("Tetraselmis_experiment/figures/all_curves_fits_skew_rsq85.png", width = 12, height = 12)
ggsave("Tetraselmis_experiment/figures/all_curves_fits_skew_rsq90.png", width = 12, height = 12)
ggsave("Tetraselmis_experiment/figures/all_curves_fits_skew_rsq90_samex_axis.png", width = 15, height = 12)
	

pos_skews <- all3 %>% 
  filter(!is.na(topt)) %>% 
  filter(mu_n > 4) %>% 
  filter(curve_code %in% pos_curves) %>% 
  mutate(skew_dir = ifelse(rel_curveskew<0, "negative skew", "positive skew")) %>% 
  filter(curvequal == "good", maxqual == "good", minqual == "good") %>% 
  filter(mu_rsqrlist > 0.85) %>% 
  filter(mu_n > 4) %>% 
  distinct(isolate_code, .keep_all = TRUE) %>% 
  select(isolate_code, source, latitude, longitude, habitat, speciesname, class, group)

write_csv(pos_skews, "Tetraselmis_experiment/data-processed/pos_skew.csv")

all3 %>% 
	filter(!is.na(topt)) %>% 
	filter(curve_code %in% pos_curves) %>% 
	mutate(skew_dir = ifelse(rel_curveskew<0, "negative skew", "positive skew")) %>% 
	# filter(curvequal == "good", maxqual == "good", minqual == "good") %>% 
	filter(mu_rsqrlist > 0.90) %>% 
	filter(mu_n > 4) %>% 
	distinct(isolate_code, .keep_all = TRUE) %>% View
	ggplot(aes(x = temperature, y = growth_rate_mu)) +
	# geom_point() +
	geom_line(aes(x = x, y = y)) + facet_wrap( ~ isolate_code, scales = "free_y") +theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=10, family = "Helvetica")) +
	theme(strip.background = element_rect(colour="white", fill="white"))
ggsave("Tetraselmis_experiment/figures/all_curves_fits_nodata_quality.png", width = 13, height = 13)
ggsave("Tetraselmis_experiment/figures/all_curves_fits_nodata_quality_nopoints.png", width = 13, height = 13)


all4 <- all3 %>% 
	filter(!is.na(topt)) %>% 
	# filter(curve_code %in% pos_curves) %>% 
	mutate(skew_dir = ifelse(rel_curveskew<0, "negative skew", "positive skew")) %>% 
	filter(curvequal == "good", maxqual == "good", minqual == "good") %>% 
	filter(mu_rsqrlist > 0.90) %>% 
	filter(mu_n > 4)


greater_90 <- unique(all4$isolate_code)

all5 <- all3 %>% 
	filter(!is.na(topt)) %>% 
	# filter(curve_code %in% pos_curves) %>% 
	mutate(skew_dir = ifelse(rel_curveskew<0, "negative skew", "positive skew")) %>% 
	filter(curvequal == "good", maxqual == "good", minqual == "good") %>% 
	filter(mu_rsqrlist > 0.95) %>% 
	filter(mu_n > 4)


greater_95 <- unique(all5$isolate_code)


