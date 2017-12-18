#### Use daily temperature data from NOAA OISST to calculate a daily temperature SD

library(tidyverse)
library(viridis)


temps1 <- read_csv("Tetraselmis_experiment/data-processed/daily_temperatures_above99.csv")
temps2 <- read_csv("Tetraselmis_experiment/data-processed/daily_temperatures_1-97.csv")
results5 <- read_csv("Tetraselmis_experiment/data-processed/results5.csv")
global_therm <- read_csv("Tetraselmis_experiment/data-processed/global_therm.csv")


temps <- bind_rows(temps1, temps2)


temps %>% 
	filter(isolate.code == 472) %>% View





temps_sd <- temps %>% 
	group_by(isolate.code, lat, lon) %>% 
	summarise_each(funs(mean, sd), sst)


all <- left_join(results5, temps_sd, by = "isolate.code")

missing_isolates <- all %>% 
	filter(is.na(sst_sd)) %>% 
	select(latitude, lat, longitude, lon, everything()) %>% 
	select(isolate.code, latitude, longitude)

write_csv(missing_isolates, "Tetraselmis_experiment/data-processed/missing_isolates_NOAA.csv")


## make a map of the original isolation locations vs. the closest NOAA ones

global_therm %>% 
	filter(Biome == "marine") %>% 
	ggplot(aes(x=Lon,y=Lat,fill=SD)) +
	geom_tile() +
	scale_fill_viridis(option = "inferno") +
	scale_color_viridis(option = "inferno") +
	# geom_contour(data = filter(coastline,y>-80),aes(x=x,y=y,z=z,fill=NULL),breaks=0, colour="white", size = 0.5) +
	geom_point(data = all, aes(x=longitude,y=latitude), color = "dodgerblue", fill = "dodgerblue", size = 1.5, shape = 1) +
	geom_point(data = all, aes(x=lon,y=lat), color = "green", fill = "green", size = 1.5, shape = 1) +
	theme(axis.line        = element_blank(),
				axis.text        = element_blank(),
				axis.ticks       = element_blank(),
				axis.title       = element_blank(),
				panel.background = element_blank(),
				panel.border     = element_blank(),
				panel.margin = unit(0,"null"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				plot.background  = element_blank(),
				panel.spacing = unit(0, "lines"),
				plot.margin = rep(unit(0,"null"),4)) 

all %>% 
	ggplot() +
	geom_point(data = all, aes(x=longitude,y=latitude), color = "dodgerblue", fill = "dodgerblue", size = 6, shape = 1) +
	geom_point(data = all, aes(x=lon,y=lat), color = "green", fill = "green", size = 6, shape = 1) + theme_classic()
	



all %>% 
	ggplot(aes(x = SD, y = sst_sd)) + geom_point() +
	geom_abline(slope = 1, intercept = 0)

all %>%
	mutate(sd_diff = sst_sd - SD) %>% 
	ggplot(aes(x = latitude, y = sd_diff)) + geom_point() +
	geom_hline(yintercept = 0) + theme_classic()



### what is the difference between SD and sst_sd (i.e. monthtly SD and daily SD)?
## looks like the biggest difference is around 40 degrees


all %>% 
filter(minqual == "good", maxqual == "good") %>% 
	# filter(mu.n > 4) %>%
	# filter(isolate.code %in% greater_90) %>% 
	# filter(rel.curveskew < 0) %>% 
	ggplot(aes(x = rel.curveskew, y = diff, color = sst_sd, label = isolate.code)) + geom_point(size = 2) +
	# geom_label()
	scale_color_viridis(option = "inferno") + theme_bw() + geom_smooth(method = "lm", color = "grey", alpha = 0.3, size = 0.4) +
	geom_hline(yintercept = 0) + ylab("Topt(cons) - Topt(var) (°C)") +
	xlab("Curve skewness") + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black")) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	geom_vline(xintercept = 0)


