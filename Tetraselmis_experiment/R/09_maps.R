library(tidyverse)
library(ggplot2)
library(viridis)
library(ggalt)
library(dplyr)
library(rgdal)
library(ggmap)
library(tools)
library(marmap)
library(dplyr)


global_therm <- read_csv("Tetraselmis_experiment/data-processed/global_therm.csv")

global_therm %>% 
	filter(Lat == 49.5, Lon == -126.5) %>% View


thomas <- read_csv("data/thermal_trait_data/Thomas_2014_traits_derived_20140606.csv")

thomas2 <- thomas %>% 
	rename(topt = mu.g.opt.list) %>% 
	rename(w = mu.wlist,
				 a = mu.alist,
				 z = mu.c.opt.list,
				 b = mu.blist) %>% 
	select(topt, w, a, b, z, everything()) 
thomas3 <- thomas2 %>% 
	rename(latitude = isolation.latitude,
				 longitude = isolation.longitude) %>% 
	filter(!is.na(latitude), !is.na(longitude)) %>% 
	filter(habitat %in% c("marine", "estuarine", "saltmarsh"))

all_thomas <- left_join(thomas3, results5, by = "isolate.code") %>% 
	filter(!is.na(SD))

## begin mapping
bathy.map <- getNOAA.bathy(lon1 = -180, lon2 = 180,lat1 = -80, lat2 = 90, resolution = 10)
coastline <- fortify(bathy.map) 

mapWorld2 <- borders("world", colour="gray50", fill="gray50", alpha = 0, size = 1) 
mapWorld <- borders("world", colour="gray50", fill="gray50", size = 0.5) 

wmap <- global_therm %>% 
	filter(Biome == "marine") 

wmap_robin <- spTransform(wmap, CRS("+proj=robin"))
wmap_df_robin <- fortify(wmap_robin)

global_therm %>% 
	filter(Biome == "marine") %>% 
ggplot(aes(x=Lon,y=Lat,fill=SD)) +
	mapWorld +
	geom_tile(alpha = 0.8) +
	# geom_contour(data = filter(coastline),aes(x=x,y=y,z=z, fill= NULL),breaks=0, colour="grey50", size = 0.75) +
	# geom_point(data=thomas3, aes(x=longitude,y=latitude), color = "dodgerblue", fill = "dodgerblue", size = 2, shape = 1) +
	scale_fill_viridis(option = "inferno") +
	scale_color_viridis(option = "inferno") +
	# mapWorld +
	theme_classic() +
	# coord_map() +
	# coord_proj("+proj=wintri +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") +
	coord_quickmap() +
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
				plot.margin = rep(unit(0,"null"),4)) +
	theme(legend.position = "none") +
	# geom_tile(alpha = 0.7) +
	geom_point(data=all_thomas, aes(x=longitude.x,y=latitude.x, color = SD), fill = "white", size = 3) +
	geom_point(data=all_thomas, aes(x=longitude.x,y=latitude.x), color = "black", fill = "black", size = 3, shape = 1, alpha = 0.5) 
ggsave("Tetraselmis_experiment/figures/SD_map.png", width = 6, height = 4)



	
	

	ggsave("Tetraselmis_experiment/figures/SD_map.png", width = 6, height = 4)
	ggsave("Tetraselmis_experiment/figures/SD_map.pdf")
 
	
	
	my_map <- get_map(source="stamen", maptype= "toner", crop=FALSE)
ggmap(my_map)

