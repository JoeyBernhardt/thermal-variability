## Extract the cell density info from the flowcam files
## V2 dataset

# load libraries ----------------------------------------------------------
library(readr)
library(tidyverse)
library(purrr)
library(lubridate)
library(stringr)
library(dplyr)

#### Step 1 #### 
## in the shell, use this command to copy the summary files from the folder of flow cam run folders to a summary-only file
## something like: cp **/*summary.csv /Users/Joey/Desktop/run-summaries

## step 1b, read in UniqueID key (if there is one)

# Unique_ID_key <- read.csv("/Users/Joey/Documents/chlamy-ktemp/k-temp/data-raw/PK-temp-UniqueID-key.csv")

#### Step 2: create a list of file names for each of the summaries ####

cell_files <- c(list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun07", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun08", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun09", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun09-eve", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun10", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun10-eve", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun11-morn", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun11-eve", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun12-morn", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun12-eve", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun13-eve", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun14-morn", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun14-eve", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun15-eve", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun16-eve", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun16-morn", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun17-eve", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun17-midday", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun17-morn", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun18-eve", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun18-midday", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun18-morn", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun19-morn", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun23", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun24", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun26", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun28", full.names = TRUE))

names(cell_files) <- cell_files %>% 
	gsub(pattern = ".csv$", replacement = "")

test <- read.csv("Tetraselmis_experiment/data-raw/flowcam-summaries-jun16-eve/5v-5-4_summary.csv")
#### Step 3: read in all the files!

all_cells <- map_df(cell_files, read_csv)



#### Step 4: pull out just the data we want, do some renaming etc.

TT_cells <- all_cells %>% 
	rename(obs_type = X1,
				 value = X2) %>% 
	filter(obs_type %in% c("List File", "Start Time", "Particles / ml", "Volume (ABD)")) %>%
	spread(obs_type, value) %>%
	separate(`List File`, into = c("temperature", "sample_day", "replicate"), sep = "[:punct:]") %>% 
	separate(temperature, into = c("temp", "variability"), sep = -2) %>%
	rename(start_time = `Start Time`,
				 cell_density = `Particles / ml`,
				 cell_volume = `Volume (ABD)`) 
write_csv(TT_cells, "Tetraselmis_experiment/data-processed/TT_cells_round3.csv")



### new 0C and 32C data

cell_files <- c(list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun23", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun24", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun26", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-jun28", full.names = TRUE))

names(cell_files) <- cell_files %>% 
	gsub(pattern = ".csv$", replacement = "")

all_cells <- map_df(cell_files, read_csv, n_max = 6, col_names = c("obs_type", "value"), .id = "file_name")


TT_extremes <- all_cells %>% 
	filter(obs_type %in% c("List File", "Start Time", "Particles / ml")) %>% 
	spread(obs_type, value) %>% 
	separate(`List File`, into = c("temperature", "sample_day", "replicate"), sep = "[:punct:]") %>% 
	separate(temperature, into = c("temp", "variability"), sep = -2) %>%
	rename(start_time = `Start Time`,
				 cell_density = `Particles / ml`) 
write_csv(TT_extremes, "Tetraselmis_experiment/data-processed/TT_cells_round3_extremes.csv")


### now bring in the 0C data from round 1
cell_files <- c(list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar18", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-may01", full.names = TRUE))

names(cell_files) <- cell_files %>% 
	gsub(pattern = ".csv$", replacement = "")

all_cells <- map_df(cell_files, read_csv, n_max = 6, col_names = c("obs_type", "value"), .id = "file_name")


TT_0 <- all_cells %>% 
	filter(obs_type %in% c("List File", "Start Time", "Particles / ml")) %>% 
	spread(obs_type, value) %>% 
	filter(grepl("^0", `List File`))%>% 
	separate(`List File`, into = c("temperature", "sample_day", "replicate"), sep = "[:punct:]") %>% 
	# separate(temperature, into = c("temp", "variability"), sep = -2) %>%
	rename(start_time = `Start Time`,
				 cell_density = `Particles / ml`)
write_csv(TT_0, "Tetraselmis_experiment/data-processed/TT_cells_round2_0.csv")

