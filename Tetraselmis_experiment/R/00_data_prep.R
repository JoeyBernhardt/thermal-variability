## Extract the cell density info from the flowcam files

# load libraries ----------------------------------------------------------

library(tidyverse)
library(purrr)
library(lubridate)
library(stringr)


#### Step 1 #### 
## in the shell, use this command to copy the summary files from the folder of flow cam run folders to a summary-only file
## something like: cp **/*summary.csv /Users/Joey/Desktop/run-summaries

## step 1b, read in UniqueID key (if there is one)

# Unique_ID_key <- read.csv("/Users/Joey/Documents/chlamy-ktemp/k-temp/data-raw/PK-temp-UniqueID-key.csv")

#### Step 2: create a list of file names for each of the summaries ####

cell_files <- c(list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-feb18", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-feb19", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-feb20", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-feb21", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-feb23", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-feb25", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar1", full.names = TRUE))

cell_files2 <- c(list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar10", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar11", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar12", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar13", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar14", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar16", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar17", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar18", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar20", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar21", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar22", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar23", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar24", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar27", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar28", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar29", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-mar30", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-apr03", full.names = TRUE),
								list.files("Tetraselmis_experiment/data-raw/flowcam-summaries-apr05", full.names = TRUE))




names(cell_files2) <- cell_files2 %>% 
	gsub(pattern = ".csv$", replacement = "")


#### Step 3: read in all the files!

all_cells <- map_df(cell_files2, read_csv, col_names = FALSE, .id = "file_name")



#### Step 4: pull out just the data we want, do some renaming etc.

TT_cells <- all_cells %>% 
	rename(obs_type = X1,
				 value = X2) %>% 
	filter(obs_type %in% c("List File", "Start Time", "Particles / ml", "Volume (ABD)")) %>%
	spread(obs_type, value) %>% 
	separate(`List File`, into = c("replicate", "sample_day"), sep = "[:punct:]") %>% 
	rename(start_time = `Start Time`,
				 cell_density = `Particles / ml`,
				 cell_volume = `Volume (ABD)`) 


TT_cells <- all_cells %>% 
	rename(obs_type = X1,
				 value = X2) %>% 
	filter(obs_type %in% c("List File", "Start Time", "Particles / ml", "Volume (ABD)")) %>%
	spread(obs_type, value) %>% 
	separate(`List File`, into = c("replicate", "sample_day"), sep = "[:punct:]") %>% 
	rename(start_time = `Start Time`,
				 cell_density = `Particles / ml`,
				 cell_volume = `Volume (ABD)`) %>% 
	filter(!grepl("in", replicate)) %>% 
	filter(!grepl("tt", replicate))
	


write_csv(TT_cells, "Tetraselmis_experiment/data-processed/TT_constant_cells.csv")
