### code for bringing in flow cam data for light project

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

cell_files <- c(list.files("data-raw/flowcam-summaries-apr30", full.names = TRUE),
								list.files("data-raw/flowcam-summaries-may01", full.names = TRUE),
								list.files("data-raw/flowcam-summaries-may02", full.names = TRUE),
								list.files("data-raw/flowcam-summaries-may03", full.names = TRUE),
								list.files("data-raw/flowcam-summaries-may04", full.names = TRUE),
								list.files("data-raw/flowcam-summaries-may05", full.names = TRUE),
								list.files("data-raw/flowcam-summaries-may06", full.names = TRUE),
								list.files("data-raw/flowcam-summaries-may08", full.names = TRUE),
								list.files("data-raw/flowcam-summaries-may12", full.names = TRUE))


names(cell_files) <- cell_files %>% 
	gsub(pattern = ".csv$", replacement = "")


#### Step 3: read in all the files!

all_cells <- map_df(cell_files, read_csv, col_names = FALSE, .id = "file_name")



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

write_csv(TT_cells, "data-processed/TT_cells.csv")
