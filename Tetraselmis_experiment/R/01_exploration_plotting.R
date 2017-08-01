### Initial plotting

# load libraries ----------------------------------------------------------

library(tidyverse)
library(purrr)
library(lubridate)
library(stringr)
library(plotrix)
library(broom)



# load data ---------------------------------------------------------------

TT_raw_variable <- read_csv("Tetraselmis_experiment/data-processed/TT_cells.csv")
TT_raw_constant <- read_csv("Tetraselmis_experiment/data-processed/TT_constant_cells.csv")
TT_current <- read_csv("Tetraselmis_experiment/data-processed/TT_constant_cells.csv")

TT_var <- TT_raw_variable %>% 
	separate(replicate, into = c("temperature", "treatment"), sep = "[I]", remove = FALSE) %>% 
	separate(replicate, into = c("temp", "rep"), sep = "[C|V]", remove = FALSE) %>% 
	separate(temperature, into = c("temp_c", "var"), sep = -2) %>% 
	select(-temp_c, -treatment) %>% 
	filter(!grepl("in", sample_day)) %>% 
	mutate(start_time = ymd_hms(start_time)) %>% 
	mutate(var = str_replace(var, "C", "constant temperature")) %>% 
	mutate(var = str_replace(var, "V", "variable temperature")) %>% 
	mutate(total_biovolume = cell_density * cell_volume) %>% 
	separate(start_time, into = c("date", "time"), sep = " ", remove = FALSE) %>% 
	mutate(date = ymd(date))

TT_var %>% 
	filter(date == "2017-02-18") %>% 
	group_by(replicate) %>% 
	summarise_each(funs(mean, std.error), cell_density) %>% View




TT_constant <- TT_raw_constant %>% 
	mutate(start_time = ymd_hms(start_time)) %>% 
	mutate(total_biovolume = cell_density * cell_volume) %>% 
	separate(start_time, into = c("date", "time"), sep = " ", remove = FALSE) %>% 
	mutate(date = ymd(date)) %>% 
	rename(temp = replicate)

TT_constant %>% 
	filter(date == "2017-03-10") %>% 
	group_by(temp) %>% 
	summarise_each(funs(mean, std.error), cell_density) %>% View



TT_constant %>% 
	# filter(date == "2017-02-18") %>% 
	ggplot(aes(x = start_time, y = cell_density)) + geom_point(size = 4) +
	facet_wrap( ~ temp) + theme_bw() + ylab("population abundance (cells/ml)")

TT %>% 
	group_by(temp, date) %>% View
	summarise_each(funs(mean, std.error), cell_density, total_biovolume) %>% 
	ggplot(aes(x = date, y = cell_density_mean, color = var)) + geom_point(size = 4) +
	geom_errorbar(aes(ymin = cell_density_mean - cell_density_std.error, ymax = cell_density_mean + cell_density_std.error), width = 0.2) + 
	facet_wrap( ~ temp) + theme_bw() + ylab("population abundance (cells/ml)")



TT %>% 
	group_by(temp, var, date) %>% 
	summarise_each(funs(mean, std.error), cell_density, total_biovolume) %>% 
	ggplot(aes(x = date, y = total_biovolume_mean, color = var)) + geom_point(size = 4) +
	geom_errorbar(aes(ymin = total_biovolume_mean - total_biovolume_std.error, ymax = total_biovolume_mean + total_biovolume_std.error), width = 0.2) + 
	facet_wrap( ~ temp) + theme_bw() + ylab("total biovolume (um3/ml)")


TT_new_var <- TT_var %>% 
	mutate(keep = NA) %>% 
	mutate(keep = ifelse(date > "2017-02-26", "drop", "keep")) %>%
	filter(keep == "keep" | temp == 10)


TT_new_var %>% 
	filter(temp == 10) %>% View

TT_new %>% 
	group_by(temp, var, date) %>% 
	summarise_each(funs(mean, std.error), cell_density, total_biovolume) %>% 
	ggplot(aes(x = date, y = cell_density_mean, color = var)) + geom_point(size = 4) +
	geom_errorbar(aes(ymin = cell_density_mean - cell_density_std.error, ymax = cell_density_mean + cell_density_std.error), width = 0.2) + 
	facet_wrap( ~ temp) + theme_bw() + ylab("population abundance (cells/ml)")


TT_new %>% 
	ggplot(aes(x = start_time, y = total_biovolume, color = var)) + geom_point(size = 4) +
	facet_wrap( ~ temp) + theme_bw() + ylab("population abundance (cells/ml)")
	
str(TT_new)




TT_new_constant <- TT_constant %>% 
	mutate(keep = NA) %>% 
	mutate(keep = ifelse(date > "2017-03-13", "drop", "keep")) %>%
	filter(keep == "keep" | temp < 29)


TT_new_var$start.time <- ymd("2017-02-18")
TT_new_var$time_since_innoc <- interval(TT_new_var$start.time, TT_new_var$date)
TT_new_constant$start.time <- ymd("2017-03-10")
TT_new_constant$time_since_innoc <- interval(TT_new_constant$start.time, TT_new_constant$date)


TT_new_constant <- TT_new_constant %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))

TT_new_var <- TT_new_var %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))

write_csv(TT_new_constant, "Tetraselmis_experiment/data-processed/TT_new_constant.csv")
write_csv(TT_new_var, "Tetraselmis_experiment/data-processed/TT_new_var.csv")



### trying to get the r's more accurate

TT_new_constant <- read_csv("Tetraselmis_experiment/data-processed/TT_new_constant.csv")
TT_new_var <- read_csv("Tetraselmis_experiment/data-processed/TT_new_var.csv")


TT_new_constant %>% 
	filter(date == "2017-03-10") %>% 
	group_by(temp) %>% 
	summarise(mean_cell_density = mean(cell_density)) %>% View



TT_constant_r <- TT_new_constant %>% 
	group_by(temp) %>% 
	do(tidy(nls(cell_density ~ 2200 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temp = as.numeric(temp)) %>% 
	mutate(var = "constant temperature")

TT_constant_15 <- TT_new_constant %>% 
	rename(temperature = temp) %>% 
	filter(temperature == 15) %>% 
	do(tidy(nls(cell_density ~ 1981.9 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	mutate(var = "constant temperature") %>% 
	mutate(temperature = "15")

TT_constant_19 <- TT_new_constant %>% 
	rename(temperature = temp) %>% 
	filter(temperature == 19) %>% 
	do(tidy(nls(cell_density ~ 1972.4* (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	mutate(var = "constant temperature") %>% 
	mutate(temperature = "19")

TT_constant_29 <- TT_new_constant %>% 
	rename(temperature = temp) %>% 
	filter(temperature == 29) %>% 
	do(tidy(nls(cell_density ~ 2494.0* (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	mutate(var = "constant temperature") %>% 
	mutate(temperature = "29")

TT_constant_35 <- TT_new_constant %>% 
	rename(temperature = temp) %>% 
	filter(temperature == 35) %>% 
	do(tidy(nls(cell_density ~ 2635.6* (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	mutate(var = "constant temperature") %>% 
	mutate(temperature = "35")


all_constant <- bind_rows(TT_constant_15, TT_constant_29, TT_constant_19, TT_constant_35)

all_constant <- all_constant %>% 
	mutate(temperature = as.numeric(temperature))


TT_var_r <- TT_new_var %>% 
	group_by(temp, var) %>% 
	do(tidy(nls(cell_density ~ 1000 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temp))


all_r <- bind_rows(TT_var_r, TT_constant_r)


all_r %>% 
	rename(variability = var) %>% 
	filter(temp > 5) %>% 
	# filter(temp != 24) %>% 
	mutate(temp = as.numeric(temp)) %>% 
	group_by(variability) %>% 
	ggplot(aes(x = temp, y = estimate, color = variability)) + geom_point(size = 2) + 
	geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) + 
	theme_bw() + xlab("temperature (C)") + ylab("intrinsic growth rate (r)") +
	theme(text = element_text(size=12)) + geom_line() + xlim(0,40)
ggsave("Tetraselmis_experiment/figures/var_constant_TPC.png", width = 8, height = 5)



all_r %>% 
	filter(temp > 5, temp < 29) %>% 
	mutate(inverse_temp = (1/(.00008617*(temp+273.15)))) %>%
	filter(var == "constant temperature") %>% 
	ggplot(aes(x = inverse_temp, y = log(estimate))) + geom_point(size = 3) + geom_smooth(method = "lm") + theme_bw()
	do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View


results <- TT_new %>% 
	group_by(temp, var) %>% 
	do(tidy(nls(cell_density ~ 2000 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temp = as.numeric(temp)) %>% 
	rename(variability = var)



# starting fresh with v3 data ---------------------------------------------

v3_27 <- TT_current %>% 
	filter(grepl("27", replicate))

v3_27 <- v3_27 %>% 
	mutate(start_time = ymd_hms(start_time)) %>% 
	mutate(total_biovolume = cell_density * cell_volume) %>% 
	separate(start_time, into = c("date", "time"), sep = " ", remove = FALSE) %>% 
	mutate(date = ymd(date)) %>% 
	rename(temp = replicate)
	

TT_current2 <- TT_current %>% 
	mutate(start_time = ymd_hms(start_time)) %>% 
	mutate(total_biovolume = cell_density * cell_volume) %>% 
	separate(start_time, into = c("date", "time"), sep = " ", remove = FALSE) %>% 
	mutate(date = ymd(date)) %>% 
	rename(temp = replicate)

TT_current2$start.time <- ymd("2017-03-10")


TT_current2$start.time[TT_current2$temp == "27"] <- ymd("2017-03-18")
TT_current2$start.time[TT_current2$temp == "27v"] <- ymd("2017-03-18")
TT_current2$start.time[TT_current2$temp == "20v"] <- ymd("2017-03-25")
TT_current2$start.time[TT_current2$temp == "5v"] <- ymd("2017-03-18")
TT_current2$start.time[TT_current2$temp == "5"] <- ymd("2017-03-10")
TT_current2$start.time[TT_current2$temp %in% c("0", "15", "19", "29", "35")] <- ymd("2017-03-25")



# TT_current3 <- TT_current2 %>%
# 	mutate(start.time = ifelse(temp == "27" | temp == "27v", ymd("2017-03-18"), start.time)) %>%
# 	mutate(start.time = ifelse(temp == "20v", ymd("2017-03-25"), start.time)) %>%
# 	mutate(start.time = ifelse(temp == "5v", ymd("2017-03-25"), start.time)) %>% 
# 	mutate(start.time = ifelse(temp %in% c("0", "5", "15", "19", "29", "35"), ymd("2017-03-10"), start.time))
	# mutate(keep = NA) %>% 
	# mutate(keep = ifelse(date > "2017-03-13", "drop", "keep")) %>%
	# filter(keep == "keep" | temp < 29)


str(TT_current3)

	
v3_27$start.time <- ymd("2017-03-18")
v3_27$time_since_innoc <- interval(v3_27$start.time, v3_27$date)


v3_27 <- v3_27 %>% 
	mutate(time_since_innoc = interval(start.time, date)) %>%
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1)) %>% 
	mutate(treatment = ifelse(grepl("v", temp), "variable temperature", "constant temperature")) %>% 
	mutate(temp = str_replace(temp, "v", ""))



TT_current4 <- TT_current2 %>% 
	mutate(time_since_innoc = interval(start.time, date)) %>%
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1)) %>% 
	mutate(treatment = ifelse(grepl("v", temp), "variable temperature", "constant temperature")) %>% 
	mutate(temp = str_replace(temp, "v", ""))


v20_fake_innoc <- read_csv("Tetraselmis_experiment/data-raw/20v-fake-innoc-data.csv") %>% 
	mutate(temp = as.character(temp)) %>% 
	mutate(treatment = "variable temperature")



TT_current4 <- bind_rows(TT_current4, v20_fake_innoc)

TT_current4 <- TT_current4 %>% 
	filter(!is.na(treatment))

TT_current4 %>% 
	filter(time_since_innoc_hours == 0) %>% 
	group_by(temp) %>% 
	summarise(mean_cell_density = mean(cell_density)) %>% View

	
TT_current4 %>% 
	filter(temp > 19) %>% 
	ggplot(aes(x = time_since_innoc_days, y = cell_density, color = treatment)) + geom_point() +
	facet_wrap( ~ temp)

TT_current4 %>% 
	filter(temp == 27) %>% 
	filter(treatment == "variable temperature") %>% View

r27 <- TT_current4 %>% 
	filter(temp == 27) %>% 
	filter(time_since_innoc_days < 5) %>% 
	group_by(temp, treatment) %>% 
	do(tidy(nls(cell_density ~ 1700 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temp = str_replace(temp, "v", "")) %>% 
	mutate(temp = as.integer(temp)) %>% 
	rename(var = treatment)

r20 <- TT_current4 %>% 
	filter(temp == 20) %>% 
	group_by(temp, treatment) %>% 
	do(tidy(nls(cell_density ~ 2000 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temp = str_replace(temp, "v", "")) %>% 
	mutate(temp = as.integer(temp)) %>% 
	rename(var = treatment)

r15 <- TT_current4 %>% 
	filter(temp == 15) %>% 
	group_by(temp, treatment) %>% 
	do(tidy(nls(cell_density ~ 1980 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temp = str_replace(temp, "v", "")) %>% 
	mutate(temp = as.integer(temp)) %>% 
	rename(var = treatment)

r19 <- TT_current4 %>% 
	filter(temp == 19) %>% 
	group_by(temp, treatment) %>% 
	do(tidy(nls(cell_density ~ 1972 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temp = str_replace(temp, "v", "")) %>% 
	mutate(temp = as.integer(temp)) %>% 
	rename(var = treatment)

r29 <- TT_current4 %>% 
	filter(temp == 29) %>% 
	group_by(temp, treatment) %>% 
	do(tidy(nls(cell_density ~ 2494 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temp = str_replace(temp, "v", "")) %>% 
	mutate(temp = as.integer(temp)) %>% 
	rename(var = treatment)

r35 <- TT_current4 %>% 
	filter(temp == 35) %>% 
	group_by(temp, treatment) %>% 
	do(tidy(nls(cell_density ~ 2635 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temp = str_replace(temp, "v", "")) %>% 
	mutate(temp = as.integer(temp)) %>% 
	rename(var = treatment)

r5 <- TT_current4 %>% 
	filter(temp == 5) %>% 
	# filter(treatment == "constant temperature") %>% View
	group_by(temp, treatment) %>% 
	do(tidy(nls(cell_density ~ 2181 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temp = str_replace(temp, "v", "")) %>% 
	mutate(temp = as.integer(temp)) %>% 
	rename(var = treatment)




v3_27_r <- v3_27 %>% 
	filter(time_since_innoc_days < 5) %>% 
	group_by(temp, treatment) %>% 
	do(tidy(nls(cell_density ~ 1700 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
	mutate(temp = str_replace(temp, "v", "")) %>% 
	mutate(temp = as.integer(temp)) %>% 
		rename(var = treatment)


TT_var_r <- TT_var_r %>% 
	mutate(temp = as.integer(temp))

str(TT_var_r)
all_r<- bind_rows(r20, r27, all_constant, r5, TT_var_r)

all_r <- all_r %>% 
	mutate(temp = ifelse(is.na(temp), temperature, temp))


### add in fake 0 growth rate data -- will fill this in later

r0 <- tribble(
	~temp, ~var, ~term, ~estimate, ~temperature,
	"0", "constant temperature", "a", "0", "0"
)

str(all_r)

r0 <- r0 %>% 
	mutate(temp = as.integer(temp)) %>% 
	mutate(estimate = as.numeric(estimate)) %>% 
	mutate(temperature = as.integer(temperature))

all_r <- bind_rows(all_r, r0)

write_csv(all_r, "Tetraselmis_experiment/data-processed/all_r.csv")


### figure of TPCs
all_r <- read_csv("Tetraselmis_experiment/data-processed/all_r.csv")

r0 <- tribble(
	~temp, ~var, ~term, ~estimate, ~temperature,
	"0", "constant temperature", "a", "0", "0"
)

r0 <- r0 %>% 
	mutate(temp = as.integer(temp)) %>% 
	mutate(estimate = as.numeric(estimate)) %>% 
	mutate(temperature = as.integer(temperature))

all_r <- bind_rows(all_r, r0)

write_csv(all_r, "Tetraselmis_experiment/data-processed/all_r_with0.csv")

all_r %>% 
	# filter(var == "constant temperature") %>% 
	mutate(var = str_replace(var, "constant temperature", "constant")) %>% 
	mutate(var = str_replace(var, "variable temperature", "fluctuating")) %>% 
	rename(`thermal environment` = var) %>% 
	mutate(temp = as.numeric(temp)) %>% 
	group_by(`thermal environment`) %>% 
	ggplot(aes(x = temp, y = estimate, color = `thermal environment`)) + geom_point(size = 4) + 
	geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) + 
	theme_bw() + xlab("temperature (C)") + ylab("intrinsic growth rate (r)") +
	theme(text = element_text(size=18)) + geom_line() + xlim(0,40) + theme_bw() + theme_minimal() +
	scale_colour_manual(values = c("black", "red"))


#### get the activation energy for growth!!
all_r_constant <- all_r %>% 
	# filter(var == "constant temperature") %>% 
	filter(temp < 29) %>% 
	filter(temp > 0) %>% 
	mutate(temp = as.numeric(temp)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temp+273.15))))
	
	
	all_r_constant %>% 
	filter(temp > 5) %>% 
	# filter(temp < 25) %>% 
	group_by(var) %>%
	do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>%
	filter(term != "(Intercept)") %>% View
		ungroup() %>% 
	mutate(var = str_replace(var, "constant temperature", "constant")) %>% 
	mutate(var = str_replace(var, "variable temperature", "fluctuating")) %>%
	rename(`thermal environment` = var) %>% 
	ggplot(aes(x = `thermal environment`, y = estimate)) + geom_point(size = 4) +
		geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) + theme_minimal()

	all_r_constant %>% 
		# filter(temp > 5) %>% 
		ggplot(aes(x = inverse_temp, y = log(estimate), group = var, color = var)) + geom_point() + geom_smooth(method = "lm") +
		scale_x_reverse()
	