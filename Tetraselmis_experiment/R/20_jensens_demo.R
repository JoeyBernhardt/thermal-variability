
library(tidyverse)
library(purrr)
library(gganimate)

set.seed(101)
temps <- data.frame(temperature = runif(50, min = 12, max = 27), time = seq_along(1:50), time_point = seq_along(1:50))
mean(temps$temperature)

ggplot(aes(x = time, y = temperature), data = temps) + geom_line() +
	theme_classic()

sd(temps$temperature)

hist(temps$temperature)

jensens <- read_csv("Tetraselmis_experiment/data-processed/jensens_demo.csv")

temps <- jensens %>% 
	select(time, temperature)


all_thermal_data <- read_csv("Tetraselmis_experiment/data-processed/all_thermal_data.csv")

curve <- all_thermal_data %>% 
	filter(isolate.code == 550)

df <- filter(temps, time == 1)
tpc1<-function(df){
	x <- df$temperature
	res<-curve$a[1]*exp(0.30*x)*(1-((x-curve$z[1])/(curve$w[1]/2))^2)
	output <- data.frame(time = df$time, predicted_r = res, temperature = df$temperature)
}

tpc_draw<-function(x){
	res<-curve$a[1]*exp(0.30*x)*(1-((x-curve$z[1])/(curve$w[1]/2))^2)
	res
}

temp_split <- temps %>% 
	split(.$time)



predicted_growth <- temp_split %>% 
	map_df(tpc1)

write_csv(predicted_growth, "Tetraselmis_experiment/data-processed/jensens_demo.csv")

predicted_growth <- read_csv("Tetraselmis_experiment/data-processed/jensens_demo.csv")


ggplot(aes(x = time, y = temperature), data = predicted_growth) + geom_line() +
	theme_classic() +
	ylab("Temperature (°C)") +
	theme(text = element_text(size=16, family = "Helvetica")) + xlab("Time") +
	geom_hline(yintercept = 19.49, color = "blue")
ggsave("Tetraselmis_experiment/figures/fake_time_series.png", width = 7, height = 3)


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + stat_function(fun = tpc_draw, color = "black", size = 1) + xlim(0, 33)


mean(predicted_growth$temperature)

q <- ggplot(predicted_growth, aes(x = temperature, y = predicted_r, frame = time)) +
	geom_point(color = "blue", size = 8) + theme_classic() +
	stat_function(fun = tpc_draw, color = "black", size = 1) + xlim(0, 33) + ylim(0, 1.5) +
	ylab("Growth rate") + xlab("Temperature (°C)") +
	theme(text = element_text(size=16, family = "Helvetica")) +
	geom_vline(xintercept = 19.49, color = "grey", linetype = "dotted", size = 2) +
	geom_hline(yintercept = 0.389, linetype = "dotted")

gganimate(q, interval = 0.5, filename = "r_time.gif", title_frame = FALSE, ani.width = 500, ani.height = 250)


time <- ggplot(temps, aes(x = time, y = temperature, frame = time)) +
	geom_point(color = "blue", size = 8) +
	# geom_vline(aes(xintercept = temperature, frame = temperature), data = temps, lty = 2, color = "red") +
	theme_classic() +
	geom_hline(yintercept = 11.53, color = "grey", linetype = "dotted", size = 2) + 
	geom_line(aes(frame = time, cumulative = TRUE), color = "blue", data = temps) +
	ylab("Temperature (°C)") + xlab("Time")
	theme(text = element_text(size=24, family = "Helvetica"))
	
gganimate(time, interval = 0.5, filename = "temp_time.gif", title_frame = FALSE, ani.width = 500, ani.height = 250)


preds <- predicted_growth %>% 
	mutate(running_avg = cummean(predicted_r)) %>% 
	mutate(running_avg_temp = cummean(temperature))

r <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
r +	geom_segment(aes())
	
average <- ggplot(data = preds) + 
geom_segment(aes(y = running_avg, yend = running_avg, x = 0, xend = 2, frame = time), size = 2, color = "blue") +
	theme_classic() +
	# geom_line(aes(frame = time, cumulative = TRUE), color = "blue", data = temps) +
	ylab("Growth rate") + xlab("") +
	scale_y_continuous(limits = c(0,1.5)) +
	geom_hline(yintercept = 0.389, linetype = "dotted") + 
theme(text = element_text(size=16, family = "Helvetica")) +
	theme(
		axis.text.x = element_blank(),
		axis.ticks = element_blank())

gganimate(average, interval = 0.5, filename = "average_time.gif", title_frame = FALSE, ani.width = 250, ani.height = 250)


average_temp <- ggplot(data = preds) + 
	geom_segment(aes(y = running_avg_temp, yend = running_avg_temp, x = 0, xend = 2, frame = time), size = 2, color = "blue") +
	theme_classic() +
	# geom_line(aes(frame = time, cumulative = TRUE), color = "blue", data = temps) +
	ylab("Temperature (°C)") + xlab("") +
	scale_y_continuous(limits = c(5,20)) +
	theme(text = element_text(size=16, family = "Helvetica")) +
	theme(
		axis.text.x = element_blank(),
		axis.ticks = element_blank()) +
	geom_hline(yintercept = 11.59355, linetype = "dotted") 
	

gganimate(average_temp, interval = 0.5, filename = "average_temp_time.gif", title_frame = FALSE, ani.width = 250, ani.height = 250)

### experiment temperatures
library(lubridate)
library(janitor)
temp24 <- read_csv("Tetraselmis_experiment/data-raw/ibuttons/24-variable-feb24.csv", skip = 14) %>% 
	clean_names()


temp24 %>% 
	mutate(date_time = mdy_hms(date_time)) %>% 
	ggplot(aes(x = date_time, y = value)) + geom_line() +
	theme_classic() +
	theme(text = element_text(size=14, family = "Helvetica")) +
	ylab("Temperature (°C)") + xlab("Time") +
	geom_hline(yintercept = 24, color = "blue") 
	ggsave("Tetraselmis_experiment/figures/temperature_ibutton_time_series.png", width = 7, height = 3)
	
