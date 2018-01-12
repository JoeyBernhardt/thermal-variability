## data exploration and visualization


library(tidyverse)
library(lubridate)
library(broom)
library(stringr)
library(simecol)
library(purrr)

ldata <- read_csv("Tetraselmis_experiment/data-processed/TT_cells_light.csv")
data5 <- read_csv("Tetraselmis_experiment/data-processed/5C_TT_cells.csv")


str(ldata)
ldata2 <- ldata %>% 
	mutate(replicate = str_replace(replicate, "01", "1")) %>% 
	mutate(replicate = str_replace(replicate, "02", "2")) %>%
	mutate(replicate = str_replace(replicate, "03", "3")) %>% 
	mutate(replicate = str_replace(replicate, "04", "4")) %>% 
	mutate(replicate = str_replace(replicate, "05", "5")) %>% 
	mutate(replicate = str_replace(replicate, "06", "6")) %>%
	mutate(replicate = str_replace(replicate, "07", "7")) %>% 
	mutate(replicate = str_replace(replicate, "08", "8")) %>% 
	mutate(replicate = str_replace(replicate, "09", "9")) %>% 
	mutate(light_level = NA) %>% 
	mutate(light_level = ifelse(replicate %in% c("11", "1", "6", "3"), "high_light", light_level)) %>%
	mutate(light_level = ifelse(replicate %in% c("2", "4", "7", "10"), "med_high_light", light_level)) %>% 
	mutate(light_level = ifelse(replicate %in% c("5", "12", "8", "9"), "med_low_light", light_level)) %>%
	mutate(light_level = ifelse(replicate %in% c("14", "13", "16", "15"), "low_light", light_level)) %>% 
	mutate(start_time = ymd_hms(start_time)) %>% 
	mutate(light = NA) %>% 
	mutate(light = ifelse(light_level == "high_light", 140, light)) %>% 
	mutate(light = ifelse(light_level == "med_high_light", 80, light)) %>% 
	mutate(light = ifelse(light_level == "med_low_light", 50, light)) %>%
	mutate(light = ifelse(light_level == "low_light", 30, light)) %>% 
	separate(start_time, into = c("date", "time"), sep = " ", remove = FALSE) %>% 
	mutate(date = ymd(date))

## get the days in

ldata2$start.time <- ymd("2017-04-30")
ldata2$time_since_innoc <- interval(ldata2$start.time, ldata2$date)


ldata3 <- ldata2 %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))

ldata4 <- ldata3 %>% 
	mutate(keep = ifelse(date > ymd("2017-05-05") & light > 45, "drop", "keep")) %>% 
	select(keep, everything())

ldata5 %>% 
	filter(keep == "keep") %>% 
	ggplot(aes(x = start_time, y = cell_density, color = light_level, group = replicate)) + geom_point(size = 3) +
	facet_wrap( ~ light) + theme_bw() + ylab("population abundance (cells/ml)") + xlab("date") + geom_line()

ldata3 %>% 
	ggplot(aes(x = start_time, y = cell_density, color = light_level, group = replicate)) + geom_point(size = 3) +
	facet_wrap( ~ light) + theme_bw() + ylab("population abundance (cells/ml)") + xlab("date") + geom_line()




growth <- ldata4 %>% 
	filter(keep == "keep") %>% 
	group_by(light_level, replicate) %>% 
	do(tidy(nls(cell_density ~ 1272.4* (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) 



growth_results <- growth %>%
	mutate(light = NA) %>% 
	mutate(light = ifelse(light_level == "high_light", 140, light)) %>% 
	mutate(light = ifelse(light_level == "med_high_light", 80, light)) %>% 
	mutate(light = ifelse(light_level == "med_low_light", 50, light)) %>%
	mutate(light = ifelse(light_level == "low_light", 30, light)) 

growth_results%>% 
	ggplot(aes(x = light, y = estimate)) + geom_point(size = 4) + theme_bw() + ylab("growth rate (r) / day") + xlab("light intensity (umols/m2/s)") +
theme(text = element_text(size=18))


library(drc)

model.drm <- drm(estimate ~ light, data = growth_results, fct = MM.2())

mml1 <- data.frame(light = seq(0, max(growth_results$light), length.out = 100))
mml1$v <- predict(model.drm, newdata = mml1)


ggplot(growth_results, aes(x = light, y = estimate)) +
	theme_bw() +
	xlab("Light [umol/m2/s]") +
	ylab("Growth rate [/day]") +
	geom_point(alpha = 0.5, size = 2) +
	geom_line(data = mml1, aes(x = light, y = v), colour = "red")

summary(model.drm)


### try something different

model.nls <- nls(estimate ~ Vm * light/(K+light), data = growth_results, 
								 start = list(K = max(growth_results$light)/2, Vm = max(growth_results$estimate)))

mml <- data.frame(light = seq(0, 150, length.out = 50))
mml$estimate <- predict(model.nls, newdata = mml)

ggplot(growth_results, aes(x = light, y = estimate)) +
	theme_bw() +
	xlab("light") +
	ylab("growth rate") +
	geom_point(alpha = 0.5, size = 3) +
	geom_line(data = mml, aes(x = light, y = estimate), colour = "blue")

summary(model.nls)

### other approaches

library(phytotools)

PI <- fitEP(growth_results$light, growth_results$estimate, normalize = FALSE, lowerlim = c(0, 0, 0), upperlim = c(100, 2000, 2000),
			fitmethod=c("Nelder-Mead"))
PI$eopt

plot(growth_results$light, growth_results$estimate, xlim=c(0,150), ylim=c(0,1), xlab="light (photons)", ylab="Growth rate")
#Add model fit
E <- seq(0,150,by=1)
with(PI,{
	P <- E/((1/(alpha[1]*eopt[1]^2))*E^2+(1/ps[1]-2/(alpha[1]*eopt[1]))*E+(1/alpha[1]))
	lines(E,P)
})


### Now fit K

Parameters <- c(r = 1, K = 10 ^ 7) ## initial parameter guesses
CRmodel <- new("odeModel",
							 main = function (time, init, parms) {
							 	with(as.list(c(init, parms)), {
							 		dp <-  r * P * (1 - (P / K))
							 		list(c(dp))
							 	})
							 },
							 parms = Parameters,
							 times = c(from = 0, to = 6, by = 0.1), # the time interval over which the model will be simulated.
							 init = c(P = 2000),
							 # init = c(P = 971410.2), # starting biovolume
							 solver = "lsoda" #lsoda will be called with tolerances of 1e-9. Default tolerances are both 1e-6. Lower is more accurate.
)

fittedparms <- c("r", "K") # for assigning fitted parameter values to fittedCRmodel

controlfit <- function(data){
	
	init(CRmodel) <- c(P = data$P[1]) # Set initial model conditions to the biovolume taken from the first measurement day
	obstime <- data$days # The X values of the observed data points we are fitting our model to
	yobs <- select(data, P) # The Y values of the observed data points we are fitting our model to
	
	
	fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
															 debuglevel = 0, fn = ssqOdeModel,
															 method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
															 control = list(trace = T)
	)
	
	r <- coef(fittedCRmodel)[1]
	K <- coef(fittedCRmodel)[2]
	ID <- data$ID[1]
	output <- data.frame(ID, r, K)
	return(output)
}

plotsinglefit <- function(data){
	
	init(CRmodel) <- c(P = data$P[1]) # Set initial model conditions to the biovolume taken from the first measurement day
	obstime <- data$days # The X values of the observed data points we are fitting our model to
	yobs <- select(data, P) # The Y values of the observed data points we are fitting our model to
	# parms(CRmodel)[TempName] <- data$temp[1] # Set the temperature parameter in CRmodel to whatever our control replicate used.
	
	# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
	# squared differences between the experimental data and our modelled data. 
	
	# The PORT algorithm is used for the model fitting, analogous to O'Connor et al.
	# "lower" is a vector containing the lower bound constraints
	# for the parameter values.
	
	fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
															 debuglevel = 0, fn = ssqOdeModel,
															 method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
															 control = list(trace = T)
	)
	
	# To display the fitted results we need to create a new OdeModel object. Here
	# we duplicate CRmodel and then alter it to use our new fitted parameters.
	plotfittedCRmodel <- CRmodel
	parms(plotfittedCRmodel)[fittedparms] <- coef(fittedCRmodel)
	
	# set model parameters to fitted values and simulate again
	times(plotfittedCRmodel) <- c(from=0, to=20, by=1)
	ysim <- out(sim(plotfittedCRmodel, rtol = 1e-9, atol = 1e-9))
	
	# Form observed data into a dataframe; the simulated data are already in a dataframe
	observeddata <- data.frame(obstime, yobs)
	observeddata$type <- "observed"
	observeddata <- rename(observeddata, time = obstime)
	simulateddata <- ysim
	simulateddata$type <- "simulated"
	output <- bind_rows(observeddata, simulateddata) ## plop everything together into one data frame
	
	
	# Alternative: Plot the results of our model fitting directly.
	# biol_plot <- ggplot() +
	# 	geom_point(data = observeddata, aes(x = obstime, y = yobs, color = "observed")) + # Observed data are points
	# 	geom_line(data = simulateddata, aes(x = time, y = P, color = "simulated")) + # Simulated data are in a continuous line
	# 	labs(x = "Time (days)", y = "Algal Biovolume")

	# Output the results as a ggplot2 object
	# output <- biol_plot
	return(output)
}


Parameters <- c(r = 0.5, K = 10 ^ 5)
# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.01, K = 10 ^ 2)
UpperBound <- c(r = 10, K = 25000) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound


l_140 <- ldata4 %>% 
	filter(keep == "keep") %>% 
	filter(light == 140) %>% 
	mutate(time_since_innoc_days = as.numeric(time_since_innoc_days)) %>% 
	mutate(time_since_innoc_hours = as.numeric(time_since_innoc_hours)) %>% 
	mutate(time_since_innoc = as.numeric(time_since_innoc)) %>% 
	rename(P = cell_density) %>% 
	rename(days = time_since_innoc_days) %>% 
	select(P, days) %>% 
	mutate(ID = 140)


controlfit(l_140)


ldata5 <- ldata3 %>% 
	mutate(keep = ifelse(date > ymd("2017-05-07") & light > 45, "drop", "keep")) %>% 
	select(keep, everything())

ldata5 %>% 
	filter(keep == "keep") %>% 
	ggplot(aes(x = start_time, y = cell_density, color = light_level, group = replicate)) + geom_point(size = 3) +
	facet_wrap( ~ light) + theme_bw() + ylab("population abundance (cells/ml)") + xlab("date") + geom_line()


growth <- ldata5 %>% 
	# filter(keep == "keep") %>% 
	mutate(time_since_innoc_days = as.numeric(time_since_innoc_days)) %>% 
	mutate(time_since_innoc_hours = as.numeric(time_since_innoc_hours)) %>% 
	mutate(time_since_innoc = as.numeric(time_since_innoc)) %>% 
	# filter(light_level != "low_light") %>% 
	rename(P = cell_density) %>% 
	rename(days = time_since_innoc_days) %>% 
	rename(ID = replicate) %>% 
	select(P, days, ID) %>% 
	filter(days < 7)

growth_split <- growth %>% 
	split(.$ID)

results <- growth_split %>% 
	map_df(controlfit)

growth %>% 
	filter(ID == 10) %>% 
plotsinglefit %>% 
	ggplot(aes(x = time, y = P, color = type)) + geom_point() + theme_bw()

## 5 c("14", "13", "16", "15")

results %>% 
	rename(replicate = ID) %>% 
	filter(replicate != "12") %>%
	filter(replicate != "2") %>% 
	mutate(light_level = NA) %>% 
mutate(light_level = ifelse(replicate %in% c("11", "1", "6", "3"), "high_light", light_level)) %>% 
	mutate(light_level = ifelse(replicate %in% c("2", "4", "7", "10"), "med_high_light", light_level)) %>% 
	mutate(light_level = ifelse(replicate %in% c("5", "12", "8", "9"), "med_low_light", light_level)) %>%
	mutate(light_level = ifelse(replicate %in% c("14", "13", "16", "15"), "low_light", light_level)) %>% 
	mutate(light = NA) %>% 
	mutate(light = ifelse(light_level == "high_light", 140, light)) %>% 
	mutate(light = ifelse(light_level == "med_high_light", 80, light)) %>% 
	mutate(light = ifelse(light_level == "med_low_light", 50, light)) %>%
	mutate(light = ifelse(light_level == "low_light", 30, light)) %>% 
	filter(K < 50000) %>% 
ggplot(aes(x = light, y = K)) + geom_point(alpha = 0.5, size = 3) + theme_bw()

