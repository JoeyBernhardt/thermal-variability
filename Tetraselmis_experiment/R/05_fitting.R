library(tidyverse)
library(simecol)
library(broom)

TTg <- read_csv("Tetraselmis_experiment/data-processed/TTg.csv")
TTg2 <- TTg %>% 
	filter(hours < 49)


tt_sub <- TTg %>% 
	filter(temp == "5", variability == "v") %>% 
	rename(P = cell_density) %>% 
	rename(days = hours) %>% 
	arrange(days) %>% 
	mutate(ID = "5v")

end_time <- max(tt_sub$days)


Parameters <- c(r = 1, K = 10 ^ 4) ## initial parameter guesses
CRmodel <- new("odeModel",
							 main = function (time, init, parms) {
							 	with(as.list(c(init, parms)), {
							 		dp <-  r * P * (1 - (P / K))
							 		list(c(dp))
							 	})
							 },
							 parms = Parameters,
							 times = c(from = 0, to = end_time, by = 0.1), # the time interval over which the model will be simulated.
							 init = c(P = 1000),
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

	fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
															 debuglevel = 0, fn = ssqOdeModel,
															 method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
															 control = list(trace = T)
	)
	
	plotfittedCRmodel <- CRmodel
	parms(plotfittedCRmodel)[fittedparms] <- coef(fittedCRmodel)
	
	# set model parameters to fitted values and simulate again
	times(plotfittedCRmodel) <- c(from=0, to=end_time, by=1)
	ysim <- out(sim(plotfittedCRmodel, rtol = 1e-9, atol = 1e-9))
	
	# Form observed data into a dataframe; the simulated data are already in a dataframe
	observeddata <- data.frame(obstime, yobs)
	observeddata$type <- "observed"
	observeddata <- rename(observeddata, time = obstime)
	simulateddata <- ysim
	simulateddata$type <- "simulated"
	output <- bind_rows(observeddata, simulateddata)
	return(output)
}
Parameters <- c(r = 0.5, K = 10 ^ 5)
LowerBound <- c(r = 0.001, K = 10.5 ^ 4)
UpperBound <- c(r = 100, K = 400000) 
ParamScaling <- 1 / UpperBound

# output_TT <- controldata_TT %>% 
	# map_df(controlfit)

output_plot <- plotsinglefit(tt_sub)

output_plot %>% 
	ggplot(aes(x = time, y = P, color = type)) + geom_point() + theme_bw()


parameters10c <- controlfit(tt_sub)
parameters10v <- controlfit(tt_sub)
parameters20c <- controlfit(tt_sub)
parameters20v <- controlfit(tt_sub)
parameters27v <- controlfit(tt_sub)
parameters27c <- controlfit(tt_sub)
parameters5c <- controlfit(tt_sub)
parameters5v <- controlfit(tt_sub)


all <- bind_rows(parameters10c, parameters10v, parameters20c, parameters20v, parameters27v, parameters27c, parameters5c, parameters5v)

all %>% 
	separate(ID, into = c("temperature", "variability"), sep = -2) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 100000) %>% 
	ggplot(aes(x = temperature, y = K, color = variability)) + geom_point()


all %>% 
	separate(ID, into = c("temperature", "variability"), sep = -2) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(K < 100000) %>% 
mutate(inverse_temp = (-1/(.00008617*(temperature+273.15)))) %>%
	do(tidy(lm(log(K) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View


all %>% 
	separate(ID, into = c("temperature", "variability"), sep = -2) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	filter(temperature < 25) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	do(tidy(lm(log(r) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View
