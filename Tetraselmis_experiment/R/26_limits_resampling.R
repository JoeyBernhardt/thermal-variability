
### get the upper and lower 95% on Tmax

fits <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample.csv")

## ok now take the fits, and make prediction curves and then take the 97.5 and 2.5% CIs

## split up the fits df by curve id
fits_split <- fits %>% 
	filter(a.list > 0, z.list > 0) %>%
	split(.$curve.id.list)

prediction_function <- function(curve1) {
	x <- seq(-15, 38, 0.01)
	predictions <- curve1$a.list[[1]]*exp(curve1$b.list[[1]]*x)*(1-((x-curve1$z.list[[1]])/(curve1$w.list[[1]]/2))^2)
	data.frame(x, predictions)
}

## make predictions for each set of parameter estimates. 
all_predictions <- fits_split %>% 
	map_df(prediction_function) 

all_predictions %>% 
	group_by(x) %>% 
	summarise(q2.5=quantile(predictions, probs=0.025),
						q97.5=quantile(predictions, probs=0.975),
						mean = mean(predictions)) %>%
	# filter(x < 0) %>% 
	ggplot(aes(x = x, y = q97.5)) + geom_point() + geom_hline(yintercept = 0) +
	xlim(-10, 0) + ylim(-0.01, 0.01)
## tmax = 31.25, 32.92
##tmin = -1.25, 




### now for variable

fits <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample_v.csv")
fits_c <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params.csv")
fits_v <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_v.csv")
## ok now take the fits, and make prediction curves and then take the 97.5 and 2.5% CIs

## split up the fits df by curve id
fits_split <- fits_v %>% 
	filter(a.list >0) %>% 
	split(.$curve.id.list)

prediction_function <- function(curve1) {
	x <- seq(-10, 38, 0.01)
	predictions <- curve1$a.list[[1]]*exp(curve1$b.list[[1]]*x)*(1-((x-curve1$z.list[[1]])/(curve1$w.list[[1]]/2))^2)
	data.frame(x, predictions)
}

## make predictions for each set of parameter estimates. 
all_predictions <- fits_split %>% 
	map_df(prediction_function) 

all_predictions %>% 
	filter(x < 0) %>% View

## 31.90, -4.47

all_predictions %>% 
	group_by(x) %>% View
	summarise(q2.5=quantile(predictions, probs=0.025),
						q97.5=quantile(predictions, probs=0.975),
						mean = mean(predictions)) %>%
	filter(x > 28) %>% View


### tmax = 28.66, 32.77


ggplot(aes(x = x, y = q2.5)) + geom_point() + geom_hline(yintercept = 0) +
	xlim(30.5, 34)



### let's just get the actual fitted 0

fits_v <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_v.csv")
fits_c <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params.csv")


