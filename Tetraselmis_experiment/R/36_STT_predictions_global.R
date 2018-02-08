
## use STT to make new TPCs for all isolates

tpc<-function(x){
  res<-all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)
  res
}

## step 2 
variable_predictions <- function(x) {
  y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(all$sst_sd^2)
}


predict_function <- function(data) {
  all <- data
  x <- seq(-2, 40, by = 0.1)
  tpc<-function(x){
    res<-all$a[1]*exp(all$b[1]*x)*(1-((x-all$z[1])/(all$w[1]/2))^2)
    res
  }
  ## step 2 
  variable_predictions <- function(x) {
    y <- tpc(x) + derivative(f = tpc, x = x, order = 2)*0.5*(all$sst_sd^2)
  }
  
  predicted_growth_variable <- sapply(x, variable_predictions)
  predicted_growth_variable2 <- data.frame(x, predicted_growth_variable) %>% 
    rename(temperature = x, 
           growth.rate = predicted_growth_variable) 
}


thomas_temps <- left_join(thomas3, sst_summary)

all_split <- thomas_temps %>% 
  split(.$isolate.code)

STT_predictions <- all_split %>% 
  map_df(predict_function, .id = "isolate.code")

write_csv(STT_predictions, "Tetraselmis_experiment/data-processed/STT_predictions.csv")
