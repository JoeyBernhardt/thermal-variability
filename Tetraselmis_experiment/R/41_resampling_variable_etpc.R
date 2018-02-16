

# Load packages -----------------------------------------------------------

library(purrr)
library(tidyverse)
library(minpack.lm)
library(broom)


# read in data ------------------------------------------------------------

cells_days_v <- read_csv("tvar/cells_days_v.csv")


avals<-seq(0,0.5,0.07)
bvals<-seq(0,0.16,0.07)
zvals<-seq(10,15,1)
wvals<-seq(33,37,1)


df <-  expand.grid(a = avals, b = bvals, z = zvals, w = wvals) %>% 
  mutate(unique_id = rownames(.)) %>% 
  mutate(sample_size = 1)


# set up resampling function ----------------------------------------------


resample_variable <- function(sample_size){
  cd <- cells_days_v %>% 
    group_by(temp, sample_group) %>% 
    sample_n(size = sample_size, replace = FALSE) 

  
  fit_growth <- function(df){
    res <- try(nlsLM(cell_density ~ 800 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
                     data= cd,  
                     start=list(z=df$z[[1]],w=df$w[[1]],a= df$a[[1]], b=df$b[[1]]),
                     lower = c(0, 0, -0.2, 0),
                     upper = c(20, 80, 0.7, 1.5),
                     control = nls.control(maxiter=1024, minFactor=1/204800000)))
    if(class(res)!="try-error"){
      out1 <- tidy(res) %>% 
        select(estimate, term) %>% 
        spread(key = term, value = estimate)
      out2 <- glance(res)
    }
    all <- bind_cols(out1, out2)
    all
  }
  
  
  df_split <- df %>% 
    split(.$unique_id)
  
  
  output <- df_split %>%
    map_df(fit_growth, .id = "run") %>% 
    filter(df.residual > 5) %>% 
    # filter(a != 1.2, w != 80, w != 0, a!= -0.2, z != 0, z != 30) %>% 
    top_n(n = -1, wt = AIC)
  
  return(output)
}


samples <- rep(1, 10)

bootstrap_time_series_fitsv <- samples %>% 
  map_df(resample_variable, .id = "replicate")

write_csv(bootstrap_time_series_fitsv, "bootstrap_time_series_fitsv.csv")


