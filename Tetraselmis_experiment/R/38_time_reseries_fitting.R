
cells_exp <- read_csv("Tetraselmis_experiment/data-processed/cells_exp.csv")
cells_v_exp <- read_csv("Tetraselmis_experiment/data-processed/cells_v_exp.csv")


avals<-seq(-0.2,1.2,0.02)
bvals<-seq(-0.2,0.3,0.02)


df <-  expand.grid(a = avals, b = bvals) %>% 
  mutate(unique_id = rownames(.)) %>% 
  mutate(sample_size = 1)


cells_days <- cells_exp %>% 
  mutate(days = time_since_innoc_hours/24) %>% 
  mutate(days = ifelse(days < 0, 0, days))

fit_growth <- function(df){
  res <- try(nlsLM(cell_density ~ 800 * (1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^(days),
                   data= cells_days  
                   start=list(z=14.4,w=35,a= df$a[[1]], b=df$b[[1]]),
                   lower = c(0, 0, -0.2, -0.2),
                   upper = c(30, 45, 1.2, 0.3),
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
  map_df(fit_growth, .id = "run") 


results <- output %>% 
  filter(df.residual > 5) %>%
  top_n(n = -1, wt = AIC)


write_csv(results, "Tetraselmis_experiment/data-processed/time_resampling_fits.csv")


tp <- results

grfunc<-function(x){
  -nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]])
}

### find limits
optinfo<-optim(c(x=tp$z[[1]]),grfunc)
opt <-optinfo$par[[1]]
maxgrowth <- -optinfo$value

library(rootSolve)
uniroot.all(function(x) nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]]),c(opt,150))
uniroot.all(function(x) nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]]),c(-15,opt))

