

library(purrr)
library(bbmle)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(minpack.lm)
library(broom)


# get data in order -------------------------------------------------------


growth_all <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling.csv")
growth_all_v <- read_csv("Tetraselmis_experiment/data-processed/growth_resampling_v.csv")

cells_exp <- read_csv("Tetraselmis_experiment/data-processed/cells_exp.csv")
cells_exp_v <- read_csv("Tetraselmis_experiment/data-processed/cells_v_exp.csv")
cells_mod <- read_csv("Tetraselmis_experiment/data-processed/cells_exp_mod.csv")
cells_days_v_mod <- read_csv("Tetraselmis_experiment/data-processed/cells_days_v_mod.csv")


cells_days_v1 <- cells_exp_v %>% 
  mutate(days = time_since_innoc_hours/24) %>% 
  mutate(days = ifelse(days < 0, 0, days)) %>% 
  mutate(time_point = trunc(time_since_innoc_hours)) %>% 
  arrange(temp, time_point)
write_csv(cells_days_v1, "Tetraselmis_experiment/data-processed/cells_days_v.csv")

fits_vals <- fits_real_constant %>% 
  summarise_each(funs(mean, sd), z, w, a, b)


cells_days <- cells_mod %>% 
  mutate(days = time_since_innoc_hours/24) %>% 
  mutate(days = ifelse(days < 0, 0, days))

cells_days_v <- cells_days_v_mod %>% 
  mutate(days = time_since_innoc_hours/24) %>% 
  mutate(days = ifelse(days < 0, 0, days))


cells_days %>% 
  filter(temp == 27) %>% 
  ggplot(aes(x = time_since_innoc_hours, y = cell_density, color = factor(temp))) + geom_point() 


# try out fitting all data together ---------------------------------------



fit_c <- nlsLM(cell_density ~ 800 * (1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^(time_since_innoc_hours),
    data= cells_days,  start=list(z=16,w=36,a=0.2,b=0.1),
    control = nls.control(maxiter=1000, minFactor=1/204800000))

avals<-seq(-0.2,1.2,0.02)
bvals<-seq(0,0.3,0.02)

cd <- cells_days %>% 
mutate(time_point = trunc(time_since_innoc_hours)) %>% 
  group_by(temp, time_point) %>% 
  sample_n(size = 1, replace = FALSE) 


df <-  expand.grid(a = avals, b = bvals) %>% 
  mutate(unique_id = rownames(.)) %>% 
  mutate(sample_size = 1)
  

(fit_d <- nlsLM(cell_density ~ 800 * (1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^(days),
               data= cells_days,  
               start=list(z=15,w=35,a=0.23, b=0.1),
               lower = c(0, 0, -0.2, -0.2),
               upper = c(30, 40, 1.2, 0.3),
               control = nls.control(maxiter=1024, minFactor=1/204800000)))


(fit_v <- nlsLM(cell_density ~ 800 * (1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^(days),
                data= cells_days_v,  
                start=list(z=15,w=35,a=0.22, b=0.1),
                lower = c(0, 0, -0.2, -0.2),
                upper = c(30, 40, 1.2, 0.3),
                control = nls.control(maxiter=1024, minFactor=1/204800000)))


# try this with a more systematic search through starting values ----------

avals<-seq(-0.2,1.2,1)
bvals<-seq(-0.2,0.3,0.1)
zvals<-seq(12,17,1)
w.guess <- 30
mod.list<-list()
AIC.list<-c()


for(ia in 1:length(avals)){
  for(ib in 1:length(bvals)){
    for(iz in 1:length(zvals)){
      a.guess<-avals[ia]
      b.guess<-bvals[ib]
      z.guess<-zvals[iz]
      res2<-nlsLM(cell_density ~ 800 * (1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^(days),
                           data= cells_days,  
                           start=list(z=z.guess,w=w.guess,a=a.guess, b=b.guess),
                           lower = c(0, 0, -0.2, -0.2),
                           upper = c(30, 40, 1.2, 0.3),
                           control = nls.control(maxiter=1024, minFactor=1/204800000))
      if(class(res2)!="try-error"){
        out1 <- tidy(res2) %>% 
          select(estimate, term) %>% 
          spread(key = term, value = estimate)
        out2 <- glance(res2)
      }
      all <- bind_cols(out1, out2)
      all <- bind_rows(all)
      }
    }
  }

mod.list[1]

if(!is.null(AIC.list)){
  bestmodind<-which(AIC.list==min(AIC.list))
  if(length(bestmodind)>1){
    bestmodind<-sample(bestmodind,1)
  }
  bestmod<-mod.list[[bestmodind]]
  cfs<-coef(bestmod)
}

class(bestmod)

mod.list[1]$getPars()

fitd1 <- tidy(fit_d)
fitv1 <- tidy(fit_v)


tidy(fit_v)
glance(fit_d)


cd <- cells_days %>% 
  mutate(time_point = trunc(time_since_innoc_hours)) %>% 
  group_by(temp, time_point) %>% 
  sample_n(size = 1, replace = FALSE) 


# fit the constant experimental data --------------------------------------

avals<-seq(-0.2,1.2,0.02)
bvals<-seq(0,0.3,0.02)
zvals<-seq(10,15,1)


df <-  expand.grid(a = avals, b = bvals, z = zvals) %>% 
  mutate(unique_id = rownames(.)) %>% 
  mutate(sample_size = 1)

resample <- function(sample_size){
  cd <- cells_days %>% 
    group_by(temp, sample_group) %>% 
    sample_n(size = sample_size, replace = FALSE) 
  

 fit_growth <- function(df){
  res <- try(nlsLM(cell_density ~ 800 * (1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^(days),
                   data= cd,  
                   start=list(z=df$z[[1]],w=35,a= df$a[[1]], b=df$b[[1]]),
                   lower = c(0, 0, -0.2, -0.2),
                   upper = c(30, 80, 1.2, 0.3),
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
  top_n(n = -1, wt = AIC)

return(output)
}


samples <- rep(1, 10) 

bootstrap_time_series_fits <- samples %>% ### this step gets us the replication
  map_df(resample, .id = "replicate")


write_csv(bootstrap_time_series_fits, "Tetraselmis_experiment/data-processed/bootstrap_time_series_fitsb.csv")




# Resample the variable ---------------------------------------------------

avals<-seq(-0.2,1.2,0.02)
bvals<-seq(0,0.3,0.02)
zvals<-seq(10,15,1)


df <-  expand.grid(a = avals, b = bvals, z = zvals) %>% 
  mutate(unique_id = rownames(.)) %>% 
  mutate(sample_size = 1)


resample_variable <- function(sample_size){
  cd <- cells_days_v_mod %>% 
    group_by(temp, sample_group) %>% 
    sample_n(size = sample_size, replace = FALSE) 
    # group_by(temp, sample_group, cycle) %>% 
    # summarise_each(funs(mean), cell_density)
  
  fit_growth <- function(df){
    res <- try(nlsLM(cell_density ~ 800 * (1+(a*exp(b*temp)*(1-((temp-z)/(w/2))^2)))^(days),
                     data= cd,  
                     start=list(z=df$z[[1]],w=35,a= df$a[[1]], b=df$b[[1]]),
                     lower = c(0, 0, -0.2, -0.2),
                     upper = c(30, 80, 5, 0.5),
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


samples <- rep(1, 100)

bootstrap_time_series_fitsv <- samples %>% 
  map_df(resample_variable, .id = "replicate")

write_csv(bootstrap_time_series_fitsv, "Tetraselmis_experiment/data-processed/bootstrap_time_series_fitsv.csv")


# other stuff and plotting ------------------------------------------------

sums <- bootstrap_time_series_fitsb %>% 
  filter(a != 1.2) %>% 
  summarise_each(funs(mean), z, w, a, b)

bind_rows(sums, fits_vals) %>% View


cd <- cells_days %>% 
  mutate(time_point = trunc(time_since_innoc_hours)) %>% 
  group_by(temp, time_point) %>% 
  sample_n(size = 1, replace = FALSE) 


  

bind_rows(tp1, tp2, tp3) %>% View

cfd <- (coef(fit_d))

tpc1<-function(x){
  res<-(cfd[["a"]]*exp(cfd[["b"]]*x)*(1-((x-cfd[["z"]])/(cfd[["w"]]/2))^2))
  res
}

tpc2<-function(x){
  res<-(tp$a[[1]]*exp(tp$b[[1]]*x)*(1-((x-tp$z[[1]])/(tp$w[[1]]/2))^2))
  res
}

k <- sums$a[[1]]/2

tpc3<-function(x){
  res<-(sums$a[[1]]*exp(sums$b[[1]]*x)*(1-((x-sums$z[[1]])/(sums$w[[1]]/2))^2))
  res
}

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + stat_function(fun = tpc3, color = "black", size = 1) +xlim(-10, 33) + 
  ylim(-2, 5) + geom_hline(yintercept = 0) +
  stat_function(fun = curve_constant_resamp, color = ic[3], size = 1.5, alpha = 0.7)
  
nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

nbcurve1<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

tp <- results

grfunc<-function(x){
  -nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]])
}


optinfo<-optim(c(x=tp$z[[1]]),grfunc)
opt <-optinfo$par[[1]]
maxgrowth <- -optinfo$value

library(rootSolve)
uniroot.all(function(x) nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]]),c(opt,150))
uniroot.all(function(x) nbcurve(x, z = tp$z[[1]],w = tp$w[[1]],a = tp$a[[1]],b = tp$b[[1]]),c(-15,opt))


tidy(fit_d) %>% View
tidy(fit_d) %>% View

days <- 2
temp <- 27
growth_fun_27 <- function(x, temp = 27){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}

growth_fun_5 <- function(x, temp = 5){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}

growth_fun_10 <- function(x, temp = 10){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}
growth_fun_24 <- function(x, temp = 24){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}

growth_fun_16 <- function(x, temp = 16){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}

growth_fun_29 <- function(x, temp = 29){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}

growth_fun_32 <- function(x, temp = 32){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}
growth_fun_0 <- function(x, temp = 0){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}
growth_fun_20 <- function(x, temp = 20){
  res <- 800 *  (1+(cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)))^(x)
  res
}


library(viridis)
library(colormap)
library(cowplot)
ic <- colormap(colormap = colormaps$viridis, nshades = 9, format = "hex",
               alpha = 1, reverse = FALSE)

cells_days %>% 
  ggplot(aes(x = days, y = cell_density, color = factor(temp))) + geom_point(size = 2) +
  geom_point(size = 2, shape= 1, color = "black") +
  stat_function(fun = growth_fun_27, color = ic[7], size = 1) +ylim(0, 30000) +
  stat_function(fun = growth_fun_5, color = ic[2], size = 1) +
  stat_function(fun = growth_fun_10, color = ic[3], size = 1) +
  stat_function(fun = growth_fun_24, color = ic[6], size = 1) +
  stat_function(fun = growth_fun_16, color = ic[4], size = 1) +
  stat_function(fun = growth_fun_29, color = ic[8], size = 1) +
  stat_function(fun = growth_fun_32, color = ic[9], size = 1) +
  stat_function(fun = growth_fun_0, color = ic[1], size = 1) +
  stat_function(fun = growth_fun_20, color = ic[5], size = 1) +
   scale_color_viridis(discrete = TRUE) +
  ylab("Cell density") + xlab("Days")



cf <- (coef(fit_c))


expected<-nbcurve(dat$temperature,cfd[[1]],cfd[[2]],cfd[[3]],cfd[[4]])
rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)


fits_c$a.list[1]*exp(fits_c$b.list[1]*x)*(1-((x-fits_c$z.list[1])/(fits_c$w.list[1]/2))^2)

tpch<-function(x){
  res<-(cf[["a"]]*exp(cf[["b"]]*x)*(1-((x-cf[["z"]])/(cf[["w"]]/2))^2))*24
  res
}



p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + stat_function(fun = tpch, color = "black", size = 1) +xlim(-2, 33) + 
  ylim(-2, 1.6)

temp <- 24
r <- cfd[["a"]]*exp(cfd[["b"]]*temp)*(1-((temp-cfd[["z"]])/(cfd[["w"]]/2))^2)

