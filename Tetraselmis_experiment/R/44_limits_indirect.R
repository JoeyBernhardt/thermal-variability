

params_raw <- read_csv("Tetraselmis_experiment/data-processed/resampling_TPC_params_exp.csv")

nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

## get Tmax
tmax_c <- uniroot.all(function(x) nbcurve(x, z = c_fits$z.list[[1]],w = c_fits$w.list[[1]],a = c_fits$a.list[[1]], b = c_fits$b.list[[1]]),c(10,150))
tmax_v <- uniroot.all(function(x) nbcurve(x, z = v_fits$z.list[[1]],w = v_fits$w.list[[1]],a = v_fits$a.list[[1]], b = v_fits$b.list[[1]]),c(10,150))
