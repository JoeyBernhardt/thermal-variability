

fits <- read_csv("Tetraselmis_experiment/data-processed/boot_fits_resample.csv")

f <- fits %>% 
  filter(curve.id.list == 677)

curvef <-function(x){
  res<-f$a.list[1]*exp(f$b.list[1]*x)*(1-((x-f$z.list[1])/(f$w.list[1]/2))^2)
  res
}

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + 
  # geom_line(aes(x = x, y = predictions, group = run), data = all_predictions) + 
  xlim(-20, 32) + ylim(-1, 2.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = -1.8) +
  stat_function(fun = curvef, color = "red")
