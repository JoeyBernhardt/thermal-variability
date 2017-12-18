
### Goal: find inflection point of the constant temperature TPC

derivative <- function(f, x, ..., order = 1, delta = 0.1, sig = 6) {
	# Numerically computes the specified order derivative of f at x
	vals <- matrix(NA, nrow = order + 1, ncol = order + 1)
	grid <- seq(x - delta/2, x + delta/2, length.out = order + 1)
	vals[1, ] <- sapply(grid, f, ...) - f(x, ...)
	for (i in 2:(order + 1)) {
		for (j in 1:(order - i + 2)) {
			stepsize <- grid[i + j - 1] - grid[i + j - 2]
			vals[i, j] <- (vals[i - 1, j + 1] - vals[i - 1, j])/stepsize
		}
	}
	return(signif(vals[order + 1, 1], sig))
}


curve_constant_resamp<-function(x){
	res<-fits_c$a.list[1]*exp(fits_c$b.list[1]*x)*(1-((x-fits_c$z.list[1])/(fits_c$w.list[1]/2))^2)
	res
}


par(mfrow = c(1, 2))
x <- seq(14.716, 14.717, by = 0.001)
for (i in 2:2) {
	plot(x, sapply(x, derivative, f = curve_constant_resamp(), order = i), type = "l", xlab = "x", 
			 ylab = "y", main = paste("Order", i, "Derivative"))
	abline(h = 0, col = "gray60")
}





derivative <- function(f, x, ..., order = i, delta = 0.1, sig = 6) {
	# Numerically computes the specified order derivative of f at x
	vals <- matrix(NA, nrow = order + 1, ncol = order + 1)
	grid <- seq(x - delta/2, x + delta/2, length.out = order + 1)
	vals[1, ] <- sapply(grid, f, ...) - f(x, ...)
	for (i in 2:(order + 1)) {
		for (j in 1:(order - i + 2)) {
			stepsize <- grid[i + j - 1] - grid[i + j - 2]
			vals[i, j] <- (vals[i - 1, j + 1] - vals[i - 1, j])/stepsize
		}
	}
	return(signif(vals[order + 1, 1], sig))
}
###

x <- seq(0, 33, by = 0.001)


deriv_2 <- function(x) {
	y <- derivative(f = curve_constant_resamp, x = x, order = 2)
}

deriv_value <- sapply(x, deriv_2)
variable_upper2 <- data.frame(x, deriv_value) %>% 
	rename(temperature = x)


variable_upper2 %>% 
	filter(deriv_value > 16.8) %>% View
	ggplot(aes(x = temperature, y = deriv_value)) + geom_line() +
	geom_hline(yintercept = 0) + xlim(16.5, 17) + ylim(-0.01, 0.01)
