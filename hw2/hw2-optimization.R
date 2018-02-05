# HW2 - Optimization

# Configurations
rootPath <- "/Users/woodie/Documents/Courses/ISyE 6416 Computational Statistics (Spring 2018)/HW/ISYE-6416/hw2"
source(paste(rootPath, "optimizer.R", sep="/"))

# Problem 1:
# The following data are an i.i.d sample from a Cauchy(theta, 1) distribution: 
data <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 3.71, 
          -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)

loglikCauchy <- function(location, scale, xs) {
  sum(dcauchy(xs, location=location, scale=scale, log=TRUE))
}

locations <- seq(-10, 10, 1/50) # Range of medians to plot
scales    <- seq(-5, 5, 1/50) # Range of (log) scales to plot
u <- as.matrix(expand.grid(locations, exp(scales))) # Points to evaluate
y <- apply(as.matrix(locations), 1, function(v) loglikCauchy(v, 1., data))
z <- matrix(apply(u, 1, function(v) loglikCauchy(v[1], v[2], data)),
            nrow=length(locations))

# a. Graph the log-likelihood function as location and scale vary
filled.contour(locations, scales, z, color.palette=heat.colors,
               xlab="Location", ylab="Log(scale)",
               main="Cauchy Log Likelihood")

#    Graph the log-likelihood function as only location vary when scale=1.
plot(locations, y, type="l", xlab="locations", ylab="log-likelihood", 
     main="Cauchy Log likelihood (scale=1.)")

#    Find the MLE for theta using the Newton-Raphson method. Try all of the 
#    following starting points: 
startPois <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
#    Discuss your results. Is the mean of the data is a good starting point?



# References:
# - https://stats.stackexchange.com/questions/54885/how-do-i-plot-loglikelihood-functions-of-the-cauchy-distribution-in-r