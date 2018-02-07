# HW2 - Optimization
# updated at Feb 5, 2018
# by Shixiang (Woody) Zhu

# Configurations
rootPath <- "/Users/woodie/Documents/Courses/ISyE 6416 Computational Statistics (Spring 2018)/HW/ISYE-6416/hw2"
source(paste(rootPath, "optimizer.R", sep="/"))
source(paste(rootPath, "plotter.R", sep="/"))

# Problem 1:
# ==============================================================================
# The following data are an i.i.d sample from a Cauchy(theta, 1) distribution: 
data <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 3.71, 
          -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)

#    Loglikelihood function of Cauchy distribution (theta, scale=1)
loglikCauchy <- function(location, scale=1., xs=data) {
  sum(dcauchy(xs, location=location, scale=scale, log=TRUE))
}
#    Derivative of loglikelihood function of Cauchy distribution (theta, 1)
derivLoglikCauchy <- function(location, xs=data) { 
  2 * sum((xs - location) / (1 + (xs - location)^2))
}

locations <- seq(-20, 60, 1/50) # Range of medians to plot
scales    <- seq(-5, 5, 1/50)   # Range of (log) scales to plot
u <- as.matrix(expand.grid(locations, exp(scales))) # Points to evaluate
y <- apply(as.matrix(locations), 1, function(v) loglikCauchy(v, 1., data))
z <- matrix(apply(u, 1, function(v) loglikCauchy(v[1], v[2], data)),
            nrow=length(locations))

# ------------------------------------------------------------------------------
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

#    Solve Maximum likelihood Estimation by Newton Raphson
newtonMLE <- newtonRaphson(f=derivLoglikCauchy, x0=mean(data))
# #    Plot solutions on their loglikelihood function and its derivative function
# plotSol(locations, derivLoglikCauchy, newtonMLE,
#         xlabel="theta", ylabel="derivative of log likelihood", 
#         subtitle="Derivative of Cauchy Log likelihood (scale=1.)",
#         title="Newton-Raphson")
# plotSol(locations, loglikCauchy, newtonMLE, 
#         xlabel="theta", ylabel="log likelihood", 
#         subtitle="Cauchy Log likelihood (scale=1.)",
#         title="Newton-Raphson",
#         zeroLine=FALSE)

# ------------------------------------------------------------------------------
# b. Apply the bisection method with starting points:
startPois <- c(-1., 1.)
#    Use additional runs to illustrate manners in which the bisection method may fail 
#    to find the global maximum. 

#    Solve Maximum likelihood Estimation by Bisection Method
bisecMLE <- bisection(f=derivLoglikCauchy, a=-1., b=1)

# ------------------------------------------------------------------------------
# c. Apply fixed-point iterations as in (2.29), starting from -1, with scaling choices
#    of alpha=1, 0.64, and 0.25. Investigate other choices of starting values and scaling 
#    factors.

#    Solve Maximum likelihood Estimation by Fixed Point Method
fixedPoiMLE <- fixedPoint(f=derivLoglikCauchy, x0=-1, alpha=0.64)

# ------------------------------------------------------------------------------
# d. From staring values of (theta^(0), theta^(1)) = (-2, -1), apply the secant method 
#    to estimate theta. What happens when (theta^(0), theta^(1)) = (-3, 3), and for 
#    other starting choices?

#    Solve Maximum likelihood Estimation by Fixed Point Method
secantMLE <- secant(f=derivLoglikCauchy, x0=-3, x1=3)

#    Plot solutions on their loglikelihood function and its derivative function
plotMLE(f=loglikCauchy, 
        newtonMLE, bisecMLE, fixedPoiMLE, secantMLE, 
        xlabel="theta", ylabel="log likelihood", 
        subtitle="Cauchy Log likelihood (scale=1.)",
        zeroLine=FALSE)
plotMLE(f=derivLoglikCauchy, 
        newtonMLE, bisecMLE, fixedPoiMLE, secantMLE, 
        xlabel="theta", ylabel="derivative of log likelihood", 
        subtitle="Derivative of Cauchy Log likelihood (scale=1.)",
        zeroLine=TRUE)

# ------------------------------------------------------------------------------
# e. Use this example to compare the speed and stability of the Newton-Raphson metho,
#    bisection, fixed-point iteration, and the secant method. Do your conclusions change
#    when you apply the methods to a random sample of size 20 from a N(theta, 1) 
#    distribution?

#    Sample data from a N(theta, 1) (theta = 0)
data <- rnorm(20, mean=0., sd=1.)
#    Loglikelihood function of Cauchy distribution (theta, scale=1)
loglikNorm <- function(theta, std=1., xs=data) {
  sum(dnorm(xs, mean=theta, sd=std, log=TRUE))
}
#    Derivative of loglikelihood function of Cauchy distribution (theta, 1)
derivLoglikNorm <- function(theta, xs=data) {
  require(numDeriv)
  genD(func=loglikNorm, x=theta)$D[1]
}

newtonMLE   <- newtonRaphson(f=derivLoglikNorm, x0=mean(data))
bisecMLE    <- bisection(f=derivLoglikNorm, a=-1., b=0)
fixedPoiMLE <- fixedPoint(f=derivLoglikNorm, x0=1, alpha=0.01)
secantMLE   <- secant(f=derivLoglikNorm, x0=-3, x1=3)

plotMLE(f=loglikNorm, 
        newtonMLE, bisecMLE, fixedPoiMLE, secantMLE, 
        zeroLine=FALSE)
plotMLE(f=derivLoglikNorm, 
        newtonMLE, bisecMLE, fixedPoiMLE, secantMLE,
        ylabel="derivative of log likelihood", 
        subtitle="Derivative of Normal Log likelihood (std=1.)",
        zeroLine=TRUE)

# Problem 2:
# ==============================================================================
# There were 46 crude oil spills of at least 1,000 barrels from tankers in U.S. waters 
# during 1974-1999. The website for this book contains the following data: the number
# of spills in the i-th year, N_i; the estimated amount of oil shipped through US waters
# as part of US import/export operations in the i-th year, adjusted for spillage in 
# international or foreign waters, b_i1; and the amount of oil shipped through U.S. 
# waters during domestic shipments in the i-th year, b_i2. The data are adapted from [11].
# Oil shipment amounts are measured in billions of barrels (Bbbl).
# 
# The volume of oil shipped is measure of exposure to spill risk. Suppose we use the 
# Possion process assumption given by N_i|b_i1, b_i2 ~ Possion(lambda_i) where 
# lambda_i = alpha_1 * b_i1 + alpha_2 * b_i2. The parameters of this model are alpha_1 and
# alpha_2, which represent the rate of spill occurence per Bbbl oil shipped during import/
# expert and domestic shipments, respectively. 

# a. 
#    



# References:
# - https://stats.stackexchange.com/questions/54885/how-do-i-plot-loglikelihood-functions-of-the-cauchy-distribution-in-r
# - https://rpubs.com/aaronsc32/newton-raphson-method
# - https://rpubs.com/aaronsc32/bisection-method-r
