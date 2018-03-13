# HW 5 MC Integration
# ---------------------------------------------------------------------
# Problem 6.7
# Consider pricing a European call option on an underlying stock with 
# current price S^(0) = 50, strike price K=52, and volatility \sigma = 
# 0.5. Suppose that there are N = 30 days to maturity and that the risk
# -free rate of return is r = 0.05.

# a. Confirm that the fair price for this option is 2.10 when the payoff
#    is based on S^(30) [i.e. a standard option with payoff as in (6.74)]
S0    <- 50
K     <- 52
sigma <- 0.5
N     <- 30
r     <- 0.05
T     <- 30
n <- 1000
m <- 100

estimatedMu <- NULL
for(j in 1:m){
  ST <- S0*exp((r-(sigma^2)/2)*T/365 + sigma*rnorm(n)*sqrt(T/365))
  C  <- NULL
  for(i in 1:n){
    C[i] <- exp(-r*T/365)*max(c(0,ST[i] - K))
  }
  estimatedMu[j] <- mean(C)
}
# mu <- mean(estimatedMu)

# b. Consider the analogous Asian option (same S^(0), K, \sigma, N, and 
#    r) with payoff based on the arithmetic mean stock price during the 
#    holding period, as in (6.77). Using simple Monte Carlo, estimate the
#    fair price for this option.
estimatedMu <- NULL
for(j in 1:m){
  #calculate MC estimate of A and theta
  A <- NULL
  for(i in 1:n){
    ST    <- NULL
    ST[1] <- S0
    for(k in 2:T){
      ST[k] <- ST[k-1]*exp(((r-(sigma^2)/2)/365) +
                             sigma*rnorm(1)/sqrt(365))
    }
    A[i] <- exp(-r*T/365)*max(c(0, mean(ST) - K))
  }
  estimatedMu[j] <- mean(A)
}
# mu <- mean(estimatedMu)

# c. Improve upon the estimate in (b) using the control variate strategy
#    described in Example 6.13
estimatedMuMC  <- NULL
estimatedTheta <- NULL
for(j in 1:m){
  # calculate MC estimate of A and theta
  A     <- NULL
  theta <- NULL
  for(i in 1:n){
    ST    <- NULL
    ST[1] <- S0
    for(k in 2:T){
      ST[k] <- ST[k-1]*exp(((r-(sigma^2)/2)/365) +
                             sigma*rnorm(1)/sqrt(365))
    }
    A[i]     <- exp(-r*T/365)*max(c(0, mean(ST) - K))
    theta[i] <- exp(-r*T/365)*max(c(0,exp(mean(log(ST))) - K))
  }
  estimatedMuMC[j]  <- mean(A)
  estimatedTheta[j] <- mean(theta)
}
# mu <- mean(estimatedMu)

# Analytic solution for theta
N  <- T
c3 <- 1 + 1/N
c2 <- sigma*((c3*T/1095)*(1 + 1/(2*N)))^.5
c1 <- (1/c2)*(log(S0/K) + (c3*T/730)*(r - (sigma^2)/2) +
             (c3*(sigma^2)*T/1095)*(1 + 1/(2*N)))
theta <- S0*pnorm(c1)*exp(-T*(r + c3*(sigma^2)/6)*(1 - 1/N)/730) -
         K*pnorm(c1-c2)*exp(-r*T/365)

# Control variate
estimatedMuCV <- estimatedMuMC - 1 * (estimatedTheta-theta)

# d. Try an antithetic approach to estimate the fair price for the option
#    described in part (b)

estimatedMu1 <- NULL
for(j in 1:m){
  #calculate MC estimate of A and theta
  A <- NULL
  for(i in 1:n/2){
    ST    <- NULL
    ST[1] <- S0
    for(k in 2:T){
      ST[k] <- ST[k-1]*exp(((r-(sigma^2)/2)/365) +
                             sigma*rnorm(1)/sqrt(365))
    }
    A[i] <- exp(-r*T/365)*max(c(0, mean(ST) - K))
  }
  estimatedMu1[j] <- mean(A)
}

estimatedMu2 <- NULL
for(j in 1:m){
  #calculate MC estimate of A and theta
  A <- NULL
  for(i in 1:n/2){
    ST    <- NULL
    ST[1] <- S0
    for(k in 2:T){
      ST[k] <- ST[k-1]*exp(((r-(sigma^2)/2)/365) +
                             sigma*rnorm(1)/sqrt(365))
    }
    A[i] <- exp(-r*T/365)*max(c(0, mean(ST) - K))
  }
  estimatedMu2[j] <- mean(A)
}

estimatedMuA <- (estimatedMu1 + estimatedMu2)/2

# Problem 6.8 
# Consider the model given by X ~ Lognormal(0, 1) and log Y = 9 + log X
# + \epsilon, where \epsilon ~ N(0, 1). We wish to estimate E{Y/X}. 
# Compare the performance of the standard Monte Carlo estimator and the 
# Rao-Blackwellized estimator.
n  <- 1000000
x  <- exp(rnorm(n, 0, 1))
y1 <- exp(9+3*log(x)+rnorm(n, 0, 1))
z1 <- y1/x
mean(z1)
var(z1)

# since E[Y/X | X] = E[exp(9+log(X))*exp(Normal(0, 1)/ X] = exp(9+log(X))*exp(1/2)/X  

x  <- exp(rnorm(n, 0, 1))
y2 <- exp(9+3*log(x))*exp(1/2)/x  
y2 <- exp(9)*x^2*exp(1/2)
z2 <- y2
mean(z2)
var(z2)
# var(z2) < var(z1)