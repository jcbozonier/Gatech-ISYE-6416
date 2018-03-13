# 		Homework 5: MC Integration

> Shixiang Zhu (GTID # 903280826)
>
> Email: [shixiang.zhu@gatech.edu](mailto:shixiang.zhu@gatech.edu)

### Problem 6.7

Consider pricing a European call option on an underlying stock withcurrent price $S^{(0)} = 50$, strike price $K=52$, and volatility $\sigma = 0.5$. Suppose that there are $N = 30$ days to maturity and that the risk-free rate of return is $r = 0.05$.

**a**. Confirm that the fair price for this option is 2.10 when the payoff is based on $S^{(30)}$ [i.e. a standard option with payoff as in (6.74)]
$$
\begin{equation}
\begin{aligned}
E[C]
&= exp(-\frac{rT}{365}) E[max(0, S^{(0)}-K)] \\
&= exp(-\frac{rT}{365})S^{(0)}exp((r-\frac{\sigma^2}{2})\frac{T}{365}) \\
&   ((1-\Phi(c-\sqrt{\frac{T}{365}}))exp(\frac{\sigma^2T}{365 \times 2})-(1-\Phi(c))\frac{K}{S^{(0)}} exp(-(\frac{r-\sigma^2}{2}\frac{T}{365})) \\
& \approx 2.101198
\end{aligned}
\end{equation}
$$

```R
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
  ST <- S0exp((r-(sigma^2)/2)T/365 + sigmarnorm(n)sqrt(T/365))
  C  <- NULL
  for(i in 1:n){
    C[i] <- exp(-rT/365)max(c(0,ST[i] - K))
  }
  estimatedMu[j] <- mean(C)
}
```

The fair price estimated by MC is `2.106494`, which is pretty close to `2.10`.



**b**. Consider the analogous Asian option (same $S^{(0)}$, $K$, $\sigma$, $N$, and $r$) with payoff based on the arithmetic mean stock price during the holding period, as in (6.77). Using simple Monte Carlo, estimate the fair price for this option.

```R
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
    A[i] <- exp(-rT/365)max(c(0, mean(ST) - K))
  }
  estimatedMu[j] <- mean(A)
}
```

The estimated fair price is `0.8562979` after average. And its variance is `0.003393549`



**c**. Improve upon the estimate in (b) using the control variate strategy described in Example 6.13

```R
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
```

The estimated fair price is `0.9440792` and its variance is `7.143709e-06`, which is way much better than the result in (b).

 

**d**. Try an antithetic approach to estimate the fair price for the option described in part (b).

```R
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
```

The estimated fair price is `0.8684389` and its variance is `0.003725535`, which is slightly better than the result in (b), but worse than (c).

**e**. Compare the sampling distributions of the estimators in (b), (c), (d).

Please see the details in (b), (c), (d).

### Problem 6.8

Consider the model given by $X ~ Lognormal(0, 1)$ and $log Y = 9 + log X + \epsilon$, where $\epsilon ~ N(0, 1)$. We wish to estimate $E\{Y/X\}$. Compare the performance of the standard Monte Carlo estimator and the Rao-Blackwellized estimator.

```R
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

# var(z2) < var(z1)
```

The estimated values are `98972.44` and `97806.38` respectively. However, the variance of Rao-Blackwellized estimator (`399722278386`) is smaller than the Monte Carlo one (`2.186901e+12`). 