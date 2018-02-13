library(MASS)

newtonRaphsonMLE <- function(loglikFunc, theta0, tol=1e-9, n=1000) {
  # Start the clock!
  ptm <- proc.time()
  
  require(numDeriv) # Package for computing f'(x) and f''(x)
  k <- matrix(0, ncol=length(theta0)) # Initialize for iteration results
  
  for (i in 1:n) {
    deriv  <- genD(func=loglikFunc, x=theta0)$D[1:length(theta0)] 
    hess   <- hessian(func=loglikFunc, x=theta0)                  
    theta1 <- theta0 - deriv %*% solve(hess)  # Calculate next value theta1 
    k      <- rbind(k, theta1)                # Store theta1

    # Once the difference between theta0 and theta1 becomes sufficiently small, 
    # output the results.
    if (sum(abs(theta1 - theta0)) < tol) {
      # Stop the clock
      dt <- proc.time() - ptm
      rootApprox <- tail(k, n=1)
      res <- list("rootApproximation" = rootApprox, 
                  "iterations" = tail(k, n=i), 
                  "time" = dt,
                  "methodName" = "Newton Raphson")
      return(res)
    }
    # If Newton-Raphson has not yet reached convergence set theta1 as theta0 
    # and continue
    theta0 <- theta1
  }
  stop("Exceeded allowed number of iterations")
}

fisherScoringMLE <- function(loglikFunc, data, theta0, 
                             tol=1e-3, n=1000) {
  # Start the clock!
  ptm <- proc.time()
  
  require(numDeriv) # Package for computing f'(x) and f''(x)
  k <- matrix(0, ncol=length(theta0)) # Initialize for iteration results
  
  for (i in 1:n) {
    # Score of loglikelihood 
    score      <- genD(func=loglikFunc, x=theta0)$D[1:length(theta0)] 
    # Fisher Information of loglikelihood
    fisherInfo <- matrix(0, ncol=length(theta0), nrow=length(theta0))
    for (j in 1:nrow(data)) {
      estFunc <- function(theta, d=data[j,]) { loglikFunc(theta, d) }
      estVec  <- genD(func=estFunc, x=theta0)$D[1:length(theta0)]
      fisherInfo <- fisherInfo + estVec%*%t(estVec)
    }
    # Update theta by fisher scoring
    theta1 <- theta0 + MASS::ginv(fisherInfo) %*% score
    k      <- rbind(k, t(theta1)) # Store theta1
    if (sum(abs(theta1 - theta0)) < tol) {
      # Stop the clock
      dt <- proc.time() - ptm
      rootApprox <- tail(k, n=1)
      res <- list("rootApproximation" = rootApprox, 
                  "iterations" = tail(k, n=i), 
                  "time" = dt,
                  "methodName" = "Fisher Scoring")
      return(res)
    }
    theta0 <- theta1
  }
  stop("Exceeded allowed number of iterations")
}

steepestAscentMLE <- function(loglikFunc, theta0, alpha0=1., 
                              tol=1e-3, n=1000, backtracking=TRUE) {
  # Start the clock!
  ptm <- proc.time()
  
  require(numDeriv) # Package for computing f'(x) and f''(x)
  k <- matrix(0, ncol=length(theta0)) # Initialize for iteration results
  
  alpha <- alpha0
  for (i in 1:n) {
    deriv  <- genD(func=loglikFunc, x=theta0)$D[1:length(theta0)] 
    # Steepest ascent
    h      <- -1 * alpha * solve(-1 * diag(length(theta0))) %*% deriv
    # Update theta
    theta1 <- theta0 + alpha * deriv
    k      <- rbind(k, t(theta1)) # Store theta1
    # Backtracking by updating alpha
    if (backtracking && (loglikFunc(theta0) - loglikFunc(theta1) > 0.)) {
      alpha <- alpha * 0.5
    }
    if (sum(abs(theta1 - theta0)) < tol) {
      # Stop the clock
      dt <- proc.time() - ptm
      rootApprox <- tail(k, n=1)
      res <- list("rootApproximation" = rootApprox, 
                  "iterations" = tail(k, n=i), 
                  "time" = dt,
                  "methodName" = "Fisher Scoring")
      return(res)
    }
    theta0 <- theta1
  }
  stop("Exceeded allowed number of iterations")
}

quasiNewtonMLE <- function(loglikFunc, theta0, alpha0=1, 
                           tol=1e-3, n=1000, backtracking=TRUE) {
  # Start the clock!
  ptm <- proc.time()
  
  k <- matrix(0, ncol=length(theta0)) # Initialize for iteration results
  # An initial matrix M0 is chosen (usually M0 = I)
  M <- diag(length(theta0))
  # An initial alpha alpha0 is chosen
  alpha <- alpha0
  for (i in 1:n) {
    # Update theta
    deriv0 <- genD(func=loglikFunc, x=theta0)$D[1:length(theta0)]
    theta1 <- theta0 - alpha * solve(M) %*% deriv0
    deriv1 <- genD(func=loglikFunc, x=theta1)$D[1:length(theta1)]
    # Update Matrix M
    z <- as.vector(theta1 - theta0)
    y <- deriv1 - deriv0
    v <- y - M %*% z
    c <- solve(t(v) %*% z)
    M <- M + c[1] * (v %*% t(v))
    # Backtracking by updating alpha
    if (backtracking && (loglikFunc(theta0) - loglikFunc(theta1) > 0.)) {
      alpha <- alpha * 0.5
    }
    # Update theta
    k <- rbind(k, t(theta1)) # Store theta1
    if (sum(abs(theta1 - theta0)) < tol) {
      # Stop the clock
      dt <- proc.time() - ptm
      rootApprox <- tail(k, n=1)
      res <- list("rootApproximation" = rootApprox, 
                  "iterations" = tail(k, n=i),
                  "time" = dt,
                  "methodName" = "Fisher Scoring")
      return(res)
    }
    theta0 <- theta1
  }
  stop("Exceeded allowed number of iterations")
}