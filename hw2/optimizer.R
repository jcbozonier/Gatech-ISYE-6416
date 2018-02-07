# Newton Raphson Method
newtonRaphson <- function(f, x0, tol=1e-3, n=1000) {
  require(numDeriv) # Package for computing f'(x)
  k <- n # Initialize for iteration results
  
  # Check the upper and lower bounds to see if approximations result in 0
  if (f(x0) == 0.) {
    res <- list("rootApproximation" = x0, "iterations" = c(x0))
  }
  
  for (i in 1:n) {
    dx   <- genD(func=f, x=x0)$D[1] # First-order derivative f'(x0)
    x1   <- x0 - (f(x0) / dx)       # Calculate next value x1
    k[i] <- x1 # Store x1
    # Once the difference between x0 and x1 becomes sufficiently small, output the results.
    if (abs(x1 - x0) < tol) {
      rootApprox <- tail(k, n=1)
      res <- list("rootApproximation" = rootApprox, 
                  "iterations" = k, 
                  "methodName" = "Newton Raphson")
      return(res)
    }
    # If Newton-Raphson has not yet reached convergence set x1 as x0 and continue
    x0 <- x1
  }
  print("Exceeded allowed number of iterations")
  rootApprox <- tail(k, n=1)
  res <- list("rootApproximation" = rootApprox, "iterations" = k)
  return(res)
}

# Bisection Method
bisection <- function(f, a, b, tol=1e-3, n=1000) {
  k <- n # Initialize for iteration results
  
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  if (!(f(a) < 0) && (f(b) > 0)) {
    stop('signs of f(a) and f(b) differ')
  } 
  else if (!(f(a) > 0) && (f(b) < 0)) {
    stop('signs of f(a) and f(b) differ')
  }
  
  for (i in 1:n) {
    c <- (a + b) / 2 # Calculate midpoint
    k[i] <- c
    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the 
    # function and return the root.
    if ((f(c) == 0) || ((b - a) / 2) < tol) {
      rootApprox <- tail(k, n=1)
      res <- list("rootApproximation" = rootApprox, 
                  "iterations" = k, 
                  "methodName" = "Bisection")
      return(res)
    }
    
    # If another iteration is required, 
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(f(c)) == sign(f(a)), 
           a <- c,
           b <- c)
  }
  # If the max number of iterations is reached and no root has been found, 
  # return message and end function.
  print("Exceeded allowed number of iterations")
}

# Fixed Point Method
fixedPoint <- function(f, x0, alpha=1, tol=1e-02, n=1000){
  # fixed-point algorithm to find x such that x + f(x) == x
  # assume that fun is a function of a single variable
  # x0 is the initial guess at the fixed point
  
  k <- n # Initialize for iteration results
  
  if (f(x0) == 0){
    res <- list("rootApproximation" = x0, "iterations" = c(x0))
  }
  for (i in 1:n) {
    x1   <- x0 + alpha * f(x0) 
    k[i] <- x1
    if ( abs((x1 - x0)) < tol ) {
      rootApprox <- tail(k, n=1)
      res <- list("rootApproximation" = rootApprox, 
                  "iterations" = k, 
                  "methodName" = "Fixed Point")
      return(res)
    }
    x0 <- x1
  } 
  stop("Exceeded allowed number of iterations")
}

secant <- function(f, x0, x1, tol=1e-03, n=1000){
  
  k <- n # Initialize for iteration results
  
  for ( i in 1:n ) {
    x2   <- x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
    k[i] <- x2
    if (abs(f(x2)) < tol) {
      rootApprox <- tail(k, n=1)
      res <- list("rootApproximation" = rootApprox, 
                  "iterations" = k,
                  "methodName" = "Secant")
      return(res)
    }
    x0 <- x1
    x1 <- x2
  }
  stop("Exceeded allowed number of iteractions")
}
