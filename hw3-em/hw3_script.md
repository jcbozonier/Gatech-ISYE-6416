# 	Homework 3: Combinatorial & EM

> Shixiang Zhu (GTID # 903280826) 
>
> Email: shixiang.zhu@gatech.edu



### Question 3.1

Implement a random starts local search algorithm for minimizing the AIC for the baseball salary regression problem. Model your algorithm after Example 3.3.

- Change the move strategy from steepest descent to immediate adoption of the first randomly selected downhill neighbor.

   

  **Code 1** *Below is a snippet of code for finding $k$-Neighborhoods of a solution.*

  ```R
  k_Neighborhood <- function(sol, p, k=1) {
    # index of positions in solution
    positions  <- 1:p
    # candidates for the positions of changes
    candidates <- combn(positions, k)
    # init neighbors with padding zeros
    neighbors  <- matrix(0, ncol=p, nrow=ncol(candidates))
    for (i in 1:ncol(candidates)){
      # init neighbor with the starting point
      neighbor <- c(sol)
      for (ind in candidates[,i]){
        # turn it to 1 if the change position is not in solution
        if (ind %in% sol){
          neighbor <- setdiff(neighbor, ind)
        }
        # turn it to 0 if the change position is in solution
        else{
          neighbor <- c(neighbor, ind)
        }
      }
      # do zero paddings for neighbor
      neighbors[i,] <- c(neighbor, rep(0, p - length(neighbor)))
    }
    return(neighbors)
  }
  ```

  **Code 2** *Below is a snippet of code for local searching via `steepest descent`*

  ```R
  localSearchForLR <- function(
      sol, x=baseballData[,2:28], y=baseballData[,1], k=1, n=1000){
    # Start the clock!
    ptm <- proc.time()
    
    p   <- ncol(x)           # number of variables
    x$y <- y                 # add y to data
    X       <- x[,sol]           # init x data according to init solution
    lr      <- lm(y ~ ., data=X) # build linear regression model on filtered data
    lastRes <- AIC(lr)           # calculate AIC for lr model
    
    iters <- c() # init iteration results of all AIC history
    for (j in 1:n){
      # get neighborhoods for the current solution
      neighbors <- k_Neighborhood(sol, p, k=k)
      
      # Rule 1: Search all neighborhood
      # calculate lr for every neighborhood of current solution
      res <- c() # init results of AIC
      for (i in 1:nrow(neighbors)){
        neighbor <- neighbors[i,]
        neighbor <- neighbor[neighbor != 0] # k-neighborhood solution
        neighborX   <- x[,neighbor]
        neighborX$y <- y
        lr  <- lm(y ~ ., data=neighborX)    # build linear regression model
        res <- c(res, AIC(lr))              # calculate AIC for lr model
      }
      # find better neighborhood in accordance with its AIC result
      minInd <- which.min(res)
      minRes <- res[minInd]
      sol    <- neighbors[minInd,]
      sol    <- sol[sol != 0]
      # Stopping criterion
      if (minRes >= lastRes){
        # Stop the clock
        dt <- proc.time() - ptm
        result <- list("solution"   = lastSol,
                       "time"       = dt,
                       "iterations" = iters)
        return(result)
      }
      
      # Logging trace history of solutions
      iters   <- c(iters, minRes)
      lastRes <- minRes
      lastSol <- sol
    }
    stop("Exceeded allowed number of iterations")
  }
  ```

  **Code 3** *Below is a snippet of code for local searching via `immediate adoption of the first randomly selected downhill neighbor`.*

  ```R
  localSearchForLR <- function(
      sol, x=baseballData[,2:28], y=baseballData[,1], k=1, n=1000){
    # Start the clock!
    ptm <- proc.time()
    
    p   <- ncol(x)           # number of variables
    x$y <- y                 # add y to data
    X       <- x[,sol]           # init x data according to init solution
    lr      <- lm(y ~ ., data=X) # build linear regression model on filtered data
    lastRes <- AIC(lr)           # calculate AIC for lr model
    
    iters <- c() # init iteration results of all AIC history
    for (j in 1:n){
      # get neighborhoods for the current solution
      neighbors <- k_Neighborhood(sol, p, k=k)
      
      # Rule 2: Immediate adoption of the first randomly selected downhill neighbor
      # randomly pick neighborhood of current solution
      counter <- 0
      for (i in sample(1:nrow(neighbors))){
        neighbor <- neighbors[i,]
        neighbor <- neighbor[neighbor != 0] # k-neighborhood solution
        neighborX   <- x[,neighbor]
        neighborX$y <- y
        lr  <- lm(y ~ ., data=neighborX) # build linear regression model on filtered data
        res <- AIC(lr)                   # calculate AIC for lr model
        # pick the first downhill neiborhood solution
        if (res < lastRes) {
          minRes <- res
          sol    <- neighbor
          break
        }
        counter <- counter + 1
      }
      # stopping criterion
      if (counter >= nrow(neighbors)){
        # Stop the clock
        dt <- proc.time() - ptm
        result <- list("solution"   = lastSol,
                       "time"       = dt,
                       "iterations" = iters)
        return(result)
      }
      
      # Logging trace history of solutions
      iters   <- c(iters, minRes)
      lastRes <- minRes
      lastSol <- sol
    }
    stop("Exceeded allowed number of iterations")
  }
  ```

  ​

- Change the algorithm to employ $2$-neighborhoods, and compare the results with those of previous runs.



### Question 3.8

Thirteen chemical measurements were carried out on each of 178 wines from three regions of Italy. These data are available from the website for this book. Using one or more heuristic search methods from this chapter, partition the wines into three groups for which the total of the within-group sum of squares is minimal. Comment on your work and the results. This is a search problem of size $3^p$ where $p=178$. If you have access to standard cluster analysis routines, check your results using a standard method like that of Hartigan and Wong.







