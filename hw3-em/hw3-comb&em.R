# HW3 - Combinatorial Optimization & EM
# updated at Feb 26, 2018
# by Shixiang (Woody) Zhu

# Question 3.1
# Implement a random starts local search algorithm for minimizing the AIC 
# for the baseball salary regression problem. Model your algorithm after 
# Example 3.3.

rootPath <- "/Users/woodie/Documents/Courses/ISyE 6416 Computational Statistics (Spring 2018)/HW/ISYE-6416"
baseballDataPath <- paste(rootPath, "datasets/baseball.dat", sep="/")
baseballData     <- read.table(baseballDataPath, header=TRUE)

# (a) Change the move strategy from steepest descent to immediate adoption 
#     of the first randomly selected downhill neighbor. 

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

localSearchForLR <- function(sol, x=baseballData[,2:28], y=baseballData[,1], k=1, n=1000){
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
      lr  <- lm(y ~ ., data=neighborX)    # build linear regression model on filtered data
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
    
    # # Rule 2: Immediate adoption of the first randomly selected downhill neighbor
    # # randomly pick neighborhood of current solution
    # counter <- 0
    # for (i in sample(1:nrow(neighbors))){
    #   neighbor <- neighbors[i,]
    #   neighbor <- neighbor[neighbor != 0] # k-neighborhood solution
    #   neighborX   <- x[,neighbor]
    #   neighborX$y <- y
    #   lr  <- lm(y ~ ., data=neighborX) # build linear regression model on filtered data
    #   res <- AIC(lr)                   # calculate AIC for lr model
    #   # pick the first downhill neiborhood solution
    #   if (res < lastRes) {
    #     minRes <- res
    #     sol    <- neighbor
    #     break
    #   }
    #   counter <- counter + 1
    # }
    # # stopping criterion
    # if (counter >= nrow(neighbors)){
    #   # Stop the clock
    #   dt <- proc.time() - ptm
    #   result <- list("solution"   = lastSol,
    #                  "time"       = dt,
    #                  "iterations" = iters)
    #   return(result)
    # }
    
    # Logging trace history of solutions
    print(sol)
    iters   <- c(iters, minRes)
    lastRes <- minRes
    lastSol <- sol
  }
  stop("Exceeded allowed number of iterations")
}

localSearchForLR(c(1,2,3,4), k=1)

# (b) Change the algorithm to employ 2-neiborhoods, and compare the 
#     results with those of previous runs.
localSearchForLR(c(1,2,3,4), k=2)

# Question 3.8
# Thirteen chemical measurements were carried out on each of 178 wines 
# from three regions of Italy. These data are available from the website 
# for this book. Using one or more heuristic search methods from this 
# chapter, partition the wines into three groups for which the total of 
# the within-group sum of squares is minimal. Comment on your work and 
# the results. This is a search problem of size 3^p where p=178. If you 
# have access to standard cluster analysis routines, check your results 
# using a standard method like that of Hartigan and Wong.

wineDataPath <- paste(rootPath, "datasets/wine.dat", sep="/")
wineData     <- read.table(wineDataPath, header=TRUE)

# Using K-Means (Hartigan and Wong) as a benchmark
set.seed(20)
wineCluster <- kmeans(wineData[, 2:14], 3, nstart = 20)
# Visualize the result via t-SNE
library(caret)  
library(Rtsne)
colorbar <- c("blue", "yellow", "red")
tsneModel <- Rtsne(
  as.matrix(wineData[, 2:14]), 
  check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
tsneEmbeddings         <- as.data.frame(tsneModel$Y)
tsneEmbeddings$cluster <- apply(as.matrix(wineCluster$cluster), 1, 
                                function(x) colorbar[x])
# Plot scatterring graph
ggplot(tsneEmbeddings, aes(x=V1, y=V2, color=cluster)) +  
  geom_point(size=2.) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  theme_light(base_size=10) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")

# Simulated Annealing
simulatedAnnealing <- function(sol, data=wineData[,2:13], n=10000, step=0.1) {
  
  # cost = group rss
  groupRSS <- function(data, cluster, K=3) {
    for (k in 1:K){
      # get indices of data entries which belongs to a same cluster
      indices <- which(cluster %in% k)
      # calculate mean value for each of attributes in a specific cluster
      mean    <- colMeans(data[indices,], dim=1)
      # get differences between original entries and mean value
      data[indices,] <- sweep(data[indices,], 2, mean)
    }
    return(sum(data^2))
  }
  
  # get neighborhood of current solution
  neighborhood <- function(sol){
    # init neighbors with padding zeros
    neighbors <- matrix(0, ncol=length(sol), nrow=length(sol)*2)
    # search all neighbors in distance 1.
    j <- 1
    for (i in 1:length(sol)){
      originVal <- sol[i]                       # original value of candidate position
      candVals  <- setdiff(c(1,2,3), originVal) # candidate values for this position
      for (val in candVals){
        neighbor      <- sol      # init as original solution
        neighbor[i]   <- val      # make change at specific position
        neighbors[j,] <- neighbor # add this neighbor to matrix
        j <- j + 1
      }
    }
    return(neighbors)
  }
  
  # Configuration
  alpha <- 0.01        # cooling rate
  beta  <- 2           # stage rate
  ptm   <- proc.time() # Start the clock!
  
  temp  <- 1    # temperature
  stage <- 1    # length of stage m   
  p    <- nrow(data) # number of variables
  cost <- groupRSS(data, sol, K=3)
  
  iters <- c() # init iteration results of cost
  for (j in 1:n){
    for (m in 1:stage){
      # get neighborhoods for the current solution
      neighbors       <- neighborhood(sol)
      res             <- c() # init results for cost
      neighborIndices <- c() # init candidates for neighbor
      for (i in 1:nrow(neighbors)){
        neighbor     <- neighbors[i,]
        neighborCost <- groupRSS(data, neighbor, K=3) # calculate cost for each neighbor
        if (neighborCost < cost){
          res             <- c(res, neighborCost)
          neighborIndices <- c(neighborIndices, i)
        }
      }
      # stop criterion
      if (length(res) <= 0){
        # Stop the clock
        dt <- proc.time() - ptm
        result <- list("solution"   = sol,
                       "time"       = dt,
                       "iterations" = iters)
        return(result)
      }
      # randomly pick a neighbor from candidates
      candInd      <- sample(1:length(res), 1)
      candNeighbor <- neighbors[neighborIndices[candInd],]
      candCost     <- res[candInd]
      # accept this solution by accept rate
      acceptRate   <- min(1, exp((cost - candCost)/temp))
      print(acceptRate)
      if (sample(c(TRUE,FALSE), size=1, replace=TRUE, 
                 prob=c(acceptRate,1-acceptRate))){
        sol   <- candNeighbor
        cost  <- candCost
        iters <- c(iters, cost)
      }
    }
    print(j)
    # update temperature and stage
    temp  <- temp/(1+alpha*temp)
    stage <- stage * beta
  }
}

initSol <- replicate(nrow(wineData), sample(c(1,2,3), 1, replace=TRUE))

# Question 4.1
# Recall the peppered moth analysis introduced in Example 4.2. In the 
# field, it is quite difficult to distinguish the insularia or typica 
# phenotypes due to variations in wing color and mottle. In addition to 
# the 662 moths mentioned in the example, suppose the sample collected 
# by the researchers actually included n_U=578 more moths that were known 
# to be insularia or typica but whose exact phenotypes could not be determined.

# (a) Derive the EM algorithm for maximum likelihood estimation of p_C, 
#     p_I, and p_I for this modified problem having observed data n_C, 
#     n_I, n_T, and n_U as given above.

# (b) Apply the algorithm to find the MLEs.

niters <- 1000
# raw data
nC <- 85
nI <- 196
nT <- 341
nU <- 578
n  <- nC + nI + nT + nU
# init value of p
pC <- 0.1
pI <- 0.1
pT <- 0.8

lastpC <- 0
lastpI <- 0
lastpT <- 1
for (i in 1:niters){
  # E (Estimation) step
  nCC <- (nC*pC^2) / (pC^2 + 2*pC*pI + 2*pC*pT)
  nCI <- (2*nC*pC*pI) / (pC^2 + 2*pC*pI + 2*pC*pT)
  nCT <- (2*nC*pC*pT) / (pC^2 + 2*pC*pI + 2*pC*pT)
  nII <- (nI*pI^2) / (pI^2 + 2*pI*pT) + (nU*pI^2) / (pI^2 + 2*pI*pT + pT^2)
  nIT <- (2*nI*pI*pT) / (pI^2 + 2*pI*pT) + (2*nU*pI*pT) / (pI^2 + 2*pI*pT + pT^2)
  nTT <- nT + (nU*pT^2) / (pI^2 + 2*pI*pT + pT^2)
  # M (Maximization) step
  pC <- (2*nCC + nCI + nCT) / (2*n)
  pI <- (2*nII + nIT + nCI) / (2*n)
  pT <- (2*nTT + nCT + nIT) / (2*n)
  # stop criterion
  if (((lastpC - pC)^2 + (lastpI - pI)^2 + (lastpT - pT)^2) < 10e-5){
    break
  } 
  lastpC <- pC
  lastpI <- pI
  lastpT <- pT
}

# (c) Estimate the standard errors and pairwise correlations for \hat{p_C}, \hat{p_I} and \hat{p_I} using the SEM algorithm.

