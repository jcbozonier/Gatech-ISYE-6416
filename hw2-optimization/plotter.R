plotSol <- function (x, yLoglikFunc, resMLE, 
                     xlabel, ylabel, subtitle, title, 
                     zeroLine=TRUE) {
  # Variables Configurations
  y           <- apply(as.matrix(x), 1, function(v) yLoglikFunc(v))
  xIters      <- resMLE$iterations 
  yIters      <- apply(as.matrix(xIters), 1, function(v) yLoglikFunc(v))
  labelIters  <- cbind(1:length(xIters), xIters)
  namebank    <- apply(as.matrix(labelIters), 1, 
                       function(v) sprintf("iter %d: %f", v[1], v[2]))
  # Plot derivative of loglikelihood function
  plot(x, y, type="l", xlab=xlabel, ylab=ylabel, sub=subtitle, main=title)
  if (zeroLine) {
    abline(h = 0, lty = 2)
  }
  points(xIters, yIters, pch=16, col="blue")
  points(resMLE$rootApproximation, # Root x
         yLoglikFunc(resMLE$rootApproximation), 
         pch=20, col="red")
  text(resMLE$rootApproximation,   # Root x
       yLoglikFunc(resMLE$rootApproximation),
       labels=sprintf("iter %d: %f\ntime: %f", length(xIters), resMLE$rootApproximation, resMLE$time[1]), 
       cex= 0.7, pos=1)
}

plotMLE <- function (f, newtonMLE, bisecMLE, fixedPoiMLE, secantMLE, 
                     step=1/50, xrange=c(-30, 30),
                     xlabel="theta", ylabel="log likelihood", 
                     subtitle="Normal Log likelihood (std=1.)",
                     zeroLine=TRUE) {
  par(mfrow=c(2,2))
  xs <- seq(xrange[1], xrange[2], step) # Range of means to plot
  plotSol(xs, f, newtonMLE, 
          xlabel=xlabel, ylabel=ylabel, subtitle=subtitle, title=newtonMLE$methodName,
          zeroLine=zeroLine)
  plotSol(xs, f, bisecMLE, 
          xlabel=xlabel, ylabel=ylabel, subtitle=subtitle, title=bisecMLE$methodName,
          zeroLine=zeroLine)
  plotSol(xs, f, fixedPoiMLE, 
          xlabel=xlabel, ylabel=ylabel, subtitle=subtitle, title=fixedPoiMLE$methodName,
          zeroLine=zeroLine)
  plotSol(xs, f, secantMLE, 
          xlabel=xlabel, ylabel=ylabel, subtitle=subtitle, title=secantMLE$methodName,
          zeroLine=zeroLine)
}

plot2DMLE <- function(f, resMLE, subtitle, xlim=c(-5, 10), ylim=c(-5, 10), step=1/50) {
  x <- seq(xlim[1], xlim[2], step)
  y <- seq(ylim[1], ylim[2], step)
  u <- as.matrix(expand.grid(x, y)) # Points to evaluate
  z <- matrix(apply(u, 1, function(v) f(v)), nrow=length(x))
  contour(x, y, z, nlevels=30, sub=subtitle)
  # points(resMLE$iterations[,1], resMLE$iterations[,2], pch=16, col="blue")
  lines(resMLE$iterations[,1], resMLE$iterations[,2], type="o", col="blue") 
  points(resMLE$rootApproximation[1], resMLE$rootApproximation[2],
         pch=20, col="red")
  text(resMLE$rootApproximation[1], resMLE$rootApproximation[2],
       labels=sprintf("iter %d: %f\ntime: %f", length(resMLE$iterations[,1]), resMLE$rootApproximation, resMLE$time[1]), 
       cex= 0.7, pos=1)
}