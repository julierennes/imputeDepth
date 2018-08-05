################################################################################
## File:             plot.one.density.R                                       ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.07.2018                                               ##
##                                                                            ##
## Contains the script for plotting imputation density (Figure 2, supplement).##
##                                                                            ##
################################################################################

library(missForest)
library(imputeDepth)
library(geometry)
library(MASS)
library(mvtnorm)
library(depth)
k <- 10000
res.M.3 <- list("")
res.z.3 <- list("")
res.h.3 <- list("")
ns <- c(50, 100, 200, 500, 1000)
d <- 2
mu <- c(1, 1)
Sigma <- matrix(c(1, 1, 1, 4), nrow = 2)
for (i in 1:length(ns)){
  res.M.3[[i]] <- rep(0, k)
  res.z.3[[i]] <- rep(0, k)
  res.h.3[[i]] <- rep(0, k)
  for (j in 1:k){
    X <- t(mu + t(rmvt(ns[i] - 1, Sigma, 1)))
    X.miss <- rbind(X, c(3, NA))
    X.imp.M <- imp.depth.Mahalanobis(X.miss)
    X.imp.z <- imp.depth.zonoid.o(X.miss)
    X.imp.h <- imp.depth.halfspace.ex.o(X.miss)
    res.M.3[[i]][j] <- X.imp.M[ns[i],2]
    res.z.3[[i]][j] <- X.imp.z[ns[i],2]
    res.h.3[[i]][j] <- X.imp.h[ns[i],2]
    cat("n =", ns[i], ": Iteration", j, "done.\n")
  }
}

cs <- c(2.5, 2.75, 2.875, 2.925, 2.955)
hs <- c(1.575, 1.25, 0.915, 0.685, 0.54)
for (i in 1:5){
  plot(density(res.h.3[[i]]), xlim = c(0, 6), ylim = c(0, 0.8), type = "l", 
       lwd = 2, main = paste("Student t1, n = ", ns[i], ", d = ", d, ", k = ", 
                             k, sep = ""), 
       xlab = "Maximum Tukey depth imputation of X[i,2] for X[i,1] = 3", 
       ylab = "Kernel density estimate of imputed values")
  grid()
  x <- seq(-2, 8, 0.01)
  hx <- dnorm(x, cs[i], hs[i])
  lines(x, hx, col = "red", lwd = 2, lty = 2)
  abline(v = cs[i], col = "red", lwd = 2, lty = 2)
  abline(v = 3, col = "green", lwd = 2, lty = 2)
}
