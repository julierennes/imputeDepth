################################################################################
## File:             plot.intro.R                                             ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.12.2016                                               ##
##                                                                            ##
## Contains the script for plotting demonstrative pictures from introduction  ##
## (Figures 1 and 2).                                                         ##
##                                                                            ##
################################################################################

# Functions ####################################################################
imputeEm <- function(X){
  s <- prelim.norm(X)
  thetahat <- em.norm(s, criterion = sqrt(.Machine$double.eps))
  params <- getparam.norm(s, thetahat)
  X.prep <- t(t(X) - params$mu)
  Inv.Sigma.tmp <- solve(params$sigma)
  miss.rowi = which(rowSums(is.na(X.prep)) > 0.5)
  X.new <- X.prep
  X.new <- X.new[miss.rowi,,drop=FALSE]
  X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllP, Inv.Sigma.tmp))
  return (t(t(X.prep) + params$mu))
}
# Functions (end) ##############################################################

library(missForest)
library(imputeDepth)
library(geometry)
library(MASS)
library(mvtnorm)
library(norm)
# Set parameters
d <- 2
mu <- c(1, 1)
Sigma <- matrix(c(1, 1, 1, 4), nrow = 2)
# Calculate regressions' coefficients
n <- 1000000
X0 <- t(mu + t(rmvt(n, Sigma, 0)))
reg.yx <- lm(X0[,2] ~ X0[,1])$coefficients
reg.xy <- lm(X0[,1] ~ X0[,2])$coefficients

# Picture 1 - Prediction
n <- 150
set.seed(1)
X1 <- t(mu + t(rmvt(n, Sigma, 0)))
X.miss1 <- X1
X.miss1[X1[,2] > 3.5,1] <- NA
X.imp.z1 <- imp.depth.zonoid.o(X.miss1)
X.imp.h1 <- imp.depth.halfspace.ex.o(X.miss1)
X.miss1.ext <- cbind(X.miss1, 1)
X.imp.f1 <- missForest(X.miss1.ext)$ximp
plot(X1[rowSums(is.na(X.miss1))<0.5,],
     main = "MAR assumption",
     ylim = c(min(X1[,2]), max(X1[,2])))
abline(h = X1[rowSums(is.na(X.miss1))>0.5,2], lty = 3)
lines(rbind(c(min(X1[,1]) * 2, min(X1[,1]) * 2 * reg.yx[2] + reg.yx[1]), 
            c(max(X1[,1]) * 2, max(X1[,1]) * 2 * reg.yx[2] + reg.yx[1])))
lines(rbind(c(min(X1[,2]) * 2 * reg.xy[2] + reg.xy[1], min(X1[,2]) * 2), 
            c(max(X1[,2]) * 2 * reg.xy[2] + reg.xy[1], max(X1[,2]) * 2)))
points(X.imp.f1[rowSums(is.na(X.miss1))>0.5,], col = "green", pch = 17)
points(X.imp.z1[rowSums(is.na(X.miss1))>0.5,], col = "red", pch = 19)

# Picture 2 - Close to data
X2 <- X1
n <- nrow(X1)
X2.miss <- X2
set.seed(1)
n.miss <- round(n * d * 0.3)
X2.miss[sample(1:(n * d), n.miss)] <- NA
X.imp.z2 <- imp.depth.zonoid.o(X2.miss)
X.imp.h2 <- imp.depth.halfspace.ex.o(X2.miss)
X.imp.M2 <- imp.depth.Mahalanobis(X2.miss)
X.imp.e2 <- imputeEm(X2.miss)
plot(X2[rowSums(is.na(X2.miss))<0.5,], main = "MCAR assumption",
     xlim = c(min(X2.miss[,1], na.rm = T), max(X2.miss[,1], na.rm = T)),
     ylim = c(min(X2.miss[,2], na.rm = T), max(X2.miss[,2], na.rm = T)))
abline(h = X2[is.na(X2.miss)[,1] & !is.na(X2.miss)[,2],2], lty = 3)
abline(v = X2[is.na(X2.miss)[,2] & !is.na(X2.miss)[,1],1], lty = 3)
lines(rbind(c(min(X2[,1]) * 2, min(X2[,1]) * 2 * reg.yx[2] + reg.yx[1]), 
            c(max(X2[,1]) * 2, max(X2[,1]) * 2 * reg.yx[2] + reg.yx[1])))
lines(rbind(c(min(X2[,2]) * 2 * reg.xy[2] + reg.xy[1], min(X2[,2]) * 2), 
            c(max(X2[,2]) * 2 * reg.xy[2] + reg.xy[1], max(X2[,2]) * 2)))
points(X.imp.e2[rowSums(is.na(X2.miss))>0.5,], col = "blue", pch = 18)
points(X.imp.z2[rowSums(is.na(X2.miss))>0.5,], col = "red", pch = 19)

# Picture 3 - Outliers
set.seed(1)
pNA <- 0.15
pOL <- 0.15
n <- 500
d <- 2
mu <- c(1, 1)
Sigma <- matrix(c(1, 1, 1, 4), nrow = 2)
X <- t(mu + t(rmvt(n * (1 - pOL), Sigma, 0)))
Xo <- t(mu + t(rmvt(n * pOL, Sigma, 1)))
X3 <- rbind(X, Xo)
n.miss <- round(n * d * pNA)
X3.miss <- t(X3)
X3.miss[sample(1:(n * (1 - pOL) * d), n.miss)] <- NA
X3.miss <- t(X3.miss)
X.imp.M3 <- imp.depth.Mahalanobis(X3.miss)
X.imp.e3 <- imputeEm(X3.miss)
X.imp.h3 <- imp.depth.halfspace.ex.o(X3.miss)
plot(X3[rowSums(is.na(X3.miss))<0.5,], main = "MCAR assumption, outliers",
     xlim = c(min(X3[1:(n * (1 - pOL)),1], na.rm = T), 
              max(X3[1:(n * (1 - pOL)),1], na.rm = T)),
     ylim = c(min(X3[1:(n * (1 - pOL)),2], na.rm = T), 
              max(X3[1:(n * (1 - pOL)),2], na.rm = T)))
abline(h = X3[is.na(X3.miss)[,1] & !is.na(X3.miss)[,2],2], lty = 3)
abline(v = X3[is.na(X3.miss)[,2] & !is.na(X3.miss)[,1],1], lty = 3)
lines(rbind(c(min(X3[,1]) * 2, min(X3[,1]) * 2 * reg.yx[2] + reg.yx[1]), 
            c(max(X3[,1]) * 2, max(X3[,1]) * 2 * reg.yx[2] + reg.yx[1])))
lines(rbind(c(min(X3[,2]) * 2 * reg.xy[2] + reg.xy[1], min(X3[,2]) * 2), 
            c(max(X3[,2]) * 2 * reg.xy[2] + reg.xy[1], max(X3[,2]) * 2)))
points(X.imp.e3[rowSums(is.na(X3.miss))>0.5,], col = "blue", pch = 18)
points(X.imp.h3[rowSums(is.na(X3.miss))>0.5,], col = "red", pch = 19)

# Picture 4 - Heavy tails
pNA <- 0.15
n <- 1000
set.seed(1000)
mu <- c(1, 1)
Sigma <- matrix(c(1, 1, 1, 4), nrow = 2)
X4 <- t(mu + t(rmvt(n, Sigma, 1)))
X4.miss <- X4
n.miss <- round(n * d * pNA)
X4.miss[sample(1:(n * d), n.miss)] <- NA
X.imp.M4 <- imp.depth.Mahalanobis(X4.miss)
X.imp.e4 <- imputeEm(X4.miss)
X.imp.h4 <- imp.depth.halfspace.ex.o(X4.miss)
plot(X4[rowSums(is.na(X4.miss))<0.5,], main = "MCAR assumption, heavy tails",
     xlim = c(-25, 25),
     ylim = c(-50, 50))
abline(h = X4[is.na(X4.miss)[,1] & !is.na(X4.miss)[,2],2], lty = 3)
abline(v = X4[is.na(X4.miss)[,2] & !is.na(X4.miss)[,1],1], lty = 3)
lines(rbind(c(min(X4[,1]) * 2, min(X4[,1]) * 2 * reg.yx[2] + reg.yx[1]), 
            c(max(X4[,1]) * 2, max(X4[,1]) * 2 * reg.yx[2] + reg.yx[1])))
lines(rbind(c(min(X4[,2]) * 2 * reg.xy[2] + reg.xy[1], min(X4[,2]) * 2), 
            c(max(X4[,2]) * 2 * reg.xy[2] + reg.xy[1], max(X4[,2]) * 2)))
points(X.imp.e4[rowSums(is.na(X4.miss))>0.5,], col = "blue", pch = 18)
points(X.imp.h4[rowSums(is.na(X4.miss))>0.5,], col = "red", pch = 19)
