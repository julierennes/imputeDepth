################################################################################
## File:             exp.impute.StudentT_MAR.R                                ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.07.2018                                               ##
##                                                                            ##
## Contains experiments on single depth-based imputation under MAR assumption ##
## (Figure 3).                                                                ##
##                                                                            ##
################################################################################

source("impute.functions.R")
library(imputeDepth)
library(geometry)
library(norm)
library(missForest)
library(missMDA)
library(rrcov)
library(VIM)
library(mvtnorm)
# Create structures
ems.depth.Tr2 <- NULL
ems.depth.zm <- NULL
ems.depth.M <- NULL
ems.depth.Mr <- NULL
ems.em <- NULL
ems.forest <- NULL
ems.knn <- NULL
ems.regPCA.1 <- NULL
ems.regPCA.2 <- NULL
ems.mean <- NULL
ems.oracle <- NULL
# Start study
d <- 3
mu <- c(1, 1, 1)
Sigma <- matrix(c(1, 1.75, 2, 1.75, 4, 4, 2, 4, 8), nrow = 3, byrow = TRUE)
n <- 100
cumPNA <- 0
k <- 1000
set.seed(1)
for (i in 1:k){
  cat("Iteration", i, "started.\n")
  # Generate data
  X <- t(mu + t(rmvt(n, Sigma, 0)))
  # Add NAs
  X.miss <- X
  missi.upper <- which(X.miss[,2] > mu[2])
  missi.lower <- which(X.miss[,2] < mu[2])
  if (length(missi.upper) > 0){
    X.miss[sample(missi.upper, length(missi.upper) * 0.48),3] <- NA
    X.miss[sample(missi.upper, length(missi.upper) * 0.08),1] <- NA
  }
  if (length(missi.lower) > 0){
    X.miss[sample(missi.lower, length(missi.lower) * 0.24),3] <- NA
    X.miss[sample(missi.lower, length(missi.lower) * 0.70),1] <- NA
  }
  cat(sum(is.na(X.miss)), " NAs produced\n")
  numNA <- sum(is.na(X.miss))
  pNA <- numNA / (n * d)
  cumPNA <- cumPNA + pNA
  # Impute
  X.imp.depth.Tr2 <- impute.depth(X.miss, depth = "halfspace",
                                  parMcd.outsiders = 0.5)
  X.imp.depth.zm <- impute.depth(X.miss, depth = "zonoid")
  X.imp.depth.M <- impute.depth(X.miss, depth = "Mahalanobis")
  X.imp.depth.Mr <- impute.depth(X.miss, depth = "Mahalanobis", 
                                 parMcd.impute = 0.85)
  X.imp.em <- imputeEm(as.matrix(X.miss))
  X.imp.forest <- missForest(X.miss)$ximp
  X.imp.knn <- imputeKnn(X.miss)
  X.imp.regPCA.1 <- imputePCA(X.miss, ncp = 1)$completeObs
  X.imp.regPCA.2 <- imputePCA(X.miss, ncp = 2)$completeObs
  X.imp.mean <- imputeMean(X.miss)
  X.imp.oracle <- imp.depth.Mahalanobis(X.miss, mu = mu, Sigma = Sigma)
  # Calculate raw statistics
  ems.depth.Tr2 <- c(ems.depth.Tr2, sqrt(sum((X.imp.depth.Tr2 - X)^2) /
                                           (n * d * pNA)))
  ems.depth.zm <- c(ems.depth.zm, sqrt(sum((X.imp.depth.zm - X)^2) /
                                         (n * d * pNA)))
  ems.depth.M <- c(ems.depth.M, sqrt(sum((X.imp.depth.M - X)^2) /
                                       (n * d * pNA)))
  ems.depth.Mr <- c(ems.depth.Mr, sqrt(sum((X.imp.depth.Mr - X)^2) /
                                         (n * d * pNA)))
  ems.em <- c(ems.em, sqrt(sum((X.imp.em - X)^2) / (n * d * pNA)))
  ems.forest <- c(ems.forest, sqrt(sum((X.imp.forest - X)^2) / (n * d * pNA)))
  ems.knn <- c(ems.knn, sqrt(sum((X.imp.knn - X)^2) / (n * d * pNA)))
  ems.regPCA.1 <- c(ems.regPCA.1, sqrt(sum((X.imp.regPCA.1 - X)^2) /
                                         (n * d * pNA)))
  ems.regPCA.2 <- c(ems.regPCA.2, sqrt(sum((X.imp.regPCA.2 - X)^2) /
                                         (n * d * pNA)))
  ems.mean <- c(ems.mean, sqrt(sum((X.imp.mean - X)^2) / (n * d * pNA)))
  ems.oracle <- c(ems.oracle, sqrt(sum((X.imp.oracle - X)^2) / (n * d * pNA)))
  # Save intermediate results
  if (i %% 10 < 1){
    save.image(paste("imp_t0-MAR-24_n100-d3-k", i, "_",
                     gsub(" ", "_", gsub(":", "_", date())), ".RData",
                     sep = ""))
  }
  # Calculate statistics
  cat("Iteration", i, "finished. Median RMSEs are:\n")
  cat("d.Tuk:", median(ems.depth.Tr2), ", c.zon:", median(ems.depth.zm),
      ", d.Mah:", median(ems.depth.M), ", d.MahR:", median(ems.depth.Mr),
      ", EM:", median(ems.em), ", rPCA1:", median(ems.regPCA.1),
      ", rPCA2:", median(ems.regPCA.2), ", kNN:", median(ems.knn), 
      ", RF:", median(ems.forest), ", mean:", median(ems.mean), 
      ", orcl", median(ems.oracle), ".\n")
}
errors <- list(TukeyR2 = ems.depth.Tr2, zonoidM = ems.depth.zm,
               Mahalanobis = ems.depth.M, MahalanobisR = ems.depth.Mr,
               em = ems.em,
               pca1 = ems.regPCA.1, pca2 = ems.regPCA.2,
               knn = ems.knn, forest = ems.forest,
               mean = ems.mean, oracle = ems.oracle)
boxplot(errors, main = paste("Large data, n = ", n, ", d = ", d, ", MCAR ", 
                             pNA, ", k = ", i, sep = ""),
        names = c("d.Tuk", "d.zon",
                  "d.Mah", "d.MahR",
                  "EM", "rPCA1", "rPCA2", "kNN", "RF",
                  "mean", "orcl"),
        ylab = "RMSE")
grid()
