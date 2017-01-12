################################################################################
## File:             exp.impute.StudentToutl_MCAR.R                           ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.12.2016                                               ##
##                                                                            ##
## Contains experiments on single depth-based imputation of outlier-          ##
## contaminated Student-t family under MCAR assumption (Table 2).             ##
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
Sigma <- matrix(c(1, 1, 1, 1, 4, 4, 1, 4, 8), nrow = 3, byrow = TRUE)
n <- 100
X.sd <- sqrt(matrix(rep(diag(Sigma), n), nrow = n, byrow = TRUE))
pNA <- 0.25
pOL <- 0.15
k <- 1000
# set.seed(1) # t1
# set.seed(2) # t2
# set.seed(3) # t3
# set.seed(5) # t5
# set.seed(10) # t10
set.seed(1000) # t0
for (i in 1:k){
  cat("Iteration", i, "started.\n")
  # Generate data with NAs
  Xn <- t(mu + t(rmvt(floor(n * (1 - pOL)), Sigma, 0)))
  Xo <- t(mu + t(rmvt(n - floor(n * (1 - pOL)), Sigma, 1)))
  X <- rbind(Xn, Xo)
  X.miss <- rbind(prodNA(Xn, pNA / (1 - pOL)), Xo)
  cat(sum(is.na(X.miss)), " NAs produced\n")
  # Impute
  X.imp.depth.Tr2 <- imp.depth.halfspace.ex.o(X.miss, 0.5)
  X.imp.depth.zm <- imp.depth.zonoid.o(X.miss, 1)
  X.imp.depth.M <- imp.depth.Mahalanobis(X.miss)
  X.imp.depth.Mr <- imp.depth.Mahalanobis(X.miss, alpha = 0.85)
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
  if (i %% 100 < 1){
    save.image(paste("imp_t0_outl-15_MCAR-25_n100-d3-k", i, "_",
                     gsub(" ", "_", gsub(":", "_", date())), ".RData",
                     sep = ""))
  }
  # Calculate statistics
  cat("Iteration", i, "finished. Median RMSEs are:\n")
  cat("Tur2:", median(ems.depth.Tr2), ", zom:", median(ems.depth.zm),
      ", Mah:", median(ems.depth.M), ", MahR:", median(ems.depth.Mr),
      ", em:", median(ems.em), ", for:", median(ems.forest),
      ", knn:", median(ems.knn), ", pc1:", median(ems.regPCA.1),
      ", pc2:", median(ems.regPCA.2), ", mean:", median(ems.mean),
      ", orcl:", median(ems.oracle), ".\n")
}
errors <- list(TukeyR2 = ems.depth.Tr2, zonoidM = ems.depth.zm,
               Mahalanobis = ems.depth.M, MahalanobisR = ems.depth.Mr,
               em = ems.em,
               forest = ems.forest, knn = ems.knn,
               pca1 = ems.regPCA.1, pca2 = ems.regPCA.2,
               mean = ems.mean, oracle = ems.oracle)
boxplot(errors, main = paste("t1 outl 100-3, MCAR ", pNA, ", k = ", i, sep = ""),
        names = c("Tur2", "zom",
                  "Mah", "MahR",
                  "em", "for", "knn", "pc1", "pc2",
                  "mean", "orcl"),
        ylab = "RMSE")
