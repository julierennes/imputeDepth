################################################################################
## File:             exp.impute.lowrank.R                                     ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     22.07.2018                                               ##
##                                                                            ##
## Contains experiments on single depth-based imputation of a low-rank model  ##
## under MCAR assumption (Table 3).                                           ##
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
library(denoiseR)
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
# Start study
n <- 50
d <- 4
rnk <- 2
snr <- 2.5
pNA <- 0.2
k <- 1000
set.seed(123)
G <- LRsim(n, d, rnk, snr)
for (i in 1:k){
  cat("Iteration", i, "started.\n")
  # Generate data
  X <- G$mu + rmvt(n, sigma = diag(d) * 0.0005, df = 1)
  # Add NAs
  X.miss <- prodNA(X, pNA, entire.rows = FALSE)
  cat(sum(is.na(X.miss)), " NAs produced\n")
  # Impute
  X.imp.depth.Tr2 <- impute.depth(X.miss, depth = "Tukey",
                                  depth.outsiders = "spatial",
                                  parMcd.outsiders = 0.5)
  X.imp.depth.zm <- impute.depth(X.miss, depth = "zonoid",
                                 depth.outsiders = "spatial",
                                 parMcd.outsiders = 1)
  X.imp.depth.M <- impute.depth(X.miss, depth = "Mahalanobis",
                                parMcd.impute = 1)
  X.imp.depth.Mr <- impute.depth(X.miss, depth = "Mahalanobis",
                                 parMcd.impute = 0.85)
  X.imp.em <- imputeEm(as.matrix(X.miss))
  X.imp.forest <- missForest(X.miss)$ximp
  X.imp.knn <- imputeKnn(X.miss)
  X.imp.regPCA.1 <- imputePCA(X.miss, ncp = 1)$completeObs
  X.imp.regPCA.2 <- imputePCA(X.miss, ncp = 2)$completeObs
  X.imp.mean <- imputeMean(X.miss)
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
  # Save intermediate results
  if (i %% 100 < 1){
    save.image(paste("imp_lowrankt1_MCAR-20_n50-d4-k", i, "_",
                     gsub(" ", "_", gsub(":", "_", date())), ".RData",
                     sep = ""))
  }
  # Calculate statistics
  cat("Iteration", i, "finished. Median RMSEs are:\n")
  cat("d.Tuk:", median(ems.depth.Tr2), ", d.zon:", median(ems.depth.zm),
      ", d.Mah:", median(ems.depth.M), ", d.MahR:", median(ems.depth.Mr),
      ", EM:", median(ems.em), ", rPCA1:", median(ems.regPCA.1),
      ", rPCA2:", median(ems.regPCA.2), ", kNN:", median(ems.knn), 
      ", RF:", median(ems.forest), ", mean:", median(ems.mean), ".\n")
}
errors <- list(TukeyR2 = ems.depth.Tr2, zonoidM = ems.depth.zm,
               Mahalanobis = ems.depth.M, MahalanobisR = ems.depth.Mr,
               em = ems.em,
               pca1 = ems.regPCA.1, pca2 = ems.regPCA.2,
               knn = ems.knn, forest = ems.forest,
               mean = ems.mean)
boxplot(errors, main = paste("t1 lowrank 50-3, MCAR ", pNA, ", k = ", i, 
                             sep = ""),
        names = c("d.Tuk", "d.zon",
                  "d.Mah", "d.MahR",
                  "EM", "rPCA1", "rPCA2", "kNN", "RF", 
                  "mean"),
        ylab = "RMSE")
grid()
