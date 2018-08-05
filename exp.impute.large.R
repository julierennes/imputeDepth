################################################################################
## File:             exp.impute.large.R                                       ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     22.07.2018                                               ##
##                                                                            ##
## Contains experiments on single depth-based imputation for the large data   ##
## set (Figure 3).                                                            ##
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
n <- 1000
d <- 6
pNA <- 0.15
k <- 500
# Create structure for times
times.depths <- matrix(NA, nrow = k, ncol = 4)
cltype <- "SOCK"
nproc <- 32
set.seed(1)
if (nproc > 1){
  # Imputing routine to be called on a tread
  imp.worker <- function(i){
    # Load packages
    source("impute.functions.R")
    library(imputeDepth)
    library(geometry)
    library(norm)
    library(missForest)
    library(missMDA)
    library(rrcov)
    library(VIM)
    library(mvtnorm)
    # Start study
    d <- 6
    n <- 1000
    mu <- rep(0, d)
    Sigma <- diag(d)
    for (i in 1:d){
      for (j in 1:d){
        Sigma[i,j] <- 2^(-abs(i - j))
      }
    }
    X.sd <- sqrt(matrix(rep(diag(Sigma), n), nrow = n, byrow = TRUE))
    pNA <- 0.15
    pOL <- 0.15
    # Generate data
    Xn <- t(mu + t(rmvt(floor(n * (1 - pOL)), Sigma, 0)))
    Xo <- t(mu + t(rmvt(n - floor(n * (1 - pOL)), Sigma, 1)))
    X <- rbind(Xn, Xo)
    # Add NAs
    X.miss <- rbind(prodNA(Xn, pNA / (1 - pOL)), Xo)
    cat(sum(is.na(X.miss)), " NAs produced\n")
    # Saving structure for times
    t.depths <- rep(-1, 4)
    # Impute
    a <- system.time(
      X.imp.depth.Tr2 <- impute.depth(X.miss, depth = "randomTukey",
                                      parMcd.outsiders = 0.85,
                                      n.proj = 1000)
    )
    t.depths[1] <- a[3]
    a <- system.time(
      X.imp.depth.zm <- impute.depth(X.miss, depth = "zonoid",
                                     parMcd.outsiders = 0.85)
    )
    t.depths[2] <- a[3]
    a <- system.time(
      X.imp.depth.M <- impute.depth(X.miss, depth = "Mahalanobis",
                                    parMcd.impute = 1)
    )
    t.depths[3] <- a[3]
    a <- system.time(
      X.imp.depth.Mr <- impute.depth(X.miss, depth = "Mahalanobis",
                                     parMcd.impute = 0.85)
    )
    t.depths[4] <- a[3]
    X.imp.em <- imputeEm(as.matrix(X.miss))
    X.imp.forest <- missForest(X.miss)$ximp
    X.imp.knn <- imputeKnn(X.miss)
    X.imp.regPCA.1 <- imputePCA(X.miss, ncp = 1)$completeObs
    X.imp.regPCA.2 <- imputePCA(X.miss, ncp = 2)$completeObs
    X.imp.mean <- imputeMean(X.miss)
    X.imp.oracle <- imp.depth.Mahalanobis(X.miss, mu = mu, Sigma = Sigma)
    # Calculate raw statistics
    ems.depth.Tr2 <- sqrt(sum((X.imp.depth.Tr2 - X)^2) / (n * d * pNA))
    ems.depth.zm <- sqrt(sum((X.imp.depth.zm - X)^2) / (n * d * pNA))
    ems.depth.M <- sqrt(sum((X.imp.depth.M - X)^2) / (n * d * pNA))
    ems.depth.Mr <- sqrt(sum((X.imp.depth.Mr - X)^2) / (n * d * pNA))
    ems.em <- sqrt(sum((X.imp.em - X)^2) / (n * d * pNA))
    ems.forest <- sqrt(sum((X.imp.forest - X)^2) / (n * d * pNA))
    ems.knn <- sqrt(sum((X.imp.knn - X)^2) / (n * d * pNA))
    ems.regPCA.1 <- sqrt(sum((X.imp.regPCA.1 - X)^2) / (n * d * pNA))
    ems.regPCA.2 <- sqrt(sum((X.imp.regPCA.2 - X)^2) / (n * d * pNA))
    ems.mean <- sqrt(sum((X.imp.mean - X)^2) / (n * d * pNA))
    ems.oracle <- sqrt(sum((X.imp.oracle - X)^2) / (n * d * pNA))
    # Return the errors and the times
    assembled <- list(errors = c(ems.depth.Tr2, ems.depth.zm, 
                                 ems.depth.M, ems.depth.Mr, ems.em, ems.forest, 
                                 ems.knn, ems.regPCA.1, ems.regPCA.2, ems.mean, 
                                 ems.oracle),
                      times = t.depths)
    return (assembled)
  }
  # Printing routine
  print.fun <-  function(outputs, B, args){
    pb <- args$mypb
    len.one.run <- args$len.one.run
    outputs <- unlist(outputs)
    outputs <- outputs[!is.null(outputs)]
    setTxtProgressBar(pb, length(outputs)/len.one.run)
  }
  # The parallel call
  library(snowFT)
  res <- performParallel(count = nproc, x = 1:k, imp.worker,
                         printfun = print.fun,
                         printargs = list(
                           mypb = txtProgressBar(min = 0,
                                                 max = k,
                                                 style = 3),
                           len.one.run = 1), printrepl = 1,
                         cltype = cltype)
  # Assemble results
  for (i in 1:k){
    ems.depth.Tr2 <- c(ems.depth.Tr2, res[[i]]$errors[1])
    ems.depth.zm <- c(ems.depth.zm, res[[i]]$errors[2])
    ems.depth.M <- c(ems.depth.M, res[[i]]$errors[3])
    ems.depth.Mr <- c(ems.depth.Mr, res[[i]]$errors[4])
    ems.em <- c(ems.em, res[[i]]$errors[5])
    ems.forest <- c(ems.forest, res[[i]]$errors[6])
    ems.knn <- c(ems.knn, res[[i]]$errors[7])
    ems.regPCA.1 <- c(ems.regPCA.1, res[[i]]$errors[8])
    ems.regPCA.2 <- c(ems.regPCA.2, res[[i]]$errors[9])
    ems.mean <- c(ems.mean, res[[i]]$errors[10])
    ems.oracle <- c(ems.oracle, res[[i]]$errors[11])
    times.depths[i,] <- res[[i]]$times
  }
  # Save results
  save.image(paste("large_k", i, "_",
                   gsub(" ", "_", gsub(":", "_", date())), ".RData",
                   sep = ""))
}else{
  for (i in 1:k){
    cat("Iteration", i, "started.\n")
    # Add NAs
    X.miss <- prodNA(X, pNA, entire.rows = FALSE)
    numNA <- sum(is.na(X.miss))
    cat((numNA), " NAs produced\n")
    # Impute
    a <- system.time(
      X.imp.depth.Tr2 <- impute.depth(X.miss, depth = "randomTukey",
                                      parMcd.outsiders = 0.85,
                                      n.proj = 1000)
    )
    times.depths[i,1] <- a[3]
    a <- system.time(
      X.imp.depth.zm <- impute.depth(X.miss, depth = "zonoid",
                                     parMcd.outsiders = 0.85)
    )
    times.depths[i,2] <- a[3]
    a <- system.time(
      X.imp.depth.M <- impute.depth(X.miss, depth = "Mahalanobis",
                                    parMcd.impute = 1)
    )
    times.depths[i,3] <- a[3]
    a <- system.time(
      X.imp.depth.Mr <- impute.depth(X.miss, depth = "Mahalanobis",
                                     parMcd.impute = 0.85)
    )
    times.depths[i,4] <- a[3]
    X.imp.em <- imputeEm(as.matrix(X.miss))
    X.imp.forest <- missForest(X.miss)$ximp
    X.imp.knn <- imputeKnn(X.miss)
    X.imp.regPCA.1 <- imputePCA(X.miss, ncp = 1)$completeObs
    X.imp.regPCA.2 <- imputePCA(X.miss, ncp = 2)$completeObs
    X.imp.mean <- imputeMean(X.miss)
    X.imp.oracle <- imp.depth.Mahalanobis(X.miss, mu = mu, Sigma = Sigma)
    # Calculate raw statistics
    ems.depth.Tr2 <- c(ems.depth.Tr2, sqrt(sum((X.imp.depth.Tr2 - X)^2,
                                               na.rm = TRUE) / numNA))
    ems.depth.zm <- c(ems.depth.zm, sqrt(sum((X.imp.depth.zm - X)^2,
                                             na.rm = TRUE) / numNA))
    ems.depth.M <- c(ems.depth.M, sqrt(sum((X.imp.depth.M - X)^2,
                                           na.rm = TRUE) / numNA))
    ems.depth.Mr <- c(ems.depth.Mr, sqrt(sum((X.imp.depth.Mr - X)^2,
                                             na.rm = TRUE) / numNA))
    ems.em <- c(ems.em, sqrt(sum((X.imp.em - X)^2, na.rm = TRUE) / numNA))
    ems.forest <- c(ems.forest, sqrt(sum((X.imp.forest - X)^2,
                                         na.rm = TRUE) / numNA))
    ems.knn <- c(ems.knn, sqrt(sum((X.imp.knn - X)^2, na.rm = TRUE) / numNA))
    ems.regPCA.1 <- c(ems.regPCA.1, sqrt(sum((X.imp.regPCA.1 - X)^2,
                                             na.rm = TRUE) / numNA))
    ems.regPCA.2 <- c(ems.regPCA.2, sqrt(sum((X.imp.regPCA.2 - X)^2,
                                             na.rm = TRUE) / numNA))
    ems.mean <- c(ems.mean, sqrt(sum((X.imp.mean - X)^2,
                                     na.rm = TRUE) / numNA))
    ems.oracle <- c(ems.oracle, sqrt(sum((X.imp.oracle - X)^2,
                                     na.rm = TRUE) / numNA))
    # Save intermediate results
    if (i %% 10 < 1){
      save.image(paste("large_k", i, "_",
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
}
ems.depth.Tr2 <- ems.depth.Tr2[ems.depth.Tr2 <= 2]
ems.depth.zm <- ems.depth.zm[ems.depth.zm <= 2]
ems.depth.M <- ems.depth.M[ems.depth.M <= 2]
ems.depth.Mr <- ems.depth.Mr[ems.depth.Mr <= 2]
ems.em <- ems.em[ems.em <= 2]
ems.forest <- ems.forest[ems.forest <= 2]
ems.knn <- ems.knn[ems.knn <= 2]
ems.regPCA.1 <- ems.regPCA.1[ems.regPCA.1 <= 2]
ems.regPCA.2 <- ems.regPCA.2[ems.regPCA.2 <= 2]
ems.mean <- ems.mean[ems.mean <= 2]
ems.oracle <- ems.oracle[ems.oracle <= 2]
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
