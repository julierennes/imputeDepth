################################################################################
## File:             exp.impute.StudentT_MCAR.R                               ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.07.2018                                               ##
##                                                                            ##
## Contains experiments on single depth-based imputation of Student-t family  ##
## contaminated with Cauchy distribution under MCAR assumption (Table 2).     ##
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
# ems.depth.TrE <- NULL
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
k <- 1000
# set.seed(1) # t1
# set.seed(2) # t2
# set.seed(3) # t3
# set.seed(5) # t5
# set.seed(10) # t10
set.seed(1000) # t0
cltype <- "SOCK"
nproc <- 32
if (nproc > 1){
  # Imputing routine to be calles on a node
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
    d <- 3
    mu <- c(1, 1, 1)
    Sigma <- matrix(c(1, 1, 1, 1, 4, 4, 1, 4, 8), nrow = 3, byrow = TRUE)
    n <- 100
    X.sd <- sqrt(matrix(rep(diag(Sigma), n), nrow = n, byrow = TRUE))
    pNA <- 0.25
    pOL <- 0.15
    # Generate data
    Xn <- t(mu + t(rmvt(floor(n * (1 - pOL)), Sigma, 0)))
    Xo <- t(mu + t(rmvt(n - floor(n * (1 - pOL)), Sigma, 1)))
    X <- rbind(Xn, Xo)
    # Add NAs
    X.miss <- rbind(prodNA(Xn, pNA / (1 - pOL)), Xo)
    cat(sum(is.na(X.miss)), " NAs produced\n")
    # Impute
    # X.imp.depth.TrE <- impute.depth(X.miss, depth = "halfspace", 
    #                                 depth.outsiders = "extremeval", 
    #                                 max.iter = 50, p.extreme = 0.1, 
    #                                 n.proj = 1000)
    X.imp.depth.Tr2 <- impute.depth(X.miss, depth = "halfspace",
                                    parMcd.outsiders = 0.5)
    X.imp.depth.zm <- impute.depth(X.miss, depth = "zonoid")
    X.imp.depth.M <- impute.depth(X.miss, depth = "Mahalanobis")
    X.imp.depth.Mr <- impute.depth(X.miss, depth = "Mahalanobis", 
                                   parMcd.impute = 0.75)
    X.imp.em <- imputeEm(as.matrix(X.miss))
    X.imp.forest <- missForest(X.miss)$ximp
    X.imp.knn <- imputeKnn(X.miss)
    X.imp.regPCA.1 <- imputePCA(X.miss, ncp = 1)$completeObs
    X.imp.regPCA.2 <- imputePCA(X.miss, ncp = 2)$completeObs
    X.imp.mean <- imputeMean(X.miss)
    X.imp.oracle <- imp.depth.Mahalanobis(X.miss, mu = mu, Sigma = Sigma)
    # Calculate raw statistics
    # ems.depth.TrE <- sqrt(sum((X.imp.depth.TrE - X)^2) / (n * d * pNA))
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
    # Return the errors
    return (c(-1, # ems.depth.TrE, 
              ems.depth.Tr2, ems.depth.zm, ems.depth.M, 
              ems.depth.Mr, ems.em, ems.forest, ems.knn, ems.regPCA.1,
              ems.regPCA.2, ems.mean, ems.oracle))
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
    # ems.depth.TrE <- c(ems.depth.TrE, res[[i]][1])
    ems.depth.Tr2 <- c(ems.depth.Tr2, res[[i]][2])
    ems.depth.zm <- c(ems.depth.zm, res[[i]][3])
    ems.depth.M <- c(ems.depth.M, res[[i]][4])
    ems.depth.Mr <- c(ems.depth.Mr, res[[i]][5])
    ems.em <- c(ems.em, res[[i]][6])
    ems.forest <- c(ems.forest, res[[i]][7])
    ems.knn <- c(ems.knn, res[[i]][8])
    ems.regPCA.1 <- c(ems.regPCA.1, res[[i]][9])
    ems.regPCA.2 <- c(ems.regPCA.2, res[[i]][10])
    ems.mean <- c(ems.mean, res[[i]][11])
    ems.oracle <- c(ems.oracle, res[[i]][12])
  }
  # Save results
  save.image(paste("imp_t0_outl-15_MCAR-25_n100-d3-k", i, "_",
                   gsub(" ", "_", gsub(":", "_", date())), ".RData",
                   sep = ""))
}else{
  # Start study
  d <- 3
  mu <- c(1, 1, 1)
  Sigma <- matrix(c(1, 1, 1, 1, 4, 4, 1, 4, 8), nrow = 3, byrow = TRUE)
  n <- 100
  pNA <- 0.25
  pOL <- 0.15
  for (i in 1:k){
    cat("Iteration", i, "started.\n")
    # Generate data
    Xn <- t(mu + t(rmvt(floor(n * (1 - pOL)), Sigma, 0)))
    Xo <- t(mu + t(rmvt(n - floor(n * (1 - pOL)), Sigma, 1)))
    X <- rbind(Xn, Xo)
    # Add NAs
    X.miss <- rbind(prodNA(Xn, pNA / (1 - pOL)), Xo)
    cat(sum(is.na(X.miss)), " NAs produced\n")
    # Impute
    # X.imp.depth.TrE <- impute.depth(X.miss, depth = "halfspace", 
    #                                 depth.outsiders = "extremeval", 
    #                                 max.iter = 50, p.extreme = 0.1, 
    #                                 n.proj = 1000)
    X.imp.depth.Tr2 <- impute.depth(X.miss, depth = "halfspace",
                                    parMcd.outsiders = 0.5)
    X.imp.depth.zm <- impute.depth(X.miss, depth = "zonoid")
    X.imp.depth.M <- impute.depth(X.miss, depth = "Mahalanobis")
    X.imp.depth.Mr <- impute.depth(X.miss, depth = "Mahalanobis", 
                                   parMcd.impute = 0.75)
    X.imp.em <- imputeEm(as.matrix(X.miss))
    X.imp.forest <- missForest(X.miss)$ximp
    X.imp.knn <- imputeKnn(X.miss)
    X.imp.regPCA.1 <- imputePCA(X.miss, ncp = 1)$completeObs
    X.imp.regPCA.2 <- imputePCA(X.miss, ncp = 2)$completeObs
    X.imp.mean <- imputeMean(X.miss)
    X.imp.oracle <- imp.depth.Mahalanobis(X.miss, mu = mu, Sigma = Sigma)
    # Calculate raw statistics
    # ems.depth.TrE <- c(ems.depth.TrE, sqrt(sum((X.imp.depth.TrE - X)^2) /
    #                                          (n * d * pNA)))
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
      save.image(paste("imp_t0_outl-15_MCAR-25_n100-d3-k", i, "_",
                       gsub(" ", "_", gsub(":", "_", date())), ".RData",
                       sep = ""))
    }
    # Calculate statistics
    cat("Iteration", i, "finished. Median RMSEs are:\n")
    cat(# "d.TukE:", median(ems.depth.TrE), 
      "d.Tuk:", median(ems.depth.Tr2), 
      ", d.zon:", median(ems.depth.zm), ", d.Mah:", median(ems.depth.M), 
      ", d.MahR:", median(ems.depth.Mr), ", EM:", median(ems.em), 
      ", rPCA1:", median(ems.regPCA.1), ", rPCA2:", median(ems.regPCA.2), 
      ", kNN:", median(ems.knn), ", RF:", median(ems.forest), 
      ", mean:", median(ems.mean), ", orcl:", median(ems.oracle), ".\n")
  }
}
errors <- list(# TukeyEE = ems.depth.TrE, 
               TukeyR2 = ems.depth.Tr2, 
               zonoidM = ems.depth.zm,
               Mahalanobis = ems.depth.M, MahalanobisR = ems.depth.Mr,
               em = ems.em,
               pca1 = ems.regPCA.1, pca2 = ems.regPCA.2,
               knn = ems.knn, forest = ems.forest,
               mean = ems.mean, oracle = ems.oracle)
boxplot(errors, main = paste("t0 outl 100-3, MCAR ", pNA, ", k = ", i, sep = ""),
        names = c(# "d.TukE", 
          "d.Tuk", "d.zon",
          "d.Mah", "d.MahR",
          "EM", "rPCA1", "rPCA2", "kNN", "RF",
          "mean", "orcl"),
        ylab = "RMSE")
grid()
