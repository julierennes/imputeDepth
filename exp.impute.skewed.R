################################################################################
## File:             exp.impute.skewed.R                                      ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     22.07.2018                                               ##
##                                                                            ##
## Contains experiments on single local-depth-based imputation of the skewed  ##
## data set (Figure 4).                                                       ##
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
library(sn)
# Create structures
ems.depth.loc1 <- NULL
ems.depth.loc2 <- NULL
ems.depth.loc3 <- NULL
ems.depth.Tr2 <- NULL
ems.depth.zm <- NULL
ems.depth.M <- NULL
ems.depth.Mr <- NULL
ems.em <- NULL
ems.forest <- NULL
ems.knn <- NULL
ems.regPCA.1 <- NULL
ems.mean <- NULL
k <- 100
pNA <- 0.15
set.seed(1)
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
    library(sn)
    # Generate the skewed sample
    d <- 2
    n <- 150
    op <- list(xi=c(1,1), Psi=matrix(c(2,2,2,3), 2, 2), lambda=c(4, -2))
    X <- rmsn(n, dp=op2dp(op,"SN"))
    # Add NAs
    pNA <- 0.15
    X.miss <- prodNA(X, pNA)
    numNA <- sum(is.na(X.miss))
    cat((numNA), " NAs produced\n")
    # Impute
    X.imp.depth.loc1 <- impute.depth.local(X.miss, par.loc = 0.8, 
                                           depth = "halfspace",
                                           parMcd.impute = 1)
    X.imp.depth.loc2 <- impute.depth.local(X.miss, par.loc = 0.8, 
                                           depth = "zonoid",
                                           parMcd.impute = 1)
    X.imp.depth.loc3 <- impute.depth.local(X.miss, par.loc = 0.8, 
                                           depth = "Mahalanobis",
                                           parMcd.impute = 1)
    X.imp.depth.Tr2 <- impute.depth(X.miss, depth = "halfspace", 
                                    depth.outsiders = "spatial", 
                                    parMcd.outsiders = 0.5)
    X.imp.depth.zm <- impute.depth(X.miss, depth = "zonoid")
    X.imp.depth.M <- impute.depth(X.miss, depth = "Mahalanobis")
    X.imp.depth.Mr <- impute.depth(X.miss, depth = "Mahalanobis", 
                                   parMcd.impute = 0.75)
    X.imp.em <- imputeEm(as.matrix(X.miss))
    X.imp.forest <- missForest(cbind(X.miss, 1))$ximp[,1:2]
    X.imp.knn <- imputeKnn(X.miss)
    X.imp.regPCA.1 <- imputePCA(X.miss, ncp = 1)$completeObs
    X.imp.mean <- imputeMean(X.miss)
    # Calculate raw statistics
    ems.depth.loc1 <- sqrt(sum((X.imp.depth.loc1 - X)^2) / (n * d * pNA))
    ems.depth.loc2 <- sqrt(sum((X.imp.depth.loc2 - X)^2) / (n * d * pNA))
    ems.depth.loc3 <- sqrt(sum((X.imp.depth.loc3 - X)^2) / (n * d * pNA))
    ems.depth.Tr2 <- sqrt(sum((X.imp.depth.Tr2 - X)^2) / (n * d * pNA))
    ems.depth.zm <- sqrt(sum((X.imp.depth.zm - X)^2) / (n * d * pNA))
    ems.depth.M <- sqrt(sum((X.imp.depth.M - X)^2) / (n * d * pNA))
    ems.depth.Mr <- sqrt(sum((X.imp.depth.Mr - X)^2) / (n * d * pNA))
    ems.em <- sqrt(sum((X.imp.em - X)^2) / (n * d * pNA))
    ems.forest <- sqrt(sum((X.imp.forest - X)^2) / (n * d * pNA))
    ems.knn <- sqrt(sum((X.imp.knn - X)^2) / (n * d * pNA))
    ems.regPCA.1 <- sqrt(sum((X.imp.regPCA.1 - X)^2) / (n * d * pNA))
    ems.mean <- sqrt(sum((X.imp.mean - X)^2) / (n * d * pNA))
    # Return the errors
    return (c(ems.depth.Tr2, ems.depth.zm, ems.depth.M, 
              ems.depth.Mr, ems.em, ems.forest, ems.knn,
              ems.regPCA.1, ems.mean,
              ems.depth.loc1, ems.depth.loc2, ems.depth.loc3))
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
    ems.depth.Tr2 <- c(ems.depth.Tr2, res[[i]][1])
    ems.depth.zm <- c(ems.depth.zm, res[[i]][2])
    ems.depth.M <- c(ems.depth.M, res[[i]][3])
    ems.depth.Mr <- c(ems.depth.Mr, res[[i]][4])
    ems.em <- c(ems.em, res[[i]][5])
    ems.forest <- c(ems.forest, res[[i]][6])
    ems.knn <- c(ems.knn, res[[i]][7])
    ems.regPCA.1 <- c(ems.regPCA.1, res[[i]][8])
    ems.mean <- c(ems.mean, res[[i]][9])
    ems.depth.loc1 <- c(ems.depth.loc1, res[[i]][10])
    ems.depth.loc2 <- c(ems.depth.loc2, res[[i]][11])
    ems.depth.loc3 <- c(ems.depth.loc3, res[[i]][12])
  }
  # Save results
  save.image(paste("imp_skewed_MCAR-15_n150-d2-k", i, "_",
                   gsub(" ", "_", gsub(":", "_", date())), ".RData",
                   sep = ""))
}else{
  # Start study
  d <- 2
  n <- 150
  pNA <- 0.15
  for (i in 1:k){
    cat("Iteration", i, "started.\n")
    # Generate the skewed sample
    op <- list(xi=c(1,1), Psi=matrix(c(2,2,2,3), 2, 2), lambda=c(4, -2))
    X <- rmsn(n, dp=op2dp(op,"SN"))
    # Add NAs
    X.miss <- prodNA(X, pNA)
    numNA <- sum(is.na(X.miss))
    cat((numNA), " NAs produced\n")
    # Impute
    X.imp.depth.loc1 <- impute.depth.local(X.miss, par.loc = 0.8, 
                                           depth = "halfspace",
                                           parMcd.impute = 1)
    X.imp.depth.loc2 <- impute.depth.local(X.miss, par.loc = 0.8, 
                                           depth = "zonoid",
                                           parMcd.impute = 1)
    X.imp.depth.loc3 <- impute.depth.local(X.miss, par.loc = 0.8, 
                                           depth = "Mahalanobis",
                                           parMcd.impute = 1)
    X.imp.depth.Tr2 <- impute.depth(X.miss, depth = "halfspace", 
                                    depth.outsiders = "spatial", 
                                    parMcd.outsiders = 0.5)
    X.imp.depth.zm <- impute.depth(X.miss, depth = "zonoid")
    X.imp.depth.M <- impute.depth(X.miss, depth = "Mahalanobis")
    X.imp.depth.Mr <- impute.depth(X.miss, depth = "Mahalanobis", 
                                   parMcd.impute = 0.75)
    X.imp.em <- imputeEm(as.matrix(X.miss))
    X.imp.forest <- missForest(cbind(X.miss, 1))$ximp[,1:2]
    X.imp.knn <- imputeKnn(X.miss)
    X.imp.regPCA.1 <- imputePCA(X.miss, ncp = 1)$completeObs
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
    ems.mean <- c(ems.mean, sqrt(sum((X.imp.mean - X)^2) / (n * d * pNA)))
    ems.depth.loc1 <- c(ems.depth.loc1, sqrt(sum((X.imp.depth.loc1 - X)^2) / 
                                               (n * d * pNA)))
    ems.depth.loc2 <- c(ems.depth.loc2, sqrt(sum((X.imp.depth.loc2 - X)^2) / 
                                               (n * d * pNA)))
    ems.depth.loc3 <- c(ems.depth.loc3, sqrt(sum((X.imp.depth.loc3 - X)^2) / 
                                               (n * d * pNA)))
    # Save intermediate results
    if (i %% 10 < 1){
      save.image(paste("imp_skewed_MCAR-15_n150-d2-k", i, "_",
                       gsub(" ", "_", gsub(":", "_", date())), ".RData",
                       sep = ""))
    }
    # Calculate statistics
    cat("Iteration", i, "finished. Median RMSEs are:\n")
    cat("d.Tuk:", median(ems.depth.Tr2), ", d.zon:", median(ems.depth.zm),
        ", d.Mah:", median(ems.depth.M), ", d.MahR:", median(ems.depth.Mr),
        ", EM:", median(ems.em), ", rPCA1:", median(ems.regPCA.1),
        ", kNN:", median(ems.knn), ", RF:", median(ems.forest), 
        ", mean:", median(ems.mean), 
        ", ld.Tuk:", median(ems.depth.loc1), 
        ", ld.zon:", median(ems.depth.loc2),
        ", ld.Mah:", median(ems.depth.loc3), ".\n")
  }
}
errors <- list(TukeyR2 = ems.depth.Tr2, zonoidM = ems.depth.zm,
               Mahalanobis = ems.depth.M, MahalanobisR = ems.depth.Mr,
               em = ems.em,
               pca1 = ems.regPCA.1,
               knn = ems.knn, forest = ems.forest,
               mean = ems.mean,
               loc1 = ems.depth.loc1, loc2 = ems.depth.loc2, 
               loc3 = ems.depth.loc3)
boxplot(errors, main = paste("Skewed, n = ", n, ", d = ", d, 
                             ", MCAR ", pNA, ", k = ", i, sep = ""),
        names = c("d.Tuk", "d.zon",
                  "d.Mah", "d.MahR",
                  "EM", "rPCA1", "kNN", "RF", 
                  "mean", "ld.Tuk", "ld.zon", "ld.Mah"),
        ylab = "RMSE")
grid()
