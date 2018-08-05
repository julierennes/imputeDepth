################################################################################
## File:             exp.impute.bloodtransfusion.R                            ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     22.07.2018                                               ##
##                                                                            ##
## Contains experiments on single depth-based imputation for the Blood        ##
## Transfusion data set (Figure 6).                                           ##
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
# Load data
bloodtransfusion.all <- read.table(file = "data//bloodtransfusion_gp.dat",
                                   header = FALSE, sep = " ")
bloodtransfusion <- bloodtransfusion.all[,1:3]
X <- bloodtransfusion
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
n <- nrow(X)
d <- ncol(X)
pNA <- 0.15
k <- 500
# Create structure for times
times.depths <- matrix(NA, nrow = k, ncol = 4)
cltype <- "SOCK"
nproc <- 32
set.seed(1)
if (nproc > 1){
  # Imputing routine to be called on a tread
  imp.worker <- function(i, ...){
    # Load packages
    source("impute.functions.R")
    library(imputeDepth)
    library(geometry)
    library(norm)
    library(missForest)
    library(missMDA)
    library(rrcov)
    library(VIM)
    # Start study
    t.depths <- rep(-1, 4)
    pNA <- 0.15
    X <- list(...)[[1]]$X
    n <- nrow(X)
    d <- ncol(X)
    # Add NAs
    X.miss <- prodNA(X, pNA, entire.rows = FALSE)
    cat(sum(is.na(X.miss)), " NAs produced\n")
    # Impute
    a <- system.time(
      X.imp.depth.Tr2 <- impute.depth(X.miss, depth = "Tukey",
                                      depth.outsiders = "spatial",
                                      parMcd.outsiders = 0.5)
    )
    t.depths[1] <- a[3]
    a <- system.time(
      X.imp.depth.zm <- impute.depth(X.miss, depth = "zonoid",
                                     depth.outsiders = "spatial",
                                     parMcd.outsiders = 1)
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
    # Return the errors and the times
    assembled <- list(errors = c(ems.depth.Tr2, ems.depth.zm, 
                                 ems.depth.M, ems.depth.Mr, ems.em, ems.forest, 
                                 ems.knn, ems.regPCA.1, ems.regPCA.2, ems.mean),
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
  parInput <- list(X = X)
  res <- performParallel(count = nproc, x = 1:k, imp.worker,
                         printfun = print.fun,
                         printargs = list(
                           mypb = txtProgressBar(min = 0,
                                                 max = k,
                                                 style = 3),
                           len.one.run = 1), printrepl = 1,
                         cltype = cltype,
                         ... = parInput)
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
    times.depths[i,] <- res[[i]]$times
  }
  # Save results
  save.image(paste("bloodtransfusion_k", i, "_",
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
      X.imp.depth.Tr2 <- impute.depth(X.miss, depth = "Tukey",
                                      depth.outsiders = "spatial",
                                      parMcd.outsiders = 0.5)
    )
    times.depths[i,1] <- a[3]
    a <- system.time(
      X.imp.depth.zm <- impute.depth(X.miss, depth = "zonoid",
                                     depth.outsiders = "spatial",
                                     parMcd.outsiders = 1)
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
    # Save intermediate results
    if (i %% 10 < 1){
      save.image(paste("bloodtransfusion_k", i, "_",
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
}
errors <- list(TukeyR2 = ems.depth.Tr2, zonoidM = ems.depth.zm,
               Mahalanobis = ems.depth.M, MahalanobisR = ems.depth.Mr,
               em = ems.em,
               pca1 = ems.regPCA.1, pca2 = ems.regPCA.2,
               knn = ems.knn, forest = ems.forest,
               mean = ems.mean)
boxplot(errors, main = paste("Blood Transfusion, n = ", n, ", d = ", d, 
                             ", MCAR ", pNA, ", k = ", i, sep = ""),
        names = c("d.Tuk", "d.zon",
                  "d.Mah", "d.MahR",
                  "EM", "rPCA1", "rPCA2", "kNN", "RF", 
                  "mean"),
        ylab = "RMSE")
grid()
