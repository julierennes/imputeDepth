################################################################################
## File:             exp.impute.pima.R                                        ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.12.2016                                               ##
##                                                                            ##
## Contains experiments on single depth-based imputation for the Pima data    ##
## set (Figure 6).                                                            ##
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
pima.all <- read.table(file = "data//pima.dat", header = FALSE, sep = " ")
pima.yes <- pima.all[pima.all[,8] == "Yes",]
X <- pima.yes[,2:5]
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
set.seed(1)
for (i in 1:k){
  cat("Iteration", i, "started.\n")
  # Add NAs
  X.miss <- prodNA(X, pNA, entire.rows = FALSE)
  numNA <- sum(is.na(X.miss))
  cat((numNA), " NAs produced\n")
  # Impute
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
  if (i %% 100 < 1){
    save.image(paste("pima_yes2-5_k", i, "_",
                     gsub(" ", "_", gsub(":", "_", date())), ".RData",
                     sep = ""))
  }
  # Calculate statistics
  cat("Iteration", i, "finished. Median RMSEs are:\n")
  cat("Tur2:", median(ems.depth.Tr2), ", zom:", median(ems.depth.zm),
      ", Mah:", median(ems.depth.M), ", MahR:", median(ems.depth.Mr),
      ", em:", median(ems.em), ", for:", median(ems.forest),
      ", knn:", median(ems.knn), ", pc1:", median(ems.regPCA.1),
      ", pc2:", median(ems.regPCA.2), ", mean:", median(ems.mean), ".\n")
}
errors <- list(TukeyR2 = ems.depth.Tr2, zonoidM = ems.depth.zm,
               Mahalanobis = ems.depth.M, MahalanobisR = ems.depth.Mr,
               em = ems.em,
               forest = ems.forest, knn = ems.knn,
               pca1 = ems.regPCA.1, pca2 = ems.regPCA.2,
               mean = ems.mean)
boxplot(errors, main = paste("Pima2-5, n = ", n, ", d = ", d, ", MCAR ",
                             pNA, ", k = ", i, sep = ""),
        names = c("Tur2", "zom",
                  "Mah", "MahR",
                  "em", "for", "knn", "pc1", "pc2",
                  "mean"),
        ylab = "RMSE")
