################################################################################
## File:             plot.impute.skewed.R                                     ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     22.07.2018                                               ##
##                                                                            ##
## Contains codes for plotting one imputation using the Tukey depth for the   ##
## skewed data set and the experiments on single local-depth-based imputation ##
## of the skewed data set (Figure 4).                                         ##
##                                                                            ##
################################################################################

library(imputeDepth)
source("impute.functions.R")
library(VIM)
library(sn)
ems.depth.loc1 <- NULL
set.seed(1)
# Generate the skewed sample
d <- 2
n <- 150
op <- list(xi=c(1,1), Psi=matrix(c(2,2,2,3), 2, 2), lambda=c(4, -2))
X <- rmsn(n, dp=op2dp(op,"SN"))
plot(X, asp = 1)
# Add NAs
pNA <- 0.15
X.miss <- prodNA(X, pNA)
numNA <- sum(is.na(X.miss))
cat((numNA), " NAs produced\n")
# Impute
X.imp.depth1 <- impute.depth(X.miss, depth = "Tukey")
X.imp.depth.loc1 <- impute.depth.local(X.miss, par.loc = 0.8, 
                                       depth = "Tukey",
                                       parMcd.impute = 1)
X.imp.depth.loc2 <- imputeKnn(X.miss)
pdf("pic-skewed-impute.pdf", width = 6, height = 6)
plot(X.miss, xlab = "", ylab = "")
points(X.imp.depth1[rowSums(is.na(X.miss)) > 0.5,], pch = "+", col = "red")
dev.off()
plot(X.miss, xlab = "", ylab = "")
points(X.imp.depth.loc1[rowSums(is.na(X.miss)) > 0.5,], pch = "+", col = "red")
load("imp_skewed_MCAR-15_n150-d2-k100_Mon_Jun_18_02_10_46_2018.RData")
pdf("pic-skewed-boxplots.pdf", width = 10, height = 7)
errors <- list(TukeyR2 = ems.depth.Tr2, 
               zonoidM = ems.depth.zm,
               Mahalanobis = ems.depth.M, MahalanobisR = ems.depth.Mr,
               em = ems.em,
               pca1 = ems.regPCA.1,
               knn = ems.knn,
               forest = ems.forest, 
               mean = ems.mean,
               loc1 = ems.depth.loc1, loc2 = ems.depth.loc2, 
               loc3 = ems.depth.loc3)
boxplot(errors, main = paste("moon, MCAR ", pNA, ", k = ", k, sep = ""),
        names = c("d.Tuk", "d.zon",
                  "d.Mah", "d.MahR",
                  "EM", 
                  "rPCA1", "kNN", "RF", "mean",
                  "ld.Tuk", "ld.zon", "ld.Mah"),
        ylab = "RMSE", cex.axis = 0.95)
dev.off()
