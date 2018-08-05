################################################################################
## File:             plot.impute.moon.R                                       ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     22.07.2018                                               ##
##                                                                            ##
## Contains codes for plotting one imputation using both local and global     ##
## Tukey depth for the moon data set and the experiments on single            ##
## local-depth-based imputation of the moon data set (Figure 5).              ##
##                                                                            ##
################################################################################


library(imputeDepth)
source("impute.functions.R")
library(VIM)
ems.depth.loc1 <- NULL
set.seed(1)
# Generate the Moon
d <- 2
n <- 150
f1 <- function(x){1.5 * (1 - x^2)}
f2 <- function(x){2 * (1 - x^2)}
X <- cbind(runif(n * 1000, -1, 1), runif(n * 1000, 0, 2))
X <- X[X[,2] >= f1(X[,1]) & X[,2] <= f2(X[,1]),]
X <- X[1:n,]
# Add NAs
pNA <- 0.15
X.miss <- X
X.miss[sample.int(nrow(X), pNA * nrow(X) * ncol(X)),2] <- NA
numNA <- sum(is.na(X.miss))
cat((numNA), " NAs produced\n")
# Impute
X.imp.depth.loc1 <- impute.depth.local(X.miss, par.loc = 0.2, 
                                       depth = "Tukey",
                                       parMcd.impute = 1)
X.imp.depth1 <- impute.depth(X.miss, depth = "Tukey",
                             parMcd.impute = 1)
pdf("pic-moon-impute.pdf", width = 6, height = 6)
plot(X.miss, xlab = "", ylab = "")
points(X.imp.depth.loc1[is.na(X.miss)[,2],], pch = "+", col = "red")
points(X.imp.depth1[is.na(X.miss)[,2],], pch = 4, col = "blue")
dev.off()
load("imp_moon_MCAR-15_n150-d2-k100_Sun_Jun_17_04_53_52_2018.RData")
pdf("pic-moon-boxplots.pdf", width = 10, height = 7)
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
