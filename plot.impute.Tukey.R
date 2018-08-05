################################################################################
## File:             plot.impute.Tukey.R                                      ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.07.2018                                               ##
##                                                                            ##
## Plots the imputation of a point using averaging for the Tukey depth        ##
## imputation (Figure 3, supplement).                                         ##
##                                                                            ##
################################################################################

library(depth)
library(MASS)
set.seed(1)
n <- 27
d <- 2
mu <- c(1, 1)
sigma <- matrix(c(1, 1, 1, 4), nrow = 2)
X <- mvrnorm(n, mu, sigma)
print(colMeans(X))
pdf("pic_Tukey_imp.pdf", width = 6, height = 6)
# 1) Draw pure data with a missing value
plot(X)
abline(v = 0.275, col = "gray", lwd = 2)
# 2) Add missing value
plot(X)
abline(v = 0.275, col = "gray", lwd = 2)
# 3) Show regions
Y <- rbind(X, c(0.275, colMeans(X)[2]))
isodepth(Y)
abline(v = 0.275, col = "gray", lwd = 2)
# 4) Show segment for imputation
isodepth(Y)
abline(v = 0.275, col = "gray", lwd = 2)
lines(rbind(c(0.275, -0.14), c(0.275, 0.775)), col = "red", lwd = 3)
# 5) Show average imputation
isodepth(Y)
abline(v = 0.275, col = "gray", lwd = 2)
lines(rbind(c(0.275, -0.14), c(0.275, 0.775)), col = "red", lwd = 3)
points(0.275, (-0.14 + 0.775) / 2, pch = 19, col = "red")
# 6) Impute
plot(X)
points(0.275, (-0.14 + 0.775) / 2, pch = 19, col = "red")
dev.off()
