################################################################################
## File:             plot3d.realdata.R                                        ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.07.2018                                               ##
##                                                                            ##
## Contains the script for plotting real data sets in 3D.                     ##
##                                                                            ##
################################################################################

pairs3d <- function(data, col = "black"){
  # Normalize data
  X.normalized <- t((t(data) - apply(data, 2, min)) /
    (apply(data, 2, max) - apply(data, 2, min)))
  # Calculate radius of a single ball and shifting constant
  radius <- 0.015
  shift <- 0.5
  # Prapare pairs
  # (1, 2, 3)
  Xplot1 <- matrix(0, nrow = nrow(data), ncol = 3)
  Xplot1[,1] <- -X.normalized[,1] - shift
  Xplot1[,2] <- X.normalized[,2] + shift
  Xplot1[,3] <- X.normalized[,3] + shift
  # (2, 3, 4)
  Xplot2 <- matrix(0, nrow = nrow(data), ncol = 3)
  Xplot2[,1] <- X.normalized[,4] + shift
  Xplot2[,2] <- X.normalized[,2] + shift
  Xplot2[,3] <- X.normalized[,3] + shift
  # (1, 3, 4)
  Xplot3 <- matrix(0, nrow = nrow(data), ncol = 3)
  Xplot3[,1] <- -X.normalized[,1] - shift
  Xplot3[,2] <- -X.normalized[,4] - shift
  Xplot3[,3] <- X.normalized[,3] + shift
  # (2, 3, 4)
  Xplot4 <- matrix(0, nrow = nrow(data), ncol = 3)
  Xplot4[,1] <- -X.normalized[,1] - shift
  Xplot4[,2] <- X.normalized[,2] + shift
  Xplot4[,3] <- -X.normalized[,4] - shift
  # Plot
  rgl.open()
  rgl.bg(color = "white")
  # Draw pairs
  rgl.spheres(Xplot1[,1], Xplot1[,2], Xplot1[,3], radius = radius, color = col,
              alpha = 1)
  rgl.spheres(Xplot2[,1], Xplot2[,2], Xplot2[,3], radius = radius, color = col,
              alpha = 1)
  rgl.spheres(Xplot3[,1], Xplot3[,2], Xplot3[,3], radius = radius, color = col,
              alpha = 1)
  rgl.spheres(Xplot4[,1], Xplot4[,2], Xplot4[,3], radius = radius, color = col,
              alpha = 1)
  # Draw boxes
  abox <- matrix(c(0, 0, 0,
                   0, 0, 1,
                   0, 0, 1,
                   0, 1, 1,
                   0, 1, 1,
                   0, 1, 0,
                   0, 1, 0,
                   0, 0, 0,
                   1, 0, 0,
                   1, 0, 1,
                   1, 0, 1,
                   1, 1, 1,
                   1, 1, 1,
                   1, 1, 0,
                   1, 1, 0,
                   1, 0, 0,
                   0, 0, 0,
                   1, 0, 0,
                   0, 0, 1,
                   1, 0, 1,
                   0, 1, 1,
                   1, 1, 1,
                   0, 1, 0,
                   1, 1, 0),
                 ncol = 3, byrow = TRUE)
  abox1 <- cbind(-abox[,1] - shift, abox[,2] + shift, abox[,3] + shift)
  rgl.lines(abox1[,1], abox1[,2], abox1[,3], color = "black")
  abox2 <- cbind(abox[,1] + shift, abox[,2] + shift, abox[,3] + shift)
  rgl.lines(abox2[,1], abox2[,2], abox2[,3], color = "black")
  abox3 <- cbind(-abox[,1] - shift, -abox[,2] - shift, abox[,3] + shift)
  rgl.lines(abox3[,1], abox3[,2], abox3[,3], color = "black")
  abox4 <- cbind(-abox[,1] - shift, abox[,2] + shift, -abox[,3] - shift)
  rgl.lines(abox4[,1], abox4[,2], abox4[,3], color = "black")
  # Indicate axis
  abox1.middles <- (abox1[1:12 * 2 - 1,] + abox1[1:12 * 2,]) / 2
  texts1 <- c("3", "2", "3", "2", "3", "2", "3", "2", "1", "1", "1", "1")
  rgl.texts(abox1.middles[,1], abox1.middles[,2], abox1.middles[,3], texts1,
            color = "black")
  abox2.middles <- (abox2[1:12 * 2 - 1,] + abox2[1:12 * 2,]) / 2
  texts2 <- texts1
  texts2 <- gsub("1", "4", texts2)
  rgl.texts(abox2.middles[,1], abox2.middles[,2], abox2.middles[,3], texts2,
            color = "black")
  abox3.middles <- (abox3[1:12 * 2 - 1,] + abox3[1:12 * 2,]) / 2
  texts3 <- texts1
  texts3 <- gsub("2", "4", texts3)
  rgl.texts(abox3.middles[,1], abox3.middles[,2], abox3.middles[,3], texts3,
            color = "black")
  abox4.middles <- (abox4[1:12 * 2 - 1,] + abox4[1:12 * 2,]) / 2
  texts4 <- texts1
  texts4 <- gsub("3", "4", texts4)
  rgl.texts(abox4.middles[,1], abox4.middles[,2], abox4.middles[,3], texts4,
            color = "black")
}

library(rgl)
library(compositions)
# Banknotes
banknotes.all <- read.table(file = "data//data_banknotes.dat", header = FALSE,
                            sep = ",")
banknotes.one <- banknotes.all[banknotes.all[,5] == 1,1:3]
X <- banknotes.one[1:100,]
pairs(X)
plot3d(X, size = 5)
# Glass
glass.all <- read.table(file = "data//glass.dat", header = FALSE, sep = " ")
glass.nf <- glass.all[glass.all[,10] == 2,1:3]
X <- glass.nf
pairs(X)
plot3d(X, size = 5)
# Blood Transfusion
bloodtransfusion.all <- read.table(file = "data//bloodtransfusion_gp.dat",
                                   header = FALSE, sep = " ")
bloodtransfusion <- bloodtransfusion.all[,1:3]
X <- bloodtransfusion
pairs(X)
plot3d(X, size = 5)
