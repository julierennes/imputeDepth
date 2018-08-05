################################################################################
## File:             plot3d.optim.R                                           ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.07.2018                                               ##
##                                                                            ##
## Contains the script for plotting 3D-pictures of the function to be         ##
## optimized for conditional depth maximization (Figure 1, supplement).       ##
##                                                                            ##
################################################################################

# Settings
xlim = c(-1.75, 4)
ylim = c(-2, 5.75)
xlim.large = c(-3.5, 7.5)
ylim.large = c(-3.5, 7.5)
z.cut = 3
# Sample points
library(mvtnorm)
d <- 3
n <- 250
mu <- c(1, 1, 1)
Sigma <- matrix(c(1, 1, 1, 1, 4, 4, 1, 4, 8), nrow = 3, byrow = TRUE)
set.seed(1)
X <- t(mu + t(rmvt(n, Sigma, 0)))
# Plot points in 3D
library(rgl)
rgl.open()
rgl.bg(color = "white")
rgl.spheres(X[,1], X[,2], X[,3], color = "red", radius = 0.1, alpha = 1)
rctg.outer <- matrix(c(xlim.large[1], ylim.large[2], 3,
                       xlim.large[1], ylim.large[1], 3,
                       xlim[1], ylim.large[1], 3,
                       xlim.large[1], ylim.large[2], 3,
                       xlim[1], ylim.large[1], 3,
                       xlim[1], ylim.large[2], 3,

                       xlim[2], ylim.large[2], 3,
                       xlim[2], ylim.large[1], 3,
                       xlim.large[2], ylim.large[1], 3,
                       xlim[2], ylim.large[2], 3,
                       xlim.large[2], ylim.large[1], 3,
                       xlim.large[2], ylim.large[2], 3,

                       xlim[1], ylim[1], 3,
                       xlim[1], ylim.large[1], 3,
                       xlim[2], ylim.large[1], 3,
                       xlim[1], ylim[1], 3,
                       xlim[2], ylim.large[1], 3,
                       xlim[2], ylim[1], 3,

                       xlim[1], ylim.large[2], 3,
                       xlim[1], ylim[2], 3,
                       xlim[2], ylim[2], 3,
                       xlim[1], ylim.large[2], 3,
                       xlim[2], ylim[2], 3,
                       xlim[2], ylim.large[2], 3), ncol = 3, byrow = TRUE)
rctg.inner <- matrix(c(xlim[1], ylim[2], 3,
                       xlim[1], ylim[1], 3,
                       xlim[2], ylim[1], 3,
                       xlim[1], ylim[2], 3,
                       xlim[2], ylim[1], 3,
                       xlim[2], ylim[2], 3), ncol = 3, byrow = TRUE)
rgl.material(color = "blue", shininess = 50.0)
rgl.triangles(rctg.outer[,1], rctg.outer[,2], rctg.outer[,3],
              col = "blue", alpha = 0.75)
rgl.triangles(rctg.inner[,1], rctg.inner[,2], rctg.inner[,3],
              col = "lightblue", alpha = 0.75)
rgl.viewpoint(theta = -120, phi = 35, fov = 60, zoom = 0.6)
# Function for 3D-surface plotting
depth.graph2 <- function(data, z,
                         depth_f = c("zonoid", "halfspace", "Mahalanobis",
                                     "projectionRandom", "spatial", "none"),
                         apoint = NULL, main = depth_f,
                         xlim = c(min(data[, 1]), max(data[, 1])),
                         ylim = c(min(data[, 2]), max(data[, 2])),
                         zlim = c(0, max(z)), xnum = 250, ynum = 250,
                         ptwd, exact, ...)
{
  x1 <- seq(xlim[1], xlim[2], length = xnum)
  x2 <- seq(ylim[1], ylim[2], length = ynum)
  x1.step <- (x1[2] - x1[1])
  x2.step <- (x2[2] - x2[1])
  all.points <- as.matrix(cbind(expand.grid(x1, x2), z))
  all.depths <- rep(0, nrow(all.points))
  df = depth_f
  if (!is.function(depth_f)) {
    depth_f = match.arg(depth_f)
    df = switch(depth_f, none = function(x, X) (0), zonoid = depth.zonoid,
                halfspace = depth.halfspace, Mahalanobis = function(x,
                                                                    X)
                  (.Mahalanobis_depth(x, colMeans(X), solve(cov(X)))),
                projectionRandom = depth.projection, spatial = depth.spatial)
    if (depth_f == "none")
      zlim = c(0, 1)
  }
  if (depth_f == "halfspace"){
    all.depths = df(all.points, data[,1:3], exact)
  }else{
    if (depth_f == "Mahalanobis"){
      all.depths = depth.Mahalanobis(all.points, data[,1:3])
    }else{
      all.depths = df(all.points, data[,1:3])
    }
  }
  z <- matrix(all.depths, ncol = ynum, nrow = xnum, byrow = FALSE)
  z.red <- NULL
  toLeft <- floor(ptwd / 2)
  toRight <- toLeft + (ptwd - 1)
  for (i in toLeft:toRight){
    for (j in toLeft:toRight){
      z.red <- c(z.red, as.integer((data[, 1] - x1[1])/x1.step + i) +
                   as.integer((data[, 2] - x2[1])/x2.step + j) * (xnum - 1))
    }
  }

  z.black <- ifelse(is.null(apoint) || !is.numeric(apoint) || length(apoint) !=
                      2, NA, as.integer((apoint[1] - x1[1])/x1.step + 1) +
                      as.integer((apoint[1] - x2[1])/x2.step + 1) * (xnum - 1))
  zfacet <- z[-1, -1] + z[-1, -ynum] + z[-xnum, -1] + z[-xnum, -ynum]
  z.indices.zero <- which(zfacet == 0)
  cols <- rep("gray", (xnum - 1) * (ynum - 1))
  cols <- replace(cols, z.indices.zero, ifelse(depth_f == "none",
                                               NA, "lightblue"))
  par(bg = "white")
  persp(x1, x2, z, xlim = xlim, ylim = ylim, zlim = zlim, r = 10,
        col = cols, main = main, ltheta = 55,
        shade = 0.55, ticktype = "detailed", xlab = "x", ylab = "y",
        zlab = "D(x|X)", border = NA, box = FALSE, ...)
}
library(ddalpha)
xlim = c(-1.75, 4)
ylim = c(-2, 5.75)
num = 500 # resolution
th = -10
dtheta = 15
phi = 35
z <- 3
depth.graph2(X, z, "halfspace", xnum = num, ynum = num, theta = th-dtheta,
             phi = phi, zlim = c(0,1.1), xlim = xlim, ylim = ylim, bold = T,
             ptwd = 9, exact = TRUE)
depth.graph2(X, z, "zonoid", xnum = num, ynum = num, theta = th-dtheta,
             phi = phi, zlim = c(0,1.1), xlim = xlim, ylim = ylim, bold = T,
             ptwd = 9)
depth.graph2(X, z, "Mahalanobis", xnum = num, ynum = num, theta = th-dtheta,
             phi = phi, zlim = c(0,1.1), xlim = xlim, ylim = ylim, bold = T,
             ptwd = 9)
