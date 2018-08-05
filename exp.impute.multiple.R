################################################################################
## File:             exp.impute.multiple.R                                    ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.07.2018                                               ##
##                                                                            ##
## Contains experiments on multiple depth-based imputation (Table 6).         ##
##                                                                            ##
################################################################################

# Functions ####################################################################
imputeEllP <- function(point, Sigma.inv){
  point.new <- point
  d <- length(point)
  index <- which(is.na(point))
  if (length(index) == 1){
    point.new[index] <- -point.new[-index] %*% Sigma.inv[index,-index] /
      Sigma.inv[index,index]
  }else{
    index <- which(is.na(point))
    A <- Sigma.inv[index,index]
    b <- -Sigma.inv[index,(1:d)[-index], drop = FALSE] %*% point[-index]
    point.new[index] <- solve(A) %*% b
  }
  return (point.new)
}

imputeEll <- function(X, num.iter = 100, X.result = NULL, P.result = NULL){
  X.prep <- X
  miss.label <- is.na(X.prep)
  miss = which(is.na(X),arr.ind=T)
  # TODO: The following line is dangerous when each observation has missingness
  X.prep[miss.label] <- matrix(rep(colMeans(X.prep, na.rm = TRUE),
                                   nrow(X.prep)), nrow = nrow(X.prep),
                               byrow = TRUE)[miss.label]
  errorX <- rep(NULL, num.iter)
  errorP <- rep(NULL, num.iter)
  for (i in 1:num.iter){
    mu.tmp <- colMeans(X.prep)
    X.prep <- t(t(X.prep) - mu.tmp)
    Sigma.tmp <- cov(X.prep)
    Inv.Sigma.tmp <- solve(Sigma.tmp)
    miss.rowi <- unique(miss[,1])
    X.new <- X.prep
    X.new[miss] <- NA
    X.new <- X.new[miss.rowi,]
    X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllP, Inv.Sigma.tmp))
    X.prep <- t(t(X.prep) + mu.tmp)
    # Diagnostics
    if (!is.null(X.result)){
      errorX[i] <- max(abs(X.prep - X.result))
    }
    if (!is.null(P.result)){
      errorP[i] <- max(abs(c(colMeans(X.prep), cov(X.prep)) - P.result))
    }
  }
  res <- list(X = X.prep)
  if (!is.null(X.result)){
    res <- c(res, list(errorX = errorX))
  }
  if (!is.null(P.result)){
    res <- c(res, list(errorP = errorP))
  }
  return (res)
}

imputeEllPRnd <- function(point, Sigma, depths){
  # Preliminary computations
  point.new <- point
  d <- length(point)
  n <- length(depths)
  index <- which(is.na(point))
  Sigma.inv <- solve(Sigma)
  # Conditional imputation
  if (length(index) == 1){
    point.new[index] <- -point.new[-index] %*% Sigma.inv[index,-index] /
      Sigma.inv[index,index]
  }else{
    index <- which(is.na(point))
    A <- Sigma.inv[index,index]
    b <- -Sigma.inv[index,(1:d)[-index], drop = FALSE] %*% point[-index]
    point.new[index] <- solve(A) %*% b
  }
  # Adding noise: draw from the conditional distribution
  depth <- as.numeric(1/(1 + point.new %*% Sigma.inv %*% point.new))
  if (depth <= depths[1]){
    rnd.depth <- depths[1]
  }else{
    if(depth >= depths[n]){
      rnd.depth <- depths[n]
    }else{
      n.new <- min(which(depth < unique(depths))) - 1
      if (n.new <= 1){
        rnd.depth <- unique(depths)[1]
      }else{
        # Create depth's conditional c.d.f.
        probs <- table(depths)[1:n.new]/sum(table(depths)[1:n.new])
        depths.new <- unique(depths)[1:n.new]
        probs <- probs / sqrt(1/depths.new - 1)^(d - 1) *
          sqrt((1/depths.new - 1) - (1/depth - 1))^(length(index) - 1) *
          (sqrt(1/depths.new - 1) / sqrt((1/depths.new - 1) - (1/depth - 1)))
        probs[n.new] <- probs[n.new - 1]
        probs <- cumsum(probs)
        # Draw the depth of the point
        rnd.prob <- runif(1, probs[1], probs[n.new])
        upperIndex <- min(which(probs > rnd.prob))
        rnd.depth <- (rnd.prob - probs[upperIndex - 1]) /
          (probs[upperIndex] - probs[upperIndex - 1]) *
          (depths[upperIndex] - depths[upperIndex - 1]) + depths[upperIndex - 1]
      }
    }
  }
  # Project the point into the sample's space
  d.na <- sum(is.na(point))
  dir <- rnorm(d.na)
  dir <- dir / sqrt(sum(dir^2))
  if (length(index) > 1){
    Sigma.cond <- Sigma[index,index,drop = FALSE] -
      Sigma[index,-index,drop = FALSE] %*%
      solve(Sigma[-index,-index,drop = FALSE]) %*%
      Sigma[-index,index,drop = FALSE]
    e <- eigen(Sigma.cond)
    L <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
    dir <- as.numeric(L %*% dir)
  }
  # Solve quadratic equation
  d.tmp <- 1/rnd.depth - 1
  r <- rep(0, d)
  r[is.na(point)] <- dir
  a.tmp <- r %*% Sigma.inv %*% r
  b.tmp <- 2 * point.new %*% Sigma.inv %*% r
  c.tmp <- point.new %*% Sigma.inv %*% point.new - d.tmp
  Dis <- b.tmp^2 - 4 * a.tmp * c.tmp
  if (Dis <= 0){
    alpha <- 0
  }else{
    x1 <- (-b.tmp + sqrt(Dis)) / (2 * a.tmp)
    x2 <- (-b.tmp - sqrt(Dis)) / (2 * a.tmp)
    if (x1 < 0 && x2 < 0){cat("Negative discriminant error!\n")}
    if (x1 > x2){alpha <- x1}else{alpha <- x2}
  }
  point.cur <- point.new + r * alpha
  # Optimize instead of solving quadratic equation
  #   min.cur <- 0
  #   max.cur <- 1000
  #   point.cur <- point.new
  #   while(abs(max.cur - min.cur) > sqrt(.Machine$double.eps)){
  #     mid.cur <- (min.cur + max.cur) / 2
  #     point.cur[is.na(point)] <- (point.new[is.na(point)] + dir * mid.cur)
  #     depth.cur <- as.numeric(1/(1 + point.cur %*% Sigma.inv %*% point.cur))
  #     if (depth.cur < rnd.depth){
  #       max.cur <- mid.cur
  #     }else{
  #       min.cur <- mid.cur
  #     }
  #   }
  return (point.cur)
}

imputeEllImproper <- function(X, num = 5, iter.burnin = 10, iter.skip = 10){
  X.prep <- X
  n <- nrow(X.prep)
  d <- ncol(X.prep)
  miss.label <- is.na(X.prep)
  miss = which(is.na(X), arr.ind=T)
  X.prep <- imputeEll(X, num.iter = 500)$X
  # means <- NULL
  # cors <- NULL
  Xs <- list("")
  for (iter in 1:(iter.burnin + (iter.skip + 1) * (num - 1) + 1)){
    mu.tmp <- colMeans(X.prep)
    X.prep <- t(t(X.prep) - mu.tmp)
    Sigma <- cov(X.prep)
    Inv.Sigma <- solve(Sigma)
    depths <- sort(1/(1 + (X.prep %*% Inv.Sigma * X.prep) %*% rep(1, d)))
    miss.rowi <- unique(miss[,1])
    X.new <- X.prep
    X.new[miss] <- NA
    X.new <- X.new[miss.rowi,]
    X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllPRnd, Sigma, depths))
    X.prep <- t(t(X.prep) + mu.tmp)
    # Diagnostics
    # means <- rbind(means, colMeans(X.prep))
    # cors <- rbind(cors, as.vector(cor(X.prep)))
    if (iter == 1){
      lines(cbind(depths, 1:length(depths)/length(depths)), col = "red")
    }else{
      if (iter == iter.burnin){
        cat("Burn-in passed ... ")
      }
      if ((iter > iter.burnin) &&
          ((iter - iter.burnin) %% (iter.skip + 1) == 1)){
        iter.sample <- floor((iter - iter.burnin) / (iter.skip + 1)) + 1
        Xs[[iter.sample]] <- X.prep
        cat(paste("Sample ", iter.sample, " taken ... ", sep = ""))
        lines(cbind(depths, 1:length(depths)/length(depths)), col = "orange")
      }else{
        lines(cbind(depths, 1:length(depths)/length(depths)), col = "green")
      }
    }
  }
  # Diagnostics
  # res <- list(Xs = Xs, means = means, cors = cors)
  cat(".\n", sep = "")
  # Diagnostics
  # return (res)
  return (Xs)
}

imputeEllProper <- function(X, m = 5, iter.burnin = 10){
  n <- nrow(X)
  d <- ncol(X)
  miss = which(is.na(X),arr.ind=T)
  Xs <- list("")
  cat("Sample", sep = "")
  for (iter in 1:m){
    X.prep <- imputeEll(X, num.iter = 500)$X
    sampleBi <- sample(x = 1:n, size = n, replace = TRUE)
    for (i in 1:iter.burnin){
      mu.tmp <- colMeans(X.prep[sampleBi,])
      Sigma <- cov(X.prep[sampleBi,])
      Inv.Sigma <- solve(Sigma)
      X.prep <- t(t(X.prep) - mu.tmp)
      depths <- sort(1/(1 + (X.prep %*% Inv.Sigma * X.prep) %*% rep(1, d)))
      miss.rowi <- unique(miss[,1])
      X.new <- X.prep
      X.new[miss] <- NA
      X.new <- X.new[miss.rowi,]
      X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllPRnd, Sigma, depths))
      X.prep <- t(t(X.prep) + mu.tmp)
      # Diagnostics
      if (i == 1){
        # lines(cbind(depths, 1:length(depths)/length(depths)), col = "red")
      }else{
        if (i == iter.burnin){
          # lines(cbind(depths, 1:length(depths)/length(depths)),
          #       col = "orange")
        }else{
          # lines(cbind(depths, 1:length(depths)/length(depths)),
          #       col = "green")
        }
      }
    }
    Xs[[iter]] <- X.prep
    cat(" ", iter, sep = "")
  }
  cat(".\n", sep = "")
  return (Xs)
}

poolResults <- function(Xs, beta.true){
  # Initialize
  betas <- NULL
  vars <- NULL
  n <- nrow(Xs[[1]])
  d <- ncol(Xs[[1]])
  m <- length(Xs)
  # Run regressions
  for (iter in 1:m){
    lmi <- lm(Xs[[iter]][,d] ~ Xs[[iter]][,1:(d - 1)])
    betas <- rbind(betas, coef(lmi))
    vars <- rbind(vars, sum(lmi$residuals^2) / lmi$df.residual *
                    diag(solve(t(cbind(rep(1, n), Xs[[iter]][,1:(d - 1)])) %*%
                                 cbind(rep(1, n), Xs[[iter]][,1:(d - 1)]))))
  }
  # Pool the estimates (Rubin's rule)
  Q <- colMeans(betas)
  U <- colMeans(vars)
  B <- colSums(t(t(betas) - Q)^2) / (m - 1)
  T <- U + (1 + 1 / m) * B
  # Calculate the coverage (adjusted Rubins' degrees of freedom)
  lambda <- (B + B / m) / T
  v.old <- (m - 1) / lambda^2
  v.com <- lmi$df.residual
  v.obs <- (v.com + 1) / (v.com + 3) * v.com * (1 - lambda)
  v <- v.old * v.obs / (v.old + v.obs)
  coverage <- (beta.true >= Q + qt(0.025, v) * sqrt(T)) &
    (beta.true <= Q + qt(0.975, v) * sqrt(T))
  iLength <- (qt(0.975, v) - qt(0.025, v)) * sqrt(T)
  return (c(Q, coverage, iLength, U, B))
}

prodNA <- function(X, noNA = 0.3, entire.rows = FALSE){
  n <- nrow(X)
  d <- ncol(X)
  X.miss <- X
  if (entire.rows){
    X[sample(1:(n * d), n * d * noNA)] <- NA
  }else{
    tmp.mat <- matrix(FALSE, nrow = n, ncol = d - 1)
    tmp.mat[sample(1:(n * (d - 1)), n * d * noNA)] <- TRUE
    tmp.nums <- rowSums(tmp.mat)
    for (i in 1:n){
      if (tmp.nums[i] > 0.5){
        X.miss[i,sample(1:d, tmp.nums[i])] <- NA
      }
    }
  }
  return (X.miss)
}
################################################################################

# Experiments - regression coefficient #########################################
library(MASS)
library(mvtnorm)
library(Amelia)
library(mice)
library(missMDA)
n <- 500
m <- 20
k <- 2000
noNA <- 0.3
# beta <- c(0.5, 1)
beta <- c(0.5, 1, 3)
# beta <- c(0.5, 1, 3, -2)
sd.reg <- 2
d <- length(beta)
Xs.T0.Amelia <- list("")
Xs.T0.mice <- list("")
Xs.T0.regPCA <- list("")
Xs.T0.Ell <- list("")
Xs.true <- list("")
Xs.miss <- list("")
# Perform experiments
betas.T0.Amelia <- NULL
betas.T0.mice <- NULL
betas.T0.regPCA <- NULL
betas.T0.Ell <- NULL
for (iter in 1:k){
  cat("Iteration ", iter, ": ... ", sep = "")
  # Generate data
  if (d == 2){
    mu <- 1
    Sigma <- 4
    X <- rnorm(n, mu, sqrt(Sigma))
  }
  if (d == 3){
    mu <- c(1, 1)
    Sigma <- matrix(c(1, 1, 1, 4), nrow = 2)
    X <- mvrnorm(n, mu, Sigma)
  }
  if (d == 4){
    mu <- c(1, 1, 1)
    Sigma <- matrix(c(1, 1, 1, 1, 4, 4, 1, 4, 10), nrow = 3, byrow = TRUE)
    X <- mvrnorm(n, mu, Sigma)
  }
  X.true <- cbind(X, cbind(rep(1, n), X) %*% beta + rnorm(n, 0, sd.reg))
  Xs.true[[iter]] <- X.true
  # Plot true depths
  X.miss <- prodNA(X.true, noNA, FALSE)
  Xs.miss[[iter]] <- X.miss
  cat(round(sum(is.na(X.miss))/(n*d) * 100, 2), " % missing ... \n", sep = "")
  # Impute
  Xs.T0.regPCA.tmp <- MIPCA(X.miss, nboot = m, ncp = 2, 
                            method.mi = "Bayes")$res.MI
  Xs.T0.regPCA[[iter]] <- list("")
  for (sub.iter in 1:m){
    Xs.T0.regPCA[[iter]][[sub.iter]] <- 
      data.matrix(Xs.T0.regPCA.tmp[[sub.iter]])
  }
  Xs.T0.Ell[[iter]] <- imputeEllProper(X.miss, m)
  Xs.T0.Amelia[[iter]] <- amelia(X.miss, m, p2s = 0)$imputations
  
  mice.iter <- mice(X.miss, m, method = "norm", maxit = 125, printFlag = FALSE)
  Xs.T0.mice[[iter]] <- list("")
  for (sub.iter in 1:m){
    Xs.T0.mice[[iter]][[sub.iter]] <- data.matrix(complete(mice.iter, sub.iter))
  }
  betas.T0.Amelia <- rbind(betas.T0.Amelia, 
                           poolResults(Xs.T0.Amelia[[iter]], beta))
  betas.T0.mice <- rbind(betas.T0.mice, 
                         poolResults(Xs.T0.mice[[iter]], beta))
  betas.T0.Ell <- rbind(betas.T0.Ell, 
                        poolResults(Xs.T0.Ell[[iter]], beta))
  betas.T0.regPCA <- rbind(betas.T0.regPCA, 
                           poolResults(Xs.T0.regPCA[[iter]], beta))
  if (iter > 1){
    cat("Amelia:     ", apply(betas.T0.Amelia[,1:d], 2, median, na.rm = TRUE),
        colMeans(betas.T0.Amelia[,(d + 1):(2*d)], na.rm = TRUE),
        apply(betas.T0.Amelia[,(2*d + 1):(3*d)], 2, median, na.rm = TRUE),
        apply(betas.T0.Amelia[,(3*d + 1):(4*d)], 2, median, na.rm = TRUE),
        apply(betas.T0.Amelia[,(4*d + 1):(5*d)], 2, median, na.rm = TRUE),
        ".\n")
    cat("mice:       ", apply(betas.T0.mice[,1:d], 2, median, na.rm = TRUE),
        colMeans(betas.T0.mice[,(d + 1):(2*d)], na.rm = TRUE),
        apply(betas.T0.mice[,(2*d + 1):(3*d)], 2, median, na.rm = TRUE),
        apply(betas.T0.mice[,(3*d + 1):(4*d)], 2, median, na.rm = TRUE),
        apply(betas.T0.mice[,(4*d + 1):(5*d)], 2, median, na.rm = TRUE), ".\n")
    cat("regPCA:     ", apply(betas.T0.regPCA[,1:d], 2, median, na.rm = TRUE),
        colMeans(betas.T0.regPCA[,(d + 1):(2*d)], na.rm = TRUE),
        apply(betas.T0.regPCA[,(2*d + 1):(3*d)], 2, median, na.rm = TRUE),
        apply(betas.T0.regPCA[,(3*d + 1):(4*d)], 2, median, na.rm = TRUE),
        apply(betas.T0.regPCA[,(4*d + 1):(5*d)], 2, median, na.rm = TRUE),".\n")
    cat("Ell:        ", apply(betas.T0.Ell[,1:d], 2, median, na.rm = TRUE),
        colMeans(betas.T0.Ell[,(d + 1):(2*d)], na.rm = TRUE),
        apply(betas.T0.Ell[,(2*d + 1):(3*d)], 2, median, na.rm = TRUE),
        apply(betas.T0.Ell[,(3*d + 1):(4*d)], 2, median, na.rm = TRUE),
        apply(betas.T0.Ell[,(4*d + 1):(5*d)], 2, median, na.rm = TRUE), ".\n")
  }
}
# Save the results
save.image("Xs_t0_d3_n500_m20_k2000_p30.RData")
# Pool the results
ests.T0.Amelia <- NULL
ests.T0.mice <- NULL
ests.T0.regPCA <- NULL
ests.T0.Ell <- NULL
for (iter in 1:k){
  ests.T0.Amelia <- rbind(ests.T0.Amelia, 
                          poolResults(Xs.T0.Amelia[[iter]], beta))
  ests.T0.mice <- rbind(ests.T0.mice, 
                        poolResults(Xs.T0.mice[[iter]], beta))
  ests.T0.regPCA <- rbind(ests.T0.regPCA, 
                          poolResults(Xs.T0.regPCA[[iter]], beta))
  ests.T0.Ell <- rbind(ests.T0.Ell, 
                       poolResults(Xs.T0.Ell[[iter]], beta))
}
boxplot(ests.T0.Amelia[,1:d])
grid()
boxplot(ests.T0.mice[,1:d])
grid()
boxplot(ests.T0.Ell[,1:d])
grid()
cat("Amelia:     ", apply(ests.T0.Amelia[,1:d], 2, median, na.rm = TRUE),
    colMeans(ests.T0.Amelia[,(d + 1):(2*d)], na.rm = TRUE),
    apply(ests.T0.Amelia[,(2*d + 1):(3*d)], 2, median, na.rm = TRUE),
    apply(ests.T0.Amelia[,(3*d + 1):(4*d)], 2, median, na.rm = TRUE),
    apply(ests.T0.Amelia[,(4*d + 1):(5*d)], 2, median, na.rm = TRUE), ".\n")
cat("mice:       ", apply(ests.T0.mice[,1:d], 2, median, na.rm = TRUE),
    colMeans(ests.T0.mice[,(d + 1):(2*d)], na.rm = TRUE),
    apply(ests.T0.mice[,(2*d + 1):(3*d)], 2, median, na.rm = TRUE),
    apply(ests.T0.mice[,(3*d + 1):(4*d)], 2, median, na.rm = TRUE),
    apply(ests.T0.mice[,(4*d + 1):(5*d)], 2, median, na.rm = TRUE), ".\n")
cat("regPCA:     ", apply(betas.T0.regPCA[,1:d], 2, median, na.rm = TRUE),
    colMeans(betas.T0.regPCA[,(d + 1):(2*d)], na.rm = TRUE),
    apply(betas.T0.regPCA[,(2*d + 1):(3*d)], 2, median, na.rm = TRUE),
    apply(betas.T0.regPCA[,(3*d + 1):(4*d)], 2, median, na.rm = TRUE),
    apply(betas.T0.regPCA[,(4*d + 1):(5*d)], 2, median, na.rm = TRUE), ".\n")
cat("Ell:        ", apply(ests.T0.Ell[,1:d], 2, median, na.rm = TRUE),
    colMeans(ests.T0.Ell[,(d + 1):(2*d)], na.rm = TRUE),
    apply(ests.T0.Ell[,(2*d + 1):(3*d)], 2, median, na.rm = TRUE),
    apply(ests.T0.Ell[,(3*d + 1):(4*d)], 2, median, na.rm = TRUE),
    apply(ests.T0.Ell[,(4*d + 1):(5*d)], 2, median, na.rm = TRUE), ".\n")
################################################################################
