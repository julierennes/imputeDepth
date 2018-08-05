################################################################################
## File:             exp.imp.quantile.R                                       ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.07.2018                                               ##
##                                                                            ##
## Contains experiments on stochastic depth-based imputation (Table 5).       ##
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
  # Diagnostics
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
      # lines(cbind(depths, 1:length(depths)/length(depths)), col = "red")
    }else{
      if (iter == iter.burnin){
        cat("Burn-in passed ... ")
      }
      if ((iter > iter.burnin) &&
          ((iter - iter.burnin) %% (iter.skip + 1) == 1)){
        iter.sample <- floor((iter - iter.burnin) / (iter.skip + 1)) + 1
        Xs[[iter.sample]] <- X.prep
        cat(paste("Sample ", iter.sample, " taken ... ", sep = ""))
        # Diagnostics
        # lines(cbind(depths, 1:length(depths)/length(depths)), col = "orange")
      }else{
        # Diagnostics
        # lines(cbind(depths, 1:length(depths)/length(depths)), col = "green")
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
      # if (i == 1){
      #   lines(cbind(depths, 1:length(depths)/length(depths)), col = "red")
      # }else{
      #   if (i == iter.burnin){
      #     lines(cbind(depths, 1:length(depths)/length(depths)),
      #           col = "orange")
      #   }else{
      #     lines(cbind(depths, 1:length(depths)/length(depths)), col = "green")
      #   }
      # }
    }
    Xs[[iter]] <- X.prep
    cat(" ", iter, sep = "")
  }
  cat(".\n", sep = "")
  return (Xs)
}

imputeEllProper.wrong1 <- function(X){
  # Initialization
  n <- nrow(X)
  d <- ncol(X)
  miss = which(is.na(X),arr.ind=T)
  # Impute bootsstrapped sample
  sampleBi <- sample(x = 1:n, size = n, replace = TRUE)
  X.Bi <- X[sampleBi,]
  X.impBi <- imputeEllImproper(X = X.Bi, num = 1, iter.burnin = 25,
                               iter.skip = 1)$Xs[[1]]
  # Get bootstrap estimates
  mu.tmp <- colMeans(X.impBi)
  X.prep <- t(t(X.impBi) - mu.tmp)
  Sigma.tmp <- cov(X.prep)
  Inv.Sigma.tmp <- solve(Sigma.tmp)
  depths <- sort(1/(1 + (X.prep %*% Inv.Sigma.tmp * X.prep) %*% rep(1, d)))
  # Impute using bootstrap estimates
  miss.rowi <- unique(miss[,1])
  X.prep <- t(t(X) - mu.tmp)
  X.new <- X.prep
  X.new <- X.new[miss.rowi,]
  X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllPRnd, Sigma.tmp, depths))
  depths <- sort(1/(1 + (X.prep %*% Inv.Sigma.tmp * X.prep) %*% rep(1, d)))
  lines(depths, 1:length(depths)/length(depths), col = "orange")
  print(depths)
  X.prep <- t(t(X.prep) + mu.tmp)
  return (X.prep)
}

imputeEmImproper <- function(X, num = 5){
  n <- nrow(X)
  d <- ncol(X)
  s <- prelim.norm(X)
  thetahat <- em.norm(s)
  Xs <- list("")
  for (iter in 1:num){
    X.prep <- imp.norm(s, thetahat, X)
    # Diagnostics
    # mu.tmp <- colMeans(X.prep)
    # Sigma <- cov(X.prep)
    # Inv.Sigma <- solve(Sigma)
    # X.ctrd <- t(t(X.prep) - mu.tmp)
    # depths <- sort(1/(1 + (X.ctrd %*% Inv.Sigma * X.ctrd) %*% rep(1, d)))
    # lines(cbind(depths, 1:length(depths)/length(depths)), col = "green")
    Xs[[iter]] <- X.prep
  }
  return (Xs)
}

imputeRegPcaImproper <- function(X, ncp = 2, scale = TRUE, 
                                 method = c("Regularized", "EM"), 
                                 threshold = 1e-04, nboot = 100, 
                                 method.mi = "Boot", Lstart = 1000, 
                                 L = 100, verbose = FALSE){
  Estim_param <- function(X, scale = F, S, method = "Regularized", 
                          threshold = 10^(-6)) {
    if (sum(is.na(X)) > 0) {
      missing <- which(is.na(X))
      impute.data <- imputePCA(X, scale = scale, ncp = S, 
                               method = method, 
                               threshold = threshold)$completeObs
      reference <- FactoMineR::PCA(impute.data, scale.unit = scale, 
                                   graph = FALSE, ncp = S)
      rec <- reconst(reference, S)
      rec.pca <- as.matrix(X)
      rec.pca[missing] <- rec[missing]
      sd_resid <- apply(impute.data, 2, sd)
      resid <- rec.pca - rec
      if (scale) {
        resid <- sweep(resid, 2, sd_resid, FUN = "/")
      }
      sigma <- sqrt((sum((resid[-missing])^2)) / 
                      (nrow(X) * ncol(X) - (length(missing) + ncol(X) + S * 
                                              (nrow(X) - 1 + ncol(X) - S))))
    }
    else {
      reference <- FactoMineR::PCA(X, scale.unit = scale, 
                                   graph = FALSE, ncp = S)
      rec <- reconst(reference, ncp = S)
      rec.pca <- as.matrix(X)
      sd_resid <- apply(X, 2, sd)
      resid <- rec.pca - rec
      if (scale) {
        resid <- sweep(resid, 2, sd_resid, FUN = "/")
      }
      sigma <- sqrt((sum((resid)^2)) / 
                      (nrow(X) * ncol(X) - (ncol(X) + S * 
                                              (nrow(X) - 1 + ncol(X) - S))))
    }
    phi <- diag(c((reference$eig[1:S, 1] - sigma^2)/reference$eig[1:S,1], 
                  rep(0, min(nrow(X) - 1, ncol(X)) - S)))
    estimation <- list(sigma = sigma, phi = phi)
    return(estimation)
  }
  BayesMIPCA <- function(X.na, nboot = 100, ncp, L = 100, verbose = T, 
                         Lstart = 1000, scale = F) {
    if (verbose) {
      cat("Multiple Imputation using Bayesian PCA using", 
          nboot, "imputed arrays", "\n")
    }
    if (scale) {
      sd.init <- apply(X.na, 2, FUN = sd, na.rm = T)
      X.na <- sweep(X.na, 2, sd.init, FUN = "/")
    }
    res.MI <- res.Over <- list()
    Mean <- colMeans(X.na, na.rm = TRUE)
    X <- sweep(X.na, 2, Mean, FUN = "-")
    X.tild <- imputePCA(X, ncp = ncp, method = "Regularized", 
                        scale = F)$fitt
    sigma <- Estim_param(X = X, S = ncp)[["sigma"]]
    for (l in 1:(Lstart + nboot * L)) {
      X[is.na(X.na)] <- X.tild[is.na(X.na)] + rnorm(sum(is.na(X.na)), 
                                                    mean = 0, sd = sigma)
      X <- sweep(X, 2, Mean, FUN = "+")
      if (l %in% c(Lstart + (1:nboot) * L)) {
        res.MI[[(l - Lstart)/L]] <- X
        res.Over[[(l - Lstart)/L]] <- sweep(X.tild + 
                                              rnorm(nrow(X.na) * ncol(X.na), 
                                                    mean = 0, sd = sigma), 2, 
                                            Mean, FUN = "+")
        if (verbose) {
          cat(paste((l - Lstart)/L, "...", sep = ""))
        }
      }
      Mean <- colMeans(X)
      X <- sweep(X, 2, Mean, FUN = "-")
      res.pca <- FactoMineR::PCA(X, ncp = ncol(X), scale.unit = F, 
                                 graph = F)
      param <- Estim_param(X = X, S = ncp)
      phi <- param[["phi"]]
      sigma <- param[["sigma"]]
      Moy <- res.pca[["svd"]][["U"]] %*% diag(res.pca[["svd"]][["vs"]]) %*% 
        phi %*% t(res.pca[["svd"]][["V"]])
      Var <- diag(rep(sigma^2/min(nrow(X) - 1, ncol(X)) * 
                        sum(phi), ncol(X)))
      X.tild <- Moy + rmvnorm(nrow(X.na), sigma = Var)
    }
    class(res.MI) <- "BayesMIPCA"
    if (scale) {
      res.MI <- lapply(res.MI, FUN = function(res) {
        sweep(res, 2, sd.init, FUN = "*")
      })
      res.Over <- lapply(res.Over, FUN = function(res) {
        sweep(res, 2, sd.init, FUN = "*")
      })
    }
    call <- list(X = X.na, nboot = nboot, ncp = ncp, L = L, 
                 verbose = verbose, Lstart = Lstart, scale = scale)
    return(list(res.MI = lapply(res.MI, as.data.frame), res.Over = res.Over, 
                call = call))
  }
  if (ncp == 0) 
    stop("No variability can be computed with ncp=0 dimension")
  res.Over <- list()
  method <- match.arg(method, c("Regularized", "regularized", 
                                "EM", "em"), several.ok = T)[1]
  method <- tolower(method)
  missing <- which(is.na(X))
  inputeData <- imputePCA(X, scale = scale, ncp = ncp, method = method, 
                          threshold = threshold)$completeObs
  if (method.mi %in% c("Bayes", "bayes")) {
    res <- BayesMIPCA(X.na = X, nboot = nboot, ncp = ncp, 
                      L = L, verbose = verbose, Lstart = Lstart, scale = scale)
    res.MI <- res$res.MI
    res.Over <- res$res.Over
  }
  else if (method.mi %in% c("Boot", "boot")) {
    if (verbose) {
      cat("Multiple Imputation using bootstrap", method, 
          "PCA using", nboot, "imputed arrays", "\n")
    }
    reference <- FactoMineR::PCA(inputeData, scale.unit = scale, 
                                 graph = FALSE, ncp = ncp)
    rec <- FactoMineR::reconst(reference, ncp)
    rec.pca <- as.matrix(X)
    rec.pca[missing] <- rec[missing]
    resid <- rec.pca - rec
    sdResid <- apply(inputeData, 2, sd)
    if (scale) 
      resid <- t(t(resid)/sdResid)
    sigma <- sqrt((sum((resid[-missing])^2)) / 
                    (nrow(X) * ncol(X) - (length(missing) + ncol(X) + ncp * 
                                            (nrow(X) - 1 + ncol(X) - ncp))))
    rownames(rec.pca) <- rownames(X)
    res.MI <- list()
    for (i in 1:nboot) {
      if (verbose) {
        cat(paste(i, "...", sep = ""))
      }
      resid.star <- matrix(rnorm(nrow(X) * ncol(X), 0, 
                                 sigma), ncol = ncol(X))
      if (scale) 
        resid.star <- t(t(resid.star) * sdResid)
      resid.star[missing] <- NA
      # Pavlo: To do it improper
      # Xstar <- rec + resid.star - matrix(mean(resid.star, na.rm = TRUE), 
      #                                    ncol = ncol(resid.star), 
      #                                    nrow = nrow(resid.star))
      Xstar <- rec
      Xstar[missing] <- NA
      acpboot <- FactoMineR::PCA(imputePCA(Xstar, scale = scale, 
                                           ncp = ncp, method = method, 
                                           threshold = threshold)$completeObs, 
                                 scale.unit = scale, ncp = ncp, graph = FALSE)
      residstar2 <- matrix(rnorm(nrow(X) * ncol(X), 0, 
                                 sigma), ncol = ncol(X))
      if (scale) 
        residstar2 <- t(t(residstar2) * sdResid)
      rec.pca[missing] <- (reconst(acpboot, ncp) + residstar2)[missing]
      res.Over[[i]] <- reconst(acpboot, ncp) + residstar2
      res.MI[[i]] <- rec.pca
      dimnames(res.MI[[i]]) = list(rownames(X), colnames(X))
      res.MI[[i]] <- as.data.frame(res.MI[[i]])
    }
  }
  else {
    stop("method.mi is misspecified")
  }
  if (verbose) {
    cat("\ndone!")
  }
  result = list(res.MI = res.MI, res.imputePCA = inputeData, 
                call = list(X = X, ncp = ncp, missing = missing, nboot = nboot, 
                            scale = scale, L = L, verbose = verbose, Lstart = Lstart, 
                            res.Over = res.Over))
  class(result) <- c("MIPCA", "list")
  return(result)
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
  return (c(Q, coverage, iLength))
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
library(norm)
library(missMDA)
library(FactoMineR)
n <- 500
m <- 1
k <- 2000
noNA <- 0.3
df <- 3
d <- 4
probs = c(0.5, 0.75, 0.85, 0.90, 0.95, 0.975, 0.99, 0.995)
Xs.T0.norm <- list("")
Xs.T0.regPCA <- list("")
Xs.T0.EllImproper <- list("")
Xs.true <- list("")
Xs.miss <- list("")
# Perform experiments
qs.T0.norm <- NULL
qs.T0.regPCA <- NULL
qs.T0.EllImproper <- NULL
# Initialize random generator
rngseed(1234567)
for (iter in 1:k){
  cat("Iteration ", iter, ": ... ", sep = "")
  # Generate data
  if (d == 2){
    mu <- c(1, 1)
    Sigma <- matrix(c(1, 1, 1, 4), nrow = 2)
  }
  if (d == 3){
    mu <- c(1, 1, 1)
    Sigma <- matrix(c(1, 1, 1, 1, 4, 4, 1, 4, 10), nrow = 3, byrow = TRUE)
  }
  if (d == 4){
    mu <- c(-1, -1, -1, -1)
    Sigma <- matrix(c(0.5, 0.5, 1, 1, 0.5, 1, 1, 1, 1, 1, 4, 4, 1, 1, 4, 10),
                    nrow = 4, byrow = TRUE)
  }
  X.true <- t(mu + t(rmvt(n, Sigma, df)))
  # Plot true depths
  Y <- X.true
  n <- nrow(Y)
  Y.prep <- t(t(Y) - colMeans(Y))
  s <- cov(Y.prep)
  s.inv <- solve(s)
  Y.depths <- 1/(1 + (Y.prep %*% s.inv * Y.prep) %*% rep(1, d))
  Y.depths <- sort(c(0, Y.depths, 1))
  plot(cbind(Y.depths, 0:(n+1)/(n+1)), type = "l", col = "black")
  grid()
  X.miss <- prodNA(X.true, noNA, FALSE)
  Xs.miss[[iter]] <- X.miss
  cat(round(sum(is.na(X.miss))/(n*d) * 100, 2), " % missing ... \n", sep = "")
  # Impute
  Xs.T0.EllImproper[[iter]] <- imputeEllImproper(X = X.miss, num = m,
                                                 iter.burnin = 10,
                                                 iter.skip = 10)
  Xs.T0.norm[[iter]] <- imputeEmImproper(X = X.miss, num = m)
  Xs.T0.regPCA[[iter]] <- imputeRegPcaImproper(X = X.miss, nboot = m)
  qs.T0.EllImproper <- rbind(qs.T0.EllImproper,
                             quantile(x = Xs.T0.EllImproper[[iter]][[1]][,1],
                                      probs = probs,
                                      names = FALSE))
  qs.T0.EllImproper <- rbind(qs.T0.EllImproper, probs)
  qs.T0.regPCA <- rbind(qs.T0.regPCA,
                      quantile(x = Xs.T0.regPCA[[iter]]$res.MI[[1]][,1],
                               probs = probs,
                               names = FALSE))
  qs.T0.norm <- rbind(qs.T0.norm,
                      quantile(x = Xs.T0.norm[[iter]][[1]][,1],
                               probs = probs,
                               names = FALSE))
  cat("done: norm:        ", apply(qs.T0.norm, 2, median, na.rm = TRUE), "\n")
  cat("done: regPCA       ", apply(qs.T0.regPCA, 2, median, na.rm = TRUE), "\n")
  cat("done: EllIpmroper: ", apply(qs.T0.EllImproper, 2, median, na.rm = TRUE),
      "\n")
}
boxplot(qs.T0.norm, main = "norm", names = paste(probs, sep = ""))
grid()
boxplot(qs.T0.regPCA, main = "regPCA", names = paste(probs, sep = ""))
grid()
boxplot(qs.T0.EllImproper, main = "EllImproper", names = paste(probs, sep = ""))
grid()
################################################################################
