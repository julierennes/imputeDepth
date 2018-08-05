################################################################################
## File:             impute.functions.R                                       ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.12.2016                                               ##
##                                                                            ##
## Contains auxiliary routines for experiments on single imputation.          ##
##                                                                            ##
################################################################################

# Functions ####################################################################
imputeMean <- function(X){
  X.prep <- X
  miss.label <- is.na(X.prep)
  miss = which(is.na(X), arr.ind=T)
  # TODO: The following line is dangerous when each observation has missingness
  X.prep[miss.label] <- matrix(rep(colMeans(X.prep, na.rm = TRUE),
                                   nrow(X.prep)), nrow = nrow(X.prep),
                               byrow = TRUE)[miss.label]
  return (X.prep)
}

imputeEm <- function(X){
  s <- prelim.norm(X)
  thetahat <- em.norm(s, criterion = sqrt(.Machine$double.eps))
  params <- getparam.norm(s, thetahat)
  X.prep <- t(t(X) - params$mu)
  Inv.Sigma.tmp <- solve(params$sigma)
  miss.rowi = which(rowSums(is.na(X.prep)) > 0.5)
  X.new <- X.prep
  X.new <- X.new[miss.rowi,,drop=FALSE]
  X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllP, Inv.Sigma.tmp))
  return (t(t(X.prep) + params$mu))
}

imputReg=function(X){
  miss = which(is.na(X),arr.ind=T)
  X[is.na(X)]=0
  colnames(X)=paste("v",1:ncol(X),sep="")
  for (i in 1:100){
    for (j in 1:ncol(X)){
      if (sum(miss[,2]==j)>0){
        aa=lm(X[,j]~X[,-j])
        X[miss[miss[,2]==j,1],j]=predict(aa)[miss[miss[,2]==j,1]]
      }
    }}
  return(X)
}

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

imputeEll <- function(X, num.iter = 500, X.result = NULL, P.result = NULL){
  X.prep <- X
  miss.label <- is.na(X.prep)
  miss = which(is.na(X), arr.ind=T)
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

prodNA <- function(X, noNA = 0.3, entire.rows = FALSE){
  n <- nrow(X)
  d <- ncol(X)
  X.miss <- X
  if (entire.rows){
    mm <- matrix(FALSE, nrow = n, ncol = d)
    mm[sample(1:(n * d), n * d * noNA)] <- TRUE
    X.miss[mm] <- NA
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

imputeKnn <- function(X, k = -1, ks = 1:15, l = 10){
  # If k specified - impute
  n <- nrow(X)
  d <- ncol(X)
  if (k %in% 1:n){
    return (kNN(data.frame(X), k = k)[,1:d])
  }
  # Initial imputation
  miss.label <- is.na(X)
  miss.n <- sum(miss.label)
  X.prep <- X
  X.prep[miss.label] <- matrix(rep(colMeans(X.prep, na.rm = TRUE),
                                   nrow(X.prep)), nrow = nrow(X.prep),
                               byrow = TRUE)[miss.label]
  # Cross-validation
  errors <- rep(0, length(ks))
  for (i in 1:l){
    X.cvi <- as.matrix(X.prep)
    X.cvi[sample(1:(n * d), miss.n)] <- NA
    for (j in 1:length(ks)){
      X.impi <- kNN(data.frame(X.cvi), k = ks[j])[,1:d]
      errors[j] <- errors[j] + sum((X.impi - X.prep)^2)
    }
  }
  k.best <- which.min(errors)
  return (kNN(data.frame(X), k = k.best)[,1:d])
}

imp.depth.Mahalanobis <- function(X, num.iter = 100, X.start = NULL, alpha = 1,
                                  mu = NULL, Sigma = NULL){
  if (!is.null(mu) && !is.null(Sigma)){
    miss <- which(is.na(X), arr.ind=T)
    X.prep <- X
    X.prep <- t(t(X.prep) - mu)
    Inv.Sigma <- solve(Sigma)
    miss.rowi <- unique(miss[,1])
    X.new <- X.prep
    X.new <- X.new[miss.rowi,,drop=FALSE]
    X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllP, Inv.Sigma))
    X.prep <- t(t(X.prep) + mu)
    return(X.prep)
  }
  if (is.null(X.start)){
    X.prep <- X
    miss.label <- is.na(X.prep)
    # TODO: The following line is dangerous when each observation has missings
    X.prep[miss.label] <- matrix(rep(colMeans(X.prep, na.rm = TRUE),
                                     nrow(X.prep)), nrow = nrow(X.prep),
                                 byrow = TRUE)[miss.label] +
      rnorm(n = sum(miss.label), mean = 0, sd = 0.001)
  }else{
    X.prep <- X.start
  }
  miss <- which(is.na(X), arr.ind=T)
  for (i in 1:num.iter){
    if (alpha > 0.99){
      mu.tmp <- colMeans(X.prep)
      Sigma.tmp <- cov(X.prep)
    }else{
      mcd.est <- covMcd(X.prep, alpha = alpha)
      mu.tmp <- mcd.est$center
      Sigma.tmp <- mcd.est$cov
    }
    X.prep <- t(t(X.prep) - mu.tmp)
    Inv.Sigma.tmp <- solve(Sigma.tmp)
    miss.rowi <- unique(miss[,1])
    X.new <- X.prep
    X.new[miss] <- NA
    X.new <- X.new[miss.rowi,,drop=FALSE]
    X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllP, Inv.Sigma.tmp))
    X.prep <- t(t(X.prep) + mu.tmp)
  }
  return (X.prep)
}
# Functions (end) ##############################################################
