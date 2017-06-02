#' @title Dunn-Smyth -residuals for GLLVM model
#' @description Calculates Dunn-Smyth -residuals for GLLVM model and plots them.
#'
#' @param object An object of class 'gllvm'
#' @param plot Logical, construct a plot of residuals. Defaults to \code{FALSE}
#' @param by.index If \code{NULL} (default), residuals are plotted by linear predictors. If \code{by.index="site"}, residuals are plotted by row index. If \code{by.index="spp"} residuals are plotted by column index.
#' @param ...	Additional graphical arguments.
#' @details
#' Computes Dunn-Smyth residuals or randomized quantile residuals (Dunn and Smyth, 1996) for GLLVM model.
#' @return A list containing \code{residuals} which is a matrix of residuals and \code{linpred} which is a matrix of linear predictors.
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @references
#' Dunn, P. K., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics, 5, 236-244.
#' @examples
#' \dontrun{
#' library(mvabund) ## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit GLLVM model
#'fit <- gllvm(y = y, family = "negative.binomial")
#'# residuals saved and plotted by linear predictors
#'res <- residuals(fit, plot = TRUE, by.index = NULL)
#'# residuals plotted by site index
#'res <- residuals(fit, plot = TRUE, by.index = "site")
#'# residuals plotted by column index
#'res <- residuals(fit, plot = TRUE, by.index = "spp")
#'}
#'@export


residuals.gllvm <- function(object, plot=TRUE, by.index=NULL, xlim=NULL, ...) {
  n <- NROW(object$y)
  p <- NCOL(object$y)

  num.lv <- object$num.lv
  X <- object$X
  y <- object$y

  num.X <- 0
  if(!is.null(object$X)) {
    num.X <- dim(object$X)[2]
    X.new <- NULL
    for (i in 1:num.X) {
      if(!is.factor(object$X[,i])) {
        if(!is.null(object$TR)){ Xi <- scale(object$X[,i]) } else { Xi <- object$X[,i] }
        X.new <- cbind(X.new,Xi); if(!is.null(colnames(object$X)[i])) colnames(X.new)[dim(X.new)[2]] <- colnames(object$X)[i]
      }
      if(is.factor(object$X[,i])) {
        dum <- model.matrix(~object$X[,i]-1)
        colnames(dum) <- paste(names(object$X)[i],levels(object$X[,i]),sep="")
        if(family=="ordinal") dum <- dum[,-1]
        X.new <- cbind(X.new,dum)
      }
    }
    X <- data.matrix(X.new); num.X <- dim(X)[2]
  }
  num.T <- 0
  if(!is.null(object$TR)) {
    num.T <- dim(object$TR)[2]
    T.new <- NULL
    for (i in 1:num.T) {
      if(!is.factor(object$TR[,i])) { T.new <- cbind(T.new,scale(object$TR[,i])); colnames(T.new)[dim(T.new)[2]] <- names(object$TR)[i] }
      if(is.factor(object$TR[,i])){
        dum <- model.matrix(~object$TR[,i]-1)
        colnames(dum) <- paste(names(object$TR)[i],levels(object$TR[,i]),sep="")
        T.new <- cbind(T.new,dum)
      }
    }
    TR <- data.matrix(T.new); num.T <- dim(TR)[2];
  }
  pzip <- function(y,mu,sigma)
  {
    pp <- NULL
      if (y > -1) {
          cdf <- ppois(y, lambda = mu, lower.tail = TRUE, log.p = FALSE)
          cdf <- sigma + (1 - sigma) * cdf
          pp <- cdf
      }
      if (y < 0) { pp <- 0 }
    pp
  }

  if (is.null(object$offset)) { offset <- matrix(0, nrow = n, ncol = p)} else {offset <- object$offset}

  eta.mat <- matrix(object$params$beta0,n,p,byrow=TRUE) + offset
  if(!is.null(object$X)) eta.mat <- eta.mat + (X %*% matrix(t(object$params$Xcoef),num.X,p))
  if(!is.null(object$TR) && object$get.fourth) eta.mat <- eta.mat + matrix( kronecker(TR,X) %*% c(t(object$params$fourth)) ,n,p)
  if(!is.null(object$TR) && object$get.trait) eta.mat <- eta.mat +  matrix(object$params$Tcoef,n,num.T,byrow = TRUE) %*% t(TR)
  if(object$row.eff) eta.mat <- eta.mat + matrix(object$params$row.params,n,p,byrow=F)
  if(num.lv > 0) eta.mat <- eta.mat  + object$lvs %*% t(object$params$theta)
  mu <- exp(eta.mat); if(any(mu==0)) mu<-mu+1e-10
  if(object$family=="binomial") mu <- binomial(link = object$link)$linkinv(eta.mat)

  ds.res <- matrix(NA, n, p)
  rownames(ds.res) <- rownames(y)
  colnames(ds.res) <- colnames(y)
  for (i in 1:n) {
    for (j in 1:p) {
      if (object$family == "poisson") {
        a <- ppois(as.vector(unlist(y[i, j])) - 1, mu[i,j])
        b <- ppois(as.vector(unlist(y[i, j])), mu[i,j])
        u <- runif(n = 1, min = a, max = b)
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "negative.binomial") {
        phis <- object$params$phi + 1e-05
        a <- pnbinom(as.vector(unlist(y[i, j])) - 1, mu = mu[i, j], size = 1/phis[j])
        b <- pnbinom(as.vector(unlist(y[i, j])), mu = mu[i, j], size = 1/phis[j])
        u <- runif(n = 1, min = a, max = b)
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "ZIP") {
        a <- pzip(as.vector(unlist(y[i, j]))-1, mu = mu[i, j], sigma = object$params$p0[j])
        b <- pzip(as.vector(unlist(y[i, j])), mu = mu[i, j], sigma = object$params$p0[j])
        u <- runif(n = 1, min = a, max = b)
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "binomial") {
        a <- pbinom(as.vector(unlist(y[i, j])) - 1, 1, mu[i, j])
        b <- pbinom(as.vector(unlist(y[i, j])), 1, mu[i, j])
        u <- runif(n = 1, min = a, max = b)
        ds.res[i, j] <- qnorm(u)
      }

      if (object$family == "tweedie") {
        phis <- object$params$phi + 1e-05
        a <- pTweedie(as.vector(unlist(y[i, j])) - 1, mu = mu[i,j], phi = phis[j], p = object$Power); if((as.vector(unlist(y[i, j]))-1)<0) a<-0
        b <- pTweedie(as.vector(unlist(y[i, j])), mu = mu[i,j], phi = phis[j], p = object$Power)
        u <- runif(n = 1, min = a, max = b)
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "ordinal") { # This doesn't work yet
        probK <- NULL
        probK[1] <- pnorm(object$params$zeta[j,1]-eta.mat[i,j],log=FALSE)
        probK[max(y[,j])] <- 1 - pnorm(object$params$zeta[j,max(y[,j]) - 1] - eta.mat[i,j])
        if(max(y[,j]) > 2) {
          j.levels <- 2:(max(y[,j])-1)
          for(k in j.levels) { probK[k] <- pnorm(object$params$zeta[j,k] - eta.mat[i,j]) - pnorm(object$params$zeta[j,k - 1] - eta.mat[i,j]) }
        }
        probK <- c(0,probK)
        cumsum.b <- sum(probK[1:(y[i,j]+1)])
        cumsum.a <- sum(probK[1:(y[i,j])])
        u <- runif(n = 1, min = cumsum.a, max = cumsum.b)
        if (abs(u - 1) < 1e-05)
          u <- 1
        if (abs(u - 0) < 1e-05)
          u <- 0
        ds.res[i, j] <- qnorm(u)
      }
    }
  }

  if(plot) {
    if(is.null(xlim)) xlim=c(-max(eta.mat),max(eta.mat))
    par(mfrow = c(1, 1), cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1.2)
    if(is.null(by.index)) {plot(eta.mat,ds.res,xlab="linear predictors",ylab="Dunn-Smyth-residuals",xlim=xlim,...)
    } else{ if(by.index=="site") plot(rep(1:n,p),ds.res,xlab="site.index",ylab="Dunn-Smyth-residuals",...)
      if(by.index=="spp") plot(rep(1:p,each=n),ds.res,xlab="spp.index",ylab="Dunn-Smyth-residuals",...)}
  }
  return(list(residuals = ds.res, linpred=eta.mat))

}
