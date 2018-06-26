#' @title Dunn-Smyth -residuals for gllvm model
#' @description Calculates Dunn-Smyth -residuals for gllvm model.
#'
#' @param object an object of class 'gllvm'.
#' @param ... not used.
#'
#' @details
#' Computes Dunn-Smyth residuals or randomized quantile residuals (Dunn and Smyth, 1996) for gllvm model.
#' For the observation \eqn{Y_{ij}} Dunn-Smyth residuals are defined as
#'
#' \deqn{r_{ij}=\Phi^{-1}(u_{ij}F_{ij}(y_{ij})  + (1-u_{ij})F_{ij}^-(y_{ij})),}
#'
#' where \eqn{\Phi(.)} and \eqn{F_{ij}(.)} are the cumulative probability functions of the standard normal
#' distribution, \eqn{F_{ij}^-(y))} is the limit as \eqn{F_{ij}(y)} is approached from the negative side, and \eqn{u_{ij}} has been
#' generated at random from the standard uniform distribution.
#'
#' @return A list containing \code{residuals} which is a matrix of residuals and \code{linpred} which is a matrix of linear predictors.
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @references
#'
#' Dunn, P. K., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics, 5, 236-244.
#'
#' Hui, F. K. C., Taskinen, S., Pledger, S., Foster, S. D., and Warton, D. I. (2015).  Model-based approaches to unconstrained ordination. Methods in Ecology and Evolution, 6:399-411.
#'
#' @examples
#' \dontrun{
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = "negative.binomial")
#'# residuals saved
#'res <- residuals(fit)
#'}
#'@export


residuals.gllvm <- function(object, ...) {
  n <- NROW(object$y)
  p <- NCOL(object$y)

  num.lv <- object$num.lv
  X <- object$X
  y <- object$y

  num.X <- ncol(object$X)
  num.T <- ncol(object$TR)

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
  if(!is.null(object$X) && is.null(object$TR)) eta.mat <- eta.mat + (X %*% matrix(t(object$params$Xcoef),num.X,p))
  if(!is.null(object$TR)) eta.mat <- eta.mat + matrix( object$X.design %*% c(object$params$B) ,n,p)
  if(object$row.eff!=FALSE) eta.mat <- eta.mat + matrix(object$params$row.params,n,p,byrow=FALSE)
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
        a <- pzip(as.vector(unlist(y[i, j]))-1, mu = mu[i, j], sigma = object$params$phi[j])
        b <- pzip(as.vector(unlist(y[i, j])), mu = mu[i, j], sigma = object$params$phi[j])
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
        a <- fishMod::pTweedie(as.vector(unlist(y[i, j])) - 1, mu = mu[i,j], phi = phis[j], p = object$Power); if((as.vector(unlist(y[i, j]))-1)<0) a<-0
        b <- fishMod::pTweedie(as.vector(unlist(y[i, j])), mu = mu[i,j], phi = phis[j], p = object$Power)
        u <- runif(n = 1, min = a, max = b)
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "ordinal") {
        
        probK <- NULL
        probK[1] <- pnorm(object$params$zeta[j,1]-eta.mat[i,j],log.p=FALSE)
        probK[max(y[,j])+1-min(y[,j])] <- 1 - pnorm(object$params$zeta[j,max(y[,j]) - min(y[,j])] - eta.mat[i,j])
        if(max(y[,j]) > 2) {
          j.levels <- 2:(max(y[,j])-min(y[,j]))#
          for(k in j.levels) { probK[k] <- pnorm(object$params$zeta[j,k] - eta.mat[i,j]) - pnorm(object$params$zeta[j,k - 1] - eta.mat[i,j]) }
        }
        probK <- c(0,probK)
        cumsum.b <- sum(probK[1:(y[i,j]+2-min(y[,j]))])
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

  return(list(residuals = ds.res, linpred=eta.mat))

}
