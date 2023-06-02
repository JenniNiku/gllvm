#' @title Dunn-Smyth residuals for gllvm model
#' @description Calculates Dunn-Smyth residuals for gllvm model.
#'
#' @param object an object of class 'gllvm'.
#' @param ... not used.
#'
#' @details
#' Computes Dunn-Smyth residuals (randomized quantile residuals, Dunn and Smyth, 1996) for gllvm model.
#' For the observation \eqn{Y_{ij}} Dunn-Smyth residuals are defined as
#'
#' \deqn{r_{ij}=\Phi^{-1}(u_{ij}F_{ij}(y_{ij})  + (1-u_{ij})F_{ij}^-(y_{ij})),}
#'
#' where \eqn{\Phi(.)} and \eqn{F_{ij}(.)} are the cumulative probability functions of the standard normal
#' distribution, \eqn{F_{ij}^-(y))} is the limit as \eqn{F_{ij}(y)} is approached from the negative side, and \eqn{u_{ij}} has been
#' generated at random from the standard uniform distribution.
#'
#' @return 
#'  \item{residuals }{matrix of residuals}
#'  \item{linpred }{matrix of linear predictors}
#'
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
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# residuals
#'res <- residuals(fit)
#'}
#'@export


residuals.gllvm <- function(object, ...) {
  n <- NROW(object$y)
  p <- NCOL(object$y)

  num.lv <- object$num.lv
  num.lv.c <- object$num.lv.c
  num.RR <- object$num.RR
  quadratic <- object$quadratic
  if(!is.null(object$X)) X <- as.matrix(object$X.design[1:n,]) else X=NULL
  y <- object$y

  if(object$row.eff != FALSE) {
    if(length(object$params$row.params) != nrow(object$y))
      object$params$row.params = c(object$TMBfn$env$data$dr0 %*% object$params$row.params)
  }
  
  num.X <- ncol(object$X.design)
  num.T <- ncol(object$TR)
  
  if (is.null(object$offset)) {
    offset <- matrix(0, nrow = n, ncol = p)
  } else {
    offset <- object$offset
  }

  eta.mat <- matrix(object$params$beta0, n, p, byrow = TRUE) + offset
  if (!is.null(object$X) && is.null(object$TR))
    eta.mat <- eta.mat + (object$X.design %*% matrix(t(object$params$Xcoef), num.X, p))
  if (!is.null(object$TR))
    eta.mat <- eta.mat + matrix(object$X.design %*% c(object$params$B) , n, p)
  if (object$row.eff != FALSE)
    eta.mat <- eta.mat + matrix(object$params$row.params, n, p, byrow = FALSE)
  if (num.lv > 0|num.lv.c>0|num.RR>0){
  lvs <- getLV(object)
  if(nrow(lvs)!=n) lvs = object$TMBfn$env$data$dr0%*%lvs # !!!
  if(num.lv>0){
  lvs[,grepl("^LV",colnames(object$lvs))] <- t(t(lvs[,grepl("^LV",colnames(object$lvs))])*object$params$sigma.lv[grepl("^LV",colnames(object$lvs))])
  }
  eta.mat <- eta.mat  + lvs %*% t(object$params$theta[,1:(num.lv+num.lv.c+num.RR),drop=F])
  }
  if(quadratic != FALSE){
   eta.mat <- eta.mat  + lvs^2 %*% t(object$params$theta[,-c(1:(num.lv+num.lv.c+num.RR)),drop=F])
  }
  
  if (!is.null(object$randomX))
    eta.mat <- eta.mat + (object$Xrandom %*% object$params$Br)
  
  mu <- exp(eta.mat)
  if (any(mu == 0))
    mu <- mu + 1e-10
  if (object$family == "binomial" || object$family == "beta")
    mu <- binomial(link = object$link)$linkinv(eta.mat)
  if (object$family == "gaussian")
    mu <- (eta.mat)
  
  ds.res <- matrix(NA, n, p)
  rownames(ds.res) <- rownames(y)
  colnames(ds.res) <- colnames(y)
  for (i in 1:n) {
    for (j in 1:p) {
      if(!is.na(object$y[i,j])){
      if (object$family == "poisson") {
        b <- ppois(as.vector(unlist(y[i, j])), mu[i, j])
        a <- min(b,ppois(as.vector(unlist(y[i, j])) - 1, mu[i, j]))
        u <- runif(n = 1, min = a, max = b)
        
        if(u==1) u=1-1e-16
        if(u==0) u=1e-16
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "negative.binomial") {
        phis <- object$params$phi + 1e-05
        b <- pnbinom(as.vector(unlist(y[i, j])), mu = mu[i, j], size = 1 / phis[j])
        a <- min(b,pnbinom(as.vector(unlist(y[i, j])) - 1, mu = mu[i, j], size = 1 / phis[j]))
        
        u <- runif(n = 1, min = a, max = b)
        if(u==1) u=1-1e-16
        if(u==0) u=1e-16
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "gaussian") {
        phis <- object$params$phi
        b <- pnorm(as.vector(unlist(y[i, j])), mu[i, j], sd = phis[j])
        a <- min(b,pnorm(as.vector(unlist(y[i, j])), mu[i, j], sd = phis[j]))
        
        u <- runif(n = 1, min = a, max = b)
        if(u==1) u=1-1e-16
        if(u==0) u=1e-16
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "gamma") {
        phis <- object$params$phi # - 1
        b <- pgamma(as.vector(unlist(y[i, j])), shape = phis[j], scale = mu[i, j]/phis[j])
        a <- min(b,pgamma(as.vector(unlist(y[i, j])), shape = phis[j], scale = mu[i, j]/phis[j]))
        
        u <- runif(n = 1, min = a, max = b)
        if(u==1) u=1-1e-16
        if(u==0) u=1e-16
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "beta") {
        b <- pbeta(as.vector(unlist(y[i, j])), shape1 = object$params$phi[j]*mu[i, j], shape2 = object$params$phi[j]*(1-mu[i, j]))
        a <- min(b,pbeta(as.vector(unlist(y[i, j])), shape1 = object$params$phi[j]*mu[i, j], shape2 = object$params$phi[j]*(1-mu[i, j])))
        u <- runif(n = 1, min = a, max = b)
        if(u==1) u=1-1e-16
        if(u==0) u=1e-16
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "exponential") {
        b <- pexp(as.vector(unlist(y[i, j])), rate = 1/mu[i, j])
        a <- min(b,pexp(as.vector(unlist(y[i, j])), rate = 1/mu[i, j]))
        u <- runif(n = 1, min = a, max = b)
        if(u==1) u=1-1e-16
        if(u==0) u=1e-16
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "ZIP") {
        b <- pzip(as.vector(unlist(y[i, j])), mu = mu[i, j], sigma = object$params$phi[j])
        a <- min(b,pzip(as.vector(unlist(y[i, j])) - 1, mu = mu[i, j], sigma = object$params$phi[j]))
        u <- runif(n = 1, min = a, max = b)
        if(u==1) u=1-1e-16
        if(u==0) u=1e-16
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "ZINB") {
        b <- pzinb(as.vector(unlist(y[i, j])), mu = mu[i, j], sigma = object$params$phi[j])
        a <- min(b,pzinb(as.vector(unlist(y[i, j])) - 1, mu = mu[i, j], sigma = object$params$phi[j]))
        u <- runif(n = 1, min = a, max = b)
        if(u==1) u=1-1e-16
        if(u==0) u=1e-16
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "binomial") {
        b <- pbinom(as.vector(unlist(y[i, j])), 1, mu[i, j])
        a <- min(b,pbinom(as.vector(unlist(y[i, j])) - 1, 1, mu[i, j]))
        u <- runif(n = 1, min = a, max = b)
        if(u==1) u=1-1e-16
        if(u==0) u=1e-16
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "tweedie") {
        phis <- object$params$phi + 1e-05
        b <- fishMod::pTweedie(as.vector(unlist(y[i, j])), mu = mu[i, j], phi = phis[j], p = object$Power)
        a <- min(b,fishMod::pTweedie(as.vector(unlist(y[i, j])) - 1, mu = mu[i, j], phi = phis[j], p = object$Power));
        if((as.vector(unlist(y[i, j])) - 1)<0)
          a<-0
        u <- runif(n = 1, min = a, max = b)
        if(u==1) u=1-1e-16
        if(u==0) u=1e-16
        ds.res[i, j] <- qnorm(u)
      }
      if (object$family == "ordinal") {
          if(object$zeta.struc == "species"){
            probK <- NULL
            probK[1] <- pnorm(object$params$zeta[j, 1] - eta.mat[i, j], log.p = FALSE)
            probK[max(y[, j]) + 1 - min(y[, j])] <- 1 - pnorm(object$params$zeta[j, max(y[, j]) - min(y[, j])] - eta.mat[i, j])
            if(length(unique(y[,j]))>2) {
              j.levels <- 2:(max(y[, j]) - min(y[, j]))#
              for (k in j.levels) {
                probK[k] <- pnorm(object$params$zeta[j, k] - eta.mat[i, j]) - pnorm(object$params$zeta[j, k - 1] - eta.mat[i, j])
              }
            }
            probK <- c(0, probK)
            cumsum.b <- sum(probK[1:(y[i,j]+ifelse(min(y[,j])==0,1,0) + 1)])
            cumsum.a <- min(cumsum.b, sum(probK[1:(y[i,j]+ifelse(min(y[,j])==0,1,0))]))
            u <- runif(n = 1, min = cumsum.a, max = cumsum.b)
            if (abs(u - 1) < 1e-05)
              u <- 1
            if (abs(u - 0) < 1e-05)
              u <- 0
            ds.res[i, j] <- qnorm(u)
          } else {
            probK <- NULL
            probK[1] <- pnorm(object$params$zeta[1] - eta.mat[i, j], log.p = FALSE)
            probK[max(y) + 1 - min(y)] <- 1 - pnorm(object$params$zeta[max(y) - min(y)] - eta.mat[i, j])
              levels <- 2:(max(y) - min(y))#
              for (k in levels) {
                probK[k] <- pnorm(object$params$zeta[k] - eta.mat[i, j]) - pnorm(object$params$zeta[k - 1] - eta.mat[i, j])
              }
            probK <- c(0, probK)
            cumsum.b <- sum(probK[1:(y[i,j]+ifelse(min(y)==0,1,0) + 1)])
            cumsum.a <- min(cumsum.b, sum(probK[1:(y[i,j]+ifelse(min(y)==0,1,0))]))
            u <- runif(n = 1, min = cumsum.a, max = cumsum.b)
            if (abs(u - 1) < 1e-05)
              u <- 1
            if (abs(u - 0) < 1e-05)
              u <- 0
            ds.res[i, j] <- qnorm(u)
          }

      }
      }
    }
  }

  return(list(residuals = ds.res, linpred = eta.mat))
}

