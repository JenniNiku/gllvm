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
#'data(antTraits, package = "mvabund")
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
  y <- object$y
  Ntrials <- object$Ntrials
  args <- list(...)
  replace = TRUE
  if(!all(c("mu", "eta.mat")%in%names(args))){
    eta.mat = predict(object, type = "link")
    mu = predict(object, type = "response")
  }else{
    eta.mat <- args$eta.mat
    mu = args$mu
    if("replace"%in%names(args))replace = args$replace
  }
  
 
  ds.res = matrix(NA, n, p)
  
      if (object$family == "poisson") {
        b <- ppois(as.vector(y), as.vector(mu))
        a <- pmin(b, ppois(as.vector(unlist(y)) - 1, as.vector(mu)))
        u = a+(b-a)*matrix(runif(n*p),n,p,byrow=TRUE)
        
        if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
        if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
        ds.res <- qnorm(u)
      }
      if (object$family == "negative.binomial") {
        phis <- object$params$phi + 1e-05
        b <- pnbinom(as.vector(y), mu = as.vector(mu), size = 1 / rep(phis, each = n))
        a <- pmin(b,pnbinom(as.vector(unlist(y)) - 1, mu = as.vector(mu), size = 1 / as.vector(phis)))
        
        u = a+(b-a)*matrix(runif(n*p),n,p,byrow=TRUE)
        
        if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
        if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
        ds.res <- qnorm(u)
      }
      if (object$family == "gaussian") {
        phis <- object$params$phi
        b <- pnorm(as.vector(y), as.vector(mu), sd = rep(phis, each = n))
        a <- pmin(b, pnorm(as.vector(y), as.vector(mu), sd = rep(phis, each = n)))
        
        u = a+(b-a)*matrix(runif(n*p),n,p,byrow=TRUE)
        
        if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
        if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
        ds.res <- qnorm(u)
      }
      if (object$family == "gamma") {
        phis <- object$params$phi # - 1
        b <- pgamma(as.vector(y), shape = rep(phis, each = n), scale = as.vector(mu)/rep(phis, each = n))
        a <- pmin(b, pgamma(as.vector(y), shape = rep(phis, each = n), scale = as.vector(mu)/rep(phis, each = n)))
        
        u = a+(b-a)*matrix(runif(n*p),n,p,byrow=TRUE)
        
        if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
        if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
        ds.res <- qnorm(u)
      }
      if (object$family == "beta") {
        b <- pbeta(as.vector(y), shape1 = rep(object$params$phi, each = n)*as.vector(mu), shape2 = rep(object$params$phi, each = n)*(1-as.vector(mu)))
        a <- pmin(b, pbeta(as.vector(y), shape1 = rep(object$params$phi, each = n)*as.vector(mu), shape2 = rep(object$params$phi, each = n)*(1-as.vector(mu))))
        
        u = a+(b-a)*matrix(runif(n*p),n,p,byrow=TRUE)
        
        if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
        if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
        ds.res <- qnorm(u)
      }
      if (object$family == "betaH") {
        for (i in 1:n) {
          for (j in 1:p) {
            # a = 0; b = 1
            if(!is.na(y[i, j])){
              if(y[i, j]==0){
                b = 1 - binomial(link = object$link)$linkinv(eta.mat[i,p+j])
                a = 0
              } else {
                b <- a <- 1 - binomial(link = object$link)$linkinv(eta.mat[i,p+j]) + binomial(link = object$link)$linkinv(eta.mat[i,p+j])*pbeta(as.vector(unlist(y[i, j])), shape1 = object$params$phi[j]*mu[i, j], shape2 = object$params$phi[j]*(1-mu[i, j]))
              }
            
              u <- runif(n = 1, min = a, max = b)
              if(u==1&&replace) u=1-1e-16
              if(u==0&&replace) u=1e-16
              ds.res[i, j] <- qnorm(u)
            }
          }
        }
      }
      if (object$family == "orderedBeta") {
        for (i in 1:n) {
          for (j in 1:p) {
            # a = 0; b = 1
            if(!is.na(y[i, j])){
              if(y[i, j]==1){
                b = 1
                a = 1 - binomial(link = object$link)$linkinv(eta.mat[i,j] - object$params$zeta[j,2])
              } else if(y[i, j]==0){
                b = 1 - binomial(link = object$link)$linkinv(eta.mat[i,j] - object$params$zeta[j,1])
                a = 0
              } else {
                b <- a <- 1 - binomial(link = object$link)$linkinv(eta.mat[i,j] - object$params$zeta[j,1]) + (binomial(link = object$link)$linkinv(eta.mat[i,j] - object$params$zeta[j,1]) - binomial(link = object$link)$linkinv(eta.mat[i,j] - object$params$zeta[j,2]))*pbeta(as.vector(unlist(y[i, j])), shape1 = object$params$phi[j]*mu[i, j], shape2 = object$params$phi[j]*(1-mu[i, j]))
              }
            
            u <- try({runif(n = 1, min = a, max = b)})
            if(u==1&&replace) u=1-1e-16
            if(u==0&&replace) u=1e-16
            ds.res[i, j] <- qnorm(u)
            }
          }
        }
      }
      if (object$family == "exponential") {
        b <- pexp(as.vector(y), rate = 1/as.vector(mu))
        a <- pmin(b, pexp(as.vector(y), rate = 1/as.vector(mu)))
        
        u = a+(b-a)*matrix(runif(n*p),n,p,byrow=TRUE)
        
        if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
        if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
        ds.res <- qnorm(u)
      }
      if (object$family == "ZIP") {
        b <- pzip(as.vector(y), mu = as.vector(mu), sigma = rep(object$params$phi, each = n))
        a <- pmin(b, pzip(as.vector(y) - 1, mu = as.vector(mu), sigma = rep(object$params$phi, each = n)))
        
        u = a+(b-a)*matrix(runif(n*p),n,p,byrow=TRUE)
        
        if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
        if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
        ds.res <- qnorm(u)
      }
      if (object$family == "ZINB") {
        b <- pzinb(as.vector(y), mu = as.vector(mu), p = rep(object$params$phi, each = n), sigma = rep(object$params$ZINB.phi, each = n))
        a <- pmin(b, pzinb(as.vector(y) - 1, mu = as.vector(mu), p = rep(object$params$phi, each = n), sigma = rep(object$params$ZINB.phi, each = n)))

        u = a+(b-a)*matrix(runif(n*p),n,p,byrow=TRUE)

        if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
        if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
        ds.res <- qnorm(u)
      }
      if (object$family == "binomial") {
        b <- pbinom(as.vector(y), Ntrials, as.vector(mu))
        a <- pmin(b, pbinom(as.vector(y) - 1, Ntrials, as.vector(mu)))

        u = a+(b-a)*matrix(runif(n*p),n,p,byrow=TRUE)

        if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
        if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
        ds.res <- qnorm(u)
      }
      if (object$family == "tweedie") {
        phis <- object$params$phi + 1e-05
        b <- fishMod::pTweedie(as.vector(y), mu = as.vector(mu), phi = rep(phis, each = n), p = object$Power)
        a <- pmin(b, fishMod::pTweedie(as.vector(y) - 1, mu = as.vector(mu), phi = rep(phis, each = n), p = object$Power));
        
        anew =  ifelse((as.vector(y) - 1)<0, 0, a)
        
        u = anew+(b-anew)*matrix(runif(n*p),n,p,byrow=TRUE)
        
        if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
        if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
        ds.res <- qnorm(u)
      }
      if (object$family == "ordinal") {
        for (i in 1:n) {
          for (j in 1:p) {
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

  rownames(ds.res) <- rownames(y)
  colnames(ds.res) <- colnames(y)

  return(list(residuals = ds.res, linpred = eta.mat))
}

