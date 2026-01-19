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
    if(any(object$family == "ordinal") & !all(object$family == "ordinal")) mu <- mu[1,,]
  }else{
    eta.mat <- args$eta.mat
    mu = args$mu
    if("replace"%in%names(args))replace = args$replace
  }
  
  if(length(object$family) != p) object$family = rep(object$family,p) [1:p]
  
  ds.res = matrix(NA, n, p)
  if(any(object$family %in% "orderedBeta") & object$zeta.struc == "common") {kz =2} else {kz=0} # For indexing
  
  
  if (any(object$family == "poisson")) {
    p_f = sum(object$family == "poisson")
    b  <- ppois(as.vector(y[,object$family == "poisson"]), as.vector(mu[,object$family == "poisson"]))
    a <- pmin(b, ppois(as.vector(unlist(y[,object$family == "poisson"])) - 1, as.vector(mu[,object$family == "poisson"])))
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "poisson"] <- matrix(qnorm(u), n, p_f)
  }
  if (any(object$family == "negative.binomial")) {
    p_f = sum(object$family == "negative.binomial")
    phis <- object$params$phi[object$family == "negative.binomial"] + 1e-05
    b <- pnbinom(as.vector(y[,object$family == "negative.binomial"]), mu = as.vector(mu[,object$family == "negative.binomial"]), size = 1 / rep(phis, each = n))
    a <- pmin(b,pnbinom(as.vector(unlist(y[,object$family == "negative.binomial"])) - 1, mu = as.vector(mu[,object$family == "negative.binomial"]), size = 1 / rep(phis, each = n)))
    
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "negative.binomial"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "negative.binomial1")) {
    p_f = sum(object$family == "negative.binomial1")
    phis <- object$params$phi[object$family == "negative.binomial1"] + 1e-05
    b <- pnbinom(as.vector(y[,object$family == "negative.binomial1"]), mu = as.vector(mu[,object$family == "negative.binomial1"]), size = as.vector(mu[,object$family == "negative.binomial1"])*rep(phis, each = n))
    a <- pmin(b,  pnbinom(as.vector(y[,object$family == "negative.binomial1"]), mu = as.vector(mu[,object$family == "negative.binomial1"]), size = as.vector(mu[,object$family == "negative.binomial1"])*rep(phis, each = n)))
    
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "negative.binomial1"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "gaussian")) {
    p_f = sum(object$family == "gaussian")
    phis <- object$params$phi[object$family == "gaussian"]
    b <- pnorm(as.vector(y[,object$family == "gaussian"]), as.vector(mu[,object$family == "gaussian"]), sd = rep(phis, each = n))
    a <- pmin(b, pnorm(as.vector(y[,object$family == "gaussian"]), as.vector(mu[,object$family == "gaussian"]), sd = rep(phis, each = n)))
    
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "gaussian"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "gamma")) {
    p_f = sum(object$family == "gamma")
    phis <- object$params$phi[object$family == "gamma"] # - 1
    b <- pgamma(as.vector(y[,object$family == "gamma"]), shape = rep(phis, each = n), scale = as.vector(mu[,object$family == "gamma"])/rep(phis, each = n))
    a <- pmin(b, pgamma(as.vector(y[,object$family == "gamma"]), shape = rep(phis, each = n), scale = as.vector(mu[,object$family == "gamma"])/rep(phis, each = n)))
    
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "gamma"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "beta")) {
    p_f = sum(object$family == "beta")
    b <- pbeta(as.vector(y[,object$family == "beta"]), shape1 = rep(object$params$phi[object$family == "beta"], each = n)*as.vector(mu[,object$family == "beta"]), shape2 = rep(object$params$phi[object$family == "beta"], each = n)*(1-as.vector(mu[,object$family == "beta"])))
    a <- pmin(b, pbeta(as.vector(y[,object$family == "beta"]), shape1 = rep(object$params$phi[object$family == "beta"], each = n)*as.vector(mu[,object$family == "beta"]), shape2 = rep(object$params$phi[object$family == "beta"], each = n)*(1-as.vector(mu[,object$family == "beta"]))))
    
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "beta"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "betaH")) {
    p_f = sum(object$family == "betaH")
    bH_ind <- c(1:p)[object$family == "betaH"]
    for (i in 1:n) {
      for (j in 1:p_f) {
        # a = 0; b = 1
        if(!is.na(y[i, bH_ind[j]])){
          if(y[i, bH_ind[j]]==0){
            b = 1 - binomial(link = object$link)$linkinv(eta.mat[i,p+j])
            a = 0
          } else {
            b <- a <- 1 - binomial(link = object$link)$linkinv(eta.mat[i,p+j]) + binomial(link = object$link)$linkinv(eta.mat[i,p+j])*pbeta(as.vector(unlist(y[i, bH_ind[j]])), shape1 = object$params$phi[bH_ind[j]]*mu[i, bH_ind[j]], shape2 = object$params$phi[bH_ind[j]]*(1-mu[i, bH_ind[j]]))
          }
          
          u <- runif(n = 1, min = a, max = b)
          if(u==1&&replace) u=1-1e-16
          if(u==0&&replace) u=1e-16
          ds.res[i, bH_ind[j]] <- qnorm(u)
        }
      }
    }
  }
  if (any(object$family == "orderedBeta")) {
    if(!is.matrix(object$params$zeta)) object$params$zeta = matrix(object$params$zeta, nrow=1)
    p_f = sum(object$family == "orderedBeta")
    oB_ind <- c(1:p)[object$family == "orderedBeta"]
    for (i in 1:n) {
      for (j in 1:p_f) {
        # a = 0; b = 1
        if(!is.na(y[i, oB_ind[j]])){
          if(y[i, oB_ind[j]]==1){
            b = 1
            a = 1 - binomial(link = object$link)$linkinv(eta.mat[i,oB_ind[j]] - object$params$zeta[min(nrow(object$params$zeta),oB_ind[j]),2])
          } else if(y[i, oB_ind[j]]==0){
            b = 1 - binomial(link = object$link)$linkinv(eta.mat[i,oB_ind[j]] - object$params$zeta[min(nrow(object$params$zeta),oB_ind[j]),1])
            a = 0
          } else {
            b <- a <- 1 - binomial(link = object$link)$linkinv(eta.mat[i,oB_ind[j]] - object$params$zeta[min(nrow(object$params$zeta),oB_ind[j]),1]) + (binomial(link = object$link)$linkinv(eta.mat[i,oB_ind[j]] - object$params$zeta[min(nrow(object$params$zeta),oB_ind[j]),1]) - binomial(link = object$link)$linkinv(eta.mat[i,oB_ind[j]] - object$params$zeta[min(nrow(object$params$zeta),oB_ind[j]),2]))*pbeta(as.vector(unlist(y[i, oB_ind[j]])), shape1 = object$params$phi[oB_ind[j]]*mu[i, oB_ind[j]], shape2 = object$params$phi[oB_ind[j]]*(1-mu[i, oB_ind[j]]))
          }
          
          u <- try({runif(n = 1, min = a, max = b)})
          if(u==1&&replace) u=1-1e-16
          if(u==0&&replace) u=1e-16
          ds.res[i, oB_ind[j]] <- qnorm(u)
        }
      }
    }
  }
  if (any(object$family == "exponential")) {
    p_f = sum(object$family == "exponential")
    b <- pexp(as.vector(y[,object$family == "exponential"]), rate = 1/as.vector(mu[,object$family == "exponential"]))
    a <- pmin(b, pexp(as.vector(y[,object$family == "exponential"]), rate = 1/as.vector(mu[,object$family == "exponential"])))
    
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "exponential"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "ZIP")) {
    p_f = sum(object$family == "ZIP")
    mu_f = exp(eta.mat[,object$family == "ZIP", drop=FALSE])
    b <- pzip(as.vector(y[,object$family == "ZIP"]), mu = as.vector(mu_f), sigma = rep(object$params$phi[object$family == "ZIP"], each = n))
    a <- pmin(b, pzip(as.vector(y[,object$family == "ZIP"]) - 1, mu = as.vector(mu_f), sigma = rep(object$params$phi[object$family == "ZIP"], each = n)))
    
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "ZIP"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "ZINB")) {
    p_f = sum(object$family == "ZINB")
    mu_f = exp(eta.mat[,object$family == "ZINB", drop=FALSE])
    b <- pzinb(as.vector(y[,object$family == "ZINB"]), mu = as.vector(mu_f), p = rep(object$params$phi[object$family == "ZINB"], each = n), sigma = rep(object$params$ZINB.phi[object$family == "ZINB"], each = n))
    a <- pmin(b, pzinb(as.vector(y[,object$family == "ZINB"]) - 1, mu = as.vector(mu_f), p = rep(object$params$phi[object$family == "ZINB"], each = n), sigma = rep(object$params$ZINB.phi[object$family == "ZINB"], each = n)))
    
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "ZINB"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "binomial")) {
    p_f = sum(object$family == "binomial")
    if(length(Ntrials)==1)Ntrials <- rep(Ntrials,p_f)
    if(length(Ntrials)==p)Ntrials <- rep(Ntrials[object$family == "binomial"], each = n)
    if(is.matrix(Ntrials))Ntrials <- c(Ntrials[,object$family == "binomial"])
    
    b <- pbinom(as.vector(y[,object$family == "binomial"]), Ntrials, as.vector(mu[,object$family == "binomial"]))
    a <- pmin(b, pbinom(as.vector(y[,object$family == "binomial"]) - 1, Ntrials, as.vector(mu[,object$family == "binomial"])))
    
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "binomial"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "ZIB")) {
    p_f = sum(object$family == "ZIB")
    if(length(Ntrials)==1)Ntrials <- rep(Ntrials,p_f)
    if(length(Ntrials)==p)Ntrials <- rep(Ntrials, each = n)
    if(is.matrix(Ntrials))Ntrials <- c(Ntrials)
    
    b <- pzib(as.vector(y[,object$family == "ZIB"]), Ntrials = Ntrials, mu = as.vector(mu[,object$family == "ZIB"]), sigma = rep(object$params$phi[object$family == "ZIB"], each = n))
    a <- pmin(b, pzib(as.vector(y[,object$family == "ZIB"]) - 1, Ntrials = Ntrials, mu = as.vector(mu[,object$family == "ZIB"]), sigma = rep(object$params$phi[object$family == "ZIB"], each = n)))
    
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "ZIB"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "ZNIB")) {
    p_f = sum(object$family == "ZNIB")
    if(length(Ntrials)==1)Ntrials <- rep(Ntrials,p_f)
    if(length(Ntrials)==p)Ntrials <- rep(Ntrials, each = n)
    if(is.matrix(Ntrials))Ntrials <- c(Ntrials)
    
    phis0 = (object$params$phi/(1+object$params$phi + object$params$ZINB.phi))[object$family == "ZNIB"];
    phisN = (object$params$ZINB.phi/(1+object$params$phi + object$params$ZINB.phi))[object$family == "ZNIB"];
    
    b <- pznib(as.vector(y[,object$family == "ZNIB"]), Ntrials = Ntrials, mu = as.vector(mu[,object$family == "ZNIB"]), p0 = rep(phis0, each = n), pN = rep(phisN, each = n))
    a <- pmin(b, pznib(as.vector(y[,object$family == "ZNIB"]) - 1, Ntrials = Ntrials, mu = as.vector(mu[,object$family == "ZNIB"]), p0 = rep(phis0, each = n), pN = rep(phisN, each = n)))
    
    u = a+(b-a)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "ZNIB"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "tweedie")) {
    p_f = sum(object$family == "tweedie")
    phis <- object$params$phi[object$family == "tweedie"] + 1e-05
    b <- fishMod::pTweedie(as.vector(y[,object$family == "tweedie"]), mu = as.vector(mu[,object$family == "tweedie"]), phi = rep(phis, each = n), p = object$Power)
    a <- pmin(b, fishMod::pTweedie(as.vector(y[,object$family == "tweedie"]) - 1, mu = as.vector(mu[,object$family == "tweedie"]), phi = rep(phis, each = n), p = object$Power));
    
    anew =  ifelse((as.vector(y) - 1)<0, 0, a)
    
    u = anew+(b-anew)*runif(n*p_f)
    
    if(any(u==1, na.rm = TRUE)&&replace)u[u==1] <- 1-1e-16
    if(any(u==0, na.rm = TRUE)&&replace)u[u==0] <- 1e-16
    ds.res[,object$family == "tweedie"] <- matrix(qnorm(u),n,p_f)
  }
  if (any(object$family == "ordinal")) {
    
    p_f = sum(object$family == "ordinal")
    o_ind <- c(1:p)[object$family == "ordinal"]
    linkfun <- switch(object$link, "probit" = pnorm, "logit" = plogis)
    for (i in 1:n) {
      for (j in 1:p_f) {
        if(object$zeta.struc == "species"){
          probK <- NULL
          probK[1] <- linkfun(object$params$zeta[o_ind[j], 1] - eta.mat[i, o_ind[j]], log.p = FALSE)
          probK[max(y[, o_ind[j]]) + 1 - min(y[, o_ind[j]])] <- 1 - linkfun(object$params$zeta[o_ind[j], max(y[, o_ind[j]]) - min(y[, o_ind[j]])] - eta.mat[i, o_ind[j]])
          if(length(unique(y[,o_ind[j]]))>2) {
            j.levels <- 2:(max(y[, o_ind[j]]) - min(y[, o_ind[j]]))#
            for (k in j.levels) {
              probK[k] <- linkfun(object$params$zeta[o_ind[j], k] - eta.mat[i, o_ind[j]]) - linkfun(object$params$zeta[o_ind[j], k - 1] - eta.mat[i, o_ind[j]])
            }
          }
          probK <- c(0, probK)
          cumsum.b <- sum(probK[1:(y[i,o_ind[j]]+ifelse(min(y[,o_ind[j]])==0,1,0) + 1)])
          cumsum.a <- min(cumsum.b, sum(probK[1:(y[i,o_ind[j]]+ifelse(min(y[,o_ind[j]])==0,1,0))]))
          u <- runif(n = 1, min = cumsum.a, max = cumsum.b)
          if (abs(u - 1) < 1e-05)
            u <- 1
          if (abs(u - 0) < 1e-05)
            u <- 0
          ds.res[i, o_ind[j]] <- qnorm(u)
        } else {
          probK <- NULL
          probK[1] <- linkfun(object$params$zeta[1+kz] - eta.mat[i, o_ind[j]], log.p = FALSE)
          probK[max(y) + 1 - min(y)] <- 1 - linkfun(object$params$zeta[max(y) - min(y)+kz] - eta.mat[i, o_ind[j]])
          levels <- 2:(max(y) - min(y))#
          for (k in levels) {
            probK[k] <- linkfun(object$params$zeta[k+kz] - eta.mat[i, o_ind[j]]) - linkfun(object$params$zeta[k - 1 +kz] - eta.mat[i, o_ind[j]])
          }
          probK <- c(0, probK)
          cumsum.b <- sum(probK[1:(y[i,o_ind[j]]+ifelse(min(y)==0,1,0) + 1)])
          cumsum.a <- min(cumsum.b, sum(probK[1:(y[i,o_ind[j]]+ifelse(min(y)==0,1,0))]))
          u <- runif(n = 1, min = cumsum.a, max = cumsum.b)
          if (abs(u - 1) < 1e-05)
            u <- 1
          if (abs(u - 0) < 1e-05)
            u <- 0
          ds.res[i, o_ind[j]] <- qnorm(u)
        }
      }
    }
  }
  
  rownames(ds.res) <- rownames(y)
  colnames(ds.res) <- colnames(y)
  
  return(list(residuals = ds.res, linpred = eta.mat))
}