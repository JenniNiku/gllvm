#' @title Extract residual covariance matrix from gllvm object
#' @description  Calculates the residual covariance matrix for gllvm model.
#'
#' @param object  an object of class 'gllvm'.
#' @param adjust  The type of adjustment used for  negative binomial, binomial and normal distribution when computing residual correlation matrix. Options are 0 (no adjustment), 1 (the default adjustment) and 2 (alternative adjustment for NB distribution), see details.
#' @param site.index A site index used used in the calculation of a GLLVM with quadratic response model, for which the residual covariance is calculated.
#'
#' @return Function returns following components:
#'  \item{cov }{residual covariance matrix}
#'  \item{trace }{trace of the residual covariance matrix, the total variance explained}
#'  \item{var.q }{trace of the residual covariance matrix per latent variable, variance explained per latent variable}
#'  \item{var.q2 }{trace of the squared term of the residual covariance matrix per latent variable, for quadratic responses. Variance explained per latent variable by the quadratic term}

#' @details 
#' Residual covariance matrix, storing information on species co-occurrence that is not explained by the environmental variables (if included), is calculated using the matrix of latent variables loadings, that is, \eqn{\Theta\Theta'}.
#' 
#' When the responses are modelled using the negative binomial distribution, the residual variances for each species must be adjusted for overdispersion. The two possible adjustement terms are \eqn{log(\phi_j + 1)} (\code{adjust = 1}) and \eqn{\psi^{(1)}(1/\phi_j)} (\code{adjust = 2}), where \eqn{\psi^{(1)}} is the trigamma function.
#' 
#' The negative binomial model can be written using different parametrizations. 
#' The residual covariance with \code{adjust = 1} can be obtained using the lognormal-Poisson parametrization, that is,
#' \deqn{Y_{ij} \sim Poisson(\mu_{ij} \lambda_j),}
#' where \eqn{\lambda_j \sim lognormal(-\sigma^2/2, \sigma^2)} and \eqn{\sigma^2 = log(\phi_j + 1)} and \eqn{log(\mu_{ij}) = \eta_{ij}}. Now \eqn{E[Y_{ij}] = \mu_{ij}} and variance \eqn{V(\mu_{ij}) = \mu_{ij} + \mu_{ij}^2 (exp(\sigma^2) - 1) = \mu_{ij} + \mu_{ij}^2 \phi_j}, which are the same as for the NB distribution.
#' Therefore, on linear predictor scale, we have the variance 
#' \deqn{V(log(\mu_{ij} \lambda_j)) = V(log\mu_{ij}) + V(log\lambda_j) = V(u_i'\theta_j) + \sigma^2 = \theta_j'\theta_j + log(\phi_j + 1).}
#' which leads to the residual covariance matrix \eqn{\Theta \Theta' + \eqn{\Psi}}, where \eqn{\Psi} is the diagonal matrix with \eqn{log(\phi_j + 1)} as diagonal elements (\code{adjust = 1}).
#' 
#' Or, for a GLLVM where species are a quadratic function of the latent variables, we instead have
#' \deqn{V(log(\mu_{ij} \lambda_j)) = V(log\mu_{ij}) + V(log\lambda_j) = V(u_i'\theta_j-u_i' D_j u_i) + \sigma^2 = \theta_j'\theta_j + 2diag(D_j)'diag(D_j)log(\phi_j + 1).}
#' which leads to the residual covariance matrix \eqn{\Theta \Theta' + 2 \Gamma_j \Gamma_j' + diag(\Phi)}, where \eqn{\Gamma_j} holds the quadratic coefficients.
#' Since the quadratic coefficients are constrained to be positive, the residual covariance in the latter case is, given the same coefficients on the linear term, equal or more positive than in the linear case.
#' 
#' The residual covariance matrix with \code{adjust = 2} can be obtained by using Poisson-Gamma parametrization
#' \deqn{Y_{ij} \sim Poisson(\mu_{ij} \lambda_j),}
#' where \eqn{\lambda_j \sim Gamma(1/\phi_j, 1/\phi_j)} and \eqn{\mu_{ij}} is as above. The mean and the variance are of similar form as above and we have that
#' \deqn{V(log(\mu_{ij} \lambda_j)) = V(log\mu_{ij}) + V(log\lambda_j) = \theta_j'\theta_j + \psi^{(1)}(1/\phi_j),}
#' where \eqn{\psi^{(1)}} is the trigamma function.
#' 
#' In the case of binomial distribution, the adjustment terms (\code{adjust = 1}) are 1 for probit link and \eqn{\pi^2/3} for logit link.
#' These are obtained by treating binomial model as latent variable model. Assume
#' \deqn{Y^*_{ij} = \eta_{ij} + e_{ij},}
#' where \eqn{e_{ij} \sim N(0, 1)} for probit model, and \eqn{e_{ij} ~ logistic(0, 1)} for logit model.
#' Then binary response is defined as \eqn{Y_{ij} = 1}, if \eqn{Y^*_{ij} > 0} and 0 otherwise.
#' Now we have that \eqn{\mu_{ij} = P(Y_{ij} = 1) = P(Y^*_{ij} > 0) = P(\eta_{ij} > -e_{ij}) = P(e_{ij} <= \eta_{ij})} which leads to probit and logit models.
#' On linear predictor scale we then have that
#' \deqn{V(\eta_{ij} + e_{ij}) = V(\eta_{ij}) + V(e_{ij}).}
#' For the probit model, the residual covariance matrix is then \eqn{\Theta\Theta' + I_m}, and for the logit model \eqn{\Theta\Theta' + \pi^2/3 I_m}.
#' Similarly as above, for a GLLVM where species are a quadratic function of the latent variables, the term \eqn{2\Gamma_j\Gamma_j'} is added to the residual covariance matrix.
#' 
#' For normal distribution, we can write
#' \deqn{Y_{ij} = \eta_{ij} + e_{ij},}
#' where \eqn{e_{ij} \sim N(0, \phi_j^2)} and thus we have that
#' \deqn{V(\eta_{ij} + e_{ij}) = V(\eta_{ij}) + V(e_{ij}).}
#' For the gaussian model, the residual covariance matrix is then \eqn{\Theta\Theta' + diag(\Phi^2)}.
#' 
#' 
#' @author Francis K.C. Hui, Jenni Niku, David I. Warton, Bert van der Veen
#'
#' @examples
#' \dontrun{
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# residual covariance:
#'rescov <- getResidualCov(fit)
#'rescov$cov
#'# Trace of the covariance matrix
#'rescov$trace
#'# Variance explained per latent variable
#'rescov$var.q
#'}
#'@aliases getResidualCov getResidualCov.gllvm
#'@method getResidualCov gllvm
#'@export
#'@export getResidualCov.gllvm
getResidualCov.gllvm = function(object, adjust = 1, site.index = NULL)
{
  if(object$quadratic!=FALSE&&is.null(site.index)&object$num.lv.c>0){
    stop("Please provide a site index vector for which the residual covariances should be calculated. \n")
  }else if(object$quadratic!=FALSE&&!is.null(site.index)){
    if(length(site.index)!=1&object$num.lv.c>0&object$quadratic!=FALSE){
      stop("Site.index should be of length 1. \n")
    }
  }
  
  if((object$num.lv+object$num.lv.c)==0){
    stop("No latent variables present in model.")
  }
  if(any(class(object)=="gllvm.quadratic")){
    ResCov <- object$params$theta[, 1:object$num.lv, drop = F] %*% t(object$params$theta[, 1:object$num.lv, drop = F]) + 2 * object$params$theta[, -c(1:object$num.lv), drop = F] %*% t(object$params$theta[, -c(1:object$num.lv), drop = F])
    if(object$num.lv.c>0)ResCov <- ResCov + Reduce("+",sapply(1:object$num.lv.c,function(q)4*c(object$lv.X[site.index,,drop=F]%*%object$params$LvXcoef[,q,drop=F]*object$lv.X[site.index,,drop=F]%*%object$params$LvXcoef[,q,drop=F])*object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F]%*%t(object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F]),simplify=F))
    ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) object$params$theta[, q] %*% t(object$params$theta[, q]), simplify = F)
    ResCov.q2 <- sapply(1:(object$num.lv+object$num.lv.c), function(q) 2*object$params$theta[, q+(object$num.lv+object$num.lv.c)] %*% t(object$params$theta[, q+(object$num.lv+object$num.lv.c)]), simplify = F)
    if(object$num.lv.c>0)ResCov <- ResCov - Reduce("+",sapply(1:object$num.lv.c,function(q)2*c(object$lv.X[site.index,,drop=F]%*%object$params$LvXcoef[,q,drop=F])*(object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F]%*%t(object$params$theta[,q,drop=F])+object$params$theta[,q,drop=F]%*%t(object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F])),simplify=F))
  }else{
    ResCov <- object$params$theta %*% t(object$params$theta)
    ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) object$params$theta[, q] %*% t(object$params$theta[, q]), simplify = F)
  }
  
  
  if(adjust > 0 && object$family %in% c("negative.binomial", "binomial", "gaussian")){
  if(object$family == "negative.binomial"){ 
    if(adjust == 1) {
      if(any(class(object)=="gllvm.quadratic")){
        ResCov <- ResCov + diag(log(object$params$phi + 1))
        ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q[[q]] + diag(log(object$params$phi + 1))/((object$num.lv+object$num.lv.c)*2), simplify = F) 
        ResCov.q2 <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q2[[q]] + diag(log(object$params$phi + 1))/((object$num.lv+object$num.lv.c)*2), simplify = F) 
      }else{
        ResCov <- ResCov + diag(log(object$params$phi + 1))
        ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q[[q]] + diag(log(object$params$phi + 1))/(object$num.lv+object$num.lv.c), simplify = F)  
      }
      }
    if(adjust == 2){
      if(any(class(object)=="gllvm.quadratic")){
        ResCov <- ResCov + diag(trigamma(1/object$params$phi))
        ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q[[q]] + diag(trigamma(1/object$params$phi))/((object$num.lv+object$num.lv.c)*2), simplify = F)   
        ResCov.q2 <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q2[[q]] + diag(trigamma(1/object$params$phi))/((object$num.lv+object$num.lv.c)*2), simplify = F)   
      }else{
     ResCov <- ResCov + diag(trigamma(1/object$params$phi))
     ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q[[q]] + diag(trigamma(1/object$params$phi))/(object$num.lv+object$num.lv.c), simplify = F)
    }
     }
    
  }
    if(object$family == "binomial"){ 
      if(object$link == "probit"){
        if(any(class(object)=="gllvm.quadratic")){
          ResCov <- ResCov + diag(ncol(object$y))
          ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q[[q]] + diag(ncol(object$y))/((object$num.lv+object$num.lv.c)*2), simplify = F)
          ResCov.q2 <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q2[[q]] + diag(ncol(object$y))/((object$num.lv+object$num.lv.c)*2), simplify = F)
        }else{
          ResCov <- ResCov + diag(ncol(object$y))
          ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q[[q]] + diag(ncol(object$y))/(object$num.lv+object$num.lv.c), simplify = F)  
        }  
        
      } 
      if(object$link == "logit"){
        ResCov <- ResCov + diag(ncol(object$y))*pi^2/3
        ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q[[q]] + (diag(ncol(object$y))*pi^2/3)/(object$num.lv+object$num.lv.c), simplify = F)
      } 
    }
    if(object$family == "gaussian"){
        if(any(class(object)=="gllvm.quadratic")){
          ResCov <- ResCov + diag((object$params$phi^2))
          ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q[[q]] + diag(ncol(object$y))/((object$num.lv+object$num.lv.c)*2), simplify = F)
          ResCov.q2 <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q2[[q]] + diag(ncol(object$y))/((object$num.lv+object$num.lv.c)*2), simplify = F)
        }else{
          ResCov <- ResCov + diag((object$params$phi^2))
          ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) ResCov.q[[q]] + diag(ncol(object$y))/(object$num.lv+object$num.lv.c), simplify = F)  
        }
      
    }
  }
  ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) sum(diag(ResCov.q[[q]])))
  if(object$num.lv.c==0)names(ResCov.q) <- paste("LV", 1:object$num.lv, sep = "")
  if(object$num.lv==0)names(ResCov.q) <- paste("CLV", 1:object$num.lv.c, sep = "")
  if(object$num.lv.c>0&object$num.lv>0)names(ResCov.q) <-c( paste("CLV", 1:object$num.lv.c, sep = ""), paste("LV", 1:object$num.lv, sep = ""))
  if(any(class(object)=="gllvm.quadratic")){
    ResCov.q2 <- sapply(1:(object$num.lv+object$num.lv.c), function(q) sum(diag(ResCov.q2[[q]])))
    if(object$num.lv.c==0)names(ResCov.q2) <- paste("LV", 1:object$num.lv, "^2",sep = "")
    if(object$num.lv==0)names(ResCov.q2) <- paste("CLV", 1:object$num.lv.c, "^2",sep = "")
    if(object$num.lv>0&object$num.lv.c>0)names(ResCov.q2) <- c(paste("CLV", 1:object$num.lv.c, "^2",sep = ""),paste("LV", 1:object$num.lv, "^2",sep = ""))
  }
  
  colnames(ResCov) <- colnames(object$y)
  rownames(ResCov) <- colnames(object$y)
  if(any(class(object)=="gllvm.quadratic")){
    out <- list(cov = ResCov, trace = sum(diag(ResCov)), var.q = ResCov.q, trace.q2 = ResCov.q2)
  }else{
    out <- list(cov = ResCov, trace = sum(diag(ResCov)), var.q = ResCov.q)  
  }
  return(out)
}

#'@export getResidualCov
getResidualCov <- function(object, adjust)
{
  UseMethod(generic = "getResidualCov")
}
