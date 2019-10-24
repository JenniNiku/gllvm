#' @title Extract residual covariance matrix from gllvm object
#' @description  Calculates the residual covariance matrix for gllvm model.
#'
#' @param object  an object of class 'gllvm'.
#' @param adjust  The type of adjustment used for  negative binomial and binomial distribution when computing residual correlation matrix. Options are 0 (no adjustment), 1 (the default adjustment) and 2 (alternative adjustment for NB distribution), see details.
#'
#' @return Function returns following components:
#'  \item{cov }{residual covariance matrix}
#'  \item{trace }{trace of the residual covariance matrix}
#'
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
#' which leads to the residual covariance matrix \eqn{\Theta \Theta' + diag(\Phi)}, where \eqn{\Psi} is the diagonal matrix with \eqn{log(\phi_j + 1)} as diagonal elements (\code{adjust = 1}).
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
#'
#' @author Francis K.C. Hui, Jenni Niku, David I. Warton
#'
#' @examples
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# residual covariance:
#'rescov <- getResidualCov(fit)
#'rescov$cov
#'# Trace of the covariance matrix
#'rescov$tr
#'
#'@aliases getResidualCov getResidualCov.gllvm
#'@method getResidualCov gllvm
#'@export
#'@export getResidualCov.gllvm
getResidualCov.gllvm = function(object, adjust = 1)
{
  ResCov <- object$params$theta %*% t(object$params$theta)
  ResCov.q <- sapply(1:object$num.lv, function(q) object$params$theta[, q] %*% t(object$params$theta[, q]), simplify = F)
  if(adjust > 0 && object$family %in% c("negative.binomial", "binomial")){
  if(object$family == "negative.binomial"){ 
    if(adjust == 1) {
      ResCov <- ResCov + diag(log(object$params$phi + 1))
      ResCov.q <- sapply(1:object$num.lv, function(q) ResCov.q[[q]] + diag(1)/object$num.lv, simplify = F)
      }
    if(adjust == 2){
     ResCov <- ResCov + diag(trigamma(1/object$params$phi))
     ResCov.q <- sapply(1:object$num.lv, function(q) ResCov.q[[q]] + diag(trigamma(1/object$params$phi))/object$num.lv, simplify = F)
     }
    
  }
    if(object$family == "binomial"){ 
      if(object$link == "probit"){
        ResCov <- ResCov + diag(ncol(object$y))
        ResCov.q <- sapply(1:object$num.lv, function(q) ResCov.q[[q]] + diag(ncol(object$y))/object$num.lv, simplify = F)
      } 
      if(object$link == "logit"){
        ResCov <- ResCov + diag(ncol(object$y))*pi^2/3
        ResCov.q <- sapply(1:object$num.lv, function(q) ResCov.q[[q]] + (diag(ncol(object$y))*pi^2/3)/object$num.lv, simplify = F)
      } 
    }
  }
  colnames(ResCov) <- colnames(object$y)
  rownames(ResCov) <- colnames(object$y)
  out <- list(cov = ResCov, trace = sum(diag(ResCov)), trace.q = ResCov.q)
  return(out)
}

#'@export getResidualCov
getResidualCov <- function(object, adjust)
{
  UseMethod(generic = "getResidualCov")
}
