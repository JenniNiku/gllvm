#' @title Extract residual covariance matrix from gllvm object
#' @description  Calculates the residual covariance matrix for gllvm model.
#'
#' @param object  an object of class 'gllvm'.
#' @param adjust  defaults to 1, when default adjustment is used for negative binomial and binomial models. Alternatives are 0, when adjustment is not used, and 2 for negative binomial, see details.
#'
#' @return Function returns following components:
#'  \item{cov }{residual covariance matrix}
#'  \item{trace }{trace of the residual covariance matrix}
#'
#' @details 
#' Residual covariance matrix, storing information on species co-occurrence that is ot explained by environmental variables, if included, is calculated based on the matrix of loadings of latent variables, denote \eqn{\Theta\Theta'}. 
#' 
#' Due the overdispersion in case of negative binomial distribution, the residual variance of each species j must be adjusted. Two alternatives for adjustment are \eqn{log(\phi_j + 1)} (\code{adjust = 1}) and trigamma term \eqn{\psi^(1)(1/\phi_j)} (\code{adjust = 2}).
#' 
#' The negative binomial model can be written using different parametrizations. 
#' The residual covariance with \code{adjust = 1} can be obtained with lognormal-Poisson distribution, assume
#' \deqn{Y_{ij} ~ Poisson(\mu_{ij} * \lambda_j),}
#' where \eqn{\lambda_j ~ lognormal(-\sigma^2/2, \sigma^2)}, \eqn{\sigma^2 = log(\phi_j + 1)} and \eqn{log(\mu_{ij}) = \eta_ij} with \eqn{\eta_ij} as usual. Now we have expectation \eqn{E[Y_{ij}] = \mu_{ij}} and variance \eqn{V(\mu_{ij}) = \mu_{ij} + \mu_{ij}^2 * (exp(\sigma^2) - 1) = \mu_{ij} + \mu_{ij}^2 * \phi_j} which are same for NB distribution.
#' Therefore, on linear predictor scale, we get variance 
#' \deqn{V(log(\mu_{ij} * \lambda_j)) = V(log\mu_{ij}) + V(log\lambda_j) = V(u_i'\theta_j) + \sigma^2 = \theta_j'\theta_j + log(\phi_j + 1).}
#' which leads to the residual covariance with argument \code{adjust = 1}, \eqn{\Theta \Theta' + diag(\Phi)}.
#' 
#' The residual covariance with \code{adjust = 2} can be obtained by assuming Poisson-Gamma distribution:
#' \deqn{Y_{ij} ~ Poisson(\mu_{ij} * \lambda_j),}
#' where \eqn{\lambda_j ~ Gamma(1/\phi_j, 1/\phi_j)} and \eqn{\mu_{ij}} as usual. We get same mean and variance for \eqn{Y_{ij}}, and we obtain
#' \deqn{V(log(\mu_{ij} * \lambda_j)) = V(log\mu_{ij}) + V(log\lambda_j) = \theta_j'\theta_j + \psi^(1) * (1/\phi_j).}
#' 
#' In case of binomial distribution, the adjustment (\code{adjust = 1}) is identity matrix for probit link and \eqn{\pi^2/3*I_m} for logit link.
#' These can be obtained by considering binomial model as a latent variable model:
#' Assume random variable
#' \deqn{Y*_{ij} = \eta_{ij} + e_{ij},}
#' where \eqn{e_{ij} ~ N(0, 1)} for probit model, and \eqn{e_{ij} ~ logistic(0, 1)} for logit model.
#' Then binary response is defined as \eqn{Y_{ij} = 1}, if \eqn{Y*_{ij} > 0} and 0 otherwise.
#' Now we have \eqn{\mu_ij = P(Y_{ij} = 1) = P(Y*_{ij} > 0) = P(\eta_{ij} > -e_{ij}) = P(e_{ij} <= \eta_{ij})} which leads to probit and logit models.
#' On linear predictor scale, we get variance
#' \deqn{V(\eta_{ij} + e_{ij}) = V(\eta_{ij}) + V(e_{ij}).}
#' For probit model we get covariance \eqn{\Theta\Theta' + I_m} and for logit \eqn{\Theta\Theta' + \pi^2/3 * I_m}.
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
  if(adjust > 0 && object$family %in% c("negative.binomial", "binomial")){
  if(object$family == "negative.binomial"){ 
    if(adjust == 1) ResCov <- ResCov + diag(log(object$params$phi + 1))
    if(adjust == 2) ResCov <- ResCov + diag(trigamma(1/object$params$phi))
  }
    if(object$family == "binomial"){ 
      if(object$link == "probit") ResCov <- ResCov + diag(ncol(object$y))
      if(object$link == "logit") ResCov <- ResCov + diag(ncol(object$y))*pi^2/3
    }
  }
  colnames(ResCov) <- colnames(object$y)
  rownames(ResCov) <- colnames(object$y)
  out <- list(cov = ResCov, trace = sum(diag(ResCov)))
  return(out)
}

#'@export getResidualCov
getResidualCov <- function(object, adjust)
{
  UseMethod(generic = "getResidualCov")
}
