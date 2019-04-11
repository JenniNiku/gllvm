#' @title Extract prediction errors for latent variables from gllvm object
#' @description  Calculates the prediction errors for latent variables for gllvm model.
#'
#' @param object   an object of class 'gllvm'.
#'
#' @details 
#' If variational approximation is used, prediction errors are based on covariances 
#' of the variational distributions, and therefore they do not take into account 
#' the uncertainty in the estimation of (fixed) parameters. 
#'
#' @return Function returns following components:
#'  \item{lvs }{prediction errors for latent variables}
#'  \item{row.effects }{prediction errors for random row effects if included}
#'
#' @author Francis K.C. Hui, Jenni Niku, David I. Warton
#'
#' @examples
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# prediction errors for latent variables:
#'getPredictErr(fit)
#'
#'
#'@aliases getPredictErr getPredictErr.gllvm
#'@method getPredictErr gllvm
#'@export
#'@export getPredictErr.gllvm
getPredictErr.gllvm = function(object)
{
  out <- list()
  if(object$method == "LA"){
    if(object$num.lv>0) out$lvs <- sqrt(apply(object$prediction.errors$lvs,1,diag))
    if(object$row.eff == "random") out$row.effects <- sqrt(abs(object$prediction.errors$row.params))
  }
  
  if(object$method == "VA"){
    if(object$num.lv>0) out$lvs <- sqrt(apply(object$A,1,diag))
    if(object$row.eff == "random") out$row.effects <- sqrt(abs(object$Ar))
  }
  if(object$num.lv > 1) out$lvs <- t(out$lvs)
  return(out)
}

#'@export getPredictErr
getPredictErr <- function(object)
{
  UseMethod(generic = "getPredictErr")
}
