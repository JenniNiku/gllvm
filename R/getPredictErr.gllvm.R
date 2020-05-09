#' @title Extract prediction errors for latent variables from gllvm object
#' @description  Calculates the prediction errors for latent variables for gllvm model.
#'
#' @param object   an object of class 'gllvm'.
#' @param CMSEP logical, if \code{TRUE} conditional mean squared errors for predictions are calculated. If \code{FALSE}, prediction errors are based on covariances of the variational distributions for \code{method ="VA"}.
#' @param ...	 not used
#'
#' @details 
#' Calculates conditional mean squared errors for predictions.
#' If variational approximation is used, prediction errors can be based on covariances 
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
getPredictErr.gllvm = function(object, CMSEP = TRUE, ...)
{
  num.lv <- object$num.lv
  out <- list()
  if(object$method == "LA"){
    if(object$num.lv>0) out$lvs <- sqrt(apply(object$prediction.errors$lvs,1,diag))
    if(object$row.eff == "random") out$row.effects <- sqrt(abs(object$prediction.errors$row.params))
  }
  
  if(object$method == "VA"){
    if(CMSEP) {
      sdb<-sdA(object)
      object$A<-sdb+object$A
    }
      r=0
      if(object$row.eff=="random"){ 
        r=1
        if(length(dim(object$A))==2){
          out$row.effects <- sqrt(object$A[,1])
        } else {
          out$row.effects <- sqrt(object$A[,1,1])
        }
        }
      if(length(dim(object$A))==2){
        out$lvs <- sqrt(object$A[,1:num.lv+r])
      } else {
        if(num.lv ==1) {
          out$lvs <- sqrt(as.matrix(object$A[,1:num.lv+r,1:num.lv+r]))
        } else {
          out$lvs <- sqrt(apply((object$A[,1:num.lv+r,1:num.lv+r]),1,diag))
        }
      }
  }
  if(object$num.lv > 1) out$lvs <- t(out$lvs)
  return(out)
}

#'@export getPredictErr
getPredictErr <- function(object, ...)
{
  UseMethod(generic = "getPredictErr")
}
