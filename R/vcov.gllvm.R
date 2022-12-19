#' @title Returns variance-covariance matrix of coefficients in a GLLVM.
#' @description Returns the variance-covariance matrix of the parameters from a GLLVM. If the variance-covariance matrix was not calculated after model fitting, this function will have to calculate the variance-covariance matrix, which may be computational intensive for a large number of species and/or sites.
#'
#'
#' @param object   an object of class 'gllvm'.
#' @param ...	not used.
#' @details
#' Calculates the variance-covariance matrix of a GLLVM object using \code{\link{se.gllvm}}, which may be computational intensive with many parameters.The parameters might have unintuitive names. Fixed-effects coefficients are labeled "b", and are ordered per species as: 1) intercepts 2) fixed-effects slopes. Coefficients of the latent variables are labled "lambda" (linear coefficients) or "lambda2".
#'
#' @author Bert van der Veen
#' 
#'@aliases vcov vcov.gllvm
#'@method vcov gllvm
#'@importFrom stats vcov
#'
#'@export
#'@export vcov.gllvm

vcov.gllvm <- function(object, ...){
  if(is.null(object$sd)){
    cat("Standard errors not present in model, calculating...")
    object$Hess<-se.gllvm(object)$Hess
  }
  V <- object$Hess$cov.mat.mod
  colnames(V) <- row.names(V) <- names(object$TMBfn$par[object$Hess$incl])
  return(V)
}

