#' @title Log-likelihood of gllvm
#' @description  Extracts Log-likelihood from 'gllvm' objects.
#'
#' @param object   an object of class 'gllvm'.
#' @param ...	not used.
#'
#' @author David I. Warton, Jenni Niku
#'
#' @examples
#' \dontrun{
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = "negative.binomial")
#'# log-Likelihood:
#'logLik(fit)
#'}
#'@export

logLik.gllvm <- function(object, ...)
{
  logL = object$logL
  if(!is.null(object$params$inv.phi)){ object$params$inv.phi<-NULL; }
  attributes(logL)$df <- length(unlist(object$params))-object$num.lv*(object$num.lv-1)/2
  attributes(logL)$nobs <- dim(object$y)[1]
  class(logL) <- "logLik"
  return( logL )
}
