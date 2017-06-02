#' @title Log-likelihood of gllvm
#' @description  Extracts Log-likelihood from 'gllvm' objects.
#'
#' @param object   An object of class 'gllvm'
#' @param ...	Not used.
#'
#' @author David Warton
#'
#' @examples
#' \dontrun{
#' library(mvabund) ## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit GLLVM model
#'fit <- gllvm(y = y, family = "negative.binomial")
#'# log-Likelihood:
#'logLik(fit)
#'}
#'@export

logLik.gllvm <- function(object, ...)
{
  logL = object$logL
  attributes(logL)$df <- length(unlist(object$params))
  attributes(logL)$nobs <- dim(object$y)[1]
  class(logL) <- "logLik"
  return( logL )
}
