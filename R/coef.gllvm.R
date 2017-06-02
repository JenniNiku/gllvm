#' @title Extract gllvm coefficients
#' @description  Extracts model coefficients from 'gllvm' objects.
#'
#' @param object   An object of class 'gllvm'
#' @param ...	Not used.
#'
#' @author David Warton, Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#' \dontrun{
#' library(mvabund) ## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit GLLVM model
#'fit <- gllvm(y = y, family = "negative.binomial")
#'# Coefficients
#'coef(fit)
#'}
#'@export


coef.gllvm <- function(object, ...)
{
  names(object$params)[names(object$params)=="beta0"]="Intercept"
  if(object$row.eff) names(object$params)[names(object$params)=="row.params"]="Row.Intercept"
  return(object$params)
}
