#' @title Extract coefficients from gllvm
#' @description  Extracts model coefficients from 'gllvm' objects.
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
#'fit <- gllvm(y, family = "negative.binomial")
#'# Coefficients
#'coef(fit)
#'}
#'@export


coef.gllvm <- function(object, ...)
{
  names(object$params)[names(object$params)=="beta0"]="Intercept"
  if(object$row.eff %in% c(TRUE,"fixed")) names(object$params)[names(object$params)=="row.params"]="Row.Intercept"
  if(object$row.eff=="random") names(object$params)[names(object$params)=="row.params"]="Random.Row.Intercept"
  return(object$params)
}
