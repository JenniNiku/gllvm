#' @title Print GLLVM
#' @description Print an object of class 'gllvm'
#'
#' @param x   An object of class 'gllvm'
#' @param ...	Not used.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#' \dontrun{
#' library(mvabund) ## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit GLLVM model
#'fit <- gllvm(y = y, family = "negative.binomial")
#'# Print:
#'print(fit)
#'}
#'@export

print.gllvm <- function(x, ...) {
  cat("Call: \n")
  print(x$call)
  cat("family: \n")
  print(x$family)
  cat("method: \n")
  print(x$method)
  cat("\n")
  cat("log-likelihood: ",x$logL,"\n")
  cat("Degrees of freedom: ",length(unlist(x$params)),"\n")
  cat("AIC: ",AIC(x),"\n")
}
#print(fit)
