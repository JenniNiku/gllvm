#' @title Print gllvm
#' @description Print an object of class 'gllvm'.
#'
#' @param x   an object of class 'gllvm'.
#' @param ...	not used.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#' \dontrun{
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
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
  if(!is.null(x$params$inv.phi)){ x$params$inv.phi<-NULL; }
  df=length(unlist(x$params))-object$num.lv*(object$num.lv-1)/2
  if(row.eff %in% c("fixed",TRUE) ) df=df-1
  cat("Degrees of freedom: ",df,"\n")
  cat("AIC: ",AIC(x),"\n")
}
#print(fit)
