#' @title Summarizing gllvm model fits
#' @description A summary of the fitted 'gllvm' object, including function call, distribution family and model parameters.
#'
#' @param object   an object of class 'gllvm'
#' @param ...	not used.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#' \dontrun{
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit GLLVM model
#'fit <- gllvm(y = y, family = "negative.binomial")
#'summary(fit)
#'}
#'@export


# object: An object of class 'gllvm'
summary.gllvm <- function(object, ...) {
  n <- NROW(object$y)
  p <- NCOL(object$y)
  nX <- dim(object$X)[2]
  nTR <- dim(object$TR)[2]
  num.lv <- object$num.lv
  family <- object$family
  
  M <- cbind(object$params$beta0, object$params$theta)
  
  newnams = c("Intercept")
  
  if (num.lv > 0)
    newnams <- c(newnams, paste("theta.LV", 1:num.lv, sep = ""))
  colnames(M) <- newnams
  rownames(M) <- colnames(object$y)
  cat("Call: \n")
  print(object$call)
  cat("family: \n")
  print(object$family)
  cat("\n")
  cat("Coefficients: \n")
  
  print(M)
  cat("\n")
  if (!is.null(object$TR)) {
    if (!is.null(object$X)) {
      cat("Covariate coefficients: \n")
      print(object$params$B)}
    cat("\n")
  } else {
    if (!is.null(object$X)) {
      cat("Covariate coefficients: \n")
      print(object$params$Xcoef)
      cat("\n")
    }
  }
  if (!is.null(object$params$row.params)) {
    cat("Row intercepts: \n")
    print(object$params$row.params)
    cat("\n")
  }
  
  if (object$row.eff=="random") {
    cat("Variance of random row intercepts: \n")
    object$params$sigma2=object$params$sigma^2;names(object$params$sigma2)="sigma^2"
    print(object$params$sigma2)
    cat("\n")
  }
  
  if(object$family=="negative.binomial"){
    cat("Dispersion parameters inv.phi: \n")
    print(object$params$inv.phi)
  }
  if(object$family=="tweedie"){
    cat("Dispersion parameters: \n")
    print(object$params$phi)
  }
  if(object$family=="ZIP"){
    cat("Zero inflation p: \n")
    print(object$params$phi)
  }
}
