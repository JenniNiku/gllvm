#' @title Summary of fitted gllvm object
#' @description A summary of the fitted gllvm object, including function call, distribution family and model parameters.
#'
#' @param object   An object of class 'gllvm'
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
    cat("X Covariate coefficients: \n")
    print(object$params$Xcoef)}
    if (object$get.trait) {
      cat("Trait covariate coefficients: \n")
      print(object$params$Tcoef)
      cat("\n")
    }
    if (object$get.trait) {
      cat("Interaction terms: \n")
      print(object$params$fourth)
      cat("\n")
    }
  } else {
    if (!is.null(object$X)) {
      cat("X Covariate coefficients: \n")
      print(object$params$Xcoef)
      cat("\n")
    }
  }
  if (!is.null(object$params$row.params)) {
    cat("Row Intercepts: \n")
    print(object$params$row.params)
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
    print(object$params$p)
  }
}
