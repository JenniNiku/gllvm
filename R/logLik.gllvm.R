#' @title Log-likelihood of gllvm
#' @description  Extracts Log-likelihood from 'gllvm' objects.
#'
#' @param object   an object of class 'gllvm'.
#' @param ...	not used.
#'
#' @author David I. Warton, Jenni Niku
#'
#' @examples
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# log-Likelihood:
#'logLik(fit)
#'
#'@export

logLik.gllvm <- function(object, ...)
{
  logL <- object$logL
  if (!is.null(object$params$inv.phi)) {
    object$params$inv.phi <- NULL

  }
    if(object$family=="ordinal"){
      if(object$zeta.struc=="species")object$params$zeta<-object$params$zeta[,-1]
      if(object$zeta.struc=="common")object$params$zeta<-object$params$zeta[-1]
    }
  if (object$row.eff %in% c("fixed", TRUE))
    object$params$row.params <- object$params$row.params[-1]
  if (object$row.eff == "random")
    object$params$row.params <- NULL
  if (!is.null(object$randomX)){
    object$params$Br <- NULL
    object$params$sigmaB <- object$params$sigmaB[lower.tri(object$params$sigmaB, diag = TRUE)]
  }

  attributes(logL)$df <- length(unlist(object$params)[!is.na(unlist(object$params))]) - object$num.lv * (object$num.lv - 1) / 2
  attributes(logL)$nobs <- dim(object$y)[1]
  class(logL) <- "logLik"
  return(logL)
}
