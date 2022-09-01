#' @title Corrected Akaike information criterion and number of observations
#' @description Calculates corrected Akaike information criterion for small sample sizes, and extracts number of observations.
#' 
#'@name AICc
#'
#' @param object an object of class 'gllvm'.
#' @param ... Not used.
#' 
#' @author Jenni Niku, Bert van der Veen
#' 
#'@aliases AICc AICc.gllvm nobs nobs.gllvm
#'
#'@export
#'@export AICc.gllvm
#'@export nobs.gllvm
#'
AICc.gllvm <- function(object, ...){
  objectlist <- list(object, ...)
  IC<-lapply(objectlist,function(x){
    abund=x$y
    n <- dim(abund)[1]
    p <- dim(abund)[2]
    k<-attributes(logLik.gllvm(x))$df
    AICc <- -2*x$logL + (k) * 2 + 2*k*(k+1)/(n*p-k-1)
    return(AICc)
  })
  return(unlist(IC))
}

#'@method AICc gllvm
#'@export AICc
AICc <- function(object, ...)
{
  UseMethod(generic = "AICc")
}

#'@rdname AICc
#'@export
nobs.gllvm <- function(object, ...){
  n <- prod(dim(object$y))
  return(n)
}