#' @title Extract loadings
#' @description  Extract loadings (species scores) from a gllvm object.
#' 
#' @param object an object of class 'gllvm'.
#' @param ... not used
#' 
#' @details
#' Function retrieves the loadings a.k.a. species scores for a GLLVM. For the optima of a quadratic response model, see \code{\link{optima.gllvm}}
#'@aliases getLoadings getLoadings.gllvm
#'@method getLoadings gllvm
#'@export
#'@export getLoadings.gllvm

getLoadings.gllvm <- function(object, ...)
{
  if(inherits(object, "gllvm.quadratic"))stop("For the quadratic response model, please use 'optima.gllvm' instead.")
  sploads <- object$params$theta
  
  if(object$num.lv>0)sploads[,(ncol(sploads)-object$num.lv+1):ncol(sploads)] <- sploads[,(ncol(sploads)-object$num.lv+1):ncol(sploads)]%*%diag(tail(object$params$sigma.lv, object$num.lv))
  # if(object$num.lv.c>0 & scale)sploads[,1:object$num.lv.c] <- sploads[,1:object$num.lv.c]%*%diag(head(object$params$sigma.lv, object$num.lv.c))
  
  return(sploads)
}

#'@export getLoadings
getLoadings <- function(object, ...)
{
  UseMethod(generic = "getLoadings")
}