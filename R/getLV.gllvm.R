#' @title Extract latent variables
#' @description  Extract latent variables from gllvm object.
#' 
#' @param object an object of class 'gllvm'.
#' 
#'@aliases getLV getLV.gllvm
#'@method getLV gllvm
#'
#'@export getLV.gllvm

getLV.gllvm <- function(object)
{
  return(object$lvs)
}

#'@export getLV

getLV <- function(object)
{
  UseMethod(generic = "getLV")
}