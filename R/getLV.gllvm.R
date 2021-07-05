#' @title Extract latent variables
#' @description  Extract latent variables from gllvm object.
#' 
#' @param object an object of class 'gllvm'.
#' @param type type of latent variable scores to retrieve from a gllvm object. For models with unconstrained latent variables, defaults to "unconstrained". For models with constrained latent variables, defaults to "constrained". A third option is "scaled", which returns latent variables multiplied with scale parameters of the loadings.
#' @param ... not used
#' 
#'@aliases getLV getLV.gllvm
#'@method getLV gllvm
#'@export
#'@export getLV.gllvm

getLV.gllvm <- function(object, type = NULL, ...)
{
  if(is.null(type)){
  if(object$num.lv.c==0&object$num.RR==0){
    type = "unconstrained"
  }else{
    type = "constrained"
  }
  }else if(type=="constrained"&object$num.lv.c==0&object$num.RR==0){
    stop("Cannot retrieve constrained latent variables for an unconstrained latent variable model.")
  }else if(type=="scaled"&object$num.lv.c==0&object$num.lv==0){
    stop("Cannot retrieve scaled latent variables for a model without scale parameter.")
  }else if(type=="unconstrained"&object$num.lv.c==0&object$num.lv==0){
    stop("Cannot retrieve unconstrained latent variables, for a model without latent effect.")
  }
  
  n <- nrow(object$y)
  if(type == "unconstrained"){
    lvs <- object$lvs
  }else if(type=="scaled"){
    lvs <- t(t(object$lvs)*object$params$sigma.lv)
  }else if(type=="constrained"){
    if(object$num.RR==0){
      lvs <- t(t(object$lvs)*object$params$sigma.lv)
    }else{
      if(object$num.lv.c>0){
        lvs<- cbind(t(t(object$lvs[,1:object$num.lv.c])*object$params$sigma.lv[1:object$num.lv.c]),matrix(0,ncol=object$num.RR,nrow=n),t(t(object$lvs[,-c(1:object$num.lv.c)])*object$params$sigma.lv[1:object$num.lv]))
      }else if(object$num.lv>0&object$num.lv.c==0){
        lvs<- cbind(matrix(0,ncol=object$num.RR,nrow=n),t(t(object$lvs)*object$params$sigma.lv))
      }else{
        lvs <- matrix(0,ncol=object$num.RR,nrow=n)
      }
    }
  }
  if((object$num.lv.c+object$num.RR)>0){
    lvs<-cbind(object$lv.X%*%object$params$LvXcoef+lvs[,1:(object$num.lv.c+object$num.RR)],lvs[,-c(1:(object$num.lv.c+object$num.RR))])
    if((object$num.lv.c+object$num.RR)>0&object$num.lv>0)colnames(lvs)<-c(paste("CLV",1:(object$num.lv.c+object$num.RR),sep=""),paste("LV",1:object$num.lv,sep=""))
    if((object$num.lv.c+object$num.RR)>0&object$num.lv==0)colnames(lvs)<-paste("CLV",1:(object$num.lv.c+object$num.RR),sep="")
  }
  
  return(lvs)
}

#'@export getLV
getLV <- function(object, ...)
{
  UseMethod(generic = "getLV")
}