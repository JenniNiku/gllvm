#' @title Extract latent variables
#' @description  Extract latent variables from gllvm object.
#' 
#' @param object an object of class 'gllvm'.
#' 
#'@aliases getLV getLV.gllvm
#'@method getLV gllvm
#'@export
#'@export getLV.gllvm

getLV.gllvm <- function(object)
{
  lvs <- t(t(object$lvs)*object$params$sigma.lv)
  
  if(object$num.lv.c>0){
    lvs<-cbind(object$lv.X%*%object$params$LvXcoef+lvs[,1:object$num.lv.c],lvs[,-c(1:object$num.lv.c)])
    if(object$num.lv.c>0&object$num.lv>0)colnames(lvs)<-c(paste("CLV",1:object$num.lv.c,sep=""),paste("LV",1:object$num.lv,sep=""))
    if(object$num.lv.c>0&object$num.lv==0)colnames(lvs)<-paste("CLV",1:object$num.lv.c,sep="")
  }
  
  return(lvs)
}

#'@export getLV
getLV <- function(object)
{
  UseMethod(generic = "getLV")
}