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
  if(object$num.RR==0){
    lvs <- object$lvs
  }else{
    if(object$num.lv.c>0){
      lvs<- cbind(t(t(object$lvs[,1:object$num.lv.c])*object$params$sigma.lv[1:num.lv.c]),matrix(0,ncol=object$num.RR,nrow=n),t(t(object$lvs[,-c(1:object$num.lv.c)])*object$params$sigma.lv[1:object$num.lv]))
    }else{
      lvs<- cbind(matrix(0,ncol=object$num.RR,nrow=n),t(t(object$lvs)*object$params$sigma.lv))
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
getLV <- function(object)
{
  UseMethod(generic = "getLV")
}