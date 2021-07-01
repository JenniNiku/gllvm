#' @title Functions to extract ecological quantities of the latent variables from a GLLVM, if species are a quadratic function of the latent variables.
#' @description Extracts species optima and tolerances, potentially with standard errors (derived with the Delta method).
#'
#'
#' @param object   an object of class 'gllvm'.
#' @param sd.errors logical. If \code{TRUE}, also returns standard errors.
#' @param ... Not used.
#'
#'@name ecoCoefs
#' 
#' @details 
#' Currently no separate method for calculating species maxima or gradient length are implemented.
#' Gradient length can be inferred from the standard deviation of the latent variables, which is reported by \code{\link{summary.gllvm}}.
#' 
#' @author Bert van der Veen
#'
#' @aliases optima optima.gllvm tolerances tolerances.gllvm
#'
#' @export
#' @export optima.gllvm 
#' @export tolerances.gllvm 
#'
optima.gllvm <- function(object,sd.errors = TRUE, ...) {
  if(!inherits(object,"gllvm.quadratic")){
    stop("Optima can only be extracted for a GLLVM where species are a quadratic function of the latent variables.")
  }
  quadratic <- object$quadratic
    num.lv <- object$num.lv
    num.lv.c <- object$num.lv.c+object$num.RR
    p <- ncol(object$y)
    opt<-object$params$theta[,1:(num.lv+num.lv.c)]/(2*abs(object$params$theta[,-c(1:(num.lv+num.lv.c))]))
    if(sd.errors==TRUE){
      if(is.null(object$sd)|all(unlist(object$sd)==FALSE)){
        cat("Standard errors not present in model, calculating...\n")
        object$Hess<-se.gllvm(object)$Hess 
      }
      V <- object$Hess$cov.mat.mod
      idx <- names(object$TMBfn$par[object$Hess$incl])%in%c("lambda","lambda2")
      V <- V[idx,idx,drop=F]
      colnames(V) <- row.names(V) <- names(object$TMBfn$par[object$Hess$incl])[idx]
      
      
      if(num.lv>0&num.lv.c==0)idx<-which(c(upper.tri(object$params$theta[,1:num.lv],diag=T)))[-1]
      if(num.lv.c>0)idx<-which(c(upper.tri(object$params$theta[,1:num.lv.c],diag=T)))[-1]
      
      #add first row and column of zeros
      V<-rbind(0,cbind(0,V))

      #add zeros where necessary
      for(q in 1:length(idx)){
        V <- rbind(V[1:(idx[q]-1),],0,V[idx[q]:ncol(V),])
        V <- cbind(V[,1:(idx[q]-1)],0,V[,idx[q]:ncol(V)])
      }
      if(num.lv>0&num.lv.c>0){
        idx<-which(c(upper.tri(object$params$theta[,(num.lv.c+1):(num.lv.c+num.lv)],diag=T)))[-1]
      #add a zero infront of second set of  parameters
      V <- cbind(V[,1:(p*num.lv.c)],0,V[,-c(1:(p*num.lv.c))])
      V <- rbind(V[1:(p*num.lv.c),],0,V[-c(1:(p*num.lv.c)),])
      
        for(q in 1:length(idx)){
          #we now have p*num.lv.c elements before the zeros need to be added
          V <- rbind(V[1:(idx[q]-1+p*num.lv.c),],0,V[(idx[q]+p*num.lv.c):ncol(V),])
          V <- cbind(V[,1:(idx[q]-1+p*num.lv.c)],0,V[,(idx[q]+p*num.lv.c):ncol(V)])
        }
      }
      colnames(V)[colnames(V)==""]<-"lambda"
      
      #re-arrange quadratic coefficients to similar order as linear: per species per latent variable (now per latent variable and then species)
      #this makes thing easier further on
      if(quadratic==TRUE){
        idx <- c(1:(p*(num.lv+num.lv.c)),(p*(num.lv+num.lv.c)+c(matrix(1:(p*(num.lv+num.lv.c)),ncol=(num.lv+num.lv.c),nrow=p,byrow=T))))  
      }else{
        idx <- c(1:(p*(num.lv+num.lv.c)),(p*(num.lv+num.lv.c))+rep(1:(num.lv+num.lv.c),each=p))
      }
      
      V<-V[idx,idx]
      du1 <- NULL #linear coefficients
      du2 <- NULL #quadratic coefficients
      #du1 needs to be species then LV, to accommmodate the order in V
      for(i in 1:(num.lv+num.lv.c)){
        for (j in 1:p) {
          du1 <- c(du1,(2 * object$params$theta[j, (num.lv+num.lv.c) + i, drop = F])^-1)
          du2 <- c(du2,object$params$theta[j, i, drop = F] / (2 * object$params$theta[j, (num.lv+num.lv.c) + i, drop = F]^2))
        }
      }
      
      #sum over parameters and their covariances
      #linear coefficients covariances
      cov.mat.lincoef<-(c(du1,du2)%*%t(c(du1,du2))*V)[1:(p*(num.lv+num.lv.c)),1:(p*(num.lv+num.lv.c))]
      #quadratic coefficients covariances
      cov.mat.quadcoef <- (c(du1,du2)%*%t(c(du1,du2))*V)[-c(1:(p*(num.lv+num.lv.c))),-c(1:(p*(num.lv+num.lv.c)))]
      #linear and quadratic coefficient covariances
      diag.cov.mat.coef <- (c(du1,du2)%*%t(c(du1,du2))*V)[1:(p*(num.lv+num.lv.c)),-c(1:(p*(num.lv+num.lv.c)))]
      
      #covariance of optima per lv per species
      #1/2 scalar here because, optima is a function of two parameters
      cov.mat.optima <- cov.mat.lincoef+cov.mat.quadcoef+2*diag.cov.mat.coef
      opt.sd <- matrix(sqrt(1/2*abs(diag(cov.mat.optima))),ncol=(num.lv+num.lv.c),nrow=p)
    }
   
      if((num.lv+num.lv.c)>1){
        if(num.lv.c==0)colnames(opt) <- paste("LV",1:num.lv,sep="") 
        if(num.lv==0)colnames(opt) <- paste("CLV",1:num.lv.c,sep="")
        if(num.lv>0&num.lv.c>0)colnames(opt) <- c(paste("CLV",1:num.lv.c,sep=""),paste("LV",1:num.lv,sep=""))
      }
        
      if((num.lv+num.lv.c)>1&sd.errors==TRUE){
        if(num.lv.c==0)colnames(opt.sd) <- paste("LV",1:num.lv,sep="")
        if(num.lv==0)colnames(opt.sd) <- paste("CLV",1:num.lv.c,sep="")
        if(num.lv>0&num.lv.c>0)colnames(opt.sd) <- c(paste("CLV",1:num.lv.c,sep=""),paste("LV",1:num.lv,sep=""))
      }
    
      if((num.lv+num.lv.c)>1){
        row.names(opt) <- colnames(object$y)
        if(sd.errors==TRUE)row.names(opt.sd) <- colnames(object$y)
      }else{
        names(opt) <- colnames(object$y)
        if(sd.errors==TRUE)names(opt.sd) <- colnames(object$y)
      }
      if(sd.errors==TRUE){
        return(list(optima=opt,sd=opt.sd)) 
      }else{
        return(optima=opt) 
      }
 
}

#'@rdname ecoCoefs
#'@export
tolerances.gllvm <- function(object,sd.errors = TRUE, ...) {
  if(!inherits(object,"gllvm.quadratic")){
    stop("Tolerances can only be extracted for a GLLVM where species are a quadratic function of the latent variables.")
  }
  num.lv <- object$num.lv
  num.lv.c <- object$num.lv.c+object$num.RR
  p <- ncol(object$y)
  quadratic <- object$quadratic
  tol<-1/sqrt(-2*object$params$theta[,-c(1:(num.lv+num.lv.c))])
  if(sd.errors==TRUE){
    if(is.null(object$sd)|all(unlist(object$sd)==FALSE)){
      cat("Standard errors not present in model, calculating...\n")
      object$Hess<-se.gllvm(object)$Hess 
    }
    V <- object$Hess$cov.mat.mod
    idx <- names(object$TMBfn$par[object$Hess$incl])%in%c("lambda2")
    V <- V[idx,idx,drop=F]
    colnames(V) <- row.names(V) <- names(object$TMBfn$par[object$Hess$incl])[idx]
    
    #re-arrange V for easier access
    if(quadratic==TRUE){
      idx <- c(matrix(1:(p*(num.lv+num.lv.c)),ncol=(num.lv+num.lv.c),nrow=p,byrow=T))
    }else{
      idx <- rep(1:(num.lv+num.lv.c),each=p)
    }
    V<-V[idx,idx]

    tol.sd<-matrix(0,nrow=p,ncol=(num.lv+num.lv.c))
    dt<-matrix(0,nrow=p,ncol=(num.lv+num.lv.c))
    for (j in 1:p) {
      for (i in 1:(num.lv+num.lv.c)) {
        dt[j,i] <- -0.5*2^-0.5*abs(object$params$theta[,-c(1:(num.lv+num.lv.c)),drop=F][j,i])^-1.5
      }
    }
    tol.sd <- sqrt(abs(dt^2*matrix(diag(V),ncol=(num.lv+num.lv.c),nrow=p)))
  }
 
  
  if((num.lv+num.lv.c)>1){
      if(num.lv.c==0)colnames(tol) <- paste("LV",1:num.lv,sep="")
      if(num.lv==0)colnames(tol) <- paste("CLV",1:num.lv.c,sep="")
      if(num.lv>0&num.lv.c>0)colnames(tol) <- c(paste("CLV",1:num.lv.c,sep=""),paste("LV",1:num.lv,sep=""))
  }
  if((num.lv+num.lv.c)>1&sd.errors==TRUE){
    if(num.lv.c==0)colnames(tol.sd) <- paste("LV",1:num.lv,sep="")
    if(num.lv==0)colnames(tol.sd) <- paste("CLV",1:num.lv.c,sep="")
    if(num.lv>0&num.lv.c>0)colnames(tol.sd) <- c(paste("CLV",1:num.lv.c,sep=""),paste("LV",1:num.lv,sep=""))
  }
  if((num.lv+num.lv.c)>1){
    row.names(tol) <- colnames(object$y)
  if(sd.errors==TRUE)row.names(tol.sd) <- colnames(object$y)
  }else{
    names(tol) <-colnames(object$y)
    if(sd.errors==TRUE)names(tol.sd) <- colnames(object$y)
  }
  if(sd.errors==TRUE){
    return(list(tolerances=tol,sd=tol.sd))   
  }else{
    return(tolerances=tol) 
  }
  
}

#'@export optima
optima <- function(object, ...)
{
  UseMethod(generic = "optima")
}
#'@method tolerances gllvm
#'@export tolerances
tolerances <- function(object, ...)
{
  UseMethod(generic = "tolerances")
}

