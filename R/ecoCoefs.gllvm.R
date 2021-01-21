#' @title Functions to extract ecological quantities of the latent variables from a GLLVM, if species are a quadratic function of the latent variables.
#' @description Extracts species optima, tolerances or gradient lengths, potentially with standard errors (derived with the Delta method).
#'
#'@name optima
#'
#' @param object   an object of class 'gllvm'.
#' @param sd.errors logical. If \code{TRUE}, also returns standard errors.
#'
#' 
#' @author Bert van der Veen
#'
#'@aliases optima optima.gllvm tolerances tolerances.gllvm gradient.length gradient.length.gllvm
#'
#'@export
#'@export optima.gllvm 
#'@export tolerances.gllvm 
#'@export gradient.length.gllvm
#'
optima.gllvm <- function(object,sd.errors = TRUE) {
  if(!inherits(object,"gllvm.quadratic")){
    stop("Optima can only be extracted for a GLLVM where species are a quadratic function of the latent variables.")
  }
  quadratic <- object$quadratic
    num.lv <- object$num.lv
    p <- ncol(object$y)
    opt<-object$params$theta[,1:num.lv]/(2*abs(object$params$theta[,-c(1:num.lv)]))
    if(sd.errors==TRUE){
      if(is.null(object$sd)|all(unlist(object$sd)==FALSE)){
        cat("Standard errors not present in model, calculating...\n")
        object$Hess<-se.gllvm(object)$Hess 
      }
      V <- object$Hess$cov.mat.mod
      idx <- names(object$TMBfn$par[object$Hess$incl])%in%c("lambda","lambda2")
      V <- V[idx,idx,drop=F]
      colnames(V) <- row.names(V) <- names(object$TMBfn$par[object$Hess$incl])[idx]
      
      
      idx<-which(c(upper.tri(object$params$theta[,1:num.lv])))
      #add zeros where necessary
      for(q in 1:length(idx)){
        V <- rbind(V[1:(idx[q]-1),],0,V[idx[q]:ncol(V),])
        V <- cbind(V[,1:(idx[q]-1)],0,V[,idx[q]:ncol(V)])
      }
      colnames(V)[colnames(V)==""]<-"lambda"
      
      #re-arrange quadratic coefficients to similar order as linear: per species per latent variable (now per latent variable and then species)
      #this makes thing easier further on
      if(quadratic==TRUE){
        idx <- c(1:(p*num.lv),(p*num.lv+c(matrix(1:(p*num.lv),ncol=num.lv,nrow=p,byrow=T))))  
      }else{
        idx <- c(1:(p*num.lv),(p*num.lv)+rep(1:num.lv,each=p))
      }
      
      V<-V[idx,idx]
      du1 <- NULL #linear coefficients
      du2 <- NULL #quadratic coefficients
      #du1 needs to be species then LV, to accommmodate the order in V
      for(i in 1:num.lv){
        for (j in 1:p) {
          du1 <- c(du1,(2 * object$params$theta[j, num.lv + i, drop = F])^-1)
          du2 <- c(du2,object$params$theta[j, i, drop = F] / (2 * object$params$theta[j, num.lv + i, drop = F]^2))
        }
      }
      
      #sum over parameters and their covariances
      #linear coefficients covariances
      cov.mat.lincoef<-(c(du1,du2)%*%t(c(du1,du2))*V)[1:(p*num.lv),1:(p*num.lv)]
      #quadratic coefficients covariances
      cov.mat.quadcoef <- (c(du1,du2)%*%t(c(du1,du2))*V)[-c(1:(p*num.lv)),-c(1:(p*num.lv))]
      #linear and quadratic coefficient covariances
      diag.cov.mat.coef <- (c(du1,du2)%*%t(c(du1,du2))*V)[1:(p*num.lv),-c(1:(p*num.lv))]
      
      #covariance of optima per lv per species
      #1/2 scalar here because, optima is a function of two parameters
      cov.mat.optima <- cov.mat.lincoef+cov.mat.quadcoef+2*diag.cov.mat.coef
      opt.sd <- matrix(sqrt(1/2*abs(diag(cov.mat.optima))),ncol=num.lv,nrow=p)
    }
   
      if(num.lv>1)colnames(opt) <- paste("LV",1:num.lv,sep="") 
      if(num.lv>1&sd.errors==TRUE)colnames(opt.sd) <- paste("LV",1:num.lv,sep="")
      if(num.lv>1){
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

#'@rdname optima
#'@export
tolerances.gllvm <- function(object,sd.errors = TRUE) {
  if(!inherits(object,"gllvm.quadratic")){
    stop("Tolerances can only be extracted for a GLLVM where species are a quadratic function of the latent variables.")
  }
  num.lv <- object$num.lv
  p <- ncol(object$y)
  quadratic <- object$quadratic
  tol<-1/sqrt(-2*object$params$theta[,-c(1:num.lv)])
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
      idx <- c(matrix(1:(p*num.lv),ncol=num.lv,nrow=p,byrow=T))
    }else{
      idx <- rep(1:num.lv,each=p)
    }
    V<-V[idx,idx]

    tol.sd<-matrix(0,nrow=p,ncol=num.lv)
    dt<-matrix(0,nrow=p,ncol=num.lv)
    for (j in 1:p) {
      for (i in 1:num.lv) {
        dt[j,i] <- -0.5*2^-0.5*mod$params$theta[,-c(1:num.lv),drop=F][j,i]
      }
    }
    tol.sd <- sqrt(abs(dt^2*matrix(diag(V),ncol=num.lv,nrow=p)))
  }
 
  
  if(num.lv>1)colnames(tol) <- paste("LV",1:num.lv,sep="")
  if(num.lv>1&sd.errors==TRUE)colnames(tol.sd) <- paste("LV",1:num.lv,sep="")
  if(num.lv>1){
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

#'@rdname optima
#'@export
gradient.length.gllvm <- function(object,sd.errors = TRUE) {
  if(!inherits(object,"gllvm.quadratic")){
    stop("Gradient length can only be extracted for a GLLVM where species are a quadratic function of the latent variables.\n")
  }
  num.lv <- object$num.lv
  p <- ncol(object$y)
  quadratic <- object$quadratic
  tol<-1/sqrt(-2*object$params$theta[,-c(1:num.lv),drop=F])
  if(sd.errors==TRUE){
if(is.null(object$sd)|all(unlist(object$sd)==FALSE)){
    cat("Standard errors not present in model, calculating...")
    object$Hess<-se.gllvm(object)$Hess
  }
  V <- object$Hess$cov.mat.mod
  idx <- names(object$TMBfn$par[object$Hess$incl])%in%c("lambda2")
  V <- V[idx,idx,drop=F]
  colnames(V) <- row.names(V) <- names(object$TMBfn$par[object$Hess$incl])[idx]
  
  #re-arrange V for easier access
  if(quadratic==TRUE){
    idx <- c(matrix(1:(p*num.lv),ncol=num.lv,nrow=p,byrow=T))
  }else{
    idx <- rep(1:num.lv)
  }
  V<-V[idx,idx]
  
  grad <- NULL
  idx<-list()
  for(q in 1:num.lv){
    if(p%%2==0&quadratic==TRUE){
      #even
      idx[[q]]<-order(1/sqrt(2*-object$params$theta[,num.lv+q]),decreasing=T)[(p/2):((p/2)+1)]
    }else{
      #odd
      idx[[q]]<-order(1/sqrt(2*-object$params$theta[,num.lv+q]),decreasing=T)[p-round(p/2)]
    }
    
    if(quadratic=="LV"|length(idx[[q]])==1){
      if(quadratic=="LV"){
        #only a function of one parameter, so this is derivative of gradient length (since there is no median)
        grad<-c(grad,4*abs(object$params$theta[1,num.lv+q])^-0.5)
      }else{
        grad<-c(grad,4*abs(2*object$params$theta[idx[[q]],num.lv+q])^-0.5)#Only needs the right idx, no transformations. Idx of 2 parameters
      }
      
    }else{
      #this one is a little more, due to the presence of the median on the tolerance scale.
      grad<-c(grad,sapply(1:2,function(i)8*(abs(object$params$theta[idx[[q]][i],num.lv+q])^1.5*sum(abs(object$params$theta[idx[[q]],num.lv+q])^-0.5)^2)^-1))
    }
    
  }
  gradSD <- NULL
    if(quadratic=="LV"){
      gradSD <- c(gradSD, sqrt(abs(diag((grad%*%t(grad)*V))))) #no scalar, only function of one parameter
    }else{
      if(p%%2==0){
        idx2 <- split(1:(num.lv*2),factor(rep(1:num.lv,each=2)))
      }else{
        idx2 <- 1:num.lv
      }
      
      
      for(q in 1:num.lv){
      gradSD <- c(gradSD, sqrt(abs(1/2*sum((grad[idx2[[q]]]%*%t(grad[idx2[[q]]])*V[(p*(q-1)+1):(p*q),(p*(q-1)+1):(p*q)][idx[[q]],idx[[q]]]))))) #scalar because function of 2 parameters, which we sum over
      }
    }
  grad.length.sd <- sqrt(abs(gradSD))
  }

  if(quadratic==TRUE){
    grad.length <- apply(abs(sapply(1:num.lv,function(q)object$params$theta[idx[[q]],q+num.lv])),2,function(x)16*sum(x^-0.5)^-1)
  }else if(quadratic=="LV"){
    grad.length <- 8*abs(object$params$theta[1,-c(1:num.lv)])^0.5
  }
  names(grad.length) <- paste("LV",1:num.lv,sep="")
  if(sd.errors==TRUE)names(grad.length.sd) <- paste("LV",1:num.lv,sep="")
  
  if(sd.errors==TRUE){
    return(list(gradient.length=grad.length,sd=grad.length.sd))   
  }else{
    return(gradient.length=grad.length)
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


#'@method gradient.length gllvm
#'@export gradient.length
gradient.length <- function(object, ...)
{
  UseMethod(generic = "gradient.length")
}
