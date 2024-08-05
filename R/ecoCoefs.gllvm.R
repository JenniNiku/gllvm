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
  
  opt<-object$params$theta[,1:(num.lv+num.lv.c),drop=FALSE]/(2*abs(object$params$theta[,-c(1:(num.lv+num.lv.c)),drop=FALSE]))
  if(num.lv>0){
    #correct for sig
    theta <- object$params$theta[,(num.lv.c+1):(num.lv.c+num.lv)]
    theta2 <- abs(object$params$theta[,(ncol(object$params$theta)-num.lv+1):ncol(object$params$theta)])
    opt[,(num.lv.c+1):ncol(opt)] <- theta%*%diag(tail(object$params$sigma.lv, num.lv), num.lv)/(2*theta2%*%diag(tail(object$params$sigma.lv, num.lv)^2, num.lv))
  }
  if(!isFALSE(object$randomB))sd.errors = FALSE
  if(sd.errors==TRUE){
    if(is.null(object$sd)|all(unlist(object$sd)==FALSE)){
      cat("Standard errors not present in model, calculating...\n")
      object$Hess<-se.gllvm(object)$Hess 
    }
    V <- object$Hess$cov.mat.mod
    colnames(V) <- row.names(V) <- names(object$TMBfn$par[object$Hess$incl])
    # switch sigmaLV tp last entries
    pars <- object$TMBfn$par[object$Hess$incl]
    pars <- pars[c((1:length(pars))[-which(names(pars)==c("sigmaLV"))],which(names(pars)==c("sigmaLV")))]
    idx <- names(pars)%in%c("lambda","lambda2")
    if(num.lv>0){
      # switch sigmaLV to last entries
      V <- V[c((1:ncol(V))[-which(names(object$TMBfn$par)==c("sigmaLV"))],which(names(object$TMBfn$par)==c("sigmaLV"))),c((1:ncol(V))[-which(names(object$TMBfn$par)==c("sigmaLV"))],which(names(object$TMBfn$par)==c("sigmaLV")))]
      # include last num.lv sigma
      idx[tail(which(names(pars)=="sigmaLV"), num.lv)] <- TRUE
    }
    V <- V[idx,idx,drop=F]
    
    if(num.lv>0&num.lv.c==0)idx<-which(c(upper.tri(object$params$theta[,1:num.lv],diag=T)))[-1]
    if(num.lv.c>0)idx<-which(c(upper.tri(object$params$theta[,1:num.lv.c],diag=T)))[-1]
    
    #add first row and column of zeros
    V<-rbind(0,cbind(0,V))
    
    if(length(idx)>0){
    #add zeros where necessary
    for(q in 1:length(idx)){
      V <- rbind(V[1:(idx[q]-1),],0,V[idx[q]:ncol(V),])
      V <- cbind(V[,1:(idx[q]-1)],0,V[,idx[q]:ncol(V)])
    }
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
    
    du1 <- NULL #linear coefficients
    du2 <- NULL #quadratic coefficients
    #du1 needs to be species then LV, to accommmodate the order in V
    if(num.lv.c>0){
      for(i in 1:(num.lv.c)){
        for (j in 1:p) {
          du1 <- c(du1,0.5*(object$params$theta[j, (num.lv+num.lv.c) + i])^-1)
          du2 <- c(du2, 0.5*object$params$theta[j, i] / (object$params$theta[j, (num.lv+num.lv.c) + i]^2))
        }
      }
    }
    ##separated due to diagonal
    du3 <- NULL # sig for num.lv
    if(num.lv>0){
      sigm <- tail(object$params$sigma.lv, num.lv)
      for(i in 1:num.lv){
        for (j in 1:p) {
          du1 <- c(du1, 0.5*(object$params$theta[j, (num.lv+2*num.lv.c) + i]*sigm[i])^-1)
          du2 <- c(du2, 0.5*(object$params$theta[j, num.lv.c+i]*sigm[i]^3) / (object$params$theta[j, (num.lv+2*num.lv.c) + i]^2*sigm[i]^4))
          du3 <- c(du3, object$params$theta[j, num.lv.c+i]*(-0.5/(object$params$theta[j, (num.lv+2*num.lv.c) + i]*sigm[i]^2) + object$params$theta[j, (num.lv+num.lv.c+num.lv.c) + i]*sigm[i]^2/(object$params$theta[j, (num.lv+num.lv.c+num.lv.c) + i]^2*sigm[i]^4))*sign(tail(object$TMBfn$par[names(object$TMBfn$par)=="sigmaLV"],num.lv))[i])
        }
      }
    }
    
    #re-arrange quadratic coefficients to similar order as linear: per species per latent variable (now per latent variable and then species)
    #this makes thing easier further on
    if(quadratic==TRUE){
      idx <- c(1:(p*(num.lv+num.lv.c)),(p*(num.lv+num.lv.c)+c(matrix(1:(p*(num.lv+num.lv.c)),ncol=(num.lv+num.lv.c),nrow=p,byrow=T))))  
    }else{
      idx <- c(1:(p*(num.lv+num.lv.c)),(p*(num.lv+num.lv.c))+rep(1:(num.lv+num.lv.c),each=p))
    }
    # Separate covariance matrix here
    if(num.lv>0){
      Vsig <- V[c(idx, tail(1:ncol(V), num.lv)), c(idx, tail(1:ncol(V), num.lv))]
      V <- V[idx,idx]
    }
    
    #sum over parameters and their covariances
    #linear coefficients covariances
    cov.mat.lincoef<-du1%*%t(du1)*V[1:(p*(num.lv+num.lv.c)),1:(p*(num.lv+num.lv.c))]
    #quadratic coefficients covariances
    cov.mat.quadcoef <- du2%*%t(du2)*V[-c(1:(p*(num.lv+num.lv.c))),-c(1:(p*(num.lv+num.lv.c)))]
    #linear and quadratic coefficient covariances
    diag.cov.mat.coef <- du1%*%t(du2)*V[1:(p*(num.lv+num.lv.c)),-c(1:(p*(num.lv+num.lv.c)))]
    #covariance of optima per lv per species
    cov.mat.optima <- cov.mat.lincoef+cov.mat.quadcoef+2*diag.cov.mat.coef
    if(num.lv>0){
      # idx vector for repeating sig
      idx <- rep(1:num.lv,each=p)
      # sig covariances
      cov.sig <- Vsig[which(colnames(Vsig)=="sigmaLV"),which(colnames(Vsig)=="sigmaLV"),drop=FALSE]
      # sig cov with lin.coef
      cov.mat.sig.lin <- Vsig[which(colnames(Vsig)=="sigmaLV"),tail(which(colnames(Vsig)=="lambda"),num.lv*p),drop=FALSE]
      # sig cov with quad.coef
      cov.mat.sig.quad <- Vsig[which(colnames(Vsig)=="sigmaLV"),tail(which(colnames(Vsig)=="lambda2"),num.lv*p),drop=FALSE]
      
      cov.mat.optima[(num.lv.c*p+1):ncol(cov.mat.optima),(num.lv.c*p+1):ncol(cov.mat.optima)] <-
        cov.mat.optima[(num.lv.c*p+1):ncol(cov.mat.optima),(num.lv.c*p+1):ncol(cov.mat.optima)] + 
        du3%*%t(du3)*cov.sig[idx,idx]+
        2*du3%*%t(tail(du1,num.lv*p))*cov.mat.sig.lin[idx,]+
        2*du3%*%t(tail(du2,num.lv*p))*cov.mat.sig.quad[idx,]
    }
    opt.sd <- matrix(sqrt(abs(diag(cov.mat.optima))),ncol=(num.lv+num.lv.c),nrow=p)
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
  theta <- object$params$theta[,-c(1:(num.lv.c+num.lv)),drop=FALSE]
  if(num.lv>0){
    theta[,(num.lv.c+1):(num.lv+num.lv.c)] <- theta[,(num.lv.c+1):(num.lv+num.lv.c)]%*%diag(tail(object$params$sigma.lv^2,num.lv), num.lv)
  }
  
  tol<-1/sqrt(-2*theta)
  
  if(!isFALSE(object$randomB))sd.errors = FALSE
  if(sd.errors==TRUE){
    if(is.null(object$sd)|all(unlist(object$sd)==FALSE)){
      cat("Standard errors not present in model, calculating...\n")
      object$Hess<-se.gllvm(object)$Hess 
    }
    thetasig = sign(matrix(object$TMBfn$par[names(object$TMBfn$par)=="lambda2"], ncol = num.lv.c+num.lv, byrow=TRUE))
    if(quadratic=="LV"){
      thetasig <- thetasig[rep(1,p),,drop=FALSE]
    }
    V <- object$Hess$cov.mat.mod
    colnames(V) <- row.names(V) <- names(object$TMBfn$par)[object$Hess$incl]
    # Place sigmaLV at the end
    V <- V[c((1:ncol(V))[-which(colnames(V)=="sigmaLV")],(1:ncol(V))[which(colnames(V)=="sigmaLV")]),c((1:ncol(V))[-which(colnames(V)=="sigmaLV")],(1:ncol(V))[which(colnames(V)=="sigmaLV")])]
    idx <- colnames(V)%in%c("lambda2")
    if(num.lv>0)idx[tail(1:ncol(V),num.lv)] <- TRUE
    
    V <- V[idx,idx,drop=F]
    
    #re-arrange V for easier access
    if(quadratic==TRUE){
      idx <- c(matrix(1:(p*(num.lv+num.lv.c)),ncol=(num.lv+num.lv.c),nrow=p,byrow=T))
    }else{
      idx <- rep(1:(num.lv+num.lv.c),each=p)
    }
    if(num.lv>0)Vsig <- V[c(idx, tail(1:ncol(V),num.lv)), c(idx, tail(1:ncol(V),num.lv)),drop=FALSE]
    V<-V[idx,idx]
    
    dt<- NULL
    if(num.lv.c>0){
    for (i in 1:(num.lv.c)) {
      for (j in 1:p) {
          dt <- c(dt, -thetasig[j,i]*0.5*2^-0.5*abs(object$params$theta[j,num.lv+num.lv.c+i])^-1.5)
        }
      }
    }
    
    if(num.lv>0){
      dsig <- NULL
      sig <- tail(object$params$sig,num.lv)
      for (i in 1:(num.lv)) {
      for (j in 1:p) {
          dt <- c(dt,  -thetasig[j,num.lv.c+i]*0.5*2^-0.5*abs(object$params$theta[j,(num.lv+2*num.lv.c)+i]*sig[i]^2)^-1.5)
          dsig <- c(dsig, -(sig[i]^2*2^0.5*abs(object$params$theta[j,(num.lv+2*num.lv.c)+i])^0.5)^-1)
        }
      }
      
      cov.mat.sig.quad <- Vsig[which(colnames(Vsig)=="sigmaLV"),tail(which(colnames(Vsig)=="lambda2"),num.lv*p),drop=FALSE][rep(1:num.lv,each=p),]
      
    }

    tol.cov <- dt%*%t(dt)*V
    #for sig
    if(num.lv>0){
      tol.cov[(num.lv.c*p+1):ncol(tol.cov),(num.lv.c*p+1):ncol(tol.cov)] = tol.cov[(num.lv.c*p+1):ncol(tol.cov),(num.lv.c*p+1):ncol(tol.cov)]+ dsig%*%t(dsig)*Vsig[colnames(Vsig)=="sigmaLV",colnames(Vsig)=="sigmaLV",drop=FALSE][rep(1:num.lv,each=p),rep(1:num.lv,each=p)] + dsig%*%t(tail(dt,num.lv*p))*cov.mat.sig.quad
    }
    tol.sd = sqrt(abs(matrix(diag(tol.cov), ncol = num.lv+num.lv.c)))
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

