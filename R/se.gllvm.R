#' @title Standard errors for gllvm model
#' @description Calculates Hessian and standard errors for gllvm model.
#'
#' @param object an object of class 'gllvm'.
#' @param ... not used.
#'
#' @details
#' Computes Hessian and standard errors for gllvm model.
#' 
#' @return 
#'  \item{sd }{ list of standard errors of parameters}
#'  \item{Hess }{ list including Hessian matrix and approximative covariance matrix of parameters}
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @references
#'
#' Dunn, P. K., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics, 5, 236-244.
#'
#' Hui, F. K. C., Taskinen, S., Pledger, S., Foster, S. D., and Warton, D. I. (2015).  Model-based approaches to unconstrained ordination. Methods in Ecology and Evolution, 6:399-411.
#'
#'@aliases se se.gllvm
#'@method se gllvm
#'@export
#'@export se.gllvm
se.gllvm <- function(object, ...){
  if(!is.finite(object$logL)) stop("Standard errors can not be calculated if log-likelihood value is not finite.")
  if(object$TMB == FALSE) stop("Function is not implemented for TMB = FALSE.")
  objrFinal <- object$TMBfn
  n <- nrow(object$y)
  p <- ncol(object$y)
  method <- object$method
  num.lv <- object$num.lv
  quadratic <- object$quadratic
  nlvr <- num.lv + (object$row.eff=="random")*1
  family = object$family
  familyn <- objrFinal$env$data$family
  out <- list()
  if (!is.null(object$TR)) {
    {
      if(object$method == "VA"){
        sdr <- objrFinal$he(objrFinal$par)
      }
      if(object$method == "LA"){
        pars <- objrFinal$par
        if(family=="ZIP") {
          p0i <- names(pars)=="lg_phi"
          p0 <- pars[p0i]
          p0 <- p0+runif(p,0,0.001)
          pars[p0i] <- p0
        }
        sdr <- optimHess(pars, objrFinal$fn, objrFinal$gr)
      }

      m <- dim(sdr)[1]; incl <- rep(TRUE,m); incld <- rep(FALSE,m)
      incl[names(objrFinal$par)=="Abb"] <- FALSE;
      if(quadratic == FALSE){incl[names(objrFinal$par)=="lambda2"]<-FALSE}
      
      incl[names(objrFinal$par)=="Au"] <- FALSE; 
      if(nlvr > 0) incld[names(objrFinal$par)=="Au"] <- TRUE
      
      if(object$beta0com){ 
        incl[names(objrFinal$par)=="b"] <- FALSE
      }
      
      if(object$row.eff=="random") {
        incl[names(objrFinal$par)=="r0"] <- FALSE; incld[names(objrFinal$par)=="r0"] <- FALSE
      } else {
        incl[names(objrFinal$par)=="log_sigma"] <- FALSE
      }
      if(object$row.eff==FALSE) incl[names(objrFinal$par)=="r0"] <- FALSE
      if(object$row.eff=="fixed") incl[1] <- FALSE
      
      
      if(is.null(object$randomX)) {
        incl[names(objrFinal$par)%in%c("Br","sigmaB","sigmaij")] <- FALSE
      } else {
        xb <- objrFinal$env$data$xb
        incl[names(objrFinal$par)=="Abb"] <- FALSE; incld[names(objrFinal$par)=="Abb"] <- TRUE
        incl[names(objrFinal$par)=="Br"] <- FALSE; incld[names(objrFinal$par)=="Br"] <- TRUE
        if(NCOL(xb)==1) incl[names(objrFinal$par) == "sigmaij"] <- FALSE
      }
      
      incl[names(objrFinal$par)=="Au"] <- FALSE; if(num.lv>0) incld[names(objrFinal$par)=="Au"] <- TRUE
      incl[names(objrFinal$par)=="u"] <- FALSE; incld[names(objrFinal$par)=="u"] <- TRUE
      
      if(familyn==0 || familyn==2 || familyn==7 || familyn==8) incl[names(objrFinal$par)=="lg_phi"] <- FALSE
      if(familyn!=7) incl[names(objrFinal$par)=="zeta"] <- FALSE
      if(familyn==7) incl[names(objrFinal$par)=="zeta"] <- TRUE
      
      if(nlvr==0){
        incl[names(objrFinal$par)=="u"] <- FALSE;
        incld[names(objrFinal$par)=="u"] <- FALSE;
        incl[names(objrFinal$par)=="lambda"] <- FALSE;
        incl[names(objrFinal$par)=="lambda2"] <- FALSE;
        incl[names(objrFinal$par)=="Au"] <- FALSE;
      }
      
      if(method=="LA" || (num.lv==0 && (object$row.eff!="random" && is.null(object$randomX)))){
        incl[names(objrFinal$par)=="Au"] <- FALSE;
        
        covM <- try(MASS::ginv(sdr[incl,incl]))
        se <- try(sqrt(diag(abs(covM))))
        if(num.lv > 0 || object$row.eff == "random" || !is.null(object$randomX)) {
          sd.random <- sdrandom(objrFinal, covM, incl)
          prediction.errors <- list()
          if(!is.null(object$randomX)){
            prediction.errors$Br  <- matrix(diag(as.matrix(sd.random))[1:(ncol(xb)*p)], ncol(xb), p);
            sd.random <- sd.random[-(1:(ncol(xb)*p)),-(1:(ncol(xb)*p))]
          }
          if(object$row.eff=="random"){
            prediction.errors$row.params <- diag(as.matrix(sd.random))[1:n];
            sd.random <- sd.random[-(1:n),-(1:n)]
          }
          if(num.lv > 0){
            cov.lvs <- array(0, dim = c(n, num.lv, num.lv))
            for (i in 1:n) {
              cov.lvs[i,,] <- as.matrix(sd.random[(0:(num.lv-1)*n+i),(0:(num.lv-1)*n+i)])
            }
            prediction.errors$lvs <- cov.lvs
          }
          out$prediction.errors <- prediction.errors
        }
        out$Hess <- list(Hess.full=sdr, incl=incl, cov.mat.mod=covM)
      } else {
        A.mat <- sdr[incl,incl] # a x a
        D.mat <- sdr[incld,incld] # d x d
        B.mat <- sdr[incl,incld] # a x d
        cov.mat.mod <- try(MASS::ginv(A.mat-B.mat%*%solve(D.mat)%*%t(B.mat)),silent=T) 
        se <- sqrt(diag(abs(cov.mat.mod)))
        
        incla<-rep(FALSE, length(incl))
        incla[names(objrFinal$par)=="u"] <- TRUE
        out$Hess <- list(Hess.full=sdr, incla = incla, incl=incl, incld=incld, cov.mat.mod=cov.mat.mod)
        
      }
      
      if(object$row.eff=="fixed") {
        se.row.params <- c(0,se[1:(n-1)]); 
        names(se.row.params)  = rownames(object$y); se <- se[-(1:(n-1))] 
      }
      if(object$beta0com){
        se.beta0 <- se[1]; se <- se[-1];
      } else {
        se.beta0 <- se[1:p]; se <- se[-(1:p)];
      }
      se.B <- se[1:length(object$params$B)]; se <- se[-(1:length(object$params$B))];
      if(num.lv>0) {
        se.theta <- matrix(0,p,num.lv); se.theta[lower.tri(se.theta, diag = TRUE)]<-se[1:(p * num.lv - sum(0:(num.lv-1)))];
        colnames(se.theta) <- paste("LV", 1:num.lv, sep="");
        rownames(se.theta) <- colnames(object$y)
        out$sd$theta <- se.theta; se <- se[-(1:(p * num.lv - sum(0:(num.lv-1))))];
        if(quadratic==TRUE){
          se.lambdas2 <- matrix(se[1:(p * num.lv)], p, num.lv, byrow = T)  
          colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
          se <- se[-(1:(num.lv*p))]
          out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
        }else if(quadratic=="LV"){
          se.lambdas2 <- matrix(se[1:num.lv], p, num.lv, byrow = T)
          colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
          se <- se[-(1:num.lv)]
          out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
        }
        # diag(out$sd$theta) <- diag(out$sd$theta)*diag(object$params$theta) !!!
      }
      out$sd$beta0 <- se.beta0; 
      if(!object$beta0com){ names(out$sd$beta0)  <-  colnames(object$y);}
      out$sd$B <- se.B; names(out$sd$B) <- colnames(objrFinal$env$data$x)
      if(object$row.eff=="fixed") {out$sd$row.params <- se.row.params}
      
      if(family %in% c("negative.binomial")) {
        se.lphis <- se[1:p];  out$sd$inv.phi <- se.lphis*object$params$inv.phi;
        out$sd$phi <- se.lphis*object$params$phi;
        names(out$sd$inv.phi) <- names(out$sd$phi) <- colnames(object$y);  se <- se[-(1:p)]
      }
      if(family %in% c("gaussian","tweedie","gamma")) {
        se.lphis <- se[1:p];
        out$sd$phi <- se.lphis*object$params$phi;
        names(out$sd$phi) <- colnames(object$y);  se <- se[-(1:p)]
      }
      if(family %in% c("ZIP")) {
        lp0 <- objrFinal$par[names(objrFinal$par)=="lg_phi"]
        se.phis <- se[1:p];
        out$sd$phi <- se.phis*exp(lp0)/(1+exp(lp0))^2;#
        names(out$sd$phi) <- colnames(object$y);  se <- se[-(1:p)]
      }
      if(!is.null(object$randomX)){
        nr <- ncol(xb)
        out$sd$sigmaB <- se[1:ncol(xb)]*c(sqrt(diag(object$params$sigmaB))); 
        names(out$sd$sigmaB) <- c(paste("sd",colnames(xb),sep = "."))
        se <- se[-(1:ncol(xb))]
        if(nr>1){
          out$sd$corrpar <- se[1:(nr*(nr-1)/2)]
          se <- se[-(1:(nr*(nr-1)/2))]
        }
      }
      if(object$row.eff=="random") { 
        out$sd$sigma <- se[1:length(object$params$sigma)]*c(object$params$sigma[1],rep(1,length(object$params$sigma)-1)); 
        names(out$sd$sigma) <- "sigma"; 
        se=se[-(1:(length(object$params$sigma)))] 
      }
      if(family %in% c("ordinal")){
        y <- object$y
        K = max(y)-min(y)
        if(min(y)==0) y <- y+1 
        se.zetanew <- se.zetas <- se;
        if(object$zeta.struc == "species"){
          se.zetanew <- matrix(NA,nrow=p,ncol=K)
          idx<-0
          for(j in 1:ncol(y)){
            k<-max(y[,j])-2
            if(k>0){
              for(l in 1:k){
                se.zetanew[j,l+1]<-se.zetas[idx+l]
              } 
            }
            idx<-idx+k
          }
          se.zetanew[,1] <- 0
          out$sd$zeta <- se.zetanew
          row.names(out$sd$zeta) <- colnames(object$y); colnames(out$sd$zeta) <- paste(min(object$y):(max(object$y)-1),"|",(min(object$y)+1):max(object$y),sep="")
          
        }else{
          se.zetanew <- c(0, se.zetanew)
          out$sd$zeta <- se.zetanew
          names(out$sd$zeta) <- paste(min(object$y):(max(object$y)-1),"|",(min(object$y)+1):max(object$y),sep="")
          
        }
      }
      
    }
  } else {

    pars <- objrFinal$par
    if(family=="ZIP") {
      p0i <- names(pars)=="lg_phi"
      p0 <- pars[p0i]
      p0 <- p0+runif(p,0,0.001)
      pars[p0i] <- p0
    }
    if(method == "VA"){
      sdr <- objrFinal$he(pars)
    }
    if(method == "LA"){
      sdr <- optimHess(pars, objrFinal$fn, objrFinal$gr)
    }
    m <- dim(sdr)[1]; incl <- rep(TRUE,m); incld <- rep(FALSE,m); inclr <- rep(FALSE,m)
    incl[names(objrFinal$par)=="B"] <- FALSE
    incl[names(objrFinal$par)%in%c("Br","sigmaB","sigmaij")] <- FALSE
    incl[names(objrFinal$par)=="Abb"]=FALSE;
    if(quadratic == FALSE){incl[names(objrFinal$par)=="lambda2"]<-FALSE}
    
    if(familyn!=7) incl[names(objrFinal$par)=="zeta"] <- FALSE
    
    if(method=="LA" || (num.lv==0 && method=="VA" && object$row.eff!="random")){
      incl[names(objrFinal$par)=="Au"] <- FALSE;
      if(object$row.eff=="random") {
        incl[names(objrFinal$par)=="r0"] <- FALSE; incld[names(objrFinal$par)=="r0"] <- FALSE
      } 
      if(object$row.eff=="fixed"){ incl[1] <- FALSE; incl[names(objrFinal$par)=="log_sigma"] <- FALSE}
      if(object$row.eff==FALSE) {incl[names(objrFinal$par)=="r0"] <- FALSE; incl[names(objrFinal$par)=="log_sigma"] <- FALSE}
      if(familyn==0 || familyn==2 || familyn==7 || familyn==8) incl[names(objrFinal$par)=="lg_phi"] <- FALSE
      if(familyn==7) incl[names(objrFinal$par)=="zeta"] <- TRUE
      if(nlvr==0){
        incl[names(objrFinal$par)=="u"] <- FALSE;
        incl[names(objrFinal$par)=="lambda"] <- FALSE;
      }
      covM <- try(MASS::ginv(sdr[incl,incl]))
      se <- try(sqrt(diag(abs(covM))))
      if(nlvr>0){
        sd.random <- sdrandom(objrFinal, covM, incl, ignore.u = FALSE)
        prediction.errors <- list()
        # if(object$row.eff=="random" && FALSE){
        #   prediction.errors$row.params <- diag(as.matrix(sd.random))[1:n];
        #   sd.random <- sd.random[-(1:n),-(1:n)]
        # }
        if(nlvr>0){
          cov.lvs <- array(0, dim = c(n, nlvr, nlvr))
          # cov.lvs <- array(0, dim = c(n, num.lv, num.lv))
          for (i in 1:n) {
            cov.lvs[i,,] <- as.matrix(sd.random[(0:(nlvr-1)*n+i),(0:(nlvr-1)*n+i)])
            # cov.lvs[i,,] <- as.matrix(sd.random[(0:(num.lv-1)*n+i),(0:(num.lv-1)*n+i)])
          }
          if(object$row.eff=="random"){
            prediction.errors$row.params <- cov.lvs[,1,1]
            if(num.lv > 0) cov.lvs <- array(cov.lvs[,-1,-1], dim = c(n, num.lv, num.lv))
          }
          
          prediction.errors$lvs <- cov.lvs
          #sd.random <- sd.random[-(1:(n*num.lv))]
        }
        out$prediction.errors <- prediction.errors
      }
      out$Hess <- list(Hess.full=sdr, incl=incl, cov.mat.mod=covM)
      
    } else {
      incl[names(objrFinal$par)=="Au"] <- FALSE;
      
      if(object$row.eff=="random") {
        inclr[names(objrFinal$par) == "r0"] <- FALSE;
        incl[names(objrFinal$par) == "r0"] <- FALSE; incld[names(objrFinal$par) == "r0"] <- FALSE
        incld[names(objrFinal$par)=="Au"] <- TRUE
      }
      if(object$row.eff=="fixed") {incl[1] <- FALSE; incl[names(objrFinal$par)=="log_sigma"] <- FALSE}
      if(object$row.eff==FALSE) {incl[names(objrFinal$par)=="r0"] <- FALSE; incl[names(objrFinal$par)=="log_sigma"] <- FALSE}
      
      if(nlvr>0){
        inclr[names(objrFinal$par)=="u"] <- TRUE;
        incl[names(objrFinal$par)=="u"] <- FALSE;
        incld[names(objrFinal$par)=="u"] <- TRUE;
        incld[names(objrFinal$par)=="Au"] <- TRUE;
      } else {
        incl[names(objrFinal$par)=="u"] <- FALSE;
        incl[names(objrFinal$par)=="lambda"] <- FALSE;
      }
      if(familyn==0 || familyn==2 || familyn==7 || familyn==8) incl[names(objrFinal$par)=="lg_phi"] <- FALSE
      
      A.mat <- sdr[incl, incl] # a x a
      D.mat <- sdr[incld, incld] # d x d
      B.mat <- sdr[incl, incld] # a x d
      cov.mat.mod <- try(MASS::ginv(A.mat-B.mat%*%solve(D.mat)%*%t(B.mat)),silent=T)
      se <- sqrt(diag(abs(cov.mat.mod)))
      
      incla<-rep(FALSE, length(incl))
      incla[names(objrFinal$par)=="u"] <- TRUE
      out$Hess <- list(Hess.full=sdr, incla = incla, incl=incl, incld=incld, cov.mat.mod=cov.mat.mod)
      
    }
    
    num.X <- 0; if(!is.null(object$X)) num.X <- dim(object$X)[2]
    if(object$row.eff == "fixed") { se.row.params <- c(0,se[1:(n-1)]); names(se.row.params) <- rownames(object$y); se <- se[-(1:(n-1))] }
    sebetaM <- matrix(se[1:((num.X+1)*p)],p,num.X+1,byrow=TRUE);  se <- se[-(1:((num.X+1)*p))]
    if(num.lv > 0) {
      se.lambdas <- matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag = TRUE)] <- se[1:(p * num.lv - sum(0:(num.lv-1)))];
      colnames(se.lambdas) <- paste("LV", 1:num.lv, sep="");
      rownames(se.lambdas) <- colnames(object$y)
      out$sd$theta <- se.lambdas; se <- se[-(1:(p * num.lv - sum(0:(num.lv-1))))];
      diag(out$sd$theta) <- diag(out$sd$theta)*diag(object$params$theta[,1:object$num.lv,drop=F]) 
      if(quadratic==TRUE){
        se.lambdas2 <- matrix(se[1:(p * num.lv)], p, num.lv, byrow = T)  
        colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
        se <- se[-(1:(num.lv*p))]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }else if(quadratic=="LV"){
        se.lambdas2 <- matrix(se[1:num.lv], p, num.lv, byrow = T)
        colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
        se <- se[-(1:num.lv)]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }
    }
    
    out$sd$beta0 <- sebetaM[,1]; names(out$sd$beta0) <- colnames(object$y);
    if(!is.null(object$X)){
      out$sd$Xcoef <- matrix(sebetaM[,-1],nrow = nrow(sebetaM));
      rownames(out$sd$Xcoef) <- colnames(object$y); colnames(out$sd$Xcoef) <- colnames(object$X);
    }
    if(object$row.eff=="fixed") {out$sd$row.params <- se.row.params}
    
    if(family %in% c("negative.binomial")) {
      se.lphis <- se[1:p];  out$sd$inv.phi <- se.lphis*object$params$inv.phi;
      out$sd$phi <- se.lphis*object$params$phi;
      names(out$sd$phi) <- colnames(object$y);  se <- se[-(1:p)]
    }
    if(family %in% c("tweedie", "gaussian", "gamma")) {
      se.lphis <- se[1:p];
      out$sd$phi <- se.lphis*object$params$phi;
      names(out$sd$phi) <- colnames(object$y);  se <- se[-(1:p)]
    }
    if(family %in% c("ZIP")) {
      lp0 <- objrFinal$par[names(objrFinal$par)=="lg_phi"]
      se.phis <- se[1:p];
      out$sd$phi <- se.phis*exp(lp0)/(1+exp(lp0))^2;#
      names(out$sd$phi) <- colnames(object$y);  se <- se[-(1:p)]
    }
    if(object$row.eff=="random") { 
      out$sd$sigma <- se[1:length(object$params$sigma)]*c(object$params$sigma[1],rep(1,length(object$params$sigma)-1)); 
      se=se[-(1:length(out$sd$sigma))] 
      names(out$sd$sigma) <- "sigma" 
    }
    if(family %in% c("ordinal")){
      y <- object$y
      K = max(y)-min(y)
      if(min(y)==0) y <- y+1 
      se.zetanew <- se.zetas <- se;
      if(object$zeta.struc == "species"){
        se.zetanew <- matrix(NA,nrow=p,ncol=K)
        idx<-0
        for(j in 1:ncol(y)){
          k<-max(y[,j])-2
          if(k>0){
            for(l in 1:k){
              se.zetanew[j,l+1]<-se.zetas[idx+l]
            } 
          }
          idx<-idx+k
        }
        se.zetanew[,1] <- 0
        out$sd$zeta <- se.zetanew
        row.names(out$sd$zeta) <- colnames(object$y); colnames(out$sd$zeta) <- paste(min(object$y):(max(object$y)-1),"|",(min(object$y)+1):max(object$y),sep="")
        
      }else{
        se.zetanew <- c(0, se.zetanew)
        out$sd$zeta <- se.zetanew
        names(out$sd$zeta) <- paste(min(object$y):(max(object$y)-1),"|",(min(object$y)+1):max(object$y),sep="")
        
      }
    }
    
  }
  return(out)
}


#' @export se
se <- function(object, ...)
{
  UseMethod(generic = "se")
}
