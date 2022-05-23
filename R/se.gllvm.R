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
  num.lv.c <- object$num.lv.c
  num.RR <- object$num.RR
  lv.X <- object$lv.X
  
  quadratic <- object$quadratic
  nlvr <- num.lv + num.lv.c 
  nlvr <- num.lv  #+ (object$row.eff=="random")*1
  cstrucn = switch(object$cstruc, "diag" = 0, "corAR1" = 1, "corExp" = 2, "corCS" = 3)
  rstruc = object$rstruc
  family = object$family
  familyn <- objrFinal$env$data$family
  disp.group <- object$disp.group
  out <- list()
  if (!is.null(object$TR)) {
    {
      if((object$method %in% c("VA", "EVA"))){
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
      
      # Variational params not included for incl

      if(quadratic == FALSE){incl[names(objrFinal$par)=="lambda2"]<-FALSE}
      
      # Terms that are not included in TMBtrait for constrained ordination
      incl[names(objrFinal$par)=="Ab_lv"] <- FALSE;
      incl[names(objrFinal$par)=="sigmab_lv"] <- FALSE;
      incl[names(objrFinal$par)=="b_lv"] <- FALSE
      
      incl[names(objrFinal$par)=="lg_Ar"] <- FALSE;
      incl[names(objrFinal$par)=="Au"] <- FALSE;
      incl[names(objrFinal$par)=="u"] <- FALSE; 

      if(quadratic == FALSE){incl[names(objrFinal$par)=="lambda2"]<-FALSE}
      if(object$beta0com){ incl[names(objrFinal$par)=="b"] <- FALSE}
      if(familyn!=7) incl[names(objrFinal$par)=="zeta"] <- FALSE
      if(familyn==0 || familyn==2 || familyn==7 || familyn==8) incl[names(objrFinal$par)=="lg_phi"] <- FALSE

      
      if(num.lv>0) {
        incld[names(objrFinal$par)=="Au"] <- TRUE
        incld[names(objrFinal$par)=="u"] <- TRUE
        incl[names(objrFinal$par)=="lambda2"] <- FALSE;
        incl[names(objrFinal$par)=="lambda"] <- FALSE;
      }
      
      
      if(object$row.eff=="random") {
        incld[names(objrFinal$par)=="lg_Ar"] <- TRUE
        incld[names(objrFinal$par)=="r0"] <- TRUE
        incl[names(objrFinal$par)=="r0"] <- FALSE; 
      } else {
        incl[names(objrFinal$par)=="log_sigma"] <- FALSE
        if(object$row.eff==FALSE) incl[names(objrFinal$par)=="r0"] <- FALSE
        if(object$row.eff=="fixed") incl[1] <- FALSE
      }
      
      
      if(is.null(object$randomX)) {
        incl[names(objrFinal$par)%in%c("Br","sigmaB","sigmaij")] <- FALSE
      } else {
        xb <- objrFinal$env$data$xb
        incl[names(objrFinal$par)=="Abb"] <- FALSE; incld[names(objrFinal$par)=="Abb"] <- TRUE
        incl[names(objrFinal$par)=="Br"] <- FALSE; incld[names(objrFinal$par)=="Br"] <- TRUE
        if(NCOL(xb)==1) incl[names(objrFinal$par) == "sigmaij"] <- FALSE
      }
 
      if(method=="LA" || (num.lv==0 && (object$row.eff!="random" && is.null(object$randomX)))){
        covM <- try(MASS::ginv(sdr[incl,incl]))
        se <- try(sqrt(diag(abs(covM))))
        if(num.lv > 0 || object$row.eff == "random" || !is.null(object$randomX)) {
          sd.random <- sdrandom(objrFinal, covM, incl)
          prediction.errors <- list()
          
          if(object$row.eff=="random"){
            prediction.errors$row.params <- sd.random$row
          }
          if(!is.null(object$randomX)){
            prediction.errors$Br  <- sd.random$Ab
          }

          if(num.lv > 0){
            prediction.errors$lvs <- sd.random$A
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
      if(num.lv > 0) {
        se.sigma.lv <- se[1:num.lv];se<-se[-c(1:num.lv)]
        se.lambdas <- matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag=FALSE)] <- se[1:(p * num.lv - sum(0:num.lv))];
        colnames(se.lambdas) <- paste("LV", 1:num.lv, sep="");
        rownames(se.lambdas) <- colnames(out$y)
        out$sd$theta <- se.lambdas; se <- se[-(1:(p * num.lv - sum(0:num.lv)))];
        
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
        out$sd$sigma.lv  <- se.sigma.lv
        names(out$sd$sigma.lv) <- colnames(object$params$theta[,1:num.lv])
      }

      out$sd$beta0 <- se.beta0; 
      if(!object$beta0com){ names(out$sd$beta0)  <-  colnames(object$y);}
      out$sd$B <- se.B; names(out$sd$B) <- colnames(objrFinal$env$data$x)
      if(object$row.eff=="fixed") {out$sd$row.params <- se.row.params}
      
      if(family %in% c("negative.binomial")) {
        se.lphis <- se[1:length(unique(disp.group))];  out$sd$inv.phi <- se.lphis*object$params$inv.phi;
        out$sd$phi <- se.lphis*object$params$phi;
        if(length(unique(disp.group))==p){
          names(out$sd$phi) <- colnames(y);
        }else if(!is.null(names(disp.group))){
          try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
        }else{
          names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
        }
        names(out$sd$inv.phi) <-  names(out$sd$phi)
        se <- se[-(1:length(unique(disp.group)))]
      }
      if(family %in% c("gaussian","tweedie","gamma", "beta")) {
        se.lphis <- se[1:length(unique(disp.group))];
        out$sd$phi <- se.lphis*object$params$phi;
        if(length(unique(disp.group))==p){
          names(out$sd$phi) <- colnames(y);
        }else if(!is.null(names(disp.group))){
          try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
        }else{
          names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
        }
        
        se <- se[-(1:length(unique(disp.group)))]
      }
      if(family %in% c("ZIP")) {
        se.phis <- se[1:length(unique(disp.group))];
        out$sd$phi <- se.phis*exp(object$lp0)/(1+exp(object$lp0))^2;#
        if(length(unique(disp.group))==p){
          names(out$sd$phi) <- colnames(y);
        }else if(!is.null(names(disp.group))){
          try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
        }else{
          names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
        }
        names(out$sd$inv.phi) <-  names(out$sd$phi)
        se <- se[-(1:length(unique(disp.group)))]
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
        if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,3))) {out$sd$rho <- se[1]*(1-object$params$rho^2)^1.5; se = se[-1]}
        if((rstruc ==2 | (rstruc == 1)) & (cstrucn ==2)) {out$sd$rho <- se[1]*object$params$rho; se = se[-1]}
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
    #Without traits#
    pars <- objrFinal$par
    if(family=="ZIP") {
      p0i <- names(pars)=="lg_phi"
      p0 <- pars[p0i]
      p0 <- p0+runif(p,0,0.001)
      pars[p0i] <- p0
    }
    if((object$method %in% c("VA", "EVA"))){
      sdr <- objrFinal$he(pars)
    }
    if(method == "LA"){
      sdr <- optimHess(pars, objrFinal$fn, objrFinal$gr)
    }
    m <- dim(sdr)[1]; incl <- rep(TRUE,m); incld <- rep(FALSE,m); inclr <- rep(FALSE,m)
    
    # Not used for this model
    incl[names(objrFinal$par)=="B"] <- FALSE
    incl[names(objrFinal$par)%in%c("Br","sigmaB","sigmaij")] <- FALSE
    
    # Variational params not included for incl
    incl[names(objrFinal$par)=="Abb"]=FALSE;
    incl[names(objrFinal$par)=="lg_Ar"] <- FALSE;
    incl[names(objrFinal$par)=="Au"] <- FALSE;
    incl[names(objrFinal$par)=="u"] <- FALSE;
    
    #slopes for reduced rank predictors
    if((num.lv.c+num.RR)==0){
      incl[names(objrFinal$par)=="sigmab_lv"] <- FALSE
      incl[names(objrFinal$par)=="b_lv"] <- FALSE
    }else if((num.lv.c+num.RR)>0&randomB==FALSE){
      incl[names(objrFinal$par)=="b_lv"] <- TRUE
      incl[names(objrFinal$par)=="sigmab_lv"] <- FALSE
    }else if((num.lv.c+num.RR>0)&randomB!=FALSE){
      incl[names(objrFinal$par)=="b_lv"] <- FALSE
      
      inclr[names(objrFinal$par)=="b_lv"] <- TRUE
      incld[names(objrFinal$par)=="b_lv"] <- TRUE
      
      incld[names(objrFinal$par)=="Ab_lv"] <- TRUE
      
      incl[names(objrFinal$par)=="sigmab_lv"] <- TRUE
    }
    
    if(quadratic == FALSE){incl[names(objrFinal$par)=="lambda2"]<-FALSE}
    if(familyn!=7) incl[names(objrFinal$par)=="zeta"] <- FALSE
    if(familyn==0 || familyn==2 || familyn==7 || familyn==8) incl[names(objrFinal$par)=="lg_phi"] <- FALSE

    if((num.lv+num.lv.c)>0){
      inclr[names(objrFinal$par)=="u"] <- TRUE;
      incld[names(objrFinal$par)=="u"] <- TRUE;
      incld[names(objrFinal$par)=="Au"] <- TRUE;
    } else {
      if(num.RR==0)incl[names(objrFinal$par)=="lambda"] <- FALSE;
      if(num.RR==0)incl[names(objrFinal$par)=="lambda2"] <- FALSE;
    }

    if(object$row.eff=="random") {
      incld[names(objrFinal$par) == "lg_Ar"] <- TRUE
      incld[names(objrFinal$par) == "r0"] <- TRUE
      inclr[names(objrFinal$par) == "r0"] <- TRUE;
      incl[names(objrFinal$par) == "r0"] <- FALSE; 
    } else {
      incl[names(objrFinal$par)=="log_sigma"] <- FALSE
      if(object$row.eff==FALSE) { incl[names(objrFinal$par)=="r0"] <- FALSE }
      if(object$row.eff=="fixed"){ incl[1] <- FALSE }
    }
    
    
    
    if(method=="LA" || ((num.lv+num.lv.c)==0 && (object$method %in% c("VA", "EVA")) && object$row.eff!="random" && object$randomB==FALSE)){
      covM <- try(MASS::ginv(sdr[incl,incl]))
      se <- try(sqrt(diag(abs(covM))))
      
      if((num.lv+num.lv.c) > 0 || object$row.eff == "random"){
        sd.random <- sdrandom(objrFinal, covM, incl, ignore.u = FALSE)
        prediction.errors <- list()
        
        if(object$row.eff=="random"){
          prediction.errors$row.params <- sd.random$row
        }
        if((num.lv+num.lv.c+num.RR)>0){
          # cov.lvs <- array(0, dim = c(n, nlvr, nlvr))
          cov.lvs <- sd.random$A
          # if(object$row.eff=="random"){
          #   prediction.errors$row.params <- cov.lvs[,1,1]
          #   if(num.lv > 0) cov.lvs <- array(cov.lvs[,-1,-1], dim = c(n, num.lv, num.lv))
          # }
          
          prediction.errors$lvs <- cov.lvs
          #sd.random <- sd.random[-(1:(n*num.lv))]
          if(randomB!=FALSE){
            prediction.errors$Ab.lv <- sd.random$Ab_lv
          }
        }
        out$prediction.errors <- prediction.errors
      }
      out$Hess <- list(Hess.full=sdr, incl=incl, cov.mat.mod=covM)
      
    } else {

      A.mat <- sdr[incl, incl] # a x a
      D.mat <- sdr[incld, incld] # d x d
      B.mat <- sdr[incl, incld] # a x d
      cov.mat.mod <- try(MASS::ginv(A.mat-B.mat%*%solve(D.mat)%*%t(B.mat)),silent=T)
      se <- sqrt(diag(abs(cov.mat.mod)))
      
      incla<-rep(FALSE, length(incl))
      incla[names(objrFinal$par)=="u"] <- TRUE
      out$Hess <- list(Hess.full=sdr, incla = incla, incl=incl, incld=incld, cov.mat.mod=cov.mat.mod)
      
    }
    
    num.X <- 0; if(!is.null(object$X)) num.X <- dim(object$X.design)[2]
    if(object$row.eff == "fixed") { se.row.params <- c(0,se[1:(n-1)]); names(se.row.params) <- rownames(object$y); se <- se[-(1:(n-1))] }
    sebetaM <- matrix(se[1:((num.X+1)*p)],p,num.X+1,byrow=TRUE);  se <- se[-(1:((num.X+1)*p))]
    
    if((num.lv.c+num.RR)>0&object$randomB==FALSE){
    se.LvXcoef <- matrix(se[1:((num.lv.c+num.RR)*ncol(lv.X))],ncol=num.lv.c+num.RR,nrow=ncol(lv.X))
    se <- se[-c(1:((num.lv.c+num.RR)*ncol(lv.X)))]
    colnames(se.LvXcoef) <- paste("CLV",1:(num.lv.c+num.RR),sep="")
    row.names(se.LvXcoef) <- colnames(lv.X)
    out$sd$LvXcoef <- se.LvXcoef
    }
    
    if((num.lv.c+num.lv)>0)se.sigma.lv <- se[1:(num.lv+num.lv.c)];se<-se[-c(1:(num.lv+num.lv.c))]
    if(num.lv > 0&(num.lv.c+num.RR)==0) {
      se.lambdas <- matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag=FALSE)] <- se[1:(p * num.lv - sum(0:num.lv))];
      colnames(se.lambdas) <- paste("LV", 1:num.lv, sep="");
      rownames(se.lambdas) <- colnames(out$y)
      out$sd$theta <- se.lambdas; se <- se[-(1:(p * num.lv - sum(0:num.lv)))];
      
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
    }else if(num.lv==0&(num.lv.c+num.RR)>0){
      se.lambdas <- matrix(0,p,(num.lv.c+num.RR)); se.lambdas[lower.tri(se.lambdas, diag=FALSE)] <- se[1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR)))];
      colnames(se.lambdas) <- paste("CLV", 1:(num.lv.c+num.RR), sep="");
      rownames(se.lambdas) <- colnames(out$y)
      out$sd$theta <- se.lambdas; se <- se[-(1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR))))];
      
      if(quadratic==TRUE){
        se.lambdas2 <- matrix(se[1:(p * (num.lv.c+num.RR))], p, (num.lv.c+num.RR), byrow = T)  
        colnames(se.lambdas2) <- paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = "")
        se <- se[-(1:((num.lv.c+num.RR)*p))]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }else if(quadratic=="LV"){
        se.lambdas2 <- matrix(se[1:(num.lv.c+num.RR)], p, (num.lv.c+num.RR), byrow = T)
        colnames(se.lambdas2) <- paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = "")
        se <- se[-(1:(num.lv.c+num.RR))]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }
      
    }else if(num.lv>0&(num.lv.c+num.RR)>0){
      se.lambdas <- matrix(0,p,num.lv+(num.lv.c+num.RR));
      se.lambdas[,1:(num.lv.c+num.RR)][lower.tri(se.lambdas[,1:(num.lv.c+num.RR),drop=F], diag=FALSE)] <- se[1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR)))];
      se <- se[-c(1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR))))];
      se.lambdas[,((num.lv.c+num.RR)+1):ncol(se.lambdas)][lower.tri(se.lambdas[,((num.lv.c+num.RR)+1):ncol(se.lambdas),drop=F], diag=FALSE)] <- se[1:(p * num.lv - sum(0:num.lv))];
      se <- se[-c(1:(p * num.lv - sum(0:num.lv)))]
      colnames(se.lambdas) <- c(paste("CLV", 1:(num.lv.c+num.RR), sep=""),paste("LV", 1:num.lv, sep=""));
      rownames(se.lambdas) <- colnames(out$y)
      out$sd$theta <- se.lambdas;
      
      if(quadratic==TRUE){
        se.lambdas2 <- matrix(se[1:(p * ((num.lv.c+num.RR)+num.lv))], p, (num.lv.c+num.RR)+num.lv, byrow = T)  
        #re-order as the order in C++ is different, num_RR is last there
        se.lambdas2 <- se.lambdas2[,c(1:num.lv.c,(num.lv.c+num.lv+1):ncol(se.lambdas2),(num.lv.c+1):(num.lv+num.lv.c))]
        colnames(se.lambdas2) <- c(paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = ""),paste("LV", 1:num.lv, "^2", sep = ""))
        se <- se[-(1:((num.lv+(num.lv.c+num.RR))*p))]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }else if(quadratic=="LV"){
        se.lambdas2 <- matrix(se[1:((num.lv.c+num.RR)+num.lv)], p, (num.lv.c+num.RR)+num.lv, byrow = T)
        se.lambdas2 <- se.lambdas2[,c(1:num.lv.c,(num.lv.c+num.lv+1):ncol(se.lambdas2),(num.lv.c+1):(num.lv+num.lv.c))]
        colnames(se.lambdas2) <- c(paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = ""),paste("LV", 1:num.lv, "^2", sep = ""))
        se <- se[-(1:((num.lv.c+num.RR)+num.lv))]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }
      
    }
    if((num.lv+num.lv.c)>0){
      out$sd$sigma.lv  <- se.sigma.lv
      names(out$sd$sigma.lv) <- colnames(object$params$theta[,1:(num.lv+num.lv.c)])
    }
    
    out$sd$beta0 <- sebetaM[,1]; names(out$sd$beta0) <- colnames(object$y);
    if(!is.null(object$X)){
      out$sd$Xcoef <- matrix(sebetaM[,-1],nrow = nrow(sebetaM));
      rownames(out$sd$Xcoef) <- colnames(object$y); colnames(out$sd$Xcoef) <- colnames(object$X.design);
    }
    if(object$row.eff=="fixed") {out$sd$row.params <- se.row.params}
    
    if(family %in% c("negative.binomial")) {
      se.lphis <- se[1:length(unique(disp.group))];  out$sd$inv.phi <- se.lphis*object$params$inv.phi;
      out$sd$phi <- se.lphis*object$params$phi;
      if(length(unique(disp.group))==p){
        names(out$sd$phi) <- colnames(y);
      }else if(!is.null(names(disp.group))){
        try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
      }else{
        names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
      }
      names(out$sd$inv.phi) <-  names(out$sd$phi)
      se <- se[-(1:length(unique(disp.group)))]
    }
    if(family %in% c("gaussian","tweedie","gamma", "beta")) {
      se.lphis <- se[1:length(unique(disp.group))];
      out$sd$phi <- se.lphis*object$params$phi;
      if(length(unique(disp.group))==p){
        names(out$sd$phi) <- colnames(y);
      }else if(!is.null(names(disp.group))){
        try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
      }else{
        names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
      }
      
      se <- se[-(1:length(unique(disp.group)))]
    }
    if(family %in% c("ZIP")) {
      se.phis <- se[1:length(unique(disp.group))];
      out$sd$phi <- se.phis*exp(object$lp0)/(1+exp(object$lp0))^2;#
      if(length(unique(disp.group))==p){
        names(out$sd$phi) <- colnames(y);
      }else if(!is.null(names(disp.group))){
        try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
      }else{
        names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
      }
      names(out$sd$inv.phi) <-  names(out$sd$phi)
      se <- se[-(1:length(unique(disp.group)))]
    }
    if((object$randomB!=FALSE)&(num.lv.c+num.RR)>0){
      se.lsigmab.lv <-  se[1:length(object$paramssigmaLvXcoef)];
      se <-  se[-c(1:length(object$params$sigmaLvXcoef))]
      out$sd$sigmaLvXcoef <- se.lsigmab.lv*object$params$sigmaLvXcoef
      if(object$randomB=="LV")names(out$sd$sigmaLvXcoef) <- paste("CLV",1:(num.lv.c+num.RR), sep="")
      if(object$randomB=="P")names(out$sd$sigmaLvXcoef) <- colnames(lv.X)
      # if(object$randomB=="all")names(out$sd$sigmaLvXcoef) <- paste(paste("CLV",1:(num.lv.c+num.RR),sep=""),rep(colnames(lv.X),each=num.RR+num.lv.c),sep=".")
      if(object$randomB=="single")names(out$sd$sigmaLvXcoef) <- NULL
    }
    if(object$row.eff=="random") { 
      out$sd$sigma <- se[1:length(object$params$sigma)]*c(object$params$sigma[1],rep(1,length(object$params$sigma)-1)); 
      se=se[-(1:length(out$sd$sigma))] 
      names(out$sd$sigma) <- "sigma" 
      if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,3))) {out$sd$rho <- se[1]*(1-object$params$rho^2)^1.5; se = se[-1]}
      if((rstruc ==2 | (rstruc == 1)) & (cstrucn ==2)) {out$sd$rho <- se[1]*object$params$rho; se = se[-1]}
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
