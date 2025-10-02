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
#'@examples
#'data(eSpider)
#'mod <- gllvm(eSpider$abund, num.lv = 2, family = "poisson", sd.errors = FALSE)
#'# Calculate standard errors after fitting
#'sdErr <- se(mod)
#'# Store the standard errors in the right place
#'mod$sd <-sdErr$sd
#'# Store the Hessian in the right place
#'mod$Hess <- sdErr$Hess
#'@aliases se se.gllvm
#'@method se gllvm
#'@export
#'@export se.gllvm
se.gllvm <- function(object, ...){
  if(!is.finite(object$logL)) stop("Standard errors can not be calculated if log-likelihood value is not finite.")
  if(object$TMB == FALSE) stop("Function is not implemented for TMB = FALSE.")
  objrFinal <- object$TMBfn
  
  if(object$family =="betaH"){
    Y01 = (object$y>0)*1; colnames(Y01) = paste("H01",colnames(object$y), sep = "_")
    object$y = cbind(object$y, Y01)
  } 
  
  # backward compatibility
  
  if(!inherits(object$row.eff, "formula") && object$row.eff == "random") object$params$row.params.random <- object$params$row.params
  if(!inherits(object$row.eff, "formula") && object$row.eff == "fixed") object$params$row.params.fixed <- object$params$row.params[-1]
  if(!is.null(object$lv.X) && is.null(object$lv.X.design))object$lv.X.design <- object$lv.X
  if(is.null(object$col.eff$col.eff))object$col.eff$col.eff <- FALSE
  
  # end backward compatibility
  n <- nrow(object$y)
  p <- ncol(object$y)
  method <- object$method
  num.lv <- object$num.lv
  num.lv.c <- object$num.lv.c
  num.lv.cor <- object$num.lvcor
  num.RR <- object$num.RR
  lv.X.design <- object$lv.X.design
  
  quadratic <- object$quadratic
  nlvr <- num.lv + num.lv.c 
  dr = object$dr
  
  cstrucn = 0
  cstruc = object$corP$cstruc
  for (i in 1:length(cstruc)) {
    cstrucn[i] = switch(cstruc[i], "ustruc" = -1, "diag" = 0, "corAR1" = 1, "corExp" = 2, "corCS" = 3, "corMatern" = 4, "propto" = 5, 
                        "proptoustruc" = 6, "corAR1ustruc" = 7, "corExpustruc" = 8, "corCSustruc" = 9, "corMaternustruc" = 10)
  }
  cstruclvn = switch(object$corP$cstruclv, "ustruc" = 0, "diag" = 0, "corAR1" = 1, "corExp" = 2, "corCS" = 3, "corMatern" = 4)
  corWithinLv <- object$corP$corWithinLV
  
  family = object$family
  familyn <- objrFinal$env$data$family
  disp.group <- object$disp.group
  out <- list()
  
  # Trait model
  if (!is.null(object$TR)) {
    {
      if((method %in% c("VA", "EVA"))){
        sdr <- objrFinal$he(objrFinal$par)
      }
      if(method == "LA"){
        pars <- objrFinal$par
        sdr <- optimHess(pars, objrFinal$fn, objrFinal$gr)
      }
      
      m <- dim(sdr)[1]; incl <- rep(TRUE,m); incld <- rep(FALSE,m)
      incl[names(objrFinal$par)=="ePower"] <- FALSE
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
      # if(object$beta0com){ incl[names(objrFinal$par)=="b"] <- FALSE}
      if(familyn!=7 & familyn!=12) incl[names(objrFinal$par)=="zeta"] <- FALSE
      if(familyn==0 || familyn==2 || familyn==7 || familyn==8) incl[names(objrFinal$par)=="lg_phi"] <- FALSE
      if(!familyn %in% c(11, 14)) incl[names(objrFinal$par)=="lg_phiZINB"] <- FALSE
      
      if(num.lv>0) {
        incld[names(objrFinal$par)=="Au"] <- TRUE
        incld[names(objrFinal$par)=="u"] <- TRUE
      } else {
        incl[names(objrFinal$par)=="lambda2"] <- FALSE;
        incl[names(objrFinal$par)=="lambda"] <- FALSE;
      }
      
      if(!is.null(object$params$row.params.random)) {
        incld[names(objrFinal$par)=="lg_Ar"] <- TRUE
        incld[names(objrFinal$par)=="r0r"] <- TRUE
        if(ncol(object$TMBfn$env$data$csR)<2)incl[names(objrFinal$par) == "sigmaijr"] <- FALSE
        incl[names(objrFinal$par)=="r0r"] <- FALSE; 
      } else {
        incl[names(objrFinal$par)=="log_sigma"] <- FALSE
        incl[names(objrFinal$par)=="r0r"] <- FALSE
        incl[names(objrFinal$par) == "sigmaijr"] <- FALSE
      }
      if(!is.null(object$params$row.params.fixed)){
        if(object$params$row.params.fixed[1]==0) incl[names(objrFinal$par)=="r0f"][1] <- FALSE
      }else{
        incl[names(objrFinal$par)=="r0f"] <- FALSE
      }
      
      if(is.null(object$randomX)) {
        incl[names(objrFinal$par)%in%c("Br","sigmaB","sigmaij")] <- FALSE
      } else {
        xb <- objrFinal$env$data$xb
        incl[names(objrFinal$par)=="Abb"] <- FALSE; incld[names(objrFinal$par)=="Abb"] <- TRUE
        incl[names(objrFinal$par)=="Br"] <- FALSE; incld[names(objrFinal$par)=="Br"] <- TRUE
        if(NCOL(xb)==1) incl[names(objrFinal$par) == "sigmaij"] <- FALSE
      }
      
      if(method=="LA" || (num.lv==0 && (is.null(object$params$row.params.random) && is.null(object$randomX)) && object$col.eff$col.eff!="random")){
        covM <- try(MASS::ginv(sdr[incl,incl]))
        if(inherits(covM, "try-error")) { stop("Standard errors for parameters could not be calculated, due to singular fit.\n") }
        se <- try(sqrt(diag(abs(covM))))
        names(se) = names(object$TMBfn$par[incl])
        
        trpred<-try({
          if(num.lv > 0 || !is.null(object$params$row.params.random) || !is.null(object$randomX) || object$col.eff$col.eff == "random") {
            sd.random <- sdrandom(objrFinal, covM, incl)
            prediction.errors <- list()
            
            if(!is.null(object$params$row.params.random)){
              prediction.errors$row.params <- sd.random$row
            }
            if(object$col.eff$col.eff == "random"){
              prediction.errors$Br <- sd.random$Ab
            }
            if(!is.null(object$randomX)){
              prediction.errors$Br  <- sd.random$Ab
            }
            
            if(num.lv > 0){
              prediction.errors$lvs <- sd.random$A
            }
            out$prediction.errors <- prediction.errors
          }
        }, silent=TRUE)
        if(inherits(trpred, "try-error")) { cat("Prediction errors for random effects could not be calculated.\n") }
        
        out$Hess <- list(Hess.full=sdr, incl=incl, cov.mat.mod=covM)
      } else {
        sds <- sqrt(abs(diag(sdr)))
        if(any(sds<1e-12))sds[sds<1e-12]<-1
        
        sdr.s <- sweep(sweep(sdr,1,sds,"/"),2,sds,"/")
        
        A.mat <- sdr.s[incl,incl] # a x a
        D.mat <- as(sdr.s[incld,incld],"TsparseMatrix") # d x d
        B.mat <- sdr.s[incl,incld] # a x d
        cov.mat.mod<- try(MASS::ginv(A.mat-B.mat%*%as.matrix(solve(D.mat, t(B.mat)))),silent=T)
        if(inherits(cov.mat.mod,"try-error")){
          # block inversion via inverse of fixed-effects block
          Ai <- try(solve(A.mat),silent=T)
          cov.mat.mod <- try(Ai+Ai%*%B.mat%*%MASS::ginv(as.matrix(D.mat-t(B.mat)%*%Ai%*%B.mat))%*%t(B.mat)%*%Ai,silent=T)
        }
        suppressWarnings(try(cov.mat.mod <- sweep(sweep(cov.mat.mod, 2, sds[incl],"/"),1,sds[incl],"/"), silent = TRUE))
        
        if(inherits(cov.mat.mod, "try-error")) { stop("Standard errors for parameters could not be calculated, due to singular fit.\n") }
        se <- sqrt(diag(abs(cov.mat.mod)))
        names(se) = names(object$TMBfn$par[incl])
        
        incla<-rep(FALSE, length(incl))
        incla[names(objrFinal$par)=="u"] <- TRUE
        out$Hess <- list(Hess.full=sdr, incla = incla, incl=incl, incld=incld, cov.mat.mod=cov.mat.mod)
        
      }
      
      
      if(any(names(se)%in%names(object$TMBfn$env$map))){
        map <- object$TMBfn$env$map[names(object$TMBfn$env$map)%in%names(se)]
        # rebuild se vector if mapped parameters
        se.new <- NULL
        for(nm in unique(names(se))){
          if(!nm%in%names(map)){
            se.new <- c(se.new,se[names(se)==nm])
          }else{
            se.new <- c(se.new, se[names(se)==nm][map[[nm]]])
          }
          
        }
        se <- se.new
      }
      # reformat SEs based on the list that went into TMB
      se <- relist_gllvm(se, object$TMBfn$env$parList())
      
      if(!is.null(object$params$row.params.fixed)) {
        se.row.params <- se$r0f;
        if(object$params$row.params.fixed[1]==0) {
          se.row.params <- c(0, se.row.params)
        }
        names(se.row.params)  = names(object$params$row.params.fixed);
      }
      se.beta0 <- se$b[1,]
      
      # if(family %in% "betaH"){
      #   out$sd$betaH <- se$bH
      # }
      se.B <- se$B
      
      if(num.lv > 0) {
        se.sigma.lv <- se$sigmaLV[1:num.lv]
        # se.sigma.lv <- se.sigma.lv*object$params$sigma.lv[1:num.lv]
        se.lambdas <- matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag=FALSE)] <- se$lambda[1:(p * num.lv - sum(0:num.lv))];
        colnames(se.lambdas) <- paste("LV", 1:num.lv, sep="");
        rownames(se.lambdas) <- colnames(object$y)
        out$sd$theta <- se.lambdas; 
        
        if(quadratic==TRUE){
          se.lambdas2 <- matrix(se$lambda2, p, num.lv, byrow = T)  
          colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
          out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
        }else if(quadratic=="LV"){
          se.lambdas2 <- matrix(se$lambda2, p, num.lv, byrow = T)
          colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
          out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
        }
        out$sd$sigma.lv  <- se.sigma.lv
        names(out$sd$sigma.lv) <- colnames(object$params$theta[,1:num.lv])
        # if(num.lv>0 & family == "betaH"){
        #   se.thetaH <- matrix(0,num.lv,p); 
        #   se.thetaH[upper.tri(se.thetaH, diag=TRUE)] <- se$thetaH;
        #   out$sd$se.thetaH <- se.thetaH
        # }
        
      }
      # For correlated LVs:
      # if(num.lv.cor>0){
      #   diag(out$sd$theta) <- c(out$sd$sigma.lv)
      # }
      out$sd$beta0 <- se.beta0; 
      # if(!object$beta0com){ names(out$sd$beta0)  <-  colnames(object$y);}
      names(out$sd$beta0)  <-  colnames(object$y);
      out$sd$B <- c(se.B); 
      names(out$sd$B) <- names(object$params$B)
      if(!is.null(object$params$row.params.fixed)) {out$sd$row.params.fixed <- se.row.params}
      
      if(family %in% c("ZINB", "ZNIB")) {
        se.ZINB.lphis <- se$lg_phiZINB[disp.group];  
        if(family == "ZINB")out$sd$ZINB.inv.phi <- se.ZINB.lphis*object$params$ZINB.inv.phi;
        out$sd$ZINB.phi <- se.ZINB.lphis*object$params$ZINB.phi;
        names(out$sd$ZINB.phi) <- colnames(object$y);
        
        if(!is.null(names(disp.group))){
          try(names(out$sd$ZINB.phi) <- names(disp.group),silent=T)
        }
        if(family == "ZINB")names(out$sd$ZINB.inv.phi) <-  names(out$sd$ZINB.phi)
      }
      
      if(family %in% c("negative.binomial", "negative.binomial1")) {
        se.lphis <- se$lg_phi[disp.group];  out$sd$inv.phi <- se.lphis*object$params$inv.phi;
        out$sd$phi <- se.lphis*object$params$phi;
        if(length(unique(disp.group))==p){
          names(out$sd$phi) <- colnames(object$y);
        }else if(!is.null(names(disp.group))){
          try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
        }else{
          names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
        }
        names(out$sd$inv.phi) <-  names(out$sd$phi)
      }
      
      if(family %in% c("gaussian","tweedie","gamma", "beta", "betaH", "orderedBeta")) {
        se.lphis <- se$lg_phi[disp.group];
        out$sd$phi <- se.lphis*object$params$phi;
        if(length(unique(disp.group))==p){
          names(out$sd$phi) <- colnames(object$y);
        }else if(!is.null(names(disp.group))){
          try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
        }else{
          names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
        }
      }
      
      if(family%in%c("ZIP","ZINB","ZIB", "ZNIB")) {
        pars <- object$TMBfn$par
        p0i <- names(pars)=="lg_phi"
        p0 <- pars[p0i]
      }
      
      if(family %in% c("ZIP","ZINB","ZIB")) {
        se.phis <- se$lg_phi[disp.group];
        out$sd$phi <- se.phis*exp(p0)/(1+exp(p0))^2;#
        if(length(unique(disp.group))==p){
          names(out$sd$phi) <- colnames(object$y);
        }else if(!is.null(names(disp.group))){
          try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
        }else{
          names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
        }
      }else if(family == "ZNIB"){
        se.phis <- se$lg_phi[disp.group];  
        out$sd$phi <- se.phis*object$params$phi;
        names(out$sd$phi) <- colnames(object$y);
        
        if(length(unique(disp.group))==p){
          names(out$sd$phi) <- colnames(object$y);
        }else if(!is.null(names(disp.group))){
          try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
        }else{
          names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
        }
      }
      
      if(!is.null(object$randomX)){
        n.xb <- ncol(xb)
        if(length(se$sigmaB)>n.xb && !is.null(object$params$rho.sp)){
          out$sd$rho.sp <- tail(se$sigmaB,ifelse(object$col.eff$colMat.rho.struct == "single", 1, n.xb))
          se$sigmaB <- head(se$sigmaB, -ifelse(object$col.eff$colMat.rho.struct == "single", 1, n.xb))
        }
        out$sd$sigmaB <- se$sigmaB*c(sqrt(diag(object$params$sigmaB))); 
        names(out$sd$sigmaB) <- c(paste("sd",colnames(xb),sep = "."))
        if(n.xb>1){
          out$sd$corrpar <- se$sigmaij
        }
      }
      
      if(!is.null(object$params$row.params.random)) { 
        iter = 1 # keep track of index
        sigma <- se$log_sigma
        if(!is.null(object$TMBfn$env$map$log_sigma)) { #clean from duplicates and NAs
          sigma = sigma[!duplicated(object$TMBfn$env$map$log_sigma) & !is.na(object$TMBfn$env$map$log_sigma)]
        }
        trmsize = object$TMBfn$env$data$trmsize
        
        for(re in 1:length(cstrucn)){
          if(cstrucn[re] %in% c(0,-1, 5, 6:10)){
            # diag, ustruc, propto, proptoustruc
            sigma[iter:(iter+trmsize[1,re]-1)] <- sigma[iter:(iter+trmsize[1,re]-1)]*object$params$sigma[iter:(iter+trmsize[1,re]-1)]
            # parse labels
            form <- parse(text = colnames(trmsize)[re])[[1]]
            trm <- terms(as.formula(bquote(~ .(substitute(foo, list(foo=form))[[2]]))))
            LHS <- labels(trm)
            if(attr(trm, "intercept"))LHS <- c("(Intercept)", LHS)
            RHS <- form[[3]]
            
            names(sigma)[iter:(iter+trmsize[1,re]-1)] <- paste0(LHS, "|", RHS)
            
            iter <- iter + trmsize[1,re]
          }
          
          if(cstrucn[re] %in% c(1,3)) {
            sigma[iter] <- sigma[iter]*object$params$sigma[iter]
            names(sigma)[iter] <- colnames(trmsize)[re]
            names(sigma)[iter+1] <- paste0(colnames(trmsize)[re],".rho")
            sigma[iter+1] <- sigma[iter+1]*(1-object$params$sigma[iter+1]^2)^1.5
            iter <- iter +2
          } else if(cstrucn[re] %in% c(2)){
            sigma[iter:(iter+1)] <- sigma[iter:(iter+1)]*object$params$sigma[iter:(iter+1)]
            names(sigma)[iter] = paste0(colnames(trmsize)[re],".Scale")
            names(sigma)[iter+1] = colnames(trmsize)[re]
            iter <- iter + 2
          } else if(cstrucn[re] %in% c(4)){
            # sigma[iter:(iter+2)] <- sigma[iter:(iter+2)]*object$params$sigma[iter:(iter+2)] # matern smoothness fixed
            sigma[iter:(iter+1)] <- sigma[iter:(iter+1)]*object$params$sigma[iter:(iter+1)]
            names(sigma)[iter] = paste0(colnames(trmsize)[re],".Scale")
            names(sigma)[iter+1] = colnames(trmsize)[re]
            iter <- iter + 2
            # Matern smoothness
            # names(sigma)[iter+1] = "Matern kappa"
            # iter <- iter +1
          }
          
          if(cstrucn[re] %in% c(7,9)){
            sigma[iter] <- sigma[iter]*(1-object$params$sigma[iter]^2)^1.5
            names(sigma)[iter] <- paste0(colnames(trmsize)[re],".rho")
            iter <- iter +1
          }else if(cstrucn[re] == 8){
            sigma[iter] <- sigma[iter]*object$params$sigma[iter]
            names(sigma)[iter] =  paste0(colnames(trmsize)[re],".Scale")
            iter <- iter + 1
          }else if(cstrucn[re] == 10){
            sigma[iter] <- sigma[iter]*object$params$sigma[iter]
            names(sigma)[iter] =  paste0(colnames(trmsize)[re],".Scale")
            iter <- iter + 1
          }
        }
        out$sd$sigma <- sigma
      }
      
      if(num.lv.cor>0 & cstruclvn>0){
        if(length(object$params$rho.lv)>0){
          if(!is.null(object$TMBfn$env$map$rho_lvc)) { #clean from duplicates and NAs
            se$rho_lvc = se$rho_lvc[!duplicated(object$TMBfn$env$map$rho_lvc) & !is.na(object$TMBfn$env$map$rho_lvc)]
          }
          if((cstruclvn %in% c(1,3))) {out$sd$rho.lv <- se$rho_lvc[1:length(object$params$rho.lv)]*(1-object$params$rho.lv^2)^1.5}
          if((cstruclvn %in% c(2,4))) {out$sd$rho.lv <- se$rho_lvc[1:length(object$params$rho.lv)]*object$params$rho.lv}
          names(out$sd$rho.lv) <- names(object$params$rho.lv)
        }
        # if((cstrucn[2] %in% c(2,4))) {out$sd$scaledc <- se[(1:length(object$params$scaledc))]*object$params$scaledc; se = se[-(1:length(out$sd$scaledc))]}
      }
      
      if(family %in% c("ordinal")){
        y <- object$y
        K = max(y)-min(y)
        if(min(y)==0) y <- y+1 
        se.zetanew <- se.zetas <- se$zeta;
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
      if(family== "orderedBeta") {
        out$sd$zeta <- matrix(se$zeta,p,2)
        colnames(out$sd$zeta) = c("cutoff0","cutoff1")
      }
      
    }
  } else {
    #Without traits#
    pars <- objrFinal$par
    if(family %in% c("ZIP","ZINB","ZIB", "ZNIB")) {
      p0i <- names(pars)=="lg_phi"
      p0 <- pars[p0i]
      p0 <- p0+runif(p,0,0.001)
      pars[p0i] <- p0
    }
    if((method %in% c("VA", "EVA"))){
      sdr <- objrFinal$he(pars)
    }
    if(method == "LA"){
      sdr <- optimHess(pars, objrFinal$fn, objrFinal$gr)
    }
    # makes a small correction to the partial derivatives of LvXcoef if fixed-effect
    # because of constrained objective function
    # assumes L(x) = f(x) + lambda*c(x) for constraint function c(x)
    # though this is not (exactly) how we are fitting the model.
    if((object$num.RR+object$num.lv.c)>1 && isFALSE(object$randomB)){
      b_lvHE <- sdr[names(pars)=="b_lv",names(pars)=="b_lv"]
      Lmult <- lambda(pars,objrFinal) #estimates  lagranian multiplier
      sdr[names(pars)=="b_lv",names(pars)=="b_lv"] = b_lvHE + b_lvHEcorrect(Lmult,K = ncol(object$lv.X.design), d = object$num.lv.c+object$num.RR)
    }
    
    m <- dim(sdr)[1]; incl <- rep(TRUE,m); incld <- rep(FALSE,m); inclr <- rep(FALSE,m)
    incl[names(objrFinal$par)=="ePower"] <- FALSE
    # Not used for this model
    incl[names(objrFinal$par)%in%c("sigmaij")] <- FALSE
    
    # Variational params not included for incl
    incl[names(objrFinal$par)=="Ab_lv"] <- FALSE;
    incl[names(objrFinal$par)=="Abb"]=FALSE;
    incl[names(objrFinal$par)=="lg_Ar"] <- FALSE;
    incl[names(objrFinal$par)=="Br"] <- FALSE;
    incl[names(objrFinal$par)=="Au"] <- FALSE;
    incl[names(objrFinal$par)=="u"] <- FALSE;
    
    #slopes for reduced rank predictors
    if((num.lv.c+num.RR)==0){
      incl[names(objrFinal$par)=="sigmab_lv"] <- FALSE
      incl[names(objrFinal$par)=="b_lv"] <- FALSE
    }else if((num.lv.c+num.RR)>0&isFALSE(object$randomB)){
      incl[names(objrFinal$par)=="b_lv"] <- TRUE
      incl[names(objrFinal$par)=="sigmab_lv"] <- FALSE
    }else if((num.lv.c+num.RR>0)&!isFALSE(object$randomB)){
      incl[names(objrFinal$par)=="b_lv"] <- FALSE
      
      inclr[names(objrFinal$par)=="b_lv"] <- TRUE
      incld[names(objrFinal$par)=="b_lv"] <- TRUE
      
      incld[names(objrFinal$par)=="Ab_lv"] <- TRUE
      
      if(object$randomB!="iid")incl[names(objrFinal$par)=="sigmab_lv"] <- TRUE
    }
    
    #loadings for quadratic models
    if(quadratic == FALSE){incl[names(objrFinal$par)=="lambda2"]<-FALSE}
    if(familyn!=7 && familyn!=12) incl[names(objrFinal$par)=="zeta"] <- FALSE
    if(familyn==0 || familyn==2 || familyn==7 || familyn==8) incl[names(objrFinal$par)=="lg_phi"] <- FALSE
    if(!familyn %in% c(11, 14)) incl[names(objrFinal$par)=="lg_phiZINB"] <- FALSE
    
    if((num.lv+num.lv.c)>0){
      inclr[names(objrFinal$par)=="u"] <- TRUE;
      incld[names(objrFinal$par)=="u"] <- TRUE;
      incld[names(objrFinal$par)=="Au"] <- TRUE;
    } else {
      if(num.RR==0)incl[names(objrFinal$par)=="lambda"] <- FALSE;
      if(num.RR==0)incl[names(objrFinal$par)=="lambda2"] <- FALSE;
      incl[names(objrFinal$par)=="sigmaLV"] <- FALSE;
    }
    
    if(!is.null(object$params$row.params.random)) {
      incld[names(objrFinal$par) == "lg_Ar"] <- TRUE
      incld[names(objrFinal$par) == "r0r"] <- TRUE
      inclr[names(objrFinal$par) == "r0r"] <- TRUE;
      if(ncol(object$TMBfn$env$data$csR)<2)incl[names(objrFinal$par) == "sigmaijr"] <- FALSE
      incl[names(objrFinal$par) == "r0r"] <- FALSE; 
    } else {
      incl[names(objrFinal$par)=="log_sigma"] <- FALSE
      incl[names(objrFinal$par)=="r0r"] <- FALSE
      incl[names(objrFinal$par)=="sigmaijr"] <- FALSE
    }
    if(!is.null(object$params$row.params.fixed)){
      if(object$params$row.params.fixed[1]==0) incl[names(objrFinal$par)=="r0f"][1] <- FALSE
    }else{
      incl[names(objrFinal$par)=="r0f"] <- FALSE
    }
    
    if(object$col.eff$col.eff=="random") {
      incld[names(objrFinal$par)=="Abb"] <- TRUE
      incld[names(objrFinal$par)=="Br"] <- TRUE
      if(is.null(object$params$B))incl[names(objrFinal$par)=="B"] <- FALSE
    } else {
      incl[names(objrFinal$par)=="sigmaB"] <- FALSE
      incl[names(objrFinal$par)=="B"] <- FALSE
    }
    
    if(method=="LA" || ((num.lv+num.lv.c)==0 && (method %in% c("VA", "EVA")) && is.null(object$params$row.params.random) && isFALSE(object$randomB)) && object$col.eff$col.eff!="random"){
      covM <- try(MASS::ginv(sdr[incl,incl]))
      if(inherits(covM, "try-error")) { stop("Standard errors for parameters could not be calculated, due to singular fit.\n") }
      se <- try(sqrt(diag(abs(covM))))
      names(se) = names(object$TMBfn$par[incl])
      
      trpred<-try({
        if((num.lv+num.lv.c) > 0 || !is.null(object$params$row.params.random) || object$col.eff$col.eff == "random"){
          sd.random <- sdrandom(objrFinal, covM, incl, ignore.u = FALSE)
          prediction.errors <- list()
          
          if(!is.null(object$params$row.params.random)){
            prediction.errors$row.params <- sd.random$row
          }
          if(object$col.eff$col.eff == "random"){
            prediction.errors$Br <- sd.random$Ab
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
            if(!isFALSE(object$randomB)){
              prediction.errors$Ab.lv <- sd.random$Ab_lv
            }
          }
          out$prediction.errors <- prediction.errors
        }
      }, silent=TRUE)
      if(inherits(trpred, "try-error")) { cat("Prediction errors for random effects could not be calculated.\n") }
      
      out$Hess <- list(Hess.full=sdr, incl=incl, cov.mat.mod=covM)
      
    } else {
      # cnrm <- apply(sdr,2,function(x)sqrt(sum(x^2)))
      # rnrm <- apply(sdr,1,function(x)sqrt(sum(x^2)))
      # sdr.s <- sweep(sweep(sdr,2,cnrm,"/"),1,rnrm,"/")
      sds <- sqrt(abs(diag(sdr)))
      if(any(sds<1e-12))sds[sds<1e-12]<-1
      
      sdr.s <- sweep(sweep(sdr,1,sds,"/"),2,sds,"/")
      
      A.mat <- sdr.s[incl, incl] # a x a
      D.mat <- as(sdr.s[incld, incld], "TsparseMatrix") # d x d
      B.mat <- sdr.s[incl, incld] # a x d
      
      cov.mat.mod<- try(MASS::ginv(A.mat-B.mat%*%as.matrix(solve(D.mat, t(B.mat)))),silent=T)
      if(inherits(cov.mat.mod,"try-error")){
        # block inversion via inverse of fixed-effects block
        Ai <- try(solve(A.mat),silent=T)
        cov.mat.mod <- try(Ai+Ai%*%B.mat%*%MASS::ginv(as.matrix(D.mat-t(B.mat)%*%Ai%*%B.mat))%*%t(B.mat)%*%Ai,silent=T)
      }
      suppressWarnings(try(cov.mat.mod <- sweep(sweep(cov.mat.mod, 2, sds[incl],"/"),1,sds[incl],"/"), silent = TRUE))
      
      if(inherits(cov.mat.mod, "try-error")) { stop("Standard errors for parameters could not be calculated, due to singular fit.\n") }
      se <- sqrt(diag(abs(cov.mat.mod)))
      names(se) = names(object$TMBfn$par[incl])
      
      incla<-rep(FALSE, length(incl))
      incla[names(objrFinal$par)=="u"] <- TRUE
      out$Hess <- list(Hess.full=sdr, incla = incla, incl=incl, incld=incld, cov.mat.mod=cov.mat.mod)
      
    }
    
    
    if(any(names(se)%in%names(object$TMBfn$env$map))){
      map <- object$TMBfn$env$map[names(object$TMBfn$env$map)%in%names(se)]
      # rebuild se vector if mapped parameters
      se.new <- NULL
      for(nm in unique(names(se))){
        if(!nm%in%names(map)){
          se.new <- c(se.new,se[names(se)==nm])
        }else{
          se.new <- c(se.new, se[names(se)==nm][map[[nm]]])
        }
      }
      if(any(is.na(se.new)))se.new[is.na(se.new)] <- 0 # Can happen if a parameter is fixed to its startig value
      se <- se.new
    }
    # reformat SEs based on the list that went into TMB
    se <- relist_gllvm(se, object$TMBfn$env$parList())
    
    num.X <- 0; if(!is.null(object$X) || !isFALSE(object$col.eff$col.eff) && !is.null(object$X.design)) num.X <- dim(object$X.design)[2]
    if(!is.null(object$params$row.params.fixed)) {
      se.row.params <- se$r0f;
      if(object$params$row.params.fixed[1]==0) {
        se.row.params <- c(0, se.row.params)
      }
      names(se.row.params)  = names(object$params$row.params.fixed);
    }
    # if(family %in% "betaH"){
    #   # out$sd$betaH <- matrix(se$bH,p,num.X+1,byrow=TRUE);
    #   # rownames(out$sd$betaH) <- rownames(object$params$betaH);
    #   # colnames(out$sd$betaH) <- colnames(object$params$betaH)
    # }
    sebetaM <- matrix(se$b,p,num.X+1,byrow=TRUE);  
    
    
    
    if((num.lv.c+num.RR)>0&isFALSE(object$randomB)){
      se.LvXcoef <- matrix(se$b_lv,ncol=num.lv.c+num.RR,nrow=ncol(lv.X.design))
      colnames(se.LvXcoef) <- paste("CLV",1:(num.lv.c+num.RR),sep="")
      row.names(se.LvXcoef) <- colnames(lv.X.design)
      out$sd$LvXcoef <- se.LvXcoef
    }
    
    if((num.lv.c+num.lv)>0){ # exp transformation
      se.sigma.lv <- se$sigmaLV;
      # se.sigma.lv <- se$sigmaLV*object$params$sigma.lv;
    }
    
    if(num.lv > 0&(num.lv.c+num.RR)==0) {
      se.lambdas <- matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag=FALSE)] <- se$lambda[1:(p * num.lv - sum(0:num.lv))];
      colnames(se.lambdas) <- paste("LV", 1:num.lv, sep="");
      rownames(se.lambdas) <- colnames(object$y)
      out$sd$theta <- se.lambdas; 
      if(((p * num.lv - sum(0:num.lv))>0) && (num.RR+num.lv.c)>0) se$lambda <- se$lambda[-(1:(p * num.lv - sum(0:num.lv)))];
      
      if(quadratic==TRUE){
        se.lambdas2 <- matrix(se$lambda2[1:(p * num.lv)], p, num.lv, byrow = T)  
        colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
        se$lambda2 <- se$lambda2[-(1:(num.lv*p))]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }else if(quadratic=="LV"){
        se.lambdas2 <- matrix(se$lambda2[1:num.lv], p, num.lv, byrow = T)
        colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
        se$lambda2 <- se$lambda2[-(1:num.lv)]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }
    }else if(num.lv==0&(num.lv.c+num.RR)>0){
      se.lambdas <- matrix(0,p,(num.lv.c+num.RR)); se.lambdas[lower.tri(se.lambdas, diag=FALSE)] <- se$lambda[1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR)))];
      colnames(se.lambdas) <- paste("CLV", 1:(num.lv.c+num.RR), sep="");
      rownames(se.lambdas) <- colnames(object$y)
      out$sd$theta <- se.lambdas; se$lambda <- se$lambda[-(1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR))))];
      
      if(quadratic==TRUE){
        se.lambdas2 <- matrix(se$lambda2[1:(p * (num.lv.c+num.RR))], p, (num.lv.c+num.RR), byrow = T)  
        colnames(se.lambdas2) <- paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = "")
        se$lambda2 <- se$lambda2[-(1:((num.lv.c+num.RR)*p))]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }else if(quadratic=="LV"){
        se.lambdas2 <- matrix(se$lambda2[1:(num.lv.c+num.RR)], p, (num.lv.c+num.RR), byrow = T)
        colnames(se.lambdas2) <- paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = "")
        se$lambda2 <- se$lambda2[-(1:(num.lv.c+num.RR))]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }
      
    }else if(num.lv>0&(num.lv.c+num.RR)>0){
      se.lambdas <- matrix(0,p,num.lv+(num.lv.c+num.RR));
      se.lambdas[,1:(num.lv.c+num.RR)][lower.tri(se.lambdas[,1:(num.lv.c+num.RR),drop=F], diag=FALSE)] <- se$lambda[1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR)))];
      se$lambda <- se$lambda[-c(1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR))))];
      se.lambdas[,((num.lv.c+num.RR)+1):ncol(se.lambdas)][lower.tri(se.lambdas[,((num.lv.c+num.RR)+1):ncol(se.lambdas),drop=F], diag=FALSE)] <- se$lambda[1:(p * num.lv - sum(0:num.lv))];
      if(((p * num.lv - sum(0:num.lv))>0) && (num.lv.c + num.RR)>0) se$lambda <- se$lambda[-c(1:(p * num.lv - sum(0:num.lv)))]
      colnames(se.lambdas) <- c(paste("CLV", 1:(num.lv.c+num.RR), sep=""),paste("LV", 1:num.lv, sep=""));
      rownames(se.lambdas) <- colnames(object$y)
      out$sd$theta <- se.lambdas;
      
      if(quadratic==TRUE){
        se.lambdas2 <- matrix(se$lambda2[1:(p * ((num.lv.c+num.RR)+num.lv))], p, (num.lv.c+num.RR)+num.lv, byrow = T)  
        colnames(se.lambdas2) <- c(paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = ""),paste("LV", 1:num.lv, "^2", sep = ""))
        se$lambda2 <- se$lambda2[-(1:((num.lv+(num.lv.c+num.RR))*p))]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }else if(quadratic=="LV"){
        se.lambdas2 <- matrix(se$lambda2[1:((num.lv.c+num.RR)+num.lv)], p, (num.lv.c+num.RR)+num.lv, byrow = T)
        colnames(se.lambdas2) <- c(paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = ""),paste("LV", 1:num.lv, "^2", sep = ""))
        se$lambda2 <- se$lambda2[-(1:((num.lv.c+num.RR)+num.lv))]
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
      }
    }
    
    # if(num.lv>0 & family == "betaH"){
    #   se.thetaH <- matrix(0,num.lv,p); 
    #   se.thetaH[upper.tri(se.thetaH, diag=TRUE)] <- se$thetaH;
    #   out$sd$se.thetaH <- se.thetaH
    # }
    
    if(family %in% c("ZINB", "ZNIB")) {
      se.ZINB.lphis <- se$lg_phiZINB[disp.group];  
      if(family == "ZINB")out$sd$ZINB.inv.phi <- se.ZINB.lphis*object$params$ZINB.inv.phi;
      out$sd$ZINB.phi <- se.ZINB.lphis*object$params$ZINB.phi;
      names(out$sd$ZINB.phi) <- colnames(object$y);
      
      if(!is.null(names(disp.group))){
        try(names(out$sd$ZINB.phi) <- names(disp.group),silent=T)
      }
      if(family == "ZINB")names(out$sd$ZINB.inv.phi) <-  names(out$sd$ZINB.phi)
    }
    
    if((num.lv+num.lv.c)>0){
      out$sd$sigma.lv  <- se.sigma.lv
      names(out$sd$sigma.lv) <- colnames(object$params$theta[,1:(num.lv+num.lv.c)])
    }
    
    # For correlated LVs:
    # if(num.lv.cor>0){
    #   diag(out$sd$theta) <- c(out$sd$sigma.lv)
    # }
    out$sd$beta0 <- sebetaM[,1]; names(out$sd$beta0) <- colnames(object$y);
    if(!is.null(object$X)){
      out$sd$Xcoef <- matrix(sebetaM[,-1],nrow = nrow(sebetaM));
      rownames(out$sd$Xcoef) <- colnames(object$y); colnames(out$sd$Xcoef) <- colnames(object$X.design);
    }
    if(!is.null(object$params$row.params.fixed)) {out$sd$row.params.fixed <- se.row.params}
    
    if(family %in% c("negative.binomial", "negative.binomial1")) {
      se.lphis <- se$lg_phi[disp.group];  out$sd$inv.phi <- se.lphis*object$params$inv.phi;
      out$sd$phi <- se.lphis*object$params$phi;
      if(length(unique(disp.group))==p){
        names(out$sd$phi) <- colnames(object$y);
      }else if(!is.null(names(disp.group))){
        try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
      }else{
        names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
      }
      names(out$sd$inv.phi) <-  names(out$sd$phi)
    }
    
    if(family %in% c("gaussian","tweedie","gamma", "beta", "betaH", "orderedBeta")) {
      se.lphis <- se$lg_phi[disp.group];
      out$sd$phi <- se.lphis*object$params$phi;
      if(length(unique(disp.group))==p){
        names(out$sd$phi) <- colnames(object$y);
      }else if(!is.null(names(disp.group))){
        try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
      }else{
        names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
      }
    }
    
    if(family %in% c("ZIP","ZINB","ZIB")) {
      se.phis <- se$lg_phi[disp.group];
      out$sd$phi <- se.phis*exp(p0)/(1+exp(p0))^2;#
      if(length(unique(disp.group))==p){
        names(out$sd$phi) <- colnames(object$y);
      }else if(!is.null(names(disp.group))){
        try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
      }else{
        names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
      }
    }else if(family == "ZNIB"){
      se.phis <- se$lg_phi[disp.group];  
      out$sd$phi <- se.phis*object$params$phi;
      names(out$sd$phi) <- colnames(object$y);
      
      if(length(unique(disp.group))==p){
        names(out$sd$phi) <- colnames(object$y);
      }else if(!is.null(names(disp.group))){
        try(names(out$sd$phi) <- unique(names(disp.group)),silent=T)
      }else{
        names(out$sd$phi) <- paste("Spp. group", as.integer(unique(disp.group)))
      }
    }
    
    if(!isFALSE(object$randomB)&(num.lv.c+num.RR)>0){
      if(object$randomB!="iid"){
        if(object$randomB=="LV"){
        se.lsigmab.lv <-  head(se$sigmab_lv, num.lv.c+num.RR)
        out$sd$sigmaLvXcoef <- se.lsigmab.lv*object$params$sigmaLvXcoef
        }else if(object$randomB=="P"){
          se.lsigmab.lv <- head(se$sigmab_lv, ncol(object$lv.X.design)+num.lv.c+num.RR-1);
          out$sd$sigmaLvXcoef <- se.lsigmab.lv*object$params$sigmaLvXcoef
        }else{
          se.lsigmab.lv <- head(se$sigmab_lv, ncol(object$lv.X.design));
          out$sd$sigmaLvXcoef <- se.lsigmab.lv*object$params$sigmaLvXcoef
        }
      }
      if(object$randomB=="LV")names(out$sd$sigmaLvXcoef) <- paste("CLV",1:(num.lv.c+num.RR), sep="")
      if(object$randomB=="P")names(out$sd$sigmaLvXcoef) <- c(colnames(lv.X.design), paste("CLV",2:(num.lv.c+num.RR), sep=""))
      # if(object$randomB=="all")names(out$sd$sigmaLvXcoef) <- paste(paste("CLV",1:(num.lv.c+num.RR),sep=""),rep(colnames(lv.X.design),each=num.RR+num.lv.c),sep=".")
      if(object$randomB=="single")names(out$sd$sigmaLvXcoef) <- NULL
    }
    if(object$col.eff$col.eff=="random"){
      n.ce <- ncol(object$col.eff$spdr)
      sigma.sp <- 2*diag(object$params$sigmaB)*se$sigmaB[1:n.ce]
      covsigma.sp <- se$sigmaB[-c(1:n.ce)]
      if(!is.null(object$params$rho.sp)){
        out$sd$rho.sp <- tail(covsigma.sp,ifelse(object$col.eff$colMat.rho.struct == "single", 1, n.ce))
        covsigma.sp <- head(covsigma.sp, -ifelse(object$col.eff$colMat.rho.struct == "single", 1, n.ce))
      }
      sigma.sp <- diag(sigma.sp, length(sigma.sp))
      if(ncol(object$TMBfn$env$data$cs)==2){
        sigma.sp[object$TMBfn$env$data$cs] <- sigma.sp[object$TMBfn$env$data$cs[,c(2,1),drop=F]] <- covsigma.sp
      }
      out$sd$sigmaB <- sigma.sp
      if(!is.null(object$params[["B"]])){
      out$sd$B <- se$B
      names(out$sd$B) <- names(out$params$B)
      }
    }

      if(!is.null(object$params$row.params.random)) { 
        iter = 1 # keep track of index
        sigma <- se$log_sigma
        if(!is.null(object$TMBfn$env$map$log_sigma)) { #clean from duplicates and NAs
          sigma = sigma[!duplicated(object$TMBfn$env$map$log_sigma) & !is.na(object$TMBfn$env$map$log_sigma)]
        }
        trmsize = object$TMBfn$env$data$trmsize
        
        for(re in 1:length(cstrucn)){
          if(cstrucn[re] %in% c(0,-1, 5, 6:10)){
            # diag, ustruc, propto, proptoustruc
            sigma[iter:(iter+trmsize[1,re]-1)] <- sigma[iter:(iter+trmsize[1,re]-1)]*object$params$sigma[iter:(iter+trmsize[1,re]-1)]
            # parse labels
            form <- parse(text = colnames(trmsize)[re])[[1]]
            trm <- terms(as.formula(bquote(~ .(substitute(foo, list(foo=form))[[2]]))))
            LHS <- labels(trm)
            if(attr(trm, "intercept"))LHS <- c("(Intercept)", LHS)
            RHS <- form[[3]]
            
            names(sigma)[iter:(iter+trmsize[1,re]-1)] <- paste0(LHS, "|", RHS)
            
            iter <- iter + trmsize[1,re]
          }
          
          if(cstrucn[re] %in% c(1,3)) {
            sigma[iter] <- sigma[iter]*object$params$sigma[iter]
            names(sigma)[iter] <- colnames(trmsize)[re]
            names(sigma)[iter+1] <- paste0(colnames(trmsize)[re],".rho")
            sigma[iter+1] <- sigma[iter+1]*(1-object$params$sigma[iter+1]^2)^1.5
            iter <- iter +2
          } else if(cstrucn[re] %in% c(2)){
            sigma[iter:(iter+1)] <- sigma[iter:(iter+1)]*object$params$sigma[iter:(iter+1)]
            names(sigma)[iter] = paste0(colnames(trmsize)[re],".Scale")
            names(sigma)[iter+1] = colnames(trmsize)[re]
            iter <- iter + 2
          } else if(cstrucn[re] %in% c(4)){
            # sigma[iter:(iter+2)] <- sigma[iter:(iter+2)]*object$params$sigma[iter:(iter+2)] # matern smoothness fixed
            sigma[iter:(iter+1)] <- sigma[iter:(iter+1)]*object$params$sigma[iter:(iter+1)]
            names(sigma)[iter] = paste0(colnames(trmsize)[re],".Scale")
            names(sigma)[iter+1] = colnames(trmsize)[re]
            iter <- iter + 2
            # Matern smoothness
            # names(sigma)[iter+1] = "Matern kappa"
            # iter <- iter +1
          }
          
          if(cstrucn[re] %in% c(7,9)){
            sigma[iter] <- sigma[iter]*(1-object$params$sigma[iter]^2)^1.5
            names(sigma)[iter] <- paste0(colnames(trmsize)[re],".rho")
            iter <- iter +1
          }else if(cstrucn[re] == 8){
            sigma[iter] <- sigma[iter]*object$params$sigma[iter]
            names(sigma)[iter] =  paste0(colnames(trmsize)[re],".Scale")
            iter <- iter + 1
          }else if(cstrucn[re] == 10){
            sigma[iter] <- sigma[iter]*object$params$sigma[iter]
            names(sigma)[iter] =  paste0(colnames(trmsize)[re],".Scale")
            iter <- iter + 1
          }
        }
        out$sd$sigma <- sigma
      }
    
    if(num.lv.cor>0 & cstruclvn>0){ 
      if(length(object$params$rho.lv)>0){
        if(!is.null(object$TMBfn$env$map$rho_lvc)) { #clean from duplicates and NAs
          se$rho_lvc = se$rho_lvc[!duplicated(object$TMBfn$env$map$rho_lvc) & !is.na(object$TMBfn$env$map$rho_lvc)]
        }
        if((cstruclvn %in% c(1,3))) {out$sd$rho.lv <- se$rho_lvc*(1-object$params$rho.lv^2)^1.5}
        if((cstruclvn %in% c(2,4))) {out$sd$rho.lv <- se$rho_lvc*object$params$rho.lv}
        names(out$sd$rho.lv) <- names(object$params$rho.lv)
      }
    }
    
    if(family %in% c("ordinal")){
      y <- object$y
      K = max(y)-min(y)
      if(min(y)==0) y <- y+1 
      se.zetanew <- se.zetas <- se$zeta;
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
    if(family== "orderedBeta") {
      out$sd$zeta <- matrix(se$zeta,p,2)
      colnames(out$sd$zeta) = c("cutoff0","cutoff1")
    }
    
  }
  return(out)
}


#' @export se
se <- function(object, ...)
{
  UseMethod(generic = "se")
}
