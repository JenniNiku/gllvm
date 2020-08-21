########################################################################################
## GLLVM, with estimation done via Variational approximation using TMB-package
## Original author: Jenni Niku
##########################################################################################

gllvm.TMB <- function(y, X = NULL, formula = NULL, num.lv = 2, family = "poisson",
      method="VA",Lambda.struc="unstructured", row.eff = FALSE, reltol = 1e-6,
      seed = NULL,maxit = 1000, start.lvs = NULL, offset=NULL, sd.errors = TRUE,
      trace=FALSE,link="logit",n.init=1,restrict=30,start.params=NULL,
      optimizer="optim",starting.val="res",Power=1.5,diag.iter=1, dependent.row = FALSE,
      Lambda.start=c(0.1,0.5), jitter.var=0, zeta.struc = "species", Xrow = NULL) {
  if(!is.null(start.params)) starting.val <- "zero"
  ignore.u=FALSE
  n <- dim(y)[1]
  p <- dim(y)[2]
  objrFinal <- optrFinal <- NULL
  
  tr <- NULL
  num.lv <- num.lv
  y <- as.matrix(y)
  formula1 <- formula
  if(method=="VA"){ link="probit"}
  jitter.var.r <- 0
  if(length(jitter.var)>1){ 
    jitter.var.r <- jitter.var[2]
    jitter.var <- jitter.var[1]
  }
  
  if (!is.numeric(y))
    stop( "y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
  if ((family %in% c("tweedie", "ZIP")) && method == "VA")
    stop("family=\"", family, "\" : family not implemented with VA method, change the method to 'LA'")
  if (is.null(rownames(y)))
    rownames(y) <- paste("Row", 1:n, sep = "")
  if (is.null(colnames(y)))
    colnames(y) <- paste("Col", 1:p, sep = "")
  if(family == "ordinal") {
     y00<-y
    if(min(y)==0){ y=y+1}
    max.levels <- apply(y,2,function(x) length(min(x):max(x)))
    if(any(max.levels == 1)&zeta.struc=="species" || all(max.levels == 2)&zeta.struc=="species")
      stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks")
    if(any(!apply(y,2,function(x)all(diff(sort(unique(x)))==1)))&zeta.struc=="species")
      stop("Can't fit ordinal model if there are species with missing classes. Please reclassify per species or use zeta.struc = `common` ")
    if(!all(min(y)==apply(y,2,min))&zeta.struc=="species")
      stop("For ordinal data and zeta.struc=`species` all species must have the same minimum category. Please reclassify per species or use zeta.struc = `common`.")
    if(any(diff(sort(unique(c(y))))!=1)&zeta.struc=="common")
      stop("Can't fit ordinal model if there are missing classes. Please reclassify.")
  }
  if(dependent.row==TRUE & !is.null(Xrow) & !all(Xrow[,1]==(1:n)))stop("Cannot fit model with dependent structured row-effects yet.")

  num.X <- 0;
  if(!is.null(X)){

    if (!is.null(formula)) {
      xb <- as.matrix(model.matrix(formula, data = data.frame(X)))
      X <- as.matrix(xb[, !(colnames(xb) %in% c("(Intercept)"))])
      colnames(X) <- colnames(xb)[!(colnames(xb) %in% c("(Intercept)"))]
      Xd <- X1 <- X

      num.X <- dim(X)[2]
    } else {
      n1 <- colnames(X)
      formula = paste("~", n1[1], sep = "")
      if (length(n1) > 1) {
        for (i1 in 2:length(n1)) {
          formula <- paste(formula, n1[i1], sep = "+")
        }
      }
      formula <- formula(formula)
      xb <- as.matrix(model.matrix(formula, data = data.frame(X)))
      X <- as.matrix(xb[, !(colnames(xb) %in% c("(Intercept)"))])
      num.X <- dim(X)[2]
      colnames(X) <- colnames(xb)[!(colnames(xb) %in% c("(Intercept)"))]
      Xd <- X1 <- X

      nxd <- colnames(Xd)
      formulab <- paste("~", nxd[1], sep = "")

      if (length(nxd) > 1) {
        for (i in 2:length(nxd))
        formulab <- paste(formulab, nxd[i], sep = "+")
      }
      formula1 <- formulab
    }}
  if (is.null(formula) && is.null(X)) {
    formula = "~ 1"
  }

  ## Set initial values for model parameters (including dispersion prm) and latent variables
  if(!is.null(seed)) {
    set.seed(seed)
  }

  n.i <- 1

  out <- list( y = y, X = X, logL = Inf, X.design = X)
  if (n.init > 1)
    seed <- sample(1:10000, n.init)

  while(n.i <= n.init){
    if(n.init > 1 && trace)
      cat("Initial run ", n.i, "\n")

    fit <- start.values.gllvm.TMB(y = y, X = X, Xrow = Xrow, TR = NULL, family = family, offset= offset, num.lv = num.lv, start.lvs = start.lvs, seed = seed[n.i], starting.val = starting.val, power = Power, jitter.var = jitter.var, row.eff = row.eff, TMB=TRUE, link=link, zeta.struc = zeta.struc)
   
    sigma <- 1
    if (is.null(start.params)) {
      beta0 <- fit$params[, 1]
      betas <- NULL
      if (!is.null(X))
        betas <- c(fit$params[, 2:(num.X + 1)])
      lambdas <- NULL

      if (num.lv > 0) {
        lambdas <- as.matrix(fit$params[, (ncol(fit$params) - num.lv + 1):ncol(fit$params)])
        if (num.lv > 0)
          lambdas[upper.tri(lambdas)] <- 0
        covM.lvs <- array(NA, dim = c(n, num.lv, num.lv))
      }
      row.params <- NULL

      if (row.eff != FALSE) {
        row.params <- fit$row.params
        if (row.eff == "random") {
          sigma <- apply(row.params,2,sd)
        }
      }#rep(0,n)
      lvs <- NULL
      if (num.lv > 0)
        lvs <- matrix(fit$index, ncol = num.lv)

    } else{
      if (all(dim(start.params$y) == dim(y)) &&
          is.null(X) == is.null(start.params$X) &&
          (row.eff == start.params$row.eff)) {
        beta0 <- start.params$params$beta0 ## column intercepts
        betas <- NULL
        if (!is.null(X))
          if(!all((dim(X) == dim(start.params$X)))) stop( "Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that predictors X are the same in both models.")
          betas <- c(start.params$params$Xcoef) ## covariates coefficients
        lambdas <- NULL
        if (num.lv > 0){
          lambdas <- start.params$params$theta
          lambdas[upper.tri(lambdas)] <- 0}
        row.params <- NULL
        if (start.params$row.eff != FALSE) {
          row.params <- start.params$params$row.params
          if(row.eff=="fixed")
            row.params[1] <- 0#needs to be changed stiill
          if(row.eff=="random")
            sigma <- start.params$params$sigma
        }## row parameters
        lvs <- NULL
        if (num.lv > 0) {
          lvs <- matrix(start.params$lvs, ncol = num.lv)
          covM.lvs <- array(NA, dim = c(n, num.lv, num.lv))
        }## LVs
      } else {
        stop( "Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that attributes y, X and row.eff match to each other.")
      }
    }
    phis <- NULL
    if (family == "negative.binomial") {
      phis <- fit$phi
      if (any(phis > 25))
        phis[phis > 25] <- 25
      if (any(phis < 0.04))
        phis[phis < 0.04] <- 0.04
      fit$phi <- phis
      phis <- 1/phis
    }
    if (family == "tweedie") {
      phis <- fit$phi
      if (any(phis > 10))
        phis[phis > 10] <- 10
      if (any(phis < 0.10))
        phis[phis < 0.10] <- 0.10
      phis = (phis)
    }
    if (family == "ZIP") {
      phis <- (colMeans(y == 0) * 0.98) + 0.01
      phis <- phis / (1 - phis)
    } # ZIP probability
    if (family %in% c("gaussian", "gamma")) {
      phis <- fit$phi
    }
    if(family=="ordinal"){
      K = max(y00)-min(y00)
      if(zeta.struc=="species"){
        zeta <- c(t(fit$zeta[,-1]))
        zeta <- zeta[!is.na(zeta)]
      }else{
        zeta <- fit$zeta[-1]
      }

    }else{
      zeta = 0
    }

    if (is.null(offset))
      offset <- matrix(0, nrow = n, ncol = p)

    current.loglik <- -1e6; iter <- 1; err <- 10;
    ## LA-likelihood
    if(!is.null(row.params)){ r0 <- row.params} else {r0 <- rep(0,n)}
    a <- c(beta0)
    b <- NULL; if(!is.null(X)) b <- matrix(betas, ncol(X), p,byrow = TRUE)
    if(num.lv > 0) {
      # diag(lambdas) <- log(diag(lambdas)) #!!!
      lambda <- lambdas[lower.tri(lambdas,diag = TRUE)]
      u <- lvs
    }
    if(!is.null(phis)) { 
      phi <- phis 
    } else { 
        phi <- rep(1, p); 
        fit$phi <- phi
    }
    
    q <- num.lv

    xb<-Br<-matrix(0); sigmaB=diag(1);sigmaij=0; Abb=0
    randoml=c(0,0)
    if(row.eff=="fixed"){xr <- matrix(1,1,p)} else {xr <- matrix(0,1,p)}
    if(row.eff=="random"&dependent.row==FALSE)xr <- matrix(1,1,p)
    if(row.eff=="random") randoml[1]=1
    if(row.eff == "random" & dependent.row==TRUE){ nlvr<-num.lv+1 } else {nlvr=num.lv}
    if(!is.null(X)){Xd <- cbind(1,X)} else {Xd <- matrix(1,n)}
    extra <- 0
    
    optr<-NULL
    timeo<-NULL
    se <- NULL
    

    if(method=="VA" && (nlvr>0)){
      if(nlvr>0){
        if(is.null(start.params) || start.params$method!="VA"){
          if(Lambda.struc=="diagonal" || diag.iter>0){
            Au <- log(rep(Lambda.start[1],nlvr*n)) #1/2, 1
          } else{
            Au <- c(log(rep(Lambda.start[1],nlvr*n)),rep(0,nlvr*(nlvr-1)/2*n)) #1/2, 1
          }
        } else {
          Au <- NULL
          for(d in 1:nlvr) {
            if(start.params$Lambda.struc=="unstructured" || length(dim(start.params$A))==3){
              Au <- c(Au,log(start.params$A[,d,d]))
            } else {
              Au <- c(Au,log(start.params$A[,d]))
              }
          }
          if(Lambda.struc!="diagonal" && diag.iter==0){
            Au <- c(Au,rep(0,nlvr*(nlvr-1)/2*n))
          }
        }} else { Au <- 0}

      extra <- 0
      if(family == "poisson") { familyn <- 0}
      if(family == "negative.binomial") { familyn <- 1}
      if(family == "binomial") { familyn <- 2}
      if(family == "gaussian") {familyn=3}
      if(family == "gamma") {familyn=4}
      if(family == "ordinal") {familyn=7}
      if(family == "exponential") {familyn=8}

      if(!is.null(Xrow)){
        if(length(Lambda.start)<2){ Ar <- matrix(1,n,ncol(Xrow))}else{Ar <- matrix(Lambda.start[2],n,ncol(Xrow))}
      }else{
        if(length(Lambda.start)<2){ Ar <- matrix(1,n,1)}else{Ar <- matrix(Lambda.start[2],n,1)}
      }
      
      
      if(!is.null(Xrow)){
        nr<-ncol(Xrow)
        Xrow<-c(Xrow)
        #create unique sequence, but take a value not present in the dataset yet
        idx<-(max(Xrow)+1):((max(Xrow)+1)+length(unique(Xrow)))
        for(i in 1:length(unique(Xrow))){
          Xrow[Xrow==unique(Xrow)[i]]<-idx[i]
        }
        #Now make sure it starts at 1
        for(i in 1:length(unique(Xrow))){
          Xrow[Xrow==unique(Xrow)[i]]<-i
        }
        #and back to a matrix
        Xrow <- matrix(Xrow,ncol=nr)
        #collect row effects based on structure in Xrow
        map = list(lg_Ar = as.factor(c(as.matrix(Xrow))), r0  = as.factor(c(as.matrix(Xrow))))#this works with 1 row effect, not with multiple, then we need to make sure there are unique levels.
        
      }else{
        #empty if there is no Xrow
        map <- list()
      }
      
      if(row.eff=="random"){
        if(dependent.row) sigma<-c(log(sigma), rep(0, num.lv))
        if(num.lv>0){
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0),deprow = as.integer(!dependent.row)), silent=TRUE,
            parameters = list(r0 = as.matrix(r0), b = rbind(a,b), B = matrix(0), Br=Br,lambda = lambda, u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=sigma, Au=Au,Abb=0,zeta=zeta, lg_Ar = log(Ar)),
            map = map,
            inner.control=list(mgcmax = 1e+200,maxit = maxit),
            DLL = "gllvm")
        } else {
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0),deprow = as.integer(!dependent.row)), silent=TRUE,
            parameters = list(r0=matrix(r0), b = rbind(a,b),B=matrix(0), Br=Br,lambda = 0, u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=sigma,Au=Au,Abb=0, zeta=zeta, lg_Ar = log(Ar)),
            map = map,
            inner.control=list(mgcmax = 1e+200,maxit = maxit),
            DLL = "gllvm")##GLLVM
        }
      } else {
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0),deprow = as.integer(!dependent.row)), silent=TRUE,
          parameters = list(r0=as.matrix(r0), b = rbind(a,b),B=matrix(0), Br=Br,lambda = lambda, u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=0,Au=Au,Abb=0, zeta=zeta, lg_Ar = log(Ar)),
          map = map,
          inner.control=list(mgcmax = 1e+200,maxit = maxit),
          DLL = "gllvm")##GLLVM
      }
      if(optimizer=="nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol)),silent = TRUE))
      }
      if(optimizer=="optim") {
        timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
      }
      if(inherits(optr,"try-error")) warning(optr[1]);
      if(diag.iter>0 && Lambda.struc=="unstructured" && nlvr>1 && !inherits(optr,"try-error")){
        objr1 <- objr
        optr1 <- optr
        param1 <- optr$par
        nam <- names(param1)
        if(!is.null(Xrow)){
          r1 <- matrix(0,n ,ncol(Xrow))  
          lg_Ar1 <- matrix(0,n ,ncol(Xrow))
          for(r in 1:ncol(Xrow)){
            r1[,r] <- param1[nam=="r0"][1:length(unique(Xrow[,r]))][Xrow[,r]]
            if(row.eff=="random")lg_Ar1[,r] <- param1[nam=="lg_Ar"][1:length(unique(Xrow))][Xrow[,r]]
          }
        }else{
          r1 <- matrix(0,n ,1)
          lg_Ar1 <- matrix(0,n ,1)
          r1[,1] <- param1[nam=="r0"]
        }
 
        b1 <- matrix(param1[nam=="b"],num.X+1,p)
        lambda1 <- param1[nam=="lambda"]
        u1 <- matrix(param1[nam=="u"],n,num.lv)
        lg_phi1 <- param1[nam=="lg_phi"]
        sigmaB1 <- param1[nam=="sigmaB"]
        sigmaij1 <- param1[nam=="sigmaij"]
        log_sigma1 <- param1[nam=="log_sigma"]
        
        Au1<- c(pmax(param1[nam=="Au"],rep(log(1e-6), nlvr*n)), rep(0,nlvr*(nlvr-1)/2*n))
        zeta <- param1[nam=="zeta"]
        

        if(row.eff == "random"){
          if(num.lv>0){
            objr <- TMB::MakeADFun(
              data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0),deprow = as.integer(!dependent.row)), silent=TRUE,
              parameters = list(r0=r1, b = b1,B=matrix(0), Br=Br,lambda = lambda1, u = u1,lg_phi=lg_phi1,sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=log_sigma1,Au=Au1,Abb=0, zeta=zeta, lg_Ar = lg_Ar1), #log(phi)
              map = map,
              inner.control=list(mgcmax = 1e+200,maxit = maxit),
              DLL = "gllvm")
          } else {
            objr <- TMB::MakeADFun(
              data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0),deprow = as.integer(!dependent.row)), silent=TRUE,
              parameters = list(r0=r1, b = b1,B=matrix(0), Br=Br,lambda = lambda1, u = u1,lg_phi=lg_phi1,sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=log_sigma1,Au=Au1,Abb=0, zeta=zeta, lg_Ar = lg_Ar1), #log(phi)
              map = map,
              inner.control=list(mgcmax = 1e+200,maxit = maxit),
              DLL = "gllvm")
          }
        } else {
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0),deprow = as.integer(!dependent.row)), silent=TRUE,
            parameters = list(r0=r1, b = b1,B=matrix(0), Br=Br,lambda = lambda1, u = u1,lg_phi=lg_phi1,sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=0,Au=Au1,Abb=0, zeta=zeta, lg_Ar = lg_Ar1), #log(phi)
            map = map,
            inner.control=list(mgcmax = 1e+200,maxit = maxit),
            DLL = "gllvm")#GLLVM#
        }
        if(optimizer=="nlminb") {
          timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol)),silent = TRUE))
        }
        if(optimizer=="optim") {
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
        }
        if(optimizer=="nlminb"){
          if(inherits(optr, "try-error") || is.nan(optr$objective) || is.na(optr$objective)|| is.infinite(optr$objective) || optr$objective < 0){optr=optr1; objr=objr1; Lambda.struc="diagonal"}
        }else if(optimizer=="optim"){
          if(inherits(optr, "try-error") || is.nan(optr$value) || is.na(optr$value)|| is.infinite(optr$value) || optr$value < 0){optr=optr1; objr=objr1; Lambda.struc="diagonal"}
        }

      }

      param<-objr$env$last.par.best
      if(family %in% c("negative.binomial", "gaussian", "gamma")) {
        phis <- exp(param[names(param)=="lg_phi"])
      }
      if(family == "ordinal"){
        zetas <- param[names(param)=="zeta"]
        if(zeta.struc=="species"){
        zetanew <- matrix(NA,nrow=p,ncol=K)
        idx<-0
        for(j in 1:ncol(y)){
          k<-max(y[,j])-2
          if(k>0){
            for(l in 1:k){
              zetanew[j,l+1]<-zetas[idx+l]
            } 
          }
          idx<-idx+k
        }
        zetanew[,1] <- 0 
        row.names(zetanew) <- colnames(y00); colnames(zetanew) <- paste(min(y00):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
        }else{
          zetanew <- c(0,zetas)
          names(zetanew) <- paste(min(y00):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
        }
        
        zetas<-zetanew
        out$y<-y00
      }
      bi <- names(param)=="b"
      li <- names(param)=="lambda"
      ui <- names(param)=="u"
      if(num.lv > 0){
        lvs<-(matrix(param[ui],n,num.lv))
        theta <- matrix(0,p,num.lv)
        if(p>1) {
          theta[lower.tri(theta,diag=TRUE)] <- param[li];
        } else {theta <- param[li]}
        #diag(theta) <- exp(diag(theta)) # !!!
      }
      
      if(row.eff!=FALSE) {
        ri <- names(param)=="r0"
          if(!is.null(Xrow)){
            row.params <- matrix(0,n ,ncol(Xrow))
            for(r in 1:ncol(Xrow)){
              row.params[,r] <- param[ri][1:length(unique(Xrow))][Xrow[,r]]
            }
          }else{
            row.params <- matrix(0,n ,1)
            row.params[,1] <- param[ri]
          }

        if(!is.null(Xrow)){
          if(row.eff=="random"){
            sigma<-exp(param["log_sigma"])[1:ncol(Xrow)]
            if(nlvr>1 && dependent.row) sigma <- c(exp(param[names(param)=="log_sigma"])[1:ncol(Xrow)],(param[names(param)=="log_sigma"])[-c(1:ncol(Xrow))])
          }else if(row.eff=="random"){
            sigma <- exp(param["log_sigma"])[1:ncol(Xrow)]
          }
        }else{
          if(row.eff=="random"){
            sigma<-exp(param["log_sigma"])
            if(nlvr>1 && dependent.row) sigma <- c(exp(param[names(param)=="log_sigma"])[1],(param[names(param)=="log_sigma"])[-1])
          }else if(row.eff=="random"){
            sigma <- exp(param["log_sigma"])[1]
          }
        }
        }
        
      betaM <- matrix(param[bi],p,num.X+1,byrow=TRUE)
      beta0 <- betaM[,1]
      if(!is.null(X)) betas <- betaM[,-1]
      
      new.loglik <- objr$env$value.best[1]

    }

    if(method=="LA" || (nlvr==0 && method=="VA" && row.eff!="random")){
      if(!is.null(X)){Xd=cbind(1,X)} else {Xd=matrix(1,n)}
      extra=0
      if(family == "poisson") {familyn=0}
      if(family == "negative.binomial") {familyn=1}
      if(family == "binomial") {
        familyn=2;
        if(link=="probit") extra=1
      }
      if(family == "gaussian") {familyn=3}
      if(family == "gamma") {familyn=4}
      if(family == "tweedie"){ familyn=5; extra=Power}
      if(family == "ZIP"){ familyn=6;}
      if(family == "ordinal"){ familyn=7}
      if(family == "exponential"){ familyn=8}
      #collect row effects based on structure in Xrow
      if(!is.null(Xrow)){
        #Make sure Xrow[,1] starts at 1
        if(min(Xrow[,1])!=1){
          Xrow[,1] <- Xrow[,1] + 1
        }
        #make sure Xrow[,1] is sequential
        if(!all(diff(sort(unique(Xrow[,1])))==1)){
          Xrow.lev <- sort(unique(Xrow[,1]))
          for(i in 1:length(Xrow.lev)){
            Xrow[Xrow[,1]==Xrow.lev[i],1] <- i
          }
        }
        #Make sure Xrow[,>1] all have unique levels
        if(ncol(Xrow)>1){
          for(j in 2:ncol(Xrow)){
            Xrow.lev <- sort(unique(Xrow[,j]))
            for(i in 1:length(Xrow.lev)){
              Xrow[Xrow[,j]==Xrow.lev[i],j]<-max(Xrow[,j-1])+i
            }
          }
        }
        #collect row effects based on structure in Xrow
        map = list(lg_Ar = as.factor(c(as.matrix(Xrow))), r0  = as.factor(c(as.matrix(Xrow))))#this works with 1 row effect, not with multiple, then we need to make sure there are unique levels.
        
      }else{
        #empty if there is no Xrow
        map <- list()
      }
      
      if(row.eff=="random"){
        if(dependent.row) sigma<-c(log(sigma), rep(0, num.lv))
        if(num.lv>0){
         # if(dependent.row)u<-cbind(r0,u)
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0),deprow = as.integer(!dependent.row)), silent=!trace,
            parameters = list(r0=as.matrix(r0), b = rbind(a,b),B=matrix(0), Br=Br,lambda = lambda, u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=c(sigma), Au=0,Abb=0, zeta=zeta, lg_Ar = log(Ar)),
            map = map,
            inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
            random = c("u","r0"), DLL = "gllvm")
        }else{
          #if(dependent.row)u<-cbind(r0)
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0),deprow = as.integer(!dependent.row)), silent=!trace,
            parameters = list(r0=as.matrix(r0), b = rbind(a,b),B=matrix(0), Br=Br,lambda = 0, u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=c(sigma),Au=0,Abb=0, zeta=zeta, lg_Ar = log(Ar)),
            map = map,
            inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
            random = "r0", DLL = "gllvm")
        }
      } else {
        if(num.lv>0){
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0),deprow = as.integer(!dependent.row)), silent=!trace,
            parameters = list(r0=as.matrix(r0), b = rbind(a,b),B=matrix(0), Br=Br,lambda = lambda, u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=0,Au=0,Abb=0, zeta=zeta, lg_Ar = log(Ar)),
            map = map,
            inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
            random = c("u"), DLL = "gllvm")
        }else{
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0),deprow = as.integer(!dependent.row)), silent=!trace,
            parameters = list(r0=as.matrix(r0), b = rbind(a,b),B=matrix(0), Br=Br,lambda = 0, u = matrix(0),lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=0,Au=0,Abb=0, zeta=zeta, lg_Ar = log(Ar)),
            map = map,
            inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
            DLL = "gllvm")
        }
      }
      if(family=="ZIP" && FALSE) {
        m <- length(objr$par)
        low <- rep(-restrict,m); upp=rep(restrict,m);
        low[names(objr$par)=="lg_phi"]=0.0; upp[names(objr$par)=="lg_phi"]=1#0.99
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol),lower = low,upper = upp),silent = TRUE))
      }
      if(optimizer=="nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,maxit=maxit)),silent = TRUE))
      }
      if(optimizer=="optim") {
        timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
      }
      if(inherits(optr,"try-error")) warning(optr[1]);


      param <- objr$env$last.par.best
      bi <- names(param)=="b"
      li <- names(param)=="lambda"
      ui <- names(param)=="u"
      
      if(nlvr > 0){
        lvs<-(matrix(param[ui],n,num.lv))
        theta <- matrix(0,p,num.lv)
        if(p>1) {
          theta[lower.tri(theta,diag=TRUE)] <- param[li];
        } else { theta <- param[li] }
        # diag(theta) <- exp(diag(theta)) # !!!
      }
      if(row.eff!=FALSE) {
        ri <- names(param)=="r0"
      if(!is.null(Xrow)){
        row.params <- matrix(0,n ,ncol(Xrow))
        for(r in 1:ncol(Xrow)){
          row.params[,r] <- param[ri][Xrow[,r]]
        }
        if(row.eff=="random"){
          sigma<-exp(param["log_sigma"])[1:ncol(Xrow)]
          if(nlvr>1 && dependent.row) sigma <- c(exp(param[names(param)=="log_sigma"])[1:ncol(Xrow)],(param[names(param)=="log_sigma"])[-c(1:ncol(Xrow))])
        }
      }else{
        row.params <- lvs[,1]; lvs<- as.matrix(lvs[,-1])
        sigma<-exp(param["log_sigma"])[1]
        if(nlvr>1 && dependent.row) sigma <- c(exp(param[names(param)=="log_sigma"])[1],(param[names(param)=="log_sigma"])[-1])
      }
      }
      betaM <- matrix(param[bi],p,num.X+1,byrow=TRUE)
      beta0 <- betaM[,1]
      if(!is.null(X)) betas=betaM[,-1]
      new.loglik <- objr$env$value.best[1]
      
      if(family %in% c("negative.binomial", "tweedie", "ZIP", "gaussian", "gamma")) {
        phis <- exp(param[names(param)=="lg_phi"])
        if(family=="ZIP") {
          lp0 <- param[names(param)=="lg_phi"]; out$lp0 <- lp0
          phis <- exp(lp0)/(1+exp(lp0));
        }
      }
    }

    if(((n.i==1 || out$logL > (new.loglik))  && is.finite(new.loglik)) && !inherits(optr, "try-error")){
      out$start <- fit
      objrFinal<-objr1 <- objr; optrFinal<-optr1<-optr;
      out$logL <- new.loglik
      if(num.lv > 0) {
        out$lvs <- lvs
        out$params$theta <- theta
        rownames(out$lvs) <- rownames(out$y);
        if(num.lv>1) {colnames(out$params$theta) <- colnames(out$lvs) <- paste("LV", 1:num.lv, sep="");
        rownames(out$params$theta) <- colnames(out$y)}
      }
      names(beta0) <- colnames(out$y); out$params$beta0 <- beta0;
      if(!is.null(X)){betas <- matrix(betas,ncol=ncol(X)); out$params$Xcoef <- betas;
      rownames(out$params$Xcoef) <- colnames(out$y); colnames(out$params$Xcoef) <- colnames(X); }

      
      if(family =="negative.binomial") {
        out$params$inv.phi <- phis; names(out$params$inv.phi) <- colnames(out$y);
        out$params$phi <- 1/phis; names(out$params$phi) <- colnames(out$y);
      }
      if(family %in% c("gaussian", "tweedie", "gamma")) {
        out$params$phi <- phis; names(out$params$phi) <- colnames(out$y);
      }
      if(family =="ZIP") {
        out$params$phi <- phis; names(out$params$phi) <- colnames(out$y);
      }
      if(row.eff!=FALSE) {
        if(row.eff=="random"){ 
          out$params$sigma=sigma; 
          names(out$params$sigma)=paste("sigma",1:length(out$params$sigma))
          if(nlvr>1 && dependent.row) names(out$params$sigma) <- paste("sigma",c("",1:num.lv), sep = "")
        }
        if(is.null(Xrow)){out$params$row.params <- row.params; names(out$params$row.params) <- rownames(out$y)}
        if(!is.null(Xrow)){out$params$row.params <- row.params; row.names(out$params$row.params) <- rownames(out$y); colnames(row.params) <- colnames(Xrow)}
      }
      if(family == "binomial") out$link <- link;
      if(family == "tweedie") out$Power <- Power;
      if(family=="ordinal"){
        out$params$zeta <- zetas
      }
      out$row.eff <- row.eff
      out$time <- timeo
      pars <- optr$par

      if(method=="VA"){
        param <- objr$env$last.par.best
        if(nlvr>0){
          param <- objr$env$last.par.best
          Au <- param[names(param)=="Au"]
          A <- array(0, dim=c(n, nlvr, nlvr))
          for (d in 1:nlvr){
            for(i in 1:n){
              A[i,d,d] <- exp(Au[(d-1)*n+i]);
            }
          }
          if(length(Au) > nlvr*n){
            k <- 0;
            for (c1 in 1:nlvr){
              r <- c1 + 1;
              while (r <= nlvr){
                for(i in 1:n){
                  A[i,r,c1] <- Au[nlvr*n+k*n+i];
                  A[i,c1,r] <- A[i,r,c1];
                }
                k <- k+1; r <- r+1;
              }
            }
          }
          out$A <- A
          if(row.eff=="random" & !is.null(Xrow)){
            Ar <- matrix(0,n ,ncol(Xrow))
            for(r in 1:ncol(Xrow)){
              Ar[,r] <- exp(param[names(param)=="lg_Ar"])[Xrow[,r]]
            }            
            out$Ar <- Ar^2
          }
        }
        }
    }
    n.i <- n.i+1;
  }
  tr<-try({
    if(sd.errors && !is.infinite(out$logL)) {
      if(trace) cat("Calculating standard errors for parameters...\n")

      if(family=="ZIP") {
        p0i <- names(pars)=="lg_phi"
        p0 <- pars[p0i]
        p0 <- p0+runif(p,0,0.001)
        pars[p0i] <- p0
      }
      sdr <- optimHess(pars, objrFinal$fn, objrFinal$gr, control = list(reltol=reltol,maxit=maxit))#maxit=maxit
      m <- dim(sdr)[1]; incl <- rep(TRUE,m); incld <- rep(FALSE,m); inclr <- rep(FALSE,m)
      incl[names(objrFinal$par)=="B"] <- FALSE
      incl[names(objrFinal$par)%in%c("Br","sigmaB","sigmaij")] <- FALSE
      incl[names(objrFinal$par)=="Abb"]=FALSE;
      
      if(familyn!=7) incl[names(objrFinal$par)=="zeta"] <- FALSE
      
      if(method=="LA" || (num.lv==0 && method=="VA" && row.eff!="random")){
        incl[names(objrFinal$par)=="Au"] <- FALSE;
        if(row.eff=="random") {
          incl[names(objrFinal$par)=="r0"] <- FALSE; incld[names(objrFinal$par)=="r0"] <- FALSE
        } 
        if(row.eff=="fixed"){ incl[1] <- FALSE; incl[names(objrFinal$par)=="log_sigma"] <- FALSE}#still need to deal with identifiability constraints..
        if(row.eff==FALSE) {incl[names(objrFinal$par)=="r0"] <- FALSE; incl[names(objrFinal$par)=="log_sigma"] <- FALSE}
        if(familyn==0 || familyn==2 || familyn==7 || familyn==8) incl[names(objrFinal$par)=="lg_phi"] <- FALSE
        if(familyn==7) incl[names(objrFinal$par)=="zeta"] <- TRUE
        if(nlvr==0){
          incl[names(objrFinal$par)=="u"] <- FALSE;
          incl[names(objrFinal$par)=="lambda"] <- FALSE;
        }
        
        covM <- try(MASS::ginv(sdr[incl,incl]))
        se <- try(sqrt(diag(abs(covM))))
        if(nlvr>0){
          sd.random <- sdrandom(objrFinal, covM, incl,ignore.u = ignore.u)
          prediction.errors <- list()
          # if(row.eff=="random" && FALSE){
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
            if(row.eff=="random"){
              prediction.errors$row.params <- cov.lvs[,1,1]
              if(num.lv > 0) cov.lvs <- array(cov.lvs[,-1,-1], dim = c(n, num.lv, num.lv))
            }
            
            prediction.errors$lvs <- cov.lvs
            #sd.random <- sd.random[-(1:(n*num.lv))]
          }
          out$prediction.errors <- prediction.errors
        }
      } else {
        incl[names(objrFinal$par)=="Au"] <- FALSE;
        incl[names(objr$par) == "lg_Ar"] <- FALSE
        if(row.eff=="random") {
          inclr[names(objrFinal$par) == "r0"] <- FALSE;
          incl[names(objrFinal$par) == "r0"] <- FALSE; incld[names(objrFinal$par) == "r0"] <- FALSE
          if(row.eff=="random"&!is.null(Xrow)){
            incld[names(objr$par) == "lg_Ar"] <- TRUE
          }else{
            incld[names(objrFinal$par)=="Au"] <- TRUE  
          }
          
        }
        if(row.eff=="fixed") {incl[1] <- FALSE; incl[names(objrFinal$par)=="log_sigma"] <- FALSE}#still need to deal with identifiability constraints
        if(row.eff==FALSE) {incl[names(objrFinal$par)=="r0"] <- FALSE; incl[names(objrFinal$par)=="log_sigma"] <- FALSE}

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
      
        if(row.eff == "fixed") {
          if(is.null(Xrow)){
            se.row.params <- c(0,se[1:(n-1)]); names(se.row.params) <- rownames(out$y); se <- se[-(1:(n-1))]
          }else{
            se.row.params <- NULL
            for(r in 1:ncol(Xrow)){
              se.row.params <- c(se.row.params, se[1:length(unique(Xrow[,r]))][Xrow[,r]]); se[-(1:length(unique(Xrow[,r])))]
            }
            names(se.row.params) <- rownames(out$y);
          }
        }
      sebetaM <- matrix(se[1:((num.X+1)*p)],p,num.X+1,byrow=TRUE);  se <- se[-(1:((num.X+1)*p))]
      if(num.lv > 0) {
        se.lambdas <- matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag = TRUE)] <- se[1:(p * num.lv - sum(0:(num.lv-1)))];
        colnames(se.lambdas) <- paste("LV", 1:num.lv, sep="");
        rownames(se.lambdas) <- colnames(out$y)
        out$sd$theta <- se.lambdas; se <- se[-(1:(p * num.lv - sum(0:(num.lv-1))))];
        diag(out$sd$theta) <- diag(out$sd$theta)*diag(out$params$theta) 
      }
      out$sd$beta0 <- sebetaM[,1]; names(out$sd$beta0) <- colnames(out$y);
      if(!is.null(X)){
        out$sd$Xcoef <- matrix(sebetaM[,-1],nrow = nrow(sebetaM));
        rownames(out$sd$Xcoef) <- colnames(y); colnames(out$sd$Xcoef) <- colnames(X);
      }
      if(row.eff=="fixed") {out$sd$row.params <- se.row.params}

      if(family %in% c("negative.binomial")) {
        se.lphis <- se[1:p];  out$sd$inv.phi <- se.lphis*out$params$inv.phi;
        out$sd$phi <- se.lphis*out$params$phi;
        names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
      }
      if(family %in% c("tweedie", "gaussian", "gamma")) {
        se.lphis <- se[1:p];
        out$sd$phi <- se.lphis*out$params$phi;
        names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
      }
      if(family %in% c("ZIP")) {
        se.phis <- se[1:p];
        out$sd$phi <- se.phis*exp(lp0)/(1+exp(lp0))^2;#
        names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
      }
      if(row.eff=="random" && dependent.row == TRUE) {
        out$sd$sigma <- se[1:length(out$params$sigma)]*c(out$params$sigma[1],rep(1,length(out$params$sigma)-1)); 
        se=se[-(1:length(out$sd$sigma))] 
        names(out$sd$sigma) <- paste("sigma",1:length(out$sd$sigma),sep="") 
      }else if(row.eff=="random"){
          out$sd$sigma <- se[1:length(out$params$sigma)]*c(out$params$sigma[1],rep(1,length(out$params$sigma)-1)); 
          se=se[-(1:length(out$sd$sigma))] 
          names(out$sd$sigma) <- paste("sigma",1:length(out$sd$sigma)) 
        }
      if(family %in% c("ordinal")){
        se.zetanew <- se.zetas <- se;
        if(zeta.struc == "species"){
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
        row.names(out$sd$zeta) <- colnames(y00); colnames(out$sd$zeta) <- paste(min(y00):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
        
        }else{
          se.zetanew <- c(0, se.zetanew)
          out$sd$zeta <- se.zetanew
          names(out$sd$zeta) <- paste(min(y00):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
          
        }
         }

    }}, silent=T)

    if(inherits(tr, "try-error")) { cat("Standard errors for parameters could not be calculated, due to singular fit.\n") }

  if(is.null(formula1)){ out$formula <- formula} else {out$formula <- formula1}

  
  # DW, 7/5/19: adding TMBfn to output:
  out$TMBfn <- objrFinal
  out$TMBfn$par <- optrFinal$par #ensure params in this fn take final values
  out$logL <- -out$logL
  
  if(method == "VA"){
    if(row.eff == "random") out$logL = out$logL + n*0.5
    if(family=="gaussian") {
      out$logL <- out$logL - n*p*log(pi)/2
    }
  }

  return(out)
}

