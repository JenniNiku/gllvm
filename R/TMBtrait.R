########################################################################################
## GLLVM fourth corner model, with estimation done via Laplace and Variational approximation using TMB-package
## Original author: Jenni Niku
##########################################################################################



trait.TMB <- function(y, X = NULL,TR=NULL,formula=NULL, num.lv = 2, family = "poisson",
      Lambda.struc="unstructured", row.eff = FALSE, reltol = 1e-6, seed = NULL,
      maxit = 1000, start.lvs = NULL, offset=NULL, sd.errors = TRUE,trace=FALSE,
      link="logit",n.init=1,start.params=NULL,start0=FALSE,optimizer="optim",
      starting.val="res",method="VA",randomX=NULL,Power=1.5,diag.iter=1,
      Lambda.start=c(0.1, 0.5), jitter.var=0, yXT = NULL, zeta.struc = zeta.struc) {
  if(is.null(X) && !is.null(TR)) stop("Unable to fit a model that includes only trait covariates")

  objrFinal <- optrFinal <- NULL

    term <- NULL
  n <- dim(y)[1]; p <- dim(y)[2];
  y <- as.data.frame(y)
  formula1 <- formula
  if(method=="VA"){ link <- "probit"}
  
  if(NCOL(X) < 1) stop("No covariates in the model, fit the model using gllvm(y,family=",family,"...)")

  # change categorical variables to dummy variables
  num.X <- 0
  X.new <- NULL
  if(!is.null(X)) {
    num.X <- dim(X)[2]
    for (i in 1:num.X) {
      if(!is.factor(X[,i])) {
        if(length(unique(X[,i]))>2){ Xi <- scale(X[,i]) } else { Xi <- X[,i] }
        X[,i] <- Xi
        X.new <- cbind(X.new,Xi); if(!is.null(colnames(X)[i])) colnames(X.new)[dim(X.new)[2]] <- colnames(X)[i]
      } else {
        dum <- model.matrix( ~ X[,i])
        dum <- dum[, !(colnames(dum) %in% c("(Intercept)"))]
        colnames(dum) <- paste(colnames(X)[i], levels(X[,i])[ - 1], sep = "")
        X.new <- cbind(X.new, dum)
      }
    }
    X.new <- data.frame(X.new);
  }
  if(!is.null(randomX)){
    method <- "LA"
    xb <- as.matrix(model.matrix(randomX,data = X.new))
    xb <- as.matrix(xb[,!(colnames(xb) %in% c("(Intercept)"))])
    Br <- matrix(0, ncol(xb), p)
    sigmaB <- diag(ncol(xb))
  }

  num.T <- 0
  T.new <- NULL
  if(!is.null(TR)) {
    num.T <- dim(TR)[2]
    T.new <- matrix(0, p, 0)
    if(num.T > 0){
    for (i in 1 : num.T) {
      if(!is.factor(TR[,i])  && length(unique(TR[,i])) > 2) {
        TR[,i] <- scale(TR[,i])
        T.new <- cbind(T.new,scale(TR[,i])); colnames(T.new)[dim(T.new)[2]] <- colnames(TR)[i]
      } else {
        dum <- model.matrix(~TR[,i]-1)
        colnames(dum) <- paste(colnames(TR)[i],levels(TR[,i]),sep="")
        T.new <- cbind(T.new,dum)
      }
    }
    T.new <- data.matrix(T.new);
    }
  }


  if(is.null(formula)){
    n1 <- colnames(X)
    n2 <- colnames(TR)
    form1 <- paste("",n1[1],sep = "")
    if(length(n1)>1){
      for(i1 in 2:length(n1)){
        form1 <- paste(form1,n1[i1],sep = "+")
      }}
    formula <- paste("y~",form1,sep = "")
    formula <- paste(formula, form1,sep = " + (")
    
    formula <- paste(formula, ") : (", sep = "")
    formula <- paste(formula, n2[1], sep = "")
    if(length(n2) > 1){
      for(i2 in 2:length(n2)){
        formula <- paste(formula, n2[i2], sep = "+")
      }}
    formula1 <- paste(formula, ")", sep = "")
    formula <- formula(formula1)
  }

  if(!is.null(X) || !is.null(TR)){
    yX <- reshape(data.frame(cbind(y, X)), direction = "long", varying = colnames(y), v.names = "y")
    TR2 <- data.frame(time = 1:p, TR)
    if(is.null(yXT)){
      yXT <- merge(yX, TR2, by = "time")
    }
    data <- yXT

    m1 <- model.frame(formula, data = data)
    term <- terms(m1)

    Xd <- as.matrix(model.matrix(formula, data = data))
    nXd <- colnames(Xd)
    Xd <- as.matrix(Xd[, !(nXd %in% c("(Intercept)"))])
    colnames(Xd) <- nXd[!(nXd %in% c("(Intercept)"))]
    if(!is.null(X.new)) fx <- apply(matrix(sapply(colnames(X.new), function(x){grepl(x, colnames(Xd))}), ncol(Xd), ncol(X.new)), 2, any)
    ft <- NULL;
    if(NCOL(T.new) > 0) {
      ft <- apply(matrix(sapply(colnames(T.new), function(x){ grepl(x, colnames(Xd)) }), ncol(Xd), ncol(T.new)), 2, any)
    }
    X1 <- as.matrix(X.new[,fx]);
    TR1 <- as.matrix(T.new[,ft]);
    colnames(X1) <- colnames(X.new)[fx]; colnames(TR1)<-colnames(T.new)[ft];
    nxd <- colnames(Xd)
    formulab <- paste("~",nxd[1],sep = "");
    for(i in 2:length(nxd)) formulab <- paste(formulab,nxd[i],sep = "+")
    formula1 <- formulab
  }


  if(!(family %in% c("poisson","negative.binomial","binomial","tweedie","ZIP", "gaussian", "ordinal")))
    stop("Selected family not permitted...sorry!")
  if(!(Lambda.struc %in% c("unstructured","diagonal")))
    stop("Lambda matrix (covariance of vartiational distribution for latent variable) not permitted...sorry!")
  if(num.lv == 1) Lambda.struc <- "diagonal" ## Prevents it going to "unstructured" loops and causing chaos
  trial.size <- 1

  y <- as.matrix(y)
  if(!is.numeric(y)) stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
  if(family == "ordinal") {
    y00<-y
    if(min(y)==0){ y=y+1}
    max.levels <- apply(y,2,function(x) length(min(x):max(x)))
    if(any(max.levels == 1) || all(max.levels == 2))
      stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks")
    
    if(any(!apply(y,2,function(x)all(diff(sort(unique(x)))==1)))&zeta.struc=="species")
      stop("Can't fit ordinal model if there are species with missing classes. Please reclassify per species or use zeta.struc = `common` ")
    
    if(any(diff(sort(unique(c(y))))!=1)&zeta.struc=="common")
      stop("Can't fit ordinal model if there are missing classes. Please reclassify.")
  }
  if(is.null(rownames(y))) rownames(y) <- paste("Row",1:n,sep="")
  if(is.null(colnames(y))) colnames(y) <- paste("Col",1:p,sep="")
  if(!is.null(X)) { if(is.null(colnames(X))) colnames(X) <- paste("x",1:ncol(X),sep="") }

  out <-  list(y = y, X = X1, TR = TR1, num.lv = num.lv, row.eff = row.eff, logL = Inf, family = family, offset=offset,randomX=randomX,X.design=Xd,terms=term)
  if(is.null(formula) && is.null(X) && is.null(TR)){formula ="~ 1"}

  n.i <- 1;
  if(n.init > 1) seed <- sample(1:10000, n.init)
  while(n.i <= n.init){

    num.X <- dim(X)[2]
    num.T <- dim(TR)[2]

    if(n.init > 1 && trace) cat("initial run ",n.i,"\n");
    res <- start.values.gllvm.TMB(y = y, X = X1, TR = TR1, family = family, offset=offset, trial.size = trial.size, num.lv = num.lv, start.lvs = start.lvs, seed = seed[n.i],starting.val=starting.val,power=Power,formula = formula, jitter.var=jitter.var,yXT=yXT, row.eff = row.eff, TMB=TRUE, link=link, zeta.struc = zeta.struc)
    if(is.null(start.params)){
      beta0 <- res$params[,1]
      # common env params or different env response for each spp
      B <- NULL
      if(!is.null(TR) && !is.null(X)) {
        B <- c(res$B)[1:ncol(Xd)]
        if(any(is.na(B))) B[is.na(B)] <- 0
      }
      row.params <- NULL;
      if(row.eff!=FALSE){
        row.params <- res$row.params
        if (row.eff == "random") {
          sigma <- sd(row.params);
        }
      }
      vameans <- theta <- lambda <- NULL

      if(num.lv > 0) {
        vameans <- res$index
        theta <- as.matrix(res$params[,(ncol(res$params) - num.lv + 1):ncol(res$params)])#fts$coef$theta#
        theta[upper.tri(theta)] <- 0
        if(Lambda.struc == "unstructured") {
          lambda <- array(NA,dim=c(n,num.lv,num.lv))
          for(i in 1:n) { lambda[i,,] <- diag(rep(1,num.lv)) }
        }
        if(Lambda.struc == "diagonal") {
          lambda <- matrix(1,n,num.lv)
        }
        zero.cons <- which(theta == 0)
      }
    } else{
      if(dim(start.params$y)==dim(y) && is.null(X)==is.null(start.params$X) && is.null(T)==is.null(start.params$TR) && row.eff == start.params$row.eff){
        beta0 <- start.params$params$beta0
        # common env params or different env response for each spp
        B <- NULL
        if(!is.null(TR) && !is.null(X)) {
          B <- start.params$params$B;
        }
        fourth <- inter <- NULL; if(!is.null(TR) ) inter <- start.params$params$fourth   # let's treat this as a vector (vec(B'))'
        vameans <- theta <- lambda <- NULL

        if(num.lv > 0)
        if(row.eff) row.params <- start.params$params$row.params ## row parameters
        if(num.lv > 0) {
          theta <- c(start.params$params$theta) ## LV coefficients
          vameans <- matrix(start.params$lvs, ncol = num.lv);
        lambda <- start.params$A}
      } else { stop("Model which is set as starting parameters isn't the suitable you are trying to fit. Check that attributes y, X, TR and row.eff match to each other.");}
    }
    if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)

    phi <- phis <- NULL;
    if(family == "negative.binomial") {
      phis <- res$phi
      if (any(phis > 10))
        phis[phis > 100] <- 100
      if (any(phis < 0.01))
        phis[phis < 0.01] <- 0.01
      res$phi <- phis
      phis <- 1/phis

      }
    if(family == "tweedie") {
      phis <- res$phi; 
      if(any(phis>10)) phis[phis>10]=10; 
      if(any(phis<0.10))phis[phis<0.10]=0.10; 
      phis= (phis)
      }
    if (family == "ZIP") {
      phis <- (colMeans(y == 0) * 0.98) + 0.01; 
      phis <- phis / (1 - phis)
      } # ZIP probability
    # if (family == "gaussian") {
    #   phis <- fit$phi
    # }
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

    q <- num.lv

    if(!is.null(row.params)){ r0 <- row.params} else {r0 <- rep(0, n)}
    a <- c(beta0)
    if(num.lv > 0) {
      # diag(theta) <- log(diag(theta)) !!!
      theta <- theta[lower.tri(theta, diag = TRUE)]
      u <- vameans
    }
    if(!is.null(phis)) {phi=(phis)} else {phi <- rep(1,p)}
    q <- num.lv
    sigma <- 1

    if(num.lv > 0){
      if(is.null(start.params) || start.params$method!="VA"){
        if(Lambda.struc=="diagonal" || diag.iter>0){
          Au <- log(rep(Lambda.start[1],num.lv*n)) #
        } else{
          Au <- c(log(rep(Lambda.start[1],num.lv*n)),rep(0,num.lv*(num.lv-1)/2*n)) #1/2, 1
        }
      } else {
        Au=NULL
        for(d in 1:num.lv) {
          if(start.params$Lambda.struc=="unstructured" || length(dim(start.params$A))==3){ Au=c(Au,log(start.params$A[,d,d]))
          } else { Au <- c(Au,log(start.params$A[,d])) }
        }
        if(Lambda.struc!="diagonal" && diag.iter==0){
          Au <- c(Au,rep(0,num.lv*(num.lv-1)/2*n))
        }
      }} else { Au <- 0}
    if(length(Lambda.start)<2){ Ar <- rep(1,n)} else {Ar <- rep(Lambda.start[2],n)}


    optr<-NULL
    timeo<-NULL
    se <- NULL


    if(row.eff==FALSE){xr <- matrix(0,1,p)} else {xr <- matrix(1,1,p)}
    extra <- 0
    if(family == "poisson") { familyn=0}
    if(family == "negative.binomial") { familyn=1}
    if(family == "binomial") {
      familyn <- 2;
      if(link=="probit") extra <- 1
    }
    if(family == "gaussian") {familyn=3}
    if(family == "tweedie"){ familyn <- 4; extra <- Power}
    if(family == "ZIP"){ familyn <- 5;}
    if(family == "ordinal") {familyn=6}
    

    if(row.eff=="random"){# || !is.null(randomX)

      if(num.lv > 0){

        if(method=="VA"){
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=1,random=1, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B),lambda = theta, u = u,lg_phi=log(phi),log_sigma=log(sigma),Au=Au,lg_Ar=log(Ar), zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = maxit),
            DLL = "gllvm")
        } else {
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=1,random=1, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B),lambda = theta, u = u,lg_phi=log(phi),log_sigma=log(sigma),Au=0,lg_Ar=0, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = maxit,tol10=0.01),
            random = c("r0","u"), DLL = "gllvm")
        }
      } else {
        if(method=="VA"){
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=1,random=1, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
            parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B),lambda = 0, u = matrix(0),lg_phi=log(phi),log_sigma=log(sigma),Au=0,lg_Ar=log(Ar), zeta=zeta), #log(phi)
            inner.control=list(mgcmax = 1e+200,maxit = maxit,tol10=0.01),
            DLL = "gllvm")
        } else {
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=1,random=1, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B),lambda = 0, u = matrix(0),lg_phi=log(phi),log_sigma=log(sigma),Au=0,lg_Ar=0, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = maxit,tol10=0.01),
            random = c("r0"), DLL = "gllvm")
        }

      }


    } else {
      if(num.lv>0) {
        if(method=="VA"){
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=1,random=0, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B),lambda = theta, u = u,lg_phi=log(phi),log_sigma=0,Au=Au,lg_Ar=log(Ar), zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = 1000),
            DLL = "gllvm")
        } else {
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=1,random=0, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B),lambda = theta, u = u,lg_phi=log(phi),log_sigma=0,Au=0,lg_Ar=0, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
            random = c("u"), DLL = "gllvm")
        }
      } else {
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=1,random=0, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
          parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B),lambda = 0, u = matrix(0),lg_phi=log(phi),log_sigma=0,Au=0,lg_Ar=0, zeta=zeta),
          inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
          DLL = "gllvm")
      }

    }

    if(optimizer=="nlminb") {
      timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,maxit=maxit)),silent = TRUE))
    }
    if(optimizer=="optim") {
      timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
    }
    if(inherits(optr,"try-error")) warning(optr[1]);

    if(diag.iter>0 && Lambda.struc=="unstructured" && method =="VA" && (num.lv>0 || row.eff=="random") && is.null(randomX)){
      objr1 <- objr
      optr1 <- optr
      param1 <- optr$par
      nam <- names(param1)
      r1 <- matrix(param1[nam=="r0"])
      b1 <- rbind(param1[nam=="b"])
      B1 <- matrix(param1[nam=="B"])
      lambda1 <- param1[nam=="lambda"]
      u1 <- matrix(param1[nam=="u"],n,num.lv)
      lg_phi1 <- param1[nam=="lg_phi"]
      lg_sigma1 <- param1[nam=="log_sigma"]

      Au1<- c(pmax(param1[nam=="Au"],rep(log(0.001), num.lv*n)), rep(0,num.lv*(num.lv-1)/2*n))
      Ar1 <- c(param1[nam=="lg_Ar"])
      zeta <- param1[nam == "zeta"]

      if(row.eff=="random"){

        if(num.lv>0){
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=1,random=1, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=r1, b = b1,B=B1,lambda = lambda1, u = u1,lg_phi=lg_phi1,log_sigma=lg_sigma1,Au=Au1,lg_Ar=Ar1, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = 1000),
            DLL = "gllvm")
        } else {
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=1,random=1, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
            parameters = list(r0=r1, b = b1,B=B1,lambda = 0, u = matrix(0),lg_phi=lg_phi1,log_sigma=lg_sigma1,Au=0,lg_Ar=Ar1, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = maxit,tol10=0.01),
            DLL = "gllvm")

        }

      } else {
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=1,random=0, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
          parameters = list(r0=r1, b = b1,B=B1,lambda = lambda1, u = u1,lg_phi=lg_phi1,log_sigma=0,Au=Au1,lg_Ar=Ar1, zeta=zeta),
          inner.control=list(mgcmax = 1e+200,maxit = 1000),
          DLL = "gllvm")
      }

      if(optimizer=="nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol)),silent = TRUE))
      }
      if(optimizer=="optim") {
        timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
      }
      if(inherits(optr, "try-error")){optr <- optr1; objr <- objr1; Lambda.struc <- "diagonal"}

    }

    param <- objr$env$last.par.best
    if(family %in% c("negative.binomial", "tweedie", "gaussian")) {
      phis=exp(param[names(param)=="lg_phi"])
    }
    if(family=="ZIP") {
      lp0 <- param[names(param)=="lg_phi"]; out$lp0=lp0
      phis <- exp(lp0)/(1+exp(lp0));#log(phis); #
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
        row.names(zetanew) <- colnames(y00); colnames(zetanew) <- paste(min(y):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
      }else{
        zetanew <- c(0,zeta)
        names(zetanew) <- paste(min(y00):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
      }
      
      zetas<-zetanew
      out$y<-y00
    }
    bi<-names(param)=="b"
    Bi<-names(param)=="B"
    li<-names(param)=="lambda"
    ui<-names(param)=="u"
    if(row.eff!=FALSE) {
      ri <- names(param)=="r0"
      if(method=="LA" || row.eff=="random"){ row.params=param[ri] } else {row.params <- param[ri]}
      if(row.eff=="random") sigma<-exp(param["log_sigma"])
    }
    if(!is.null(randomX)){
      Bri <- names(param)=="Br"
      Br <- matrix(param[Bri],ncol(xb),p)
      Sri <- names(param)=="sigmaB"
      L <- diag(ncol(xb))
      if(ncol(xb)>1){
        sigmaB <- diag(exp(param[Sri]))
        Srij <- names(param)=="sigmaij"
        Sr <- param[Srij]
        L[upper.tri(L)] <- Sr
        D <- diag(diag(t(L)%*%L))
      } else{
        D <- 1
        sigmaB <- (exp(param[Sri]))
      }
      sigmaB_ <- solve(sqrt(D))%*%(t(L)%*%L)%*%solve(sqrt(D))
      sigmaB <- sigmaB%*%sigmaB_%*%t(sigmaB)

    }
    beta0 <- param[bi]
    B <- param[Bi]
    if(num.lv > 0){
      lvs <- (matrix(param[ui],n,q))
      theta <- matrix(0,p,num.lv)
      if(p>1) {
        theta[lower.tri(theta,diag=TRUE)] <- param[li];
      } else {theta <- param[li]}
      # diag(theta) <- exp(diag(theta)) !!!
    }
    new.loglik<-objr$env$value.best[1]


    if((n.i==1 || out$logL > (new.loglik)) && is.finite(new.loglik) && !inherits(optr, "try-error")){
      objrFinal<-objr1 <- objr; optrFinal<-optr1 <- optr;
      out$logL <- new.loglik
      if(num.lv > 0) {
        out$lvs <- lvs
        out$params$theta <- theta
        rownames(out$lvs) <- rownames(out$y);
        colnames(out$params$theta) <- colnames(out$lvs) <- paste("LV", 1:num.lv, sep="");
        rownames(out$params$theta) <- colnames(out$y)
      }
      names(beta0) <- colnames(out$y); out$params$beta0 <- beta0;
      out$params$B <- B; names(out$params$B)=colnames(Xd)

      if(row.eff!=FALSE) {
        if(row.eff=="random"){ out$params$sigma=sigma; names(out$params$sigma)="sigma"}
        out$params$row.params <- row.params; names(out$params$row.params) <- rownames(out$y)
      }
      if(family %in% c("negative.binomial")) {
        out$params$phi <- 1/phis; names(out$params$phi) <- colnames(out$y);
        out$params$inv.phi <- phis; names(out$params$inv.phi) <- colnames(out$y);
      }
      if(family %in% c("gaussian","tweedie")) {
        out$params$phi <- phis; names(out$params$phi) <- colnames(out$y);
      }
      if(family =="ZIP") {
        out$params$phi <- phis; names(out$params$phi) <- colnames(out$y);
      }
      if (family == "ordinal") {
        out$params$zeta <- zetas
      }
      if(!is.null(randomX)){
        out$params$Br <- Br
        out$params$sigmaB <- sigmaB
        out$corr <- sigmaB_
      }
      if(family == "binomial") out$link <- link;
      out$row.eff <- row.eff
      out$time <- timeo
      out$start<-res
      out$Power <- Power
      pars <- optr$par

      if(method=="VA" && num.lv>0){
        param <- objr$env$last.par.best
        Au <- param[names(param)=="Au"]
        A <- array(0,dim=c(n,num.lv,num.lv))
        for (d in 1:num.lv){
          for(i in 1:n){
            A[i,d,d] <- exp(Au[(d-1)*n+i]);
          }
        }
        if(length(Au)>num.lv*n){
          k <- 0;
          for (c1 in 1:num.lv){
            r <- c1 + 1;
            while (r <= num.lv){
              for(i in 1:n){
                A[i,r,c1] <- Au[num.lv*n+k*n+i];
                A[i,c1,r] <- A[i,r,c1];
              }
              k <- k+1; r <- r+1;
            }
          }
        }
        out$A <- A
      }

      if(method=="VA" && row.eff=="random"){
        Ar <- exp(param[names(param)=="lg_Ar"])
        out$Ar <- (Ar)^2
      }
    }



    tr<-try({
      if(sd.errors && !is.infinite(out$logL)) {
        if(trace) cat("Calculating standard errors for parameters...\n")

        sdr <- optimHess(optrFinal$par, objrFinal$fn, objrFinal$gr,control = list(reltol=reltol,maxit=maxit))

        m <- dim(sdr)[1]; incl <- rep(TRUE,m); incld <- rep(FALSE,m)
        incl[names(objrFinal$par)=="lg_Ar"] <- FALSE;
        if(row.eff!="random") {incl[names(objrFinal$par)=="log_sigma"] <- FALSE}
        if(familyn!=6) incl[names(objrFinal$par)=="zeta"] <- FALSE
        
        incl[names(objrFinal$par)=="Au"] <- FALSE; if(num.lv>0) incld[names(objrFinal$par)=="Au"] <- TRUE

        if(row.eff=="random") {
          incl[names(objrFinal$par)=="lg_Ar"] <- FALSE; incld[names(objrFinal$par)=="lg_Ar"] <- TRUE
          incl[names(objrFinal$par)=="r0"] <- FALSE; incld[names(objrFinal$par)=="r0"] <- TRUE
        }
        if(row.eff==FALSE) incl[names(objrFinal$par)=="r0"] <- FALSE
        if(row.eff=="fixed") incl[1] <- FALSE
        incl[names(objrFinal$par)=="u"] <- FALSE; incld[names(objrFinal$par)=="u"] <- TRUE
        if(familyn==0 || familyn==2 || familyn==6) incl[names(objrFinal$par)=="lg_phi"] <- FALSE
        if(familyn==6) incl[names(objrFinal$par)=="zeta"] <- TRUE
        if(num.lv==0){
          incl[names(objrFinal$par)=="u"] <- FALSE;
          incl[names(objrFinal$par)=="lambda"] <- FALSE;
          incl[names(objrFinal$par)=="Au"] <- FALSE;
        }

        if(method=="LA" || (num.lv==0 && row.eff!="random")){
          incl[names(objrFinal$par)=="Au"] <- FALSE;

          covM <- try(MASS::ginv(sdr[incl,incl]))
          se <- try(sqrt(diag(abs(covM))))
          if(num.lv>0 || row.eff=="random" || !is.null(randomX)) {
            sd.random <- sdrandom(objrFinal, covM, incl)
            prediction.errors <- list()
            if(row.eff=="random"){
              prediction.errors$row.params <- diag(as.matrix(sd.random))[1:n];
              sd.random <- sd.random[-(1:n),-(1:n)]
            }
            if(num.lv>0){
              cov.lvs <- array(0, dim = c(n, num.lv, num.lv))
              for (i in 1:n) {
                cov.lvs[i,,] <- as.matrix(sd.random[(0:(num.lv-1)*n+i),(0:(num.lv-1)*n+i)])
              }
              prediction.errors$lvs <- cov.lvs
            }
            out$prediction.errors <- prediction.errors
          }
        } else {
          A.mat <- -sdr[incl,incl] # a x a
          D.mat <- -sdr[incld,incld] # d x d
          B.mat <- -sdr[incl,incld] # a x d
          cov.mat.mod <- try(MASS::ginv(A.mat-B.mat%*%solve(D.mat)%*%t(B.mat)))
          se <- sqrt(diag(abs(cov.mat.mod)))
          
          incla<-rep(FALSE, length(incl))
          incla[names(objrFinal$par)=="u"] <- TRUE
          out$Hess <- list(Hess.full=sdr, incla = incla, incl=incl, incld=incld, cov.mat.mod=cov.mat.mod)
          
        }

        if(row.eff=="fixed") { se.row.params <- c(0,se[1:(n-1)]); names(se.row.params)  = rownames(out$y); se <- se[-(1:(n-1))] }
        se.beta0 <- se[1:p]; se <- se[-(1:p)];
        se.B <- se[1:length(B)]; se <- se[-(1:length(B))];
        if(num.lv>0) {
          se.theta <- matrix(0,p,num.lv); se.theta[lower.tri(se.theta, diag = TRUE)]<-se[1:(p * num.lv - sum(0:(num.lv-1)))];
          colnames(se.theta) <- paste("LV", 1:num.lv, sep="");
          rownames(se.theta) <- colnames(out$y)
          out$sd$theta <- se.theta; se <- se[-(1:(p * num.lv - sum(0:(num.lv-1))))];
          # diag(out$sd$theta) <- diag(out$sd$theta)*diag(out$params$theta) !!!
        }
        out$sd$beta0 <- se.beta0; names(out$sd$beta0)  <-  colnames(out$y);
        out$sd$B <- se.B; names(out$sd$B) <- colnames(Xd)
        if(row.eff=="fixed") {out$sd$row.params <- se.row.params}

        if(family %in% c("negative.binomial")) {
          se.lphis <- se[1:p];  out$sd$inv.phi <- se.lphis*out$params$inv.phi;
          out$sd$phi <- se.lphis*out$params$phi;
          names(out$sd$inv.phi) <- names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
        }
        if(family %in% c("gaussian","tweedie")) {
          se.lphis <- se[1:p];
          out$sd$phi <- se.lphis*out$params$phi;
          names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
        }
        if(family %in% c("ZIP")) {
          se.phis <- se[1:p];
          out$sd$phi <- se.phis*exp(lp0)/(1+exp(lp0))^2;#
          names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
        }
        if(row.eff=="random") { out$sd$sigma <- se[1]*out$params$sigma; names(out$sd$sigma) <- "sigma"; se=se[-1] }
        if(!is.null(randomX)){
          nr <- ncol(xb)
          out$sd$sigmaB <- se*c(diag(out$params$sigmaB), rep(1,nr*(nr-1)/2 ) )
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
        
      }})

    n.i <- n.i+1;
  }
  if(inherits(tr, "try-error")) { cat("Standard errors for parameters could not be calculated.\n") }

  if(is.null(formula1)){ out$formula <- formula} else {out$formula <- formula1}

  out$D <- Xd
  out$TMBfn <- objrFinal
  out$TMBfn$par <- optrFinal$par #ensure params in this fn take final values
  out$logL <- -out$logL
  
  if(method == "VA"){
    #if(num.lv > 0) out$logL = out$logL + n*0.5*num.lv
    if(row.eff == "random") out$logL = out$logL + n*0.5
    #if(!is.null(randomX)) out$logL = out$logL + p*0.5*ncol(xb)
    if(family=="gaussian") {
      out$logL <- out$logL - n*p*log(pi)/2
    }
  }

  return(out)
}

