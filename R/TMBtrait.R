########################################################################################
## GLLVM fourth corner model, with estimation done via Laplace and Variational approximation using TMB-package
## Original author: Jenni Niku
##########################################################################################



trait.TMB <- function(
      y, X = NULL,TR=NULL,formula=NULL, num.lv = 2, family = "poisson",
      Lambda.struc = "unstructured", Ab.struct = "unstructured", row.eff = FALSE, reltol = 1e-6, seed = NULL,
      maxit = 1000, start.lvs = NULL, offset=NULL, sd.errors = TRUE,trace=FALSE,
      link="logit",n.init=1,start.params=NULL,start0=FALSE,optimizer="optim",
      starting.val="res",method="VA",randomX=NULL,Power=1.5,diag.iter=1, Ab.diag.iter = 0, dependent.row = TRUE,
      Lambda.start=c(0.1, 0.5), jitter.var=0, yXT = NULL, scale.X = FALSE, randomX.start = "zero", beta0com = FALSE
      ,zeta.struc = "species") {
  if(is.null(X) && !is.null(TR)) stop("Unable to fit a model that includes only trait covariates")
  if(!is.null(start.params)) starting.val <- "zero"
  
  objrFinal <- optrFinal <- NULL

    term <- NULL
  n <- dim(y)[1]; p <- dim(y)[2];
  y <- as.data.frame(y)
  formula1 <- formula
  beta0com0 = beta0com
  if(method=="VA"){ link <- "probit"}
  jitter.var.r <- 0
  if(length(jitter.var)>1){ 
    jitter.var.r <- jitter.var[2]
    jitter.var <- jitter.var[1]
  }
  
  if(NCOL(X) < 1) stop("No covariates in the model, fit the model using gllvm(y,family=",family,"...)")

  # change categorical variables to dummy variables
  num.X <- 0
  X.new <- NULL
  if(!is.null(X)) {
    num.X <- dim(X)[2]
    for (i in 1:num.X) {
      if(!is.factor(X[,i])) {
        if(length(unique(X[,i]))>2){ Xi <- scale(X[,i], scale = scale.X, center = scale.X) } else { Xi <- X[,i] }
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
  
  randomXb <- NULL
  
  if(!is.null(randomX)){
    #
    if(num.lv>0 && randomX.start == "res" && starting.val == "res") {randomXb <- randomX}
    #
    xb <- as.matrix(model.matrix(randomX, data = data.frame(X)))
    rnam <- colnames(xb)[!(colnames(xb) %in% c("(Intercept)"))]
    xb <- as.matrix(xb[, rnam]); #as.matrix(X.new[, rnam])
    if(NCOL(xb) == 1) colnames(xb) <- rnam
    bstart <- start.values.randomX(y, xb, family, starting.val = starting.val, power = Power)
    Br <- bstart$Br
    sigmaB <- bstart$sigmaB
    sigmaij <- rep(0,(ncol(xb)-1)*ncol(xb)/2)

    # method <- "LA"
    # xb <- as.matrix(model.matrix(randomX,data = X.new))
    # xb <- as.matrix(xb[,!(colnames(xb) %in% c("(Intercept)"))])
    # Br <- matrix(0, ncol(xb), p)
    # sigmaB <- diag(ncol(xb))
  } else {
    xb <- Br <- matrix(0); sigmaB <- diag(1); sigmaij <- 0; Abb <- 0 
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
        T.new <- cbind(T.new,scale(TR[,i], scale = scale.X, center = scale.X)); colnames(T.new)[dim(T.new)[2]] <- colnames(TR)[i]
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
    phi<-phis <- NULL
    sigma <- 1
    phi <- phis <- NULL;
    
    if(n.init > 1 && trace) cat("initial run ",n.i,"\n");

    res <- start.values.gllvm.TMB(y = y, X = X1, TR = TR1, family = family, offset=offset, trial.size = trial.size, num.lv = num.lv, start.lvs = start.lvs, seed = seed[n.i],starting.val=starting.val,power=Power,formula = formula, jitter.var=0,yXT=yXT, row.eff = row.eff, TMB=TRUE, link=link, randomX=randomXb, beta0com = beta0com0, zeta.struc = zeta.struc)

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
        if(!is.null(randomXb)){
          Br <- res$Br
          sigmaB <- (res$sigmaB)
          if(length(sigmaB)>1) sigmaij <- rep(0,length(res$sigmaij))
        }
        
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
        
        if(n.init > 1 && !is.null(res$mu) && starting.val == "res" && family != "tweedie") {
          if(family=="ZIP") {
            lastart <- FAstart(res$mu, family="poisson", y=y, num.lv = num.lv, jitter.var = jitter.var[1])
          } else {
            lastart <- FAstart(res$mu, family=family, y=y, num.lv = num.lv, phis = res$phi, jitter.var = jitter.var[1])
          }
          theta <- lastart$gamma#/lastart$gamma
          vameans<-lastart$index#/max(lastart$index)
        }
      }
      
    } else{
      if(all(dim(start.params$y)==dim(y)) && is.null(X)==is.null(start.params$X) && is.null(T)==is.null(start.params$TR) && row.eff == start.params$row.eff){
        beta0 <- start.params$params$beta0
        # common env params or different env response for each spp
        B <- NULL
        if(!is.null(TR) && !is.null(X)) {
          B <- start.params$params$B;
        }
        fourth <- inter <- NULL; if(!is.null(TR) ) inter <- start.params$params$fourth   # let's treat this as a vector (vec(B'))'
        vameans <- theta <- lambda <- NULL

        row.params <- NULL
        if(row.eff %in% c("fixed","random",TRUE)) {
          if(row.eff == start.params$row.eff){ 
            res$row.params <- row.params <- start.params$params$row.params
            if(row.eff %in% c("random")) res$sigma <- sigma <- start.params$params$sigma
          } else {
            row.params <- res$row.params
          }
          
        }     
        if(num.lv > 0) {
          theta <- (start.params$params$theta) ## LV coefficients
          vameans <- matrix(start.params$lvs, ncol = num.lv);
          lambda <- start.params$A
        }
        if(family == "negative.binomial" && start.params$family == "negative.binomial" && !is.null(start.params$params$phi)) {res$phi<-start.params$params$phi}
      } else { stop("Model which is set as starting parameters isn't the suitable you are trying to fit. Check that attributes y, X, TR and row.eff match to each other.");}
    }
    if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)

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
    if(jitter.var.r>0){ 
      if(row.eff == "random") row.params <- row.params + rnorm(n, 0, sd = sqrt(jitter.var.r));
      if(!is.null(randomX)) Br <- Br + t(mvtnorm::rmvnorm(p, rep(0, nrow(Br)),diag(nrow(Br))*jitter.var.r));
    }
    
    q <- num.lv
    
    a <- c(beta0)
    if(num.lv > 0) {
      # diag(theta) <- log(diag(theta)) # !!!
      theta <- theta[lower.tri(theta, diag = TRUE)]
      u <- vameans
    }
    if(!is.null(phis)) {phi=(phis)} else {phi <- rep(1,p)}
    q <- num.lv

    if(!is.null(row.params)){ r0 <- row.params} else {r0 <- rep(0, n)}
    if(row.eff == "random"){ nlvr<-num.lv+1 } else {nlvr=num.lv}
    if(row.eff=="fixed"){xr <- matrix(1,1,p)} else {xr <- matrix(0,1,p)}

    # set starting values for variational distribution covariances
    if(nlvr > 0){
      if(Lambda.struc=="diagonal" || diag.iter>0){
        Au <- log(rep(Lambda.start[1],nlvr*n)) #
      } else{
        Au <- c(log(rep(Lambda.start[1],nlvr*n)),rep(0,nlvr*(nlvr-1)/2*n))
      }
    } else { Au <- 0}
    if(length(Lambda.start)<2){ Ar <- rep(1,n)} else {Ar <- rep(Lambda.start[2],n)}

    if(!is.null(randomX)){ 
      if(length(Lambda.start)>2) { 
        a.var <- Lambda.start[3];
      } else {a.var <- 1;}
      if(Ab.struct == "diagonal" || Ab.diag.iter>0){
        Abb <- c(log(rep(a.var, ncol(xb) * p)))
      } else {
        Abb <- c(log(rep(a.var, ncol(xb) * p)), rep(0, ncol(xb) * (ncol(xb) - 1) / 2 * p))
      } 
    } else { Abb <- 0 }

    optr<-NULL
    timeo<-NULL
    se <- NULL


    randoms=c("u","Br")
    randoml=c(0,0)
    
    # family settings
    extra <- c(0,1)
    if(family == "poisson") { familyn=0}
    if(family == "negative.binomial") { familyn=1}
    if(family == "binomial") {
      familyn <- 2;
      if(link=="probit") extra[1] <- 1
    }
    if(family == "gaussian") {familyn=3}
    if(family == "tweedie"){ familyn <- 4; extra[1] <- Power}
    if(family == "ZIP"){ familyn <- 5;}
    if(family == "ordinal") {familyn=6}
    if(beta0com){
      extra[2] <- 0
      Xd<-cbind(1,Xd)
      a <- a*0
      B<-c(mean(a),B)
    }

    if(row.eff=="random" || !is.null(randomX)){#

      if(!is.null(randomX)){
        randoml[2]=1
        res$Br <- Br
        res$sigmaB <- sigmaB
      }
      if(row.eff=="random"){
        randoml[1] <- 1
        if(dependent.row) sigma<-c(sigma[1], rep(0, num.lv))
        if(num.lv>0){
          u<-cbind(r0,u)
        }else {
          u<-cbind(r0)
        }
      }
      if(row.eff!="random" && num.lv==0){
        u=matrix(0)
      }
      if(num.lv > 0){

        if(method == "VA"){
          objr <- TMB::MakeADFun(

            data = list(y = y, x = Xd, xr=xr, xb=xb, offset=offset, num_lv = num.lv, family=familyn, extra=extra, method=0, model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a), B=matrix(B), Br=Br, lambda = theta, u = u, lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij, log_sigma=c(sigma), Au=Au, Abb=Abb, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = maxit),
            DLL = "gllvm")
        } else {
          objr <- TMB::MakeADFun(

            data = list(y = y, x = Xd,xr=xr, xb=xb, offset=offset, num_lv = num.lv, family=familyn,extra=extra,method=1,model=1,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B), Br=Br,lambda = theta, u = u, lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij, log_sigma=c(sigma), Au=0, Abb=0, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = maxit,tol10=0.01),
            random = c(("Br")[randoml[2]>0],"u"), DLL = "gllvm")
        }
      } else {
        if(method=="VA"){
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd, xr=xr, xb=xb, offset=offset, num_lv = num.lv, family=familyn, extra=extra, method=0, model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a), B=matrix(B), Br=Br,lambda = 0, u = u, lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij, log_sigma=c(sigma), Au=Au, Abb=Abb, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = maxit,tol10=0.01),
            DLL = "gllvm")
        } else {
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd, xr=xr, xb=xb, offset=offset, num_lv = num.lv, family=familyn, extra=extra, method=1, model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a), B=matrix(B), Br=Br, lambda = 0, u = u, lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij, log_sigma=c(sigma), Au=0, Abb=0, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = maxit,tol10=0.01),
            random = c(randoms[randoml>0]), DLL = "gllvm")
        }

      }


    } else {
      if(num.lv>0) {
        if(method=="VA"){
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd, xr=xr, xb=xb, offset=offset, num_lv = num.lv, family=familyn,extra=extra,method=0,model=1,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B), Br=Br,lambda = theta, u = u,lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij,log_sigma=0,Au=Au, Abb=0, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = 1000),
            DLL = "gllvm")
        } else {
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd, xr=xr, xb=xb, offset=offset, num_lv = num.lv, family=familyn,extra=extra,method=1,model=1,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B), Br=Br,lambda = theta, u = u,lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij,log_sigma=0,Au=0, Abb=0, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
            random = c("u"), DLL = "gllvm")
        }
      } else {
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd, xr=xr, xb=xb, offset=offset, num_lv = num.lv, family=familyn,extra=extra,method=1,model=1,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
          parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B), Br=Br,lambda = 0, u = matrix(0),lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij,log_sigma=0,Au=0, Abb=0, zeta=zeta),
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

    if(diag.iter>0 && Lambda.struc=="unstructured" && method =="VA" && (nlvr>0 || !is.null(randomX)) && !inherits(optr,"try-error")){
      objr1 <- objr
      optr1 <- optr
      param1 <- optr$par
      nam <- names(param1)
      r1 <- matrix(param1[nam=="r0"])
      b1 <- rbind(param1[nam=="b"])
      B1 <- matrix(param1[nam=="B"])
      Br1 <- matrix(param1[nam=="Br"], ncol(xb), p)
      lambda1 <- param1[nam=="lambda"]
      u1 <- matrix(param1[nam=="u"],n,nlvr)
      lg_phi1 <- param1[nam=="lg_phi"]
      sigmaB1 <- param1[nam=="sigmaB"]
      sigmaij1 <- param1[nam=="sigmaij"]
      lg_sigma1 <- param1[nam=="log_sigma"]
      
      Au1<- c(pmax(param1[nam=="Au"],rep(log(1e-4), nlvr*n)), rep(0,nlvr*(nlvr-1)/2*n))
      Abb1 <- param1[nam=="Abb"]
      if(Ab.diag.iter>0 && Ab.struct == "unstructured")
        Abb1 <- c(Abb1 + 0.1, rep(0,ncol(xb)*(ncol(xb)-1)/2*p))
      zeta <- param1[nam == "zeta"]

      if(row.eff=="random" || !is.null(randomX)){

        if(nlvr>0){
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd, xr=xr, xb=xb, offset=offset, num_lv = num.lv, family=familyn, extra=extra, method=0, model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
            parameters = list(r0=r1, b = b1, B=B1, Br=Br1, lambda = lambda1, u = u1, lg_phi=lg_phi1, sigmaB=sigmaB1, sigmaij=sigmaij1, log_sigma=lg_sigma1, Au=Au1, Abb=Abb1, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = 1000),
            DLL = "gllvm")
        } else {
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=1,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
            parameters = list(r0=r1, b = b1,B=B1,Br=Br1,lambda = 0, u = matrix(0),lg_phi=lg_phi1,sigmaB=sigmaB1,sigmaij=sigmaij1,log_sigma=lg_sigma1,Au=0, Abb=Abb1, zeta=zeta),
            inner.control=list(mgcmax = 1e+200,maxit = maxit,tol10=0.01),
            DLL = "gllvm")

        }

      } else {
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=1,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0)), silent=!trace,
          parameters = list(r0=r1, b = b1,B=B1,Br=Br1,lambda = lambda1, u = u1,lg_phi=lg_phi1, sigmaB=sigmaB1, sigmaij=sigmaij1, log_sigma=0, Au=Au1, Abb=Abb1, zeta=zeta),
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
        zetanew <- c(0,zetas)
        names(zetanew) <- paste(min(y00):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
      }
      
      zetas<-zetanew
      out$y<-y00
    }
    bi<-names(param)=="b"
    Bi<-names(param)=="B"
    li<-names(param)=="lambda"
    ui<-names(param)=="u"
    if(nlvr > 0){
      lvs <- (matrix(param[ui],n,nlvr))
      theta <- matrix(0,p,num.lv)
      if(p>1) {
        theta[lower.tri(theta,diag=TRUE)] <- param[li];
      } else {theta <- param[li]}
      # diag(theta) <- exp(diag(theta))#!!!
    }
    if(row.eff!=FALSE) {
      ri <- names(param)=="r0"
      if(method=="LA" || row.eff=="random"){ row.params=param[ri] } else {row.params <- param[ri]}
      if(row.eff=="random") {
        row.params <- lvs[,1]; lvs<- as.matrix(lvs[,-1])
        sigma<-exp(param["log_sigma"])[1]
        if(nlvr>1 && dependent.row) sigma <- c(exp(param[names(param)=="log_sigma"])[1],(param[names(param)=="log_sigma"])[-1])
      }
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
    if(beta0com){
      beta0=B[1]
      B = B[-1]
      cn<-colnames(Xd)
      Xd<-as.matrix(Xd[,-1])
      colnames(Xd)<-cn[-1]
    }
    new.loglik<-objr$env$value.best[1]


    if((n.i==1 || out$logL > abs(new.loglik)) && is.finite(new.loglik) && !inherits(optr, "try-error")){
      objrFinal<-objr1 <- objr; optrFinal<-optr1 <- optr;
      out$logL <- new.loglik
      if(num.lv > 0) {
        out$lvs <- lvs
        out$params$theta <- theta
        rownames(out$lvs) <- rownames(out$y);
        colnames(out$params$theta) <- colnames(out$lvs) <- paste("LV", 1:num.lv, sep="");
        rownames(out$params$theta) <- colnames(out$y)
      }
      if(!beta0com) names(beta0) <- colnames(out$y); 
      names(beta0) <- colnames(out$y); out$params$beta0 <- beta0;
      out$params$B <- B; names(out$params$B)=colnames(Xd)

      if(row.eff!=FALSE) {
        if(row.eff=="random"){  
          out$params$sigma <- sigma; 
          names(out$params$sigma) <- "sigma"
          if(num.lv>0 && dependent.row) names(out$params$sigma) <- paste("sigma",c("",1:num.lv), sep = "")
          }
        out$params$row.params <- row.params; 
        names(out$params$row.params) <- rownames(out$y)
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
        rownames(out$params$Br) <- rownames(out$params$sigmaB) <- colnames(out$params$sigmaB) <- colnames(xb)
      }
      if(family == "binomial") out$link <- link;
      out$row.eff <- row.eff
      out$time <- timeo
      out$start <- res
      out$Power <- Power
      pars <- optr$par

      if(method=="VA" && num.lv>0){
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
      }

      if(method == "VA" && !is.null(randomX)){
        Abb <- param[names(param) == "Abb"]
        dr <- ncol(xb)
        Ab <- array(0,dim=c(p,dr,dr))
        for (d in 1:dr){
          for(j in 1:p){
            Ab[j,d,d] <- exp(Abb[(d-1)*p + j]);
          }
        }
        if(length(Abb)>dr*p){
          k <- 0;
          for (c1 in 1:dr){
            r <- c1+1;
            while (r <= dr){
              for(j in 1:p){
                Ab[j,r,c1] <- Abb[dr*p+k*p+j];
                Ab[j,c1,r] <- Ab[j,r,c1];
              }
              k <- k+1; r <- r+1;
            }
          }
        }
        out$Ab <- Ab
      }
    }

    n.i <- n.i+1;
  }

    tr<-try({
      if(sd.errors && is.finite(out$logL)) {
        if(trace) cat("Calculating standard errors for parameters...\n")

        sdr <- optimHess(optrFinal$par, objrFinal$fn, objrFinal$gr,control = list(reltol=reltol,maxit=maxit))

        m <- dim(sdr)[1]; incl <- rep(TRUE,m); incld <- rep(FALSE,m)
        incl[names(objrFinal$par)=="Abb"] <- FALSE;
        incl[names(objrFinal$par)=="Au"] <- FALSE; 
        if(nlvr > 0) incld[names(objrFinal$par)=="Au"] <- TRUE
        
        if(beta0com){ 
          incl[names(objrFinal$par)=="b"] <- FALSE
        }
        
        if(row.eff=="random") {
          incl[names(objrFinal$par)=="r0"] <- FALSE; incld[names(objrFinal$par)=="r0"] <- FALSE
        } else {
          incl[names(objrFinal$par)=="log_sigma"] <- FALSE
        }
        if(row.eff==FALSE) incl[names(objrFinal$par)=="r0"] <- FALSE
        if(row.eff=="fixed") incl[1] <- FALSE
        
        
        if(is.null(randomX)) {
          incl[names(objrFinal$par)%in%c("Br","sigmaB","sigmaij")] <- FALSE
        } else {
          incl[names(objrFinal$par)=="Abb"] <- FALSE; incld[names(objrFinal$par)=="Abb"] <- TRUE
          incl[names(objrFinal$par)=="Br"] <- FALSE; incld[names(objrFinal$par)=="Br"] <- TRUE
          if(NCOL(xb)==1) incl[names(objrFinal$par) == "sigmaij"] <- FALSE
        }
        
        incl[names(objrFinal$par)=="Au"] <- FALSE; if(num.lv>0) incld[names(objrFinal$par)=="Au"] <- TRUE
        incl[names(objrFinal$par)=="u"] <- FALSE; incld[names(objrFinal$par)=="u"] <- TRUE
        
        if(familyn==0 || familyn==2 || familyn==6) incl[names(objrFinal$par)=="lg_phi"] <- FALSE
        if(familyn!=6) incl[names(objrFinal$par)=="zeta"] <- FALSE
        if(familyn==6) incl[names(objrFinal$par)=="zeta"] <- TRUE
        
        if(nlvr==0){
          incl[names(objrFinal$par)=="u"] <- FALSE;
          incld[names(objrFinal$par)=="u"] <- FALSE;
          incl[names(objrFinal$par)=="lambda"] <- FALSE;
          incl[names(objrFinal$par)=="Au"] <- FALSE;
        }

        if(method=="LA" || (num.lv==0 && (row.eff!="random" && is.null(randomX)))){
          incl[names(objrFinal$par)=="Au"] <- FALSE;

          covM <- try(MASS::ginv(sdr[incl,incl]))
          se <- try(sqrt(diag(abs(covM))))
          if(num.lv > 0 || row.eff == "random" || !is.null(randomX)) {
            sd.random <- sdrandom(objrFinal, covM, incl)
            prediction.errors <- list()
            if(!is.null(randomX)){
              prediction.errors$Br  <- matrix(diag(as.matrix(sd.random))[1:(ncol(xb)*ncol(y))], ncol(xb), ncol(y));
              sd.random <- sd.random[-(1:(ncol(xb)*ncol(y))),-(1:(ncol(xb)*ncol(y)))]
            }
            if(row.eff=="random"){
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
        } else {
          A.mat <- -sdr[incl,incl] # a x a
          D.mat <- -sdr[incld,incld] # d x d
          B.mat <- -sdr[incl,incld] # a x d
          cov.mat.mod <- try(MASS::ginv(A.mat-B.mat%*%solve(D.mat)%*%t(B.mat)),silent=T) 
          se <- sqrt(diag(abs(cov.mat.mod)))
          
          incla<-rep(FALSE, length(incl))
          incla[names(objrFinal$par)=="u"] <- TRUE
          out$Hess <- list(Hess.full=sdr, incla = incla, incl=incl, incld=incld, cov.mat.mod=cov.mat.mod)
          
        }

        if(row.eff=="fixed") {
          se.row.params <- c(0,se[1:(n-1)]); 
          names(se.row.params)  = rownames(out$y); se <- se[-(1:(n-1))] 
        }
        if(beta0com){
          se.beta0 <- se[1]; se <- se[-1];
        } else {
          se.beta0 <- se[1:p]; se <- se[-(1:p)];
        }
        se.B <- se[1:length(B)]; se <- se[-(1:length(B))];
        if(num.lv>0) {
          se.theta <- matrix(0,p,num.lv); se.theta[lower.tri(se.theta, diag = TRUE)]<-se[1:(p * num.lv - sum(0:(num.lv-1)))];
          colnames(se.theta) <- paste("LV", 1:num.lv, sep="");
          rownames(se.theta) <- colnames(out$y)
          out$sd$theta <- se.theta; se <- se[-(1:(p * num.lv - sum(0:(num.lv-1))))];
          # diag(out$sd$theta) <- diag(out$sd$theta)*diag(out$params$theta) !!!
        }
        out$sd$beta0 <- se.beta0; 
        if(!beta0com){ names(out$sd$beta0)  <-  colnames(out$y);}
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
        if(row.eff=="random") { 
          out$sd$sigma <- se[1:length(out$params$sigma)]*c(out$params$sigma[1],rep(1,length(out$params$sigma)-1)); 
          names(out$sd$sigma) <- "sigma"; 
          se=se[-(1:(1+num.lv))] 
          }
        if(!is.null(randomX)){
          nr <- ncol(xb)
          #out$sd$sigmaB <- se[1:ncol(xb)]*c(diag(out$params$sigmaB), rep(1,nr*(nr-1)/2 ) )
          out$sd$sigmaB <- se[1:ncol(xb)]*c(sqrt(diag(out$params$sigmaB))); 
          names(out$sd$sigmaB) <- c(paste("sd",colnames(xb),sep = "."))
          if(length(se)>1){
            se <- se[-(1:ncol(xb))]
            out$sd$corrpar <- se
          }
          
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
  if(inherits(tr, "try-error")) { cat("Standard errors for parameters could not be calculated.\n") }

  if(is.null(formula1)){ out$formula <- formula} else {out$formula <- formula1}

  out$Xrandom <- xb
  out$D <- Xd
  out$TMBfn <- objrFinal
  out$TMBfn$par <- optrFinal$par #ensure params in this fn take final values
  out$logL <- -out$logL
  
  if(method == "VA"){
    if(num.lv > 0) out$logL = out$logL + n*0.5*num.lv
    if(row.eff == "random") out$logL = out$logL + n*0.5
    if(!is.null(randomX)) out$logL = out$logL + p*0.5*ncol(xb)
    if(family=="gaussian") {
      out$logL <- out$logL - n*p*log(pi)/2
    }
  }

  return(out)
}

