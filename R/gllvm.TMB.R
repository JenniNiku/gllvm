########################################################################################
## GLLVM, with estimation done via Variational approximation using TMB-package
## Original author: Jenni Niku
##########################################################################################


gllvm.TMB <- function(y, X = NULL, formula = NULL, num.lv = 2, family = "poisson",method="VA",Lambda.struc="unstructured", row.eff = FALSE, reltol = 1e-6, seed = NULL,maxit = 1000, start.lvs = NULL, offset=NULL, sd.errors = TRUE,trace=TRUE,link="logit",n.init=1,restrict=30,start.params=NULL, optimizer="optim",starting.val="res",Power=1.5,diag.iter=1,Lambda.start=0.1, jitter.var=0) {
  n <- dim(y)[1]; p <- dim(y)[2];
  num.lv <- num.lv
  y <- as.matrix(y)
  if(!is.numeric(y)) stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
  if((family %in% c("tweedie","ZIP")) && method=="VA") stop("family=\"",family,"\" : family not implemented with VA method, change the method to 'LA'")
  if(is.null(rownames(y))) rownames(y) <- paste("Row",1:n,sep="")
  if(is.null(colnames(y))) colnames(y) <- paste("Col",1:p,sep="")

  num.X=0;
  if(!is.null(X)){
  if(!is.null(formula)){
    xb=as.matrix(model.matrix(formula,data = data.frame(X)))
    X=as.matrix(xb[,!(colnames(xb) %in% c("(Intercept)"))])
    num.X <- dim(X)[2]
  } else {
    X.new <- NULL
    num.X <- dim(X)[2];
    if(is.null(colnames(X)))	colnames(X)=paste("x",1:num.X,sep = "")
    for (i in 1:num.X) {
      if(is.factor(X[,i])) {
        dum <- model.matrix(~X[,i])#-1
        dum <- dum[,!(colnames(dum) %in% c("(Intercept)"))]
        colnames(dum)<-paste(names(X)[i],levels(X[,i])[-1],sep="")
        X.new <- cbind(X.new,dum)
      } else {
        X.new <- cbind(X.new,X[,i]); if(!is.null(colnames(X)[i])) colnames(X.new)[dim(X.new)[2]] <- colnames(X)[i]
        }
    }
    X <- data.matrix(X.new);
    num.X <- dim(X)[2];
  }}
  ## Set initial values for model parameters (including dispersion prm) and latent variables
  if(!is.null(seed)) {set.seed(seed);}

  n.i<-1;
  out <- list(y=y,X=X, logL = Inf)
  if(n.init>1) seed<-sample(1:10000,n.init)
  while(n.i<=n.init){
    if(n.init>1 && trace) cat("Initial run ",n.i,"\n");

    fit <- start.values.gllvm.TMB(y = y, X = X, TR = NULL, family = family, offset= offset, num.lv = num.lv, start.lvs = start.lvs, seed = seed[n.i],starting.val=starting.val,power=Power, jitter.var=jitter.var)
    if(is.null(start.params)){
      beta0 <- fit$params[,1]
      betas <- NULL; if(!is.null(X)) betas <- c(fit$params[,2:(num.X + 1)])
      lambdas <- NULL;
      if(num.lv > 0) {
        lambdas <- as.matrix(fit$params[,(ncol(fit$params) - num.lv + 1):ncol(fit$params)])
        if(num.lv > 0) lambdas[upper.tri(lambdas)] <- 0
        covM.lvs<- array(NA,dim=c(n,num.lv,num.lv))
      }
      row.params <- NULL;
      if(row.eff!=FALSE) {
        row.params <- rep(0,n);
        #if(row.eff=="random") row.params <-  log(rowMeans(y))-log(mean(y));
        fit$row.params<-row.params }#rep(0,n)
      lvs <- NULL; if(num.lv > 0) lvs <- matrix(fit$index, ncol = num.lv);
    } else{
      if(dim(start.params$y)==dim(y) && is.null(X)==is.null(start.params$X) && (row.eff==start.params$row.eff)){
        beta0 <- start.params$params$beta0 ## column intercepts
        betas <- NULL; if(!is.null(X)) betas <- c(start.params$params$Xcoef) ## covariates coefficients
        lambdas <- NULL; if(num.lv > 0) lambdas <- (start.params$params$theta); if(num.lv > 0) lambdas[upper.tri(lambdas)] <- 0#lambdas=matrix(1,p,num.lv) ## LV coefficients
        row.params <- NULL; if(start.params$row.eff!=FALSE){ row.params <- start.params$params$row.params; row.params[1]=0 }## row parameters
        lvs <- NULL; if(num.lv > 0){ lvs <- matrix(start.params$lvs, ncol = num.lv);
        covM.lvs<- array(NA,dim=c(n,num.lv,num.lv))}## LVs
      } else { stop("Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that attributes y, X and row.eff match to each other.");}
    }
    phis <- NULL; if(family == "negative.binomial") {phis <- 1/fit$phi; if(any(phis>10))phis[phis>10]=10; if(any(phis<0.10))phis[phis<0.10]=0.10;}# && phi.upd=="inv.phi"
    if(family == "tweedie") {phis <- fit$phi; if(any(phis>10))phis[phis>10]=10; if(any(phis<0.10))phis[phis<0.10]=0.10; phis= (phis)}
    if (family == "ZIP") {phis <- (colMeans(y==0)*0.98)+0.01; phis<-phis/(1-phis)} # ZIP probability

    if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)

    current.loglik <- -1e6; iter <- 1; err <- 10;
    ## LA-likelihood
    if(!is.null(row.params)){ r0=row.params} else {r0=rep(0,n)}
    a=c(beta0)
    b=NULL; if(!is.null(X)) b=matrix(betas,ncol(X),p,byrow =TRUE)
    if(num.lv > 0) lambda=lambdas[lower.tri(lambdas,diag = TRUE)]
    if(num.lv > 0) u=lvs
    if(!is.null(phis)) {phi=(phis)} else {phi=rep(1,p)}
    q = num.lv
    sigma=1

    optr<-NULL
    timeo<-NULL
    se=NULL


  if(method=="VA" && num.lv>0){
    if(is.null(start.params) || start.params$method!="VA"){
      if(Lambda.struc=="diagonal" || diag.iter>0){
        Au=log(rep(Lambda.start,num.lv*n)) #1/2, 1
      } else{
        Au=c(log(rep(Lambda.start,num.lv*n)),rep(0,num.lv*(num.lv-1)/2*n)) #1/2, 1
      }
    } else {
      Au=NULL
      for(d in 1:num.lv) {
        if(start.params$Lambda.struc=="unstructured" || length(dim(start.params$Lambda))==3){ Au=c(Au,log(start.params$Lambda[,d,d]))
        } else { Au=c(Au,log(start.params$Lambda[,d])) }
      }
      if(Lambda.struc!="diagonal" && diag.iter==0){
        Au=c(Au,rep(0,num.lv*(num.lv-1)/2*n))
      }
    }
    Ar=rep(1,n)

    if(row.eff==FALSE){xr=matrix(0,1,p)} else {xr=matrix(1,1,p)}
    if(!is.null(X)){Xd=cbind(1,X)} else {Xd=matrix(1,n)}
    extra=0
    if(family == "poisson") { familyn=0}
    if(family == "negative.binomial") { familyn=1}
    if(family == "binomial") { familyn=2}

    if(row.eff=="random"){
      if(num.lv>0){
        #dyn.load(dynlib("VArandom"))
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn), silent=TRUE,
          parameters = list(r0=matrix(r0), b = rbind(a,b),lambda = lambda, u = u,lg_phi=log(phi),log_sigma=log(sigma),Au=Au,lg_Ar=log(Ar)),
          inner.control=list(mgcmax = 1e+200,maxit = 1000),
          DLL = "VArandom")
      } else {
        #dyn.load(dynlib("VArandom0"))
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd,xr=xr,offset=offset,family=familyn), silent=TRUE,
          parameters = list(r0=matrix(r0), b = rbind(a,b),lg_phi=log(phi),log_sigma=log(sigma),lg_Ar=log(Ar)),
          inner.control=list(mgcmax = 1e+200,maxit = 1000),
          DLL = "VArandom0")
      }
    } else {
      #dyn.load(dynlib("VAfixed2"))
      objr <- TMB::MakeADFun(
        data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=0), silent=TRUE,
        parameters = list(r0=matrix(r0), b = rbind(a,b),B=matrix(0),lambda = lambda, u = u,lg_phi=log(phi),Au=Au),
        inner.control=list(mgcmax = 1e+200,maxit = maxit),
        DLL = "gllvm")##GLLVM
    }
    if(optimizer=="nlminb") {
      timeo<-system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol)),silent = TRUE))
    }
    if(optimizer=="optim") {
      timeo<-system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
    }
    if(diag.iter>0 && Lambda.struc=="unstructured" && num.lv>1){
      param1=optr$par
      nam=names(param1)
      r1=matrix(param1[nam=="r0"])
      b1=matrix(param1[nam=="b"]+rnorm(p,0,0.02),num.X+1,p)
      lambda1=param1[nam=="lambda"]
      u1=matrix(param1[nam=="u"],n,num.lv)
      lg_phi1=param1[nam=="lg_phi"]
      log_sigma1=param1[nam=="log_sigma"]
      Au1=c(param1[nam=="Au"],rep(0,num.lv*(num.lv-1)/2*n))
      lg_Ar1=matrix(param1[nam=="lg_Ar"])
      if(row.eff=="random"){
        if(num.lv>0){
          #dyn.load(dynlib("VArandom"))
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn), silent=TRUE,
            parameters = list(r0=r1, b = b1,lambda = lambda1, u = u1,lg_phi=lg_phi1,log_sigma=log_sigma1,Au=Au1,lg_Ar=lg_Ar1),
            inner.control=list(mgcmax = 1e+200,maxit = 1000),
            DLL = "VArandom")
        } else {
          #dyn.load(dynlib("VArandom0"))
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset,family=familyn), silent=TRUE,
            parameters = list(r0=r1, b = b1,lg_phi=lg_phi1,log_sigma=log_sigma1,lg_Ar=lg_Ar1),
            inner.control=list(mgcmax = 1e+200,maxit = 1000),
            DLL = "VArandom0")
        }
      } else {
        #dyn.load(dynlib("VAfixed2"))
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=0), silent=TRUE,
          parameters = list(r0=r1, b = b1,B=matrix(0),lambda = lambda1, u = u1,lg_phi=lg_phi1,Au=Au1), #log(phi)
          inner.control=list(mgcmax = 1e+200,maxit = maxit),
          DLL = "gllvm")#GLLVM#
      }
      if(optimizer=="nlminb") {
        timeo<-system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol)),silent = TRUE))
      }
      if(optimizer=="optim") {
        timeo<-system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
      }
    }

    param<-objr$env$last.par.best
    if(family =="negative.binomial") {
      phis=exp(param[names(param)=="lg_phi"])
    }
    bi<-names(param)=="b"
    li<-names(param)=="lambda"
    ui<-names(param)=="u"
    if(row.eff!=FALSE) {
      ri=names(param)=="r0"
      if(row.eff=="fixed") row.params=param[ri]#c(0,param[ri])
      if(row.eff=="random"){ sigma<-exp(param["log_sigma"]); row.params=param[ri]}
    }
    betaM<-matrix(param[bi],p,num.X+1,byrow=TRUE)
    beta0=betaM[,1]
    if(!is.null(X)) betas=betaM[,-1]
    if(num.lv > 0){
      lvs<-(matrix(param[ui],n,q))
      theta=matrix(0,p,num.lv)
      if(p>1) {
        theta[lower.tri(theta,diag=TRUE)] <- param[li];
      } else {theta=param[li]}

    }
    new.loglik<-objr$env$value.best[1]

    }

    if(method=="LA" || num.lv==0){
      if(row.eff==FALSE){xr=matrix(0,1,p)} else {xr=matrix(1,1,p)}
      if(!is.null(X)){Xd=cbind(1,X)} else {Xd=matrix(1,n)}
      extra=0
      if(family == "poisson") {familyn=0}
      if(family == "negative.binomial") {familyn=1}
      if(family == "binomial") {
        familyn=2;
        if(link=="probit") extra=1
      }
      if(family=="tweedie"){ familyn=3; extra=Power}
      if(family=="ZIP"){ familyn=4;}

      if(row.eff=="random"){
        if(num.lv>0){
          #dyn.load(dynlib("LArandom"))
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,offset=offset, num_lv = num.lv,family=familyn,extra=extra), silent=!trace,
            parameters = list(r0=r0, b = rbind(a,b),lambda = lambda, u = u,lg_phi=log(phi),log_sigma=log(sigma)),
            inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
            random = c("r0","u"), DLL = "LArandom")
        }else{
          #dyn.load(dynlib("LArandom0"))
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,offset=offset,family=familyn,extra=extra), silent=!trace,
            parameters = list(r0=r0, b = rbind(a,b),lg_phi=log(phi),log_sigma=log(sigma)),
            inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
            random = c("r0"), DLL = "LArandom0")
        }
      } else {
        if(num.lv>0){
          #dyn.load(dynlib("gllvm")) # r0(1) asetettu 0:ksi
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=0), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a,b),B=matrix(0),lambda = lambda, u = u,lg_phi=log(phi),Au=0),
            inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
            random = c("u"), DLL = "gllvm")
        }else{
          #dyn.load(dynlib("LAfixed0"))
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=0), silent=!trace,
            parameters = list(r0=matrix(r0), b = rbind(a,b),B=matrix(0),lambda = 0, u = matrix(0),lg_phi=log(phi),Au=0),
            inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
            DLL = "gllvm")
        }
      }
      if(family=="ZIP" && FALSE) {
        m<-length(objr$par)
        low=rep(-restrict,m); upp=rep(restrict,m);
        low[names(objr$par)=="lg_phi"]=0.0; upp[names(objr$par)=="lg_phi"]=1#0.99
        timeo<-system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol),lower = low,upper = upp),silent = TRUE))
      }
        if(optimizer=="nlminb") {
          timeo<-system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,maxit=maxit)),silent = TRUE))
        }
        if(optimizer=="optim") {
          timeo<-system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
        }


      param<-objr$env$last.par.best
      bi<-names(param)=="b"
      li<-names(param)=="lambda"
      ui<-names(param)=="u"
      if(row.eff!=FALSE) {
        ri=names(param)=="r0"
        row.params=param[ri]
        if(row.eff=="random") sigma<-exp(param["log_sigma"])
      }
      betaM<-matrix(param[bi],p,num.X+1,byrow=TRUE)
      beta0=betaM[,1]
      if(!is.null(X)) betas=betaM[,-1]
      if(num.lv > 0){
        lvs<-(matrix(param[ui],n,q))
        theta=matrix(0,p,num.lv)
        if(p>1) {
          theta[lower.tri(theta,diag=TRUE)] <- param[li];
        } else {theta=param[li]}
      }
      new.loglik<-objr$env$value.best[1]
      if(family %in% c("negative.binomial","tweedie","ZIP")) {
        phis=exp(param[names(param)=="lg_phi"])
        if(family=="ZIP") {
          lp0=param[names(param)=="lg_phi"]; out$lp0=lp0
          phis=exp(lp0)/(1+exp(lp0));#log(phis); #
          }
      }

    }


    if(((n.i==1 || out$logL > abs(new.loglik)) && new.loglik>0) && !inherits(optr, "try-error")){
      out$start<-fit
      objr1=objr; optr1=optr;
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
      if(family =="tweedie") {
        out$params$phi <- phis; names(out$params$phi) <- colnames(out$y);
      }
      if(family =="ZIP") {
        out$params$phi <- phis; names(out$params$phi) <- colnames(out$y);
      }
      if(row.eff!=FALSE) {
        if(row.eff=="random"){ out$params$sigma=sigma; names(out$params$sigma)="sigma"}
        out$params$row.params <- row.params; names(out$params$row.params) <- rownames(out$y)
      }
      if(family == "binomial") out$link <- link;
      if(family == "tweedie") out$Power <- Power;

      out$row.eff=row.eff
      out$time=timeo
      pars=optr$par

      if(method=="VA" && num.lv>0){
      param<-objr$env$last.par.best
      if(num.lv>0){
        Au=param[names(param)=="Au"]
        A=array(0,dim=c(n,num.lv,num.lv))
        for (d in 1:num.lv){
          for(i in 1:n){
            A[i,d,d]=exp(Au[(d-1)*n+i]);
          }
        }
        if(length(Au)>num.lv*n){
          k=0;
          for (c1 in 1:num.lv){
            r=c1+1;
            while (r <=num.lv){
              for(i in 1:n){
                A[i,r,c1]=Au[num.lv*n+k*n+i];
                A[i,c1,r]=A[i,r,c1];
              }
              k=k+1; r=r+1;
            }
          }
        }
        out$A=A
      }
      if(row.eff=="random"){
        Ar=exp(param[names(param)=="lg_Ar"])
        out$Ar=Ar
      }}
    }

    n.i <- n.i+1;
  }
  tr<-try({
    if(sd.errors && !is.infinite(out$logL)) {
      if(trace) cat("Calculating standard errors for parameters...\n")

      if(family=="ZIP") {
        p0i<-names(pars)=="lg_phi"
        p0=pars[p0i]
        p0 = p0+runif(p,0,0.001)
        pars[p0i] = p0
      }
      sdr <- optimHess(pars, objr$fn, objr$gr, control = list(reltol=reltol,maxit=maxit))#maxit=maxit
      m<-dim(sdr)[1]; incl=rep(TRUE,m); incld=rep(FALSE,m)
      incl[names(objr$par)=="B"]=FALSE
      if(method=="LA" || num.lv==0){
        incl[names(objr$par)=="Au"]=FALSE;
        if(row.eff=="fixed") incl[1]=FALSE
        if(row.eff==FALSE) incl[names(objr$par)=="r0"]=FALSE
        if(familyn==0 || familyn==2) incl[names(objr$par)=="lg_phi"]=FALSE
        if(num.lv==0){
          incl[names(objr$par)=="u"]=FALSE;
          incl[names(objr$par)=="lambda"]=FALSE;
        }
        se<- try(sqrt(diag(abs(MASS::ginv(-sdr[incl,incl])))))
      } else {
        incl[names(objr$par)=="Au"]=FALSE; incld[names(objr$par)=="Au"]=TRUE
        if(row.eff=="random") {
          incl[names(objr$par)=="lg_Ar"]=FALSE; incld[names(objr$par)=="lg_Ar"]=TRUE
          incl[names(objr$par)=="r0"]=FALSE; incld[names(objr$par)=="r0"]=TRUE
        }
        incl[names(objr$par)=="u"]=FALSE; incld[names(objr$par)=="u"]=TRUE
        if(row.eff=="fixed") incl[1]=FALSE
        if(row.eff==FALSE) incl[names(objr$par)=="r0"]=FALSE
        if(familyn==0 || familyn==2) incl[names(objr$par)=="lg_phi"]=FALSE
        A.mat=-sdr[incl,incl] # a x a
        D.mat=-sdr[incld,incld] # d x d
        B.mat=-sdr[incl,incld] # a x d
        cov.mat.mod <- MASS::ginv(A.mat-B.mat%*%solve(D.mat)%*%t(B.mat))
        #cov.mat.mod <- MASS::ginv(A.mat)
        se<- sqrt(diag(abs(cov.mat.mod)))
      }

      if(row.eff=="fixed") { se.row.params <- c(0,se[1:(n-1)]); names(se.row.params)  = rownames(out$y); se <- se[-(1:(n-1))] }
      sebetaM<-matrix(se[1:((num.X+1)*p)],p,num.X+1,byrow=TRUE);  se<-se[-(1:((num.X+1)*p))]
      if(num.lv>0) {
        se.lambdas<-matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag = TRUE)]<-se[1:(p * num.lv - sum(0:(num.lv-1)))];
        colnames(se.lambdas) = paste("LV", 1:num.lv, sep="");
        rownames(se.lambdas) = colnames(out$y)
        out$sd$theta<-se.lambdas; se <- se[-(1:(p * num.lv - sum(0:(num.lv-1))))];
      }
      out$sd$beta0 <- sebetaM[,1]; names(out$sd$beta0)  = colnames(out$y);
      if(!is.null(X)){
        out$sd$Xcoef <- matrix(sebetaM[,-1],nrow = nrow(sebetaM));
        rownames(out$sd$Xcoef) <- colnames(y); colnames(out$sd$Xcoef) <- colnames(X);
      }
      if(row.eff=="fixed") {out$sd$row.params=se.row.params}

      if(family %in% c("negative.binomial")) {
        se.lphis=se[1:p];  out$sd$inv.phi=se.lphis*out$params$inv.phi;
        #out$sd$lphi=se.lphis; out$params$lphi=log(out$params$inv.phi)
        out$sd$phi=se.lphis*out$params$phi;
        names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
      }
      if(family %in% c("tweedie")) {
        se.lphis=se[1:p];
        out$sd$phi=se.lphis*out$params$phi;
        names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
      }
      if(family %in% c("ZIP")) {
        se.phis=se[1:p];
        out$sd$phi=se.phis*exp(lp0)/(1+exp(lp0))^2;#
        names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
      }
      if(row.eff=="random") { out$sd$sigma <- se*out$params$sigma; names(out$sd$sigma) <- "sigma" }

    }})
  if(inherits(tr, "try-error")) { cat("Standard errors for parameters could not be calculated.\n") }


  out$logL=-out$logL
  return(out)
}

