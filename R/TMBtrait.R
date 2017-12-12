########################################################################################
## GLLVM fourth corner model, with estimation done via Laplace and Variational approximation using TMB-package
## Original author: Jenni Niku
##########################################################################################



trait.TMB <- function(y, X = NULL,TR=NULL,formula=NULL, num.lv = 2, family = "poisson",Lambda.struc="unstructured", row.eff = FALSE, reltol = 1e-5, seed = NULL,maxit = 1000, start.lvs = NULL, offset=NULL, sd.errors = TRUE,trace=TRUE,link="logit",n.init=1,start.params=NULL,start0=FALSE,optimizer="optim",starting.val="res",method="VA",randomX=NULL,Power=1.5,diag.iter=1,Lambda.start=0.1, jitter.var=0) {
  if(is.null(X) && !is.null(TR)) stop("Unable to fit a model that includes only trait covariates")

  n<-dim(y)[1]; p<-dim(y)[2];
  y <- as.data.frame(y)

  # change categorical variables to dummy variables
  num.X <- 0
  if(!is.null(X)) {
    num.X <- dim(X)[2]
    X.new <- NULL
      for (i in 1:num.X) {
      if(!is.factor(X[,i])) {
        if(!is.null(TR) && length(unique(X[,i]))>2){ Xi <- scale(X[,i]) } else { Xi <- X[,i] }
        X[,i] <- Xi
        X.new <- cbind(X.new,Xi); if(!is.null(colnames(X)[i])) colnames(X.new)[dim(X.new)[2]] <- colnames(X)[i]
      } else {
        dum <- model.matrix(~X[,i])#-1
        dum <- dum[,!(colnames(dum) %in% c("(Intercept)"))]
        #dum <- model.matrix(~X[,i]-1)
          colnames(dum)<-paste(colnames(X)[i],levels(X[,i])[-1],sep="")
          X.new <- cbind(X.new,dum)

      }
    }
    X.new <- data.frame(X.new);
  }
  if(!is.null(randomX)){
    method <- "LA"
    xb <- as.matrix(model.matrix(randomX,data = X.new))
    xb <- as.matrix(xb[,!(colnames(xb) %in% c("(Intercept)"))])
    Br <- matrix(0,ncol(xb),p)
    sigmaB <- diag(ncol(xb))
  }


  num.T <- 0
  if(!is.null(TR)) {
    num.T <- dim(TR)[2]
    T.new <- NULL
    for (i in 1:num.T) {
      if(!is.factor(TR[,i])  && length(unique(TR[,i]))>2) {
        TR[,i]=scale(TR[,i])
        T.new <- cbind(T.new,scale(TR[,i])); colnames(T.new)[dim(T.new)[2]] <- colnames(TR)[i]
      } else {
        dum <- model.matrix(~TR[,i]-1)
        colnames(dum) <- paste(colnames(TR)[i],levels(TR[,i]),sep="")
        T.new <- cbind(T.new,dum)
      }
    }
    T.new <- data.matrix(T.new);
  }


  if(is.null(formula)){
    n1 <- colnames(X)
    n2 <- colnames(TR)
    formula=paste("y~(",n1[1],sep = "")
  if(length(n1)>1){
    for(i1 in 2:length(n1)){
      formula <- paste(formula,n1[i1],sep = "+")
    }}
  formula=paste(formula,")*(",sep = "")
  formula=paste(formula,n2[1],sep = "")

  if(length(n2)>1){
    for(i2 in 1:length(n2)){
      formula <- paste(formula,n2[i2],sep = "+")
    }}
  formula=paste(formula,")",sep = "")
  formula=formula(formula)
  }

  if(is.null(formula)){
    X <- X.new; TR <- T.new;
    yX <- reshape(data.frame(cbind(y,X)),direction = "long", varying = colnames(y),v.names = "y")
    TR2<-data.frame(time=1:p,TR)
    yXT <- merge(yX,TR2,by="time")

    data <- yXT[,!(colnames(yXT) %in% c("time","y","id"))]
    TX <- kronecker(as.matrix(TR),as.matrix(X))
    Xd <- as.matrix(cbind(data,TX))
  } else {
    yX <- reshape(data.frame(cbind(y,X)),direction = "long", varying = colnames(y),v.names = "y")
    TR2<-data.frame(time=1:p,TR)
    yXT <- merge(yX,TR2,by="time")
    data <- yXT

    Xd <- as.matrix(model.matrix(formula,data = data))
    Xd <- as.matrix(Xd[,!(colnames(Xd) %in% c("(Intercept)"))])
    X <- as.matrix(X.new[,colnames(X.new) %in% colnames(Xd)]); TR=as.matrix(T.new[,colnames(T.new) %in% colnames(Xd)]);
    colnames(X) <- colnames(X.new)[colnames(X.new) %in% colnames(Xd)]; colnames(TR)=colnames(T.new)[colnames(T.new) %in% colnames(Xd)];
    nxd <- colnames(Xd)
    formulab <- paste("~",nxd[1],sep = "");
    for(i in 2:length(nxd)) formulab <- paste(formulab,nxd[i],sep = "+")
    formula <- formulab
  }
  #num.X <- dim(X)[2]
  #num.T <- dim(TR)[2];


  if(!(family %in% c("poisson","negative.binomial","binomial","ordinal")))
    stop("Selected family not permitted...sorry!")
  if(!(Lambda.struc %in% c("unstructured","diagonal")))
    stop("Lambda matrix (covariance of vartiational distribution for latent variable) not permitted...sorry!")
  if(num.lv == 1) Lambda.struc <- "diagonal" ## Prevents it going to "unstructured" loops and causing chaos
  trial.size <- 1

  y <- as.matrix(y)
  if(!is.numeric(y)) stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")

  if(is.null(rownames(y))) rownames(y) <- paste("Row",1:n,sep="")
  if(is.null(colnames(y))) colnames(y) <- paste("Col",1:p,sep="")
  if(!is.null(X)) { if(is.null(colnames(X))) colnames(X) <- paste("x",1:ncol(X),sep="") }

  out <-  list(y = y, X = X, TR = TR, num.lv = num.lv, row.eff = row.eff, logL = Inf, family = family, offset=offset,randomX=randomX,X.design=Xd)

  n.i<-1;
  if(n.init>1) seed <- sample(1:10000,n.init)
  while(n.i<=n.init){

num.X <- dim(X)[2]
num.T <- dim(TR)[2]

if(n.init>1 && trace) cat("initial run ",n.i,"\n");
#offset=NULL
res <- start.values.gllvm.TMB(y = y, X = X, TR = TR, family = family, offset=offset, trial.size = trial.size, num.lv = num.lv, start.lvs = start.lvs, seed = seed[n.i],starting.val=starting.val,formula = formula, jitter.var=jitter.var)
if(is.null(start.params)){
  beta0 <- res$params[,1]
  # common env params or different env response for each spp
  B <- NULL
  if(!is.null(TR) && !is.null(X)) {
    B <- c(res$B)[1:ncol(Xd)]
  }
  row.params <- NULL;
  if(row.eff!=FALSE){
    rmeany <- rowMeans(y)
    if(family=="binomial"){
      rmeany=1e-3+0.99*rmeany
      row.params <-  binomial(link = link)$linkfun(rmeany)-binomial(link = link)$linkfun(rmeany[1])
    } else {
      row.params <- rep(0,n);#log(rmeany)-log(rmeany[1])
    }
    res$row.params<-row.params
    }#rep(0,n)
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
      B=start.params$params$B;
    }
    fourth<-inter <- NULL; if(!is.null(TR) ) inter <- start.params$params$fourth   # let's treat this as a vector (vec(B'))'
    vameans <- theta <- lambda <- NULL

    if(num.lv > 0) theta<- c(start.params$params$theta) ## LV coefficients
    if(row.eff) row.params <- start.params$params$row.params ## row parameters
    if(num.lv > 0) { vameans <- matrix(start.params$lvs, ncol = num.lv);
    lambda <- start.params$Lambda}
  } else { stop("Model which is set as starting parameters isn't the suitable you are trying to fit. Check that attributes y, X, TR and row.eff match to each other.");}
}
if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)

phi <- NULL; if(family == "negative.binomial") { phis <- res$phi; if(any(phis>10))phis[phis>10]=10; if(any(phis<0.10))phis[phis<0.10]=0.10; phi <- 1/phis }
q=num.lv

if(!is.null(row.params)){ r0=row.params} else {r0=rep(0,n)}
a=c(beta0)
if(num.lv > 0) theta=theta[lower.tri(theta,diag = TRUE)]
if(num.lv > 0) u=vameans
if(!is.null(phi)) {phi=(phi)} else {phi=rep(1,p)}
q = num.lv
sigma=1

if(num.lv>0){
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
    #if(start.params$Lambda.struc=="diagonal") Au=c(Au,log(start.params$Lambda[,d]))
  }
  if(Lambda.struc!="diagonal" && diag.iter==0){
    Au=c(Au,rep(0,num.lv*(num.lv-1)/2*n))
  }
}}
Ar=rep(1,n)


optr<-NULL
timeo<-NULL
se=NULL


if(row.eff==FALSE){xr=matrix(0,1,p)} else {xr=matrix(1,1,p); xr[1,1]=0}
#if(!is.null(X)){Xd=cbind(1,X)} else {Xd=matrix(1,n)}
extra=0
if(family == "poisson") { familyn=0}
if(family == "negative.binomial") { familyn=1}
if(family == "binomial") { familyn=2}


if(row.eff=="random" || !is.null(randomX)){
  if(num.lv>0){
    #compile("LAtraitsrandomb.cpp")
    #dyn.load(dynlib("LAtraitsrandomb"))
    objr <- TMB::MakeADFun(
      data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, num_lv = num.lv,family=familyn,extra=extra), silent=!trace,
      parameters = list(r0=matrix(r0),Br=Br, b = a,B=matrix(B),lambda = theta, u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaB[upper.tri(sigmaB)]),#sigmaB=c(log(diag(sigmaB)),sigmaB[lower.tri(sigmaB)])),
      inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
      random = c("Br","u"), DLL = "LAtraitsrandomb")
  } else {
    #dyn.load(dynlib("LAtraitsrandomb0"))
    objr <- TMB::MakeADFun(
      data = list(y = y, x = Xd,xr=xr,xb=xb,offset=offset, family=familyn,extra=extra), silent=!trace,
      parameters = list(r0=matrix(r0),Br=Br, b = a,B=matrix(B),lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaB[upper.tri(sigmaB)]),#sigmaB=c(log(diag(sigmaB)),sigmaB[lower.tri(sigmaB)])),
      inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
      random = c("Br"), DLL = "LAtraitsrandomb0")
  }

} else {
  if(num.lv>0) {
  if(method=="VA"){
  #compile("VAtraits.cpp")
  #dyn.load(dynlib("VAtraits"))
  objr <- TMB::MakeADFun(
    data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=1), silent=!trace,
    parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B),lambda = theta, u = u,lg_phi=log(phi),Au=Au),
    inner.control=list(mgcmax = 1e+200,maxit = 1000),
    DLL = "gllvm")
  } else {
  #compile("LAtraits.cpp")
  #dyn.load(dynlib("LAtraits"))
  objr <- TMB::MakeADFun(
    data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=1), silent=!trace,
    parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B),lambda = theta, u = u,lg_phi=log(phi),Au=0),
    inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
    random = c("u"), DLL = "gllvm")
  }
  } else {
    #dyn.load(dynlib("LAtraits0"))
    objr <- TMB::MakeADFun(
      data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=1,model=1), silent=!trace,
      parameters = list(r0=matrix(r0), b = rbind(a),B=matrix(B),lambda = 0, u = matrix(0),lg_phi=log(phi),Au=0),
      inner.control=list(mgcmax = 1e+200,maxit = 1000,tol10=0.01),
      DLL = "gllvm")
  }

}

if(optimizer=="nlminb") {
  timeo<-system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,maxit=maxit)),silent = TRUE))
}
if(optimizer=="optim") {
  timeo<-system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
}

if(diag.iter>0 && Lambda.struc=="unstructured" && method =="VA" && num.lv>0 && row.eff!="random" && is.null(randomX)){
  param1=optr$par
  nam=names(param1)
  r1=matrix(param1[nam=="r0"])
  b1=rbind(param1[nam=="b"]+rnorm(p,0,0.01))
  B1=matrix(param1[nam=="B"])
  lambda1=param1[nam=="lambda"]
  u1=matrix(param1[nam=="u"],n,num.lv)
  lg_phi1=param1[nam=="lg_phi"]
  Au1=c(param1[nam=="Au"],rep(0,num.lv*(num.lv-1)/2*n))
  
        #dyn.load(dynlib("VAtraits"))
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,method=0,model=1), silent=!trace,
          parameters = list(r0=r1, b = b1,B=B1,lambda = lambda1, u = u1,lg_phi=lg_phi1,Au=Au1),
          inner.control=list(mgcmax = 1e+200,maxit = 1000),
          DLL = "gllvm")

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
if(family=="ZIP") {
  lp0=param[names(param)=="lg_phi"]; out$lp0=lp0
  phis=exp(lp0)/(1+exp(lp0));#log(phis); #
}
bi<-names(param)=="b"
Bi<-names(param)=="B"
li<-names(param)=="lambda"
ui<-names(param)=="u"
if(row.eff!=FALSE) {
  ri=names(param)=="r0"
  if(method=="LA" || row.eff=="random"){ row.params=param[ri] } else {row.params=param[ri]}#c(0,param[ri])}
  if(row.eff=="random") sigma<-exp(param["log_sigma"])
}
if(!is.null(randomX)){
  Bri=names(param)=="Br"
  Br=matrix(param[Bri],ncol(xb),p)
  Sri=names(param)=="sigmaB"
    L=diag(ncol(xb))
    if(ncol(xb)>1){
      sigmaB=diag(exp(param[Sri]))
      Srij=names(param)=="sigmaij"
      Sr=param[Srij]
      L[upper.tri(L)]=Sr
      D=diag(diag(t(L)%*%L))
    } else{
      D=1
      sigmaB=(exp(param[Sri]))
    }
    sigmaB_=solve(sqrt(D))%*%(t(L)%*%L)%*%solve(sqrt(D))
    sigmaB=sigmaB%*%sigmaB_%*%t(sigmaB)

}
beta0=param[bi]
B=param[Bi]
if(num.lv > 0){
  lvs=(matrix(param[ui],n,q))
  theta=matrix(0,p,num.lv)
  theta[lower.tri(theta,diag = TRUE)] <- param[li];
}
new.loglik<-objr$env$value.best[1]


if((n.i==1 || out$logL > abs(new.loglik)) && !inherits(optr, "try-error") && new.loglik>0){
  objr1=objr;optr1=optr;
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
  if(family %in% c("negative.binomial","tweedie")) {
    out$params$phi <- 1/phis; names(out$params$phi) <- colnames(out$y);
    out$params$inv.phi <- phis; names(out$params$inv.phi) <- colnames(out$y);
  }
  if(family =="ZIP") {
    out$params$phi <- phis; names(out$params$phi) <- colnames(out$y);
  }
  if(!is.null(randomX)){
    out$params$Br=Br
    out$params$sigmaB=sigmaB
    out$corr=sigmaB_
  }
  if(family == "binomial") out$link <- link;
  out$row.eff=row.eff
  out$time=timeo
  out$start<-res
  pars=optr$par

if(method=="VA" && (num.lv>0 || row.eff=="random")){
  param<-objr$env$last.par.best
  Au=param[names(param)=="Au"]
  A=array(0,dim=c(n,num.lv,num.lv))
  #A=array(0,dim=c(num.lv,num.lv,n))
  for (d in 1:num.lv){
    for(i in 1:n){
      A[i,d,d]=exp(Au[(d-1)*n+i]);
      #A[d,d,i]=Au[(d-1)*n+i];
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
          #A[r,c1,i]=Au[num.lv*n+k*n+i];
          #A[c1,r,i]=A[r,c1,i];
        }
        k=k+1; r=r+1;
      }
    }
  }
  out$A=A
  if(row.eff=="random"){
    Ar=exp(param[names(param)=="lg_Ar"])
    out$Ar
  }
  }
}

n.i <- n.i+1;
  }

  tr<-try({
    if(sd.errors) {
      if(trace) cat("Calculating standard errors for parameters...\n")

      sdr <- optimHess(pars, objr$fn, objr$gr,control = list(reltol=reltol,maxit=maxit))

      m<-dim(sdr)[1]; incl=rep(TRUE,m); incld=rep(FALSE,m)
      incl[names(objr$par)=="Au"]=FALSE; incld[names(objr$par)=="Au"]=TRUE
      if(row.eff=="random") {
        incl[names(objr$par)=="lg_Ar"]=FALSE; incld[names(objr$par)=="lg_Ar"]=TRUE
        incl[names(objr$par)=="r0"]=FALSE; incld[names(objr$par)=="r0"]=TRUE
      }
      if(row.eff==FALSE) incl[names(objr$par)=="r0"]=FALSE
      if(row.eff=="fixed") incl[1]=FALSE
      incl[names(objr$par)=="u"]=FALSE; incld[names(objr$par)=="u"]=TRUE
      if(familyn==0 || familyn==2) incl[names(objr$par)=="lg_phi"]=FALSE
      if(num.lv==0){
        incl[names(objr$par)=="u"]=FALSE;
        incl[names(objr$par)=="lambda"]=FALSE;
      }
      if(method=="LA" || num.lv==0){
        incl[names(objr$par)=="Au"]=FALSE;
        se<- try(sqrt(diag(abs(MASS::ginv(-sdr[incl,incl])))))
      } else {
        A.mat=-sdr[incl,incl] # a x a
        D.mat=-sdr[incld,incld] # d x d
        B.mat=-sdr[incl,incld] # a x d
        cov.mat.mod <- MASS::ginv(A.mat-B.mat%*%solve(D.mat)%*%t(B.mat))
        se<- sqrt(diag(abs(cov.mat.mod)))
      }

      if(row.eff=="fixed") { se.row.params <- c(0,se[1:(n-1)]); names(se.row.params)  = rownames(out$y); se <- se[-(1:(n-1))] }
      se.beta0 <- se[1:p]; se<-se[-(1:p)];
      #if(!is.null(X)){ out$sd$Xcoef <- as.matrix(sebetaM[,-1]);
      #  rownames(out$sd$Xcoef) <- colnames(y); colnames(out$sd$Xcoef) <- colnames(X);}
      se.B=se[1:length(B)]; se<-se[-(1:length(B))];
      if(num.lv>0) {
        se.theta<-matrix(0,p,num.lv); se.theta[lower.tri(se.theta, diag = TRUE)]<-se[1:(p * num.lv - sum(0:(num.lv-1)))];
        colnames(se.theta) = paste("LV", 1:num.lv, sep="");
        rownames(se.theta) = colnames(out$y)
        out$sd$theta<-se.theta; se <- se[-(1:(p * num.lv - sum(0:(num.lv-1))))];
      }
      out$sd$beta0=se.beta0; names(out$sd$beta0)  = colnames(out$y);
      out$sd$B=se.B; names(out$sd$B)  = colnames(Xd)
      if(row.eff=="fixed") {out$sd$row.params=se.row.params}

      if(family %in% c("negative.binomial","tweedie")) {
        se.lphis=se[1:p];  out$sd$inv.phi=se.lphis*out$params$inv.phi;
        #out$sd$lphi=se.lphis; out$params$lphi=log(out$params$phi)
        out$sd$phi=se.lphis*out$params$phi;
        names(out$sd$inv.phi) <- names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
      }
      if(family %in% c("ZIP")) {
        se.phis=se[1:p];
        out$sd$phi=se.phis*exp(lp0)/(1+exp(lp0))^2;#
        names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
      }
      if(row.eff=="random") { out$sd$sigma <- se[1]*out$params$sigma; names(out$sd$sigma) <- "sigma"; se=se[-1] }
      if(!is.null(randomX)){
        nr=ncol(xb)
        out$sd$sigmaB <- se*c(diag(out$params$sigmaB), rep(1,nr*(nr-1)/2 ) )
      }

    }})
  if(inherits(tr, "try-error")) { cat("Standard errors for parameters could not be calculated.\n") }

  out$D=Xd
  out$logL=-out$logL
  return(out)
}

