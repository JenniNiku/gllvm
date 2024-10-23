##############################################################
## Fits GLLVM model using VA method.
## Authors: Hui, Niku, Taskinen
##
## If environmental covariates X are included in the model,
## the regression coeffs are estimated for each spp separately.
## When additional trait covariates T are included, we'll fit
## the "fourth corner lv model".
##
##
## 27.10.2015
##############################################################


gllvm.VA <- function(y, X = NULL, TR = NULL, formula=NULL, family = "poisson",
      num.lv = 2, max.iter = 200, eps = 1e-6, row.eff = FALSE,
      Lambda.struc = "unstructured", trace = TRUE, plot = FALSE, sd.errors = FALSE,
      start.lvs = NULL, offset = NULL, maxit = 100, diag.iter = 5, seed = NULL,
      get.fourth = TRUE, get.trait = TRUE, n.init = 1, constrOpt = FALSE,
      restrict = 30, start.params = NULL, starting.val = "res", Lambda.start = 0.1,
      jitter.var = 0, yXT = NULL) {

  if(is.null(X) && !is.null(TR)) stop("Unable to fit a model that includes only trait covariates")

  term=NULL
  n<-dim(y)[1]; p<-dim(y)[2];
  y=as.data.frame(y)
  X1=X; TR1=TR;
  formula1=formula

  if(!is.null(TR) && NCOL(X)<1) stop("No covariates in the model, fit the model using gllvm(y,family=",family,"...)")

  # change categorical variables to dummy variables
  num.X <- 0
  if(!is.null(X)) {
    num.X <- dim(X)[2]
    X.new <- NULL
    for (i in 1:num.X) {
      if(!is.factor(X[,i])) {
        if(!is.null(TR) && length(unique(X[,i]))>2){ Xi <- scale(X[,i]) } else { Xi <- X[,i] }
        X[,i]=Xi
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

  num.T <- 0
  T.new <- NULL
  if(!is.null(TR)) {
    num.T <- dim(TR)[2]
    T.new <- matrix(0,p,0)
    if(num.T>0){
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
  }


  if(is.null(TR)){
    if(!is.null(formula) && !is.null(X)){
      xb=as.matrix(model.matrix(formula,data = data.frame(X)))
      X=as.matrix(xb[,!(colnames(xb) %in% c("(Intercept)"))])
      colnames(X)<- colnames(xb)[!(colnames(xb) %in% c("(Intercept)"))]
      num.X <- dim(X)[2]
    }
    if(is.null(formula) && !is.null(X)){
      n1 <- colnames(X)
      formula=paste("~",n1[1],sep = "")
      if(length(n1)>1){
        for(i1 in 2:length(n1)){
          formula <- paste(formula,n1[i1],sep = "+")
        }}
      formula1=formula
      formula=formula(formula)
      xb=as.matrix(model.matrix(formula,data = data.frame(X)))
      X<-as.matrix(xb[,!(colnames(xb) %in% c("(Intercept)"))])
      colnames(X)<- colnames(xb)[!(colnames(xb) %in% c("(Intercept)"))]

      num.X=ncol(X);
    }
    Xd<-X1<-X
  } else {
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
      if(length(n2)>1){
        for(i2 in 2:length(n2)){
          formula <- paste(formula,n2[i2],sep = "+")
        }}
      formula1 <- paste(formula,")",sep = "")
      formula <- formula(formula1)
    }

    yX <- reshape(data.frame(cbind(y,X)),direction = "long", varying = colnames(y),v.names = "y")
    TR2<-data.frame(time=1:p,TR)
    if(is.null(yXT)){
      yXT=merge(yX,TR2,by="time")
    }
    data <- yXT

    m1<-model.frame(formula,data = data)
    term<-terms(m1)

    Xd <- as.matrix(model.matrix(formula,data = data))
    nXd <- colnames(Xd)
    Xd <- as.matrix(Xd[,!(nXd %in% c("(Intercept)"))])
    colnames(Xd) = nXd[!(nXd %in% c("(Intercept)"))]

    if(!is.null(X.new)) fx <- apply(matrix(sapply(colnames(X.new),function(x){grepl(x, colnames(Xd))}),ncol(Xd),ncol(X.new)),2,any)
    ft <- NULL;
    if(NCOL(T.new)>0) ft <- apply(matrix(sapply(colnames(T.new),function(x){grepl(x, colnames(Xd))}),ncol(Xd),ncol(T.new)),2,any)

    X1=as.matrix(X.new[,fx]);
    TR1=as.matrix(T.new[,ft]);
    colnames(X1) <- colnames(X.new)[fx]; colnames(TR1)=colnames(T.new)[ft];
    nxd <- colnames(Xd)
    formulab <- paste("~",nxd[1],sep = "");
    for(i in 2:length(nxd)) formulab <- paste(formulab,nxd[i],sep = "+")
    formula1 <- formulab
    nd=ncol(Xd)

  }

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
  if(is.null(formula) && is.null(X) && is.null(TR)){formula1 =  "~ 1"}

  y00=y
  if(family == "ordinal") {
    if(min(y)==0){ y=y+1}
    max.levels <- apply(y,2,function(x) length(min(x):max(x)))
    if(any(max.levels == 1) || all(max.levels == 2))
      stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks")
                        
    if(any(!apply(y,2,function(x)all(diff(sort(unique(x)))==1))))
      stop("Can't fit ordinal model if there are species with missing classes. Please reclassify per species.")                    
  }

  out.list <-  list(y = y, X = X1, TR = TR1, num.lv = num.lv, row.eff = row.eff, logLik = -Inf, family = family, offset=offset,X.design=Xd, terms=term)

  tstart<-Sys.time()
  n.i<-1;
  if(n.init>1){
    if(!is.null(seed))set.seed(seed)
    seed<-sample(1:10000,n.init)
  } 
  while(n.i<=n.init){
    if(n.init>1 && trace) cat("initial run ",n.i,"\n");
    set.seed(seed[n.i])
    xr = dr = matrix(0)
    if(row.eff == "random")dr <- diag(n)
    if(row.eff == "fixed") xr <- model.matrix(~site, data.frame(site = factor(1:n)))[,-1,drop=FALSE]
    res <- start_values_gllvm_TMB(y = y, X = X1, TR = TR1, family = family, formula = formula, offset=offset, trial.size = trial.size, num.lv = num.lv, start.lvs = start.lvs, starting.val=starting.val, jitter.var=jitter.var, yXT=yXT, TMB=FALSE, link="probit",zeta.struc="species")
    if(is.null(start.params)){
      new.beta0 <- beta0 <- res$params[,1]
      # common env params or different env response for each spp
      new.env <- env <- NULL
      new.trait <- trait <- NULL
      if(is.null(TR)) {
        if(!is.null(X)) new.env <- env <- res$params[,2:(num.X + 1)]
      }
      new.row.params <- row.params <- NULL;
      if(row.eff){
        new.row.params <- row.params <- res$row.params
        }
      new.vameans <- vameans <- new.theta <- theta <- new.lambda <- lambda <- NULL

      if(num.lv > 0) {
        new.vameans <- vameans <- res$index
        new.theta <- theta <- as.matrix(res$params[,(ncol(res$params) - num.lv + 1):ncol(res$params)])#fts$coef$theta#
        new.theta[upper.tri(new.theta)] <- theta[upper.tri(theta)] <- 0
        if(Lambda.struc == "unstructured") {
          new.lambda <- lambda <- array(NA,dim=c(n,num.lv,num.lv))
          for(i in 1:n) { new.lambda[i,,] <- lambda[i,,] <- diag(rep(0.1,num.lv)) }
        }
        if(Lambda.struc == "diagonal") {
          new.lambda <- lambda <- matrix(Lambda.start,n,num.lv)
        }
        zero.cons <- which(new.theta == 0)
      }
      new.B=B=NULL; if(!is.null(TR)) B=c(res$B)[1:nd]
      if(any(is.na(B))) B[is.na(B)]=0
    } else{
      if(dim(start.params$y)==dim(y) && is.null(X)==is.null(start.params$X) && is.null(TR)==is.null(start.params$TR) && row.eff == start.params$row.eff){
        new.beta0 <- beta0 <- start.params$params$beta0
        # common env params or different env response for each spp
        new.env <- env <- NULL
        new.trait <- trait <- NULL
        if(is.null(TR)) {
          if(!is.null(X)) new.env <- env <- start.params$params$Xcoef
        }
        new.B=B=NULL; if(!is.null(TR)) B=c(start.params$params$B)[1:nd]

        new.vameans <- vameans <- new.theta <- theta <- new.lambda <- lambda <- NULL

        if(num.lv > 0) new.theta <- theta<- c(start.params$params$theta) ## LV coefficients
        if(row.eff) new.row.params <- row.params <- start.params$params$row.params ## row parameters
        if(num.lv > 0) { new.vameans <- vameans <- matrix(start.params$lvs, ncol = num.lv);
        new.lambda <- lambda <- start.params$A}
      } else { stop("Model which is set as starting parameters isn't the suitable you are trying to fit. Check that attributes y, X, TR and row.eff match to each other.");}
    }
    if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)

    new.zeta <- zeta <- NULL; if(family == "ordinal") { new.zeta <- zeta <- res$zeta }
    new.phi <- phi <- NULL; if(family == "negative.binomial") { phis <- res$phi; if(any(phis>10))phis[phis>100]=100; if(any(phis<0.10))phis[phis<0.01]=0.01; new.phi <- phi <- 1/phis; res$phi=phis }

    if(constrOpt ){
      const=min(restrict,15)/5
      new.beta0 <- beta0<-beta0*const/max(c(abs(beta0),1))
      if(!is.null(X)) new.env <- env<-env*const/max(c(abs(env),1))
      if(num.lv>0) new.theta <- theta<-theta*const/max(c(abs(theta),1))
    }
    q=num.lv


    current.loglik <- -1e6; iter <- 1; err <- 10; div=1e5
    while((div> eps*(abs(current.loglik)+eps)) && iter <= max.iter) {
      #while((err > (1 + eps) || err < (1 - eps)) && iter <= max.iter) {

      if(trace) cat("Iteration:", iter, "\n")

      if(Lambda.struc == "unstructured" & iter <= diag.iter && num.lv>0) {
        tmp.Lambda.struc <- "diagonal"
        if(iter == 1)  new.lambda <- lambda <- matrix(Lambda.start,n,num.lv)#0.1
        if(iter == diag.iter) { tmp.Lambda.struc <- "unstructured"; new.lambda <- lambda <- lambda.convert(lambda,type=2) }
      } else { tmp.Lambda.struc <- Lambda.struc }


      ## log-likelihood
      ll0 <- function(x,v=NULL,lambda=NULL,phi,zeta=NULL) {
        x2 <- x
        new.phi <- phi
        new.lambda <- lambda
        new.vameans <- new.theta <- NULL
        if(num.lv > 0) {
          new.vameans <- matrix(v,n,num.lv);
          new.theta <- matrix(c(x2[1:(p * num.lv)]),p,num.lv); x2 <- x2[-(1:(p * num.lv))]
          new.theta[upper.tri(new.theta)] <- 0
        }
        new.zeta <- zeta
        new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
        new.env <- NULL
        if(!is.null(X)) {
          if(is.null(TR)) { new.env <- matrix(x2[1:(p * num.X)],p,num.X); x2 <- x2[-(1:(p * num.X))]
          } else { B=x2[1:nd]; x2 <- x2[-(1:nd)]}
        }
        new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }

        mu.mat <- matrix(new.beta0,n,p,byrow=TRUE) + offset
        if(!is.null(X) && is.null(TR)) mu.mat <- mu.mat + X %*% t(new.env)
        if(!is.null(TR)) mu.mat <- mu.mat + matrix((Xd %*% B),n,p)
        if(row.eff) mu.mat <- mu.mat + matrix(new.row.params,n,p,byrow=FALSE)
        if(num.lv > 0) mu.mat <- mu.mat  + new.vameans %*% t(new.theta)
        eta.mat <- mu.mat
        if(num.lv > 0) eta.mat <- eta.mat + calc.quad(new.lambda,new.theta,tmp.Lambda.struc)$mat

        if(family=="poisson") { out1 <- sum(y * mu.mat - lfactorial(y) - exp(eta.mat)) }

        if(family=="negative.binomial") {
          phi.mat <- matrix(new.phi,n,p,byrow=TRUE)
          if(num.lv > 0) eta.mat <- mu.mat - calc.quad(new.lambda,new.theta,tmp.Lambda.struc)$mat
          out1 <- sum(-phi.mat * mu.mat - lfactorial(y) + phi.mat * log(phi.mat) - lgamma(phi.mat) - (y + phi.mat) * log(1 + phi.mat / exp(eta.mat)) + lgamma(y + phi.mat))
        }

        if(family=="binomial") {
          probs<-pnorm(mu.mat);
          out1 <- dbinom(as.matrix(y),size=trial.size,prob=probs,log=TRUE)
          out1 <- sum(out1[is.finite(out1)])
          if(num.lv > 0) out1 <- out1 - calc.quad(new.lambda,new.theta,tmp.Lambda.struc)$mat.sum
        }

        if(family == "ordinal") {
          out1 <- matrix(NA,n,p)
          for(j in 1:p) {
            out1[y[,j] == 1,j] <- pnorm(new.zeta[j,1]-mu.mat[y[,j] == 1,j],log.p=TRUE)
            out1[y[,j] == max(y[,j]),j] <- log(1 - pnorm(new.zeta[j,max(y[,j]) - 1] - mu.mat[y[,j] == max(y[,j]),j]))
            if(max(y[,j]) > 2) {
              j.levels <- 2:(max(y[,j])-1)
              for(k in j.levels) { out1[y[,j] == k,j] <- log(pnorm(new.zeta[j,k] - mu.mat[y[,j] == k,j]) - pnorm(new.zeta[j,k - 1] - mu.mat[y[,j] == k,j])) }
            }
          }
          out1 <- sum(out1[is.finite(out1)])
          if(num.lv > 0) out1 <- out1 - calc.quad(new.lambda,new.theta,tmp.Lambda.struc)$mat.sum
        }

        if(num.lv > 0) {
          if(tmp.Lambda.struc == "unstructured") {
            foo2 <- function(i) { 0.5 * (log(det(new.lambda[i,,])) - sum(diag(new.lambda[i,,])) - sum(new.vameans[i,]^2)) }
          }

          if(tmp.Lambda.struc == "diagonal") {
            foo2 <- function(i) { 0.5 * (sum(log(t(new.lambda)[,i])) - sum(t(new.lambda)[,i]) - sum(new.vameans[i,]^2)) }
          }

          out1 <- out1 + sum(sapply(1:n,foo2))
        }

        return(out1)
      }

      ## Update model params by maximizing lower bound for loglik
      grad.mod <- function(x,v=NULL,lambda=NULL,phi,zeta=NULL) {
        x2 <- x
        new.phi <- phi
        new.lambda <- lambda
        new.vameans <- new.theta <- NULL
        if(num.lv > 0) {
          new.vameans <- matrix(v,n,num.lv);
          new.theta <- matrix(c(x2[1:(p * num.lv)]),p,num.lv); x2 <- x2[-(1:(p * num.lv))]
          new.theta[upper.tri(new.theta)] <- 0
        }
        new.zeta <- zeta
        new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
        new.env <- NULL
        if(!is.null(X)) {
          if(is.null(TR)) { new.env <- matrix(x2[1:(p * num.X)],p,num.X); x2 <- x2[-(1:(p * num.X))]
          } else { B=x2[1:nd]; x2 <- x2[-(1:nd)]}
        }
        new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }

        grad.theta <- grad.beta0 <- grad.env <- grad.B <- grad.row.params <- NULL

        if(tmp.Lambda.struc == "unstructured" & num.lv > 0) {
          Lambda.theta <- array(NA,dim=c(p,n,num.lv))
          for(i in 1:n) { Lambda.theta[,i,] <- new.theta %*% new.lambda[i,,] }
        }


        eta.mat <- matrix(new.beta0,n,p,byrow=TRUE) + offset
        if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
        if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% B),n,p)
        if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
        if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta)

        if(family=="poisson") {
          if(num.lv > 0) eta.mat <- eta.mat + calc.quad(new.lambda,new.theta,tmp.Lambda.struc)$mat
          grad.beta0 <- colSums(y - exp(eta.mat))
          if(num.lv > 0) {
            for(l in 1:num.lv) {
              if(tmp.Lambda.struc == "unstructured") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep(t(Lambda.theta[,,l]),1,new.vameans[,l],"+") * exp(eta.mat)
              if(tmp.Lambda.struc == "diagonal") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep((new.lambda[,l]%*%t(new.theta[,l])),1,new.vameans[,l],"+") * exp(eta.mat)
              grad.theta <- c(grad.theta,colSums(sum1))
            }
          }
          if(!is.null(X) && is.null(TR)) {
            for (l in 1:num.X) {
              sum1 <- sweep(y-exp(eta.mat),1,X[,l],"*")
              grad.env <- c(grad.env,colSums(sum1)) # else{grad.env <- c(grad.env,sum(sum1)) }
            }
          }
          if(!is.null(TR)) {
            sum1 <- c(y-exp(eta.mat)) * Xd
            grad.B <- c(grad.B,colSums(sum1))
          }
          if(row.eff)  grad.row.params <- rowSums(y - exp(eta.mat))
        }

        if(family=="negative.binomial") {
          phi.mat <- matrix(new.phi,n,p,byrow=TRUE)
          if(num.lv > 0) eta.mat <- eta.mat - calc.quad(new.lambda,new.theta,tmp.Lambda.struc)$mat
          grad.beta0 <- colSums(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(eta.mat) + phi.mat)))
          if(num.lv > 0) {
            for(l in 1:num.lv) {
              if(tmp.Lambda.struc == "unstructured") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep(t(Lambda.theta[,,l]),1,-new.vameans[,l],"+") * (-(y + phi.mat) * (phi.mat / (exp(eta.mat) + phi.mat)))
              if(tmp.Lambda.struc == "diagonal") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep((new.lambda[,l] %*% t(new.theta[,l])),1,-new.vameans[,l],"+") * (-(y + phi.mat) * (phi.mat / (exp(eta.mat) + phi.mat)))
              grad.theta <- c(grad.theta,colSums(sum1))
            }
          }
          if(!is.null(X) && is.null(TR)) {
            for(l in 1:num.X) {
              sum1 <- sweep(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(eta.mat) + phi.mat)),1,X[,l],"*")
              grad.env <- c(grad.env,colSums(sum1))
            }
          }
          if(!is.null(TR)) {
            sum1 <- c(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(eta.mat) + phi.mat))) * Xd
            grad.B <- c(grad.B,colSums(sum1))
          }
          if(row.eff)  grad.row.params <- rowSums(-phi.mat - (y + phi.mat) * (-phi.mat/(exp(eta.mat) + phi.mat)))
        }

        if(family=="binomial") {
          probs<-pnorm(eta.mat);

          grad.beta0 <- colSums(dnorm(eta.mat)*(y-trial.size*probs)/(probs*(1-probs)+1e-5),na.rm=TRUE)
          if(num.lv > 0) {
            for(l in 1:num.lv) {
              if(tmp.Lambda.struc == "unstructured") sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * probs) / (probs * (1 - probs) + 1e-5),1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
              if(tmp.Lambda.struc == "diagonal") sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * probs) / (probs * (1 - probs) + 1e-5),1,new.vameans[,l],"*") - (new.lambda[,l] %*% t(new.theta[,l]))
              grad.theta <- c(grad.theta,colSums(sum1,na.rm=TRUE))
            }
          }
          if(!is.null(X) && is.null(TR)) {
            for (l in 1:num.X) {
              sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * probs) / (probs * (1 - probs) + 1e-5),1,X[,l],"*")
              grad.env <- c(grad.env,colSums(sum1))
            }
          }
          if(!is.null(TR)) {
            sum1 <- c(dnorm(eta.mat) * (y - trial.size * probs) / (probs * (1 - probs) + 1e-5)) * Xd
            grad.B <- c(grad.B,colSums(sum1))
          }
          if(row.eff)  grad.row.params <- rowSums(dnorm(eta.mat) * (y - trial.size * probs) / (probs * (1 - probs) + 1e-5),na.rm=TRUE)
        }

        if(family=="ordinal") {
          eta.mat <- matrix(new.beta0,n,p,byrow=TRUE)  + offset
          if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
          if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
          if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd),n,p)
          if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta)
          deriv.trunnorm <- matrix(0,n,p)

          for(j in 1:p) {
            deriv.trunnorm[y[,j] == 1,j] <- -dnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j]) / pnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j])
            deriv.trunnorm[y[,j] == max(y[,j]),j] <- dnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]) / (1 - pnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]))
            if(max(y[,j]) > 2) {
              j.levels <- 2:(max(y[,j]) - 1)
              for(k in j.levels) { deriv.trunnorm[y[,j] == k,j] <- (-dnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) + dnorm(new.zeta[j,k - 1] - eta.mat[y[,j] == k,j])) / (pnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) - pnorm(new.zeta[j,k - 1] - eta.mat[y[,j] == k,j])) }                                 }
          }
          deriv.trunnorm[!is.finite(deriv.trunnorm)] <- 0
          grad.beta0 <- colSums(deriv.trunnorm)

          if(num.lv > 0) {
            for(l in 1:num.lv) {
              if(tmp.Lambda.struc == "unstructured") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
              if(tmp.Lambda.struc == "diagonal") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - (new.lambda[,l]%*%t(new.theta[,l]))
              grad.theta <- c(grad.theta,colSums(sum1,na.rm=TRUE))
            }
          }
          if(!is.null(X) && is.null(TR)) {
            for(l in 1:num.X) { sum1 <- sweep(deriv.trunnorm,1,X[,l],"*")
            grad.env <- c(grad.env,colSums(sum1,na.rm=TRUE))
            }
          }
          if(!is.null(TR)) {
            sum1 <- c(deriv.trunnorm) * Xd
            grad.B <- c(grad.B,colSums(sum1))
          }
          if(row.eff) { grad.row.params <- rowSums(deriv.trunnorm,na.rm=TRUE) }
        }

        if(num.lv>0){
          grad.theta=matrix(grad.theta,nrow=p,ncol=num.lv)
          grad.theta[upper.tri(grad.theta)] <- 0}
        if(row.eff) grad.row.params[1]=0;
        return(c(grad.theta, grad.beta0, grad.env, grad.B, grad.row.params))
      }


      if(constrOpt) {
        A=NULL;
        if(row.eff){for(i in 1:p){
          A=rbind(A,diag(n))
        }}

        B0<-matrix(1,n,1)
        B1=B0
        for(i in 1:(p-1)){
          B1=as.matrix(Matrix::bdiag(B1,B0))
        }

        B2=NULL;
        if(!is.null(X) && is.null(TR)){
          for(s in 1:num.X){
            B2.=X[,s]
            B0=X[,s]
            for(i in 1:(p-1)){
              B2.=as.matrix(Matrix::bdiag(B2.,B0))
            }
            B2=cbind(B2,B2.)
          }
        }

        G=NULL;
        if(num.lv>0){for(s in 1:num.lv){
          G0=new.vameans[,s]
          Gi=new.vameans[,s]
          for(i in 1:(p-1)){
            Gi=as.matrix(Matrix::bdiag(Gi,G0))
          }
          G=cbind(G,Gi)
        }}

        Bf=NULL;
        if(!is.null(TR)){
          Bf=Xd
        }

        U=cbind(G,B1,B2,Bf,A)

        q <- try(constrOptim(c(theta,beta0,env,B,row.params), v = c(new.vameans), lambda = new.lambda, phi = phi, zeta = zeta, method = "BFGS", f = ll0, grad = grad.mod, control = list(trace = 0, fnscale = -1, maxit = maxit),ui=rbind(U,-U),ci=rep(-restrict,NROW(U)*2)), silent = TRUE)
      } else {
        q <- try(optim(c(theta,beta0,env,B,row.params), v = c(new.vameans), lambda = new.lambda, phi = phi, zeta = zeta, method = "BFGS", fn = ll0, gr = grad.mod, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      }


      if(!inherits(q, "try-error")) {
        if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when trying to update model parameters on iteration step ", iter,"\n") }
        if(iter > 1 && current.loglik > q$value) { if(trace) cat("Optimization did not improve estimates of model parameters on iteration step ",iter,"\n"); new.theta <- theta; new.beta0 <- beta0; new.env <- env; new.row.params <- row.params}
        else {
          if(trace) cat("Model parameters updated", "\n")
          x2 <- q$par; #plot(q$par)
          if(num.lv > 0) {
            new.theta <- matrix(x2[1:(p * num.lv)],p,num.lv); x2 <- x2[-(1:(p * num.lv))];
            new.theta[upper.tri(new.theta)] <- 0
            if(starting.val=="zero" && iter==1)
              new.theta[!upper.tri(new.theta)]=1
          }
          new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]; #plot(new.beta0)
          new.env <- NULL
          if(!is.null(X)) {
            if(is.null(TR)) { new.env <- matrix(x2[1:(p * num.X)],p,num.X); x2 <- x2[-(1:(p * num.X))] }
            else {
              new.B <- x2[1:nd]; x2 <- x2[-(1:nd)]
            }
          }
          new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
        }
      }
      else { new.theta <- theta; new.beta0 <- beta0; new.env <- env; new.B <- B; new.row.params <- row.params }



      # Update dispersion parameters for NB distribution
      if(family=="negative.binomial") {
        grad.phi <- function(x,v=NULL,lambda=NULL,phi,zeta = NULL) {
          x2 <- x; new.phi <- phi; new.lambda <- lambda; new.vameans <- new.theta <- NULL
          if(num.lv > 0) {
            new.vameans <- matrix(v,n,num.lv)
            new.theta <- matrix(c(x2[1:(p * num.lv)]),p,num.lv); x2 <- x2[-(1:(p * num.lv))]
          }
          new.zeta <- zeta; new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]; new.env <- NULL
          if(!is.null(X)) {
            if(is.null(TR)) { new.env <- matrix(x2[1:(p * num.X)],p,num.X); x2 <- x2[-(1:(p * num.X))]
            } else { B=x2[1:nd]; x2 <- x2[-(1:nd)]}
          }
          new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }

          phi.mat <- matrix(new.phi,n,p,byrow=TRUE)

          #grad.theta <- grad.beta0 <- grad.env <- grad.row.params <- NULL

          if(tmp.Lambda.struc == "unstructured" & num.lv > 0) {
            Lambda.theta <- array(NA,dim=c(p,n,num.lv));
            for(i in 1:n) { Lambda.theta[,i,] <- new.theta%*%new.lambda[i,,] }
          }

          eta.mat <- matrix(new.beta0,n,p,byrow=TRUE)  + offset
          if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
          if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% B),n,p)
          if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
          if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta)
          mu.mat <- eta.mat
          if(num.lv > 0) mu.mat <- eta.mat - calc.quad(new.lambda,new.theta,tmp.Lambda.struc)$mat

          out <- colSums(-eta.mat + 1 + log(phi.mat) - digamma(phi.mat) - log(1 + phi.mat / exp(mu.mat)) - (y + phi.mat) / (phi.mat + exp(mu.mat)) + digamma(y + phi.mat))
          return(out)
        }

        q <- try(optim(phi, x = c(new.theta,new.beta0,new.env,new.B,new.row.params), v = c(new.vameans), lambda = new.lambda, zeta = zeta, method = "L-BFGS-B", lower = 1e-4, upper = 1e4, fn = ll0, gr = grad.phi, control = list(trace=0, fnscale=-1)), silent = TRUE) #, factr = 1e-3,maxit=maxit

        if(!inherits(q, "try-error")) {
          if(iter > 1 && current.loglik > q$value) { if(trace) cat("Optimization did not improve estimates of dispersion parameters on iteration step ",iter,"\n"); new.phi <- phi }
          else {
            if(trace) cat("Dispersion parameters updated", "\n")
            new.phi <- q$par[1:p]
          }
        }
        else { new.phi <- phi }
      }

      ## Update zeta (cutoffs for ordinal data)
      if(family == "ordinal") {
        eta.mat <- matrix(new.beta0,n,p,byrow=TRUE) + offset
        if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
        if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% new.B),n,p)
        if(!is.null(row.params)) eta.mat <- eta.mat + matrix(new.row.params,n,p)
        if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta)

        func.zetaj <- function(cw.zeta, j) { ## Exclude the first column of zeta, which is fixed at 0
          zeta0 <- c(0,cw.zeta); out <- 0
          out <- out + sum(pnorm(zeta0[1]-eta.mat[which(y[,j] == 1),j],log.p=TRUE))
          out <- out + sum(log(1 - pnorm(zeta0[max(y[,j]) - 1]-eta.mat[which(y[,j] == max(y[,j])),j])))
          if(max(y[,j])>2) {
            j.levels <- 2:(max(y[,j]) - 1)
            for(k in j.levels) { out <- out + sum(log(pnorm(zeta0[k] - eta.mat[y[,j] == k,j]) - pnorm(zeta0[k - 1] - eta.mat[y[,j] == k,j]))) }
          }
          out
        }

        grad.zetaj <- function(cw.zeta, j) { ## Exclude the first column of zeta, which is fixed at 0
          zeta0 <- c(0,cw.zeta); deriv.trunnorm <- numeric(length(zeta0)); ## L-1 length
          deriv.trunnorm[length(deriv.trunnorm)] <- deriv.trunnorm[length(deriv.trunnorm)] - sum(dnorm(zeta0[max(y[,j]) - 1] - eta.mat[which(y[,j] == max(y[,j])),j]) / (1 - pnorm(zeta0[max(y[,j]) - 1] - eta.mat[which(y[,j] == max(y[,j])),j])))
          if(max(y[,j])>2) {
            j.levels <- 2:(max(y[,j]) - 1)
            for(k in j.levels) {
              deriv.trunnorm[k] <- deriv.trunnorm[k] + sum(dnorm(zeta0[k] - eta.mat[y[,j] == k,j]) / (pnorm(zeta0[k] - eta.mat[y[,j] == k,j]) - pnorm(zeta0[k - 1] - eta.mat[y[,j] == k,j])),na.rm=TRUE)
              deriv.trunnorm[k - 1] <- deriv.trunnorm[k - 1] - sum(dnorm(zeta0[k - 1]-eta.mat[y[,j] == k,j]) / (pnorm(zeta0[k] - eta.mat[y[,j] == k,j]) - pnorm(zeta0[k - 1] - eta.mat[y[,j] == k,j])),na.rm = TRUE)
            }
          }
          deriv.trunnorm[-1]
        }

        for(j in 1:p) {
          if(max(y[,j]) == 2) { new.zeta[j,] <- zeta[j,] }
          if(max(y[,j]) > 2) {
            constraint.mat <- matrix(0,max(y[,j])-2,max(y[,j])-2); constraint.mat[1,1] <- 1 ## Constrains zeta2 > 0
            if(nrow(constraint.mat) > 1) {
              for(k in 2:nrow(constraint.mat)) constraint.mat[k,(k-1):k] <- c(-1,1)
            }

            update.zeta <- constrOptim(theta = zeta[j,2:(max(y[,j])-1)], f = func.zetaj, grad = grad.zetaj, ui = constraint.mat, ci = rep(0,max(y[,j])-2), j = j, outer.eps = 1e-3, control = list(trace=0, fnscale = -1))
            if(!inherits(update.zeta, "try-error")) new.zeta[j,2:(max(y[,j])-1)] <- update.zeta$par
            if(inherits(update.zeta, "try-error")) new.zeta[j,] <- zeta[j,]
          }
        }
        if(trace) cat("Cutoffs updated \n")
      }


      if(num.lv > 0) {
        ## Update vameans by maximizing lower bound for loglik
        grad.var <- function(x,v=NULL,lambda=NULL,phi,zeta=NULL) {
          x2 <- x
          new.phi <- phi
          new.lambda <- lambda
          new.vameans <- new.theta <- NULL
          if(num.lv > 0) {
            new.vameans <- matrix(v,n,num.lv);
            new.theta <- matrix(c(x2[1:(p * num.lv)]),p,num.lv); x2 <- x2[-(1:(p * num.lv))]
          }
          new.zeta <- zeta
          new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
          new.env <- NULL
          if(!is.null(X)) {
            if(is.null(TR)) { new.env <- matrix(x2[1:(p * num.X)],p,num.X); x2 <- x2[-(1:(p * num.X))]
            } else { B=x2[1:nd]; x2 <- x2[-(1:nd)]}
          }
          new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }

          grad.vameans <- grad.lambda <- NULL

          eta.mat <- matrix(new.beta0,n,p,byrow=TRUE)  + offset
          if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
          if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% B),n,p)
          if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
          if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta)

          if(family=="poisson") {
            eta.mat <- eta.mat + calc.quad(new.lambda,new.theta,tmp.Lambda.struc)$mat
            sum1 <- (y-exp(eta.mat))
            for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,rowSums(sweep(sum1,2,new.theta[,l],"*"))-new.vameans[,l]) }
          }

          if(family=="negative.binomial") {
            phi.mat <- matrix(new.phi,n,p,byrow=TRUE)
            eta.mat <- eta.mat - calc.quad(new.lambda,new.theta,tmp.Lambda.struc)$mat
            sum1 <- -phi.mat - (y + phi.mat) * (-phi.mat / (exp(eta.mat) + phi.mat))
            for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,rowSums(sweep(sum1,2,new.theta[,l],"*"))-new.vameans[,l]) }
          }

          if(family=="binomial") {
            probs<-pnorm(eta.mat);
            sum1 <- dnorm(eta.mat) * (y - trial.size * probs) / (probs * (1-probs) + 1e-5)
            for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,rowSums(sweep(sum1,2,new.theta[,l],"*"),na.rm=TRUE)-new.vameans[,l]) }
          }

          if(family=="ordinal") {
            deriv.trunnorm <- matrix(NA,n,p)
            for(j in 1:p) {
              deriv.trunnorm[y[,j] == 1,j] <- -dnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j]) / pnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j])
              deriv.trunnorm[y[,j] == max(y[,j]),j] <- dnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]) / (1 - pnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]))
              if(max(y[,j]) > 2) {
                j.levels <- 2:(max(y[,j])-1)
                for(k in j.levels) { deriv.trunnorm[y[,j] == k,j] <- (-dnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) + dnorm(new.zeta[j,k-1] - eta.mat[y[,j] == k,j])) / (pnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) - pnorm(new.zeta[j,k-1] - eta.mat[y[,j] == k,j])) }
              }
            }
            deriv.trunnorm[!is.finite(deriv.trunnorm)] <- 0

            for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,rowSums(sweep(deriv.trunnorm,2,new.theta[,l],"*")) - new.vameans[,l]) }
          }

          return(c(grad.vameans))
        }

        if(constrOpt)  {
          U=NULL
          for(d in 1:num.lv){
            U=cbind(U,kronecker(new.theta[,d],diag(n)))
          }
          etas <- matrix(new.beta0,n,p,byrow=TRUE)  + offset
          if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
          if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% new.B),n,p)
          if(row.eff) etas <- etas + matrix(new.row.params,n,p)
          ci=c(-restrict-c(etas),-restrict+c(etas))
          q <- constrOptim(c(vameans), method = "BFGS", f = ll0, grad = grad.var, x = c(new.theta,new.beta0,new.env,new.B,new.row.params), lambda = lambda, phi = new.phi, zeta = zeta, control = list(trace = 0, fnscale = -1, maxit = maxit),ui=rbind(U,-U),ci=ci)
        } else {
          q <- try(optim(c(vameans), method = "BFGS", fn = ll0, gr = grad.var, x = c(new.theta,new.beta0,new.env,new.B,new.row.params), lambda = lambda, phi = new.phi, zeta = zeta, control = list(trace = 0, fnscale = -1, maxit = maxit,reltol=1e-6)), silent = TRUE)
        }
        if(!inherits(q, "try-error")) {
          if(q$convergence != 0) { if(trace) cat("Optimization algorithm did not converge when it tried to update variational parameters on step ", iter,"\n") }
          if(iter > 1 && current.loglik > q$value) { if(trace) cat("Optimization did not improve estimates of variational parameters.\n"); new.vameans <- vameans }
          else {
            if(trace) cat("Variational parameters updated", "\n")
            new.vameans <- q$par[1:(num.lv*n)]; new.vameans <- matrix(new.vameans,n,num.lv)
          }
        }
        else { new.vameans <- vameans }

        ## Update Lambda_i via fixed-point algorithm (covariance matrix for VA distribution)
        eta.mat <- matrix(new.beta0,n,p,byrow=TRUE) + new.vameans %*% t(new.theta)   + offset
        if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
        if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% new.B),n,p)
        if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
        if(family == c("poisson")) mu.mat <- exp(eta.mat+calc.quad(lambda,new.theta,tmp.Lambda.struc)$mat)
        if(family == c("negative.binomial")) {
          phi.mat <- matrix(new.phi,n,p,byrow=TRUE)
          eta.mat <- eta.mat -  calc.quad(lambda,new.theta,tmp.Lambda.struc)$mat
          mu.mat <- (y + phi.mat) * (phi.mat / (exp(eta.mat) + phi.mat))
        }

        for(i in 1:n) {
          error <- 1; lambda.iter <- 0
          if(tmp.Lambda.struc == "unstructured") new.lambda.mat <- lambda[i,,]
          if(tmp.Lambda.struc == "diagonal") new.lambda.mat <- diag(x=lambda[i,],nrow=num.lv)
          while(error > 1e-2 && lambda.iter < 100) {

            cw.lambda.mat <- new.lambda.mat

            if(tmp.Lambda.struc == "unstructured") {
              theta2 <- sapply(1:p,function(j,theta) theta[j,] %*% t(theta[j,]), theta = new.theta)
              theta2 <- t(theta2)
              if(family %in% c("poisson","negative.binomial"))  new.lambda.mat <- solve(diag(rep(1,num.lv)) + matrix(apply(mu.mat[i,]*theta2,2,sum),nrow=num.lv))
              if(family %in% c("binomial","ordinal"))  new.lambda.mat <- solve(diag(rep(1,num.lv)) + matrix(apply(theta2,2,sum),nrow=num.lv))
            }

            if(tmp.Lambda.struc == "diagonal") {
              theta2 <- new.theta^2
              if(family %in% c("poisson","negative.binomial"))  new.lambda.mat <- solve(diag(rep(1,num.lv)) + diag(apply(mu.mat[i,]*theta2,2,sum),num.lv,num.lv))
              if(family %in% c("binomial","ordinal"))  new.lambda.mat <- solve(diag(rep(1,num.lv)) + diag(apply(theta2,2,sum),num.lv,num.lv))
            }

            error <-  sum((new.lambda.mat - cw.lambda.mat)^2)
            if(tmp.Lambda.struc == "unstructured") lambda[i,,] <- new.lambda.mat
            if(tmp.Lambda.struc == "diagonal") lambda[i,] <- diag(new.lambda.mat)

            eta.mat <- matrix(new.beta0,n,p,byrow=TRUE) + offset + new.vameans %*% t(new.theta)
            if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
            if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% new.B),n,p)
            if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
            if(family == c("poisson")) mu.mat <- exp(eta.mat+calc.quad(lambda,new.theta,tmp.Lambda.struc)$mat)
            if(family == c("negative.binomial")) {
              phi.mat <- matrix(new.phi,n,p,byrow=TRUE)
              eta.mat <- eta.mat -  calc.quad(lambda,new.theta,tmp.Lambda.struc)$mat
              mu.mat <- (y + phi.mat) * (phi.mat / (exp(eta.mat) + phi.mat))
            }

            lambda.iter <- lambda.iter + 1
          }
          if(tmp.Lambda.struc == "unstructured") new.lambda[i,,] <- new.lambda.mat
          if(tmp.Lambda.struc == "diagonal") new.lambda[i,] <- diag(new.lambda.mat)

          if((family %in% c("binomial","ordinal")) & i == 1) break; ## Since same Lambda matrix for all i, then do one then finish
        }

        if(family %in% c("binomial","ordinal")) {
          if(tmp.Lambda.struc == "diagonal") for(i2 in 2:n) new.lambda[i2,] <- new.lambda[1,]
          if(tmp.Lambda.struc == "unstructured") for(i2 in 2:n) new.lambda[i2,,] <- new.lambda[1,,]
        }
      }


      ## Take values of loglik from optim -function to define stopping rule
      q <- list(value = ll0(c(new.theta,new.beta0,new.env,new.B,new.row.params), v = c(new.vameans), lambda = new.lambda, phi = new.phi, zeta = new.zeta))
      new.loglik <- q$value
      div=abs(new.loglik-current.loglik)
      err <- abs(new.loglik/current.loglik);
      if(trace) cat("New Loglik:", new.loglik,"Current Loglik:", current.loglik, "Ratio", err, ". Difference in log-likelihoods:",div,"\n")

      ## Plot new ordination points for spp's and sites
      if( num.lv <= 2 && num.lv > 0 && plot == TRUE) {
        par(mfrow=c(1,2))
        plot(new.vameans[!duplicated(new.vameans),], type = "n", xlab = ifelse(num.lv == 1, "Row index (unique elements only)", "LV1"),ylab = "LV2", main = "Ordination of rows")
        text(new.vameans[!duplicated(new.vameans),],labels = which(!duplicated(new.vameans) == 1))#,col=color)

        plot(new.theta, type="n", xlab=ifelse(num.lv==1, "Column index", "Coefs for LV1"), ylab="Coefs for LV2", main="Ordination of columns")
        text(new.theta,labels = seq(1,dim(y)[2]))
      }

      current.loglik <- new.loglik
      beta0 <- new.beta0; env <- new.env; theta <- new.theta;  phi <- new.phi; B<-new.B
      vameans <- new.vameans; lambda <- new.lambda; row.params <- new.row.params; zeta <- new.zeta;
      iter <- iter + 1
    }

    tstop <- Sys.time(); time <- difftime(tstop,tstart)




    ## Bling up the output  ##
    if(n.i==1 || out.list$logL < current.loglik){
      out.list$logLik<-out.list$logL<-current.loglik
      out.list$iter=iter-1
      out.list$Lambda.struc = tmp.Lambda.struc

      if(num.lv > 0) {
        rownames(theta) <- colnames(y); colnames(theta) <- paste("LV",1:num.lv,sep = "")
        theta[upper.tri(theta)] <- 0
        rownames(vameans) <- rownames(y); colnames(vameans) <- 1:num.lv
        out.list$Lambda <- lambda; out.list$lvs <- vameans; out.list$coef$theta <- theta;
        out.list$lvs[!is.finite(out.list$lvs)] <- 0
        out.list$coef$theta[!is.finite(out.list$coef$theta)] <- 0
      }

      env.coefs <- NULL
      spp.coefs <- beta0;
      out.list$coef$beta0=spp.coefs
      if(!is.null(X)) {
        if(is.null(TR)) {
          spp.coefs <- cbind(env)
          colnames(spp.coefs) <- colnames(X); rownames(spp.coefs) <- colnames(y)
          out.list$coef$Xcoef=spp.coefs
        } else {
          names(B) <- colnames(Xd)
          out.list$coef$B=B
        }
      }
      out.list$start=res

      if(row.eff) {
        names(row.params) <- rownames(y); out.list$coef$row.params <- row.params
      }
      if(family == "negative.binomial") { out.list$coef$phi <- 1/phi; names(out.list$coef$phi) <- colnames(y); } ## Flip it back so that it is V = mu + phi * mu^2

      if(family == "ordinal") {
        out.list$coef$zeta <- zeta; rownames(out.list$coef$zeta) <- colnames(y);
        colnames(out.list$coef$zeta) <- paste(min(y00):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
      }}
    n.i=n.i+1;
  }

  if(sd.errors) {
    if(trace) cat("Calculating standard errors for parameters...\n")
    tr <- try({get.sds <- calc.infomat(theta=out.list$coef$theta, beta0=out.list$coef$beta0, env=out.list$coef$Xcoef, row.params=out.list$coef$row.params, vameans=out.list$lvs, lambda=out.list$Lambda, phi=1/out.list$coef$phi, zeta=out.list$coef$zeta, num.lv = num.lv, family = family, Lambda.struc = tmp.Lambda.struc, row.eff = row.eff, y = y, X = X, TR = TR, Xd=Xd,offset=offset, B = out.list$coef$B)#, formula=formula
    out.list$sd <- get.sds})
    if(inherits(tr, "try-error")) { cat("Standard errors for parameters could not be calculated.\n") }
    if(family == "ordinal") {
      colnames(out.list$sd$zeta) <- paste(min(y00):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
    }
  }
  if(is.null(formula1)){ out.list$formula=formula} else {out.list$formula=formula1}

  return(out.list)
}
