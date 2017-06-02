########################################################################################
## Tweedie GLLVM
## Original Author: Jenni Niku
## 1<power<2
# E(Y)=mu
# Var(Y)=phi*mu^power
########################################################################################

#library(fishMod) # for tweedie-distribution function
#library(statmod) # for initial values fit

gllvm.tweedie<-function(y,X=NULL,num.lv=2,row.eff = TRUE,max.iter = 100, eps = 1e-4, seed = NULL,maxit = 1000, start.lvs = NULL,offset=NULL,fixed.power=FALSE,Power=1.5,n.init=1,plot=TRUE,trace=TRUE,info=TRUE){
  n<-dim(y)[1]
  p<-dim(y)[2]
  num.X=0
  if(!is.null(X)) {X<-as.matrix(X); num.X=dim(X)[2]}

  num.lv <- num.lv
  y <- as.matrix(y)
  if(is.null(rownames(y))) rownames(y) <- paste("Row",1:n,sep="")
  if(is.null(colnames(y))) colnames(y) <- paste("Col",1:p,sep="")
  ## Set initial values for model parameters (including dispersion prm) and latent variables
  if(!is.null(seed)) set.seed(seed)

  n.i<-1;
  out <- list(y=y,X=X, logL = -Inf)
  if(n.init>1) seed<-sample(1:10000,n.init)
  while(n.i<=n.init){
    if(n.init>1 && trace) cat("initial run ",n.i,"\n");

    fit <- start.values.gllvm(y, X=X, num.lv=num.lv, family = "tweedie", start.lvs=start.lvs,power = Power)

    beta0 <- fit$params[,1] ## column intercepts
    betas <- NULL; if(!is.null(X)) betas <- c(fit$params[,-(1:(num.lv + 1))]) ## covariates coefficients
    lambdas <- NULL;
    if(num.lv > 0) {
      lambdas <- as.matrix(fit$params[,2:(num.lv+1)]) ## LV coefficients
      if(num.lv > 0) lambdas[upper.tri(lambdas)] <- 0
      covM.lvs<- array(NA,dim=c(n,num.lv,num.lv))
    }
    row.params <- NULL; if(row.eff) row.params <- rep(0,n) ## row parameters
    lvs <- NULL; if(num.lv > 0) lvs <- matrix(fit$index, ncol = num.lv); ## LVs

    if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)

    pp=Power # Power parameter
    phis=fit$phi

    new.params <- params <- c(row.params, beta0, betas, lambdas)
    if(any(abs(params)>100)) params[abs(params)>100]=sign(params[abs(params)>100])*100

    current.loglik <- -1e6; iter <- 1; err <- 10;
    while((err > (1 + eps) || err < (1 - eps)) && iter <= max.iter) {

      ll <- function(params, lvs, y, phis, pp,lw) {
        if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
        if(!row.eff) { row.params <- 0; x2 <- params }
        beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
        betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
        if(num.lv > 0) {
          lambdas.mat <- matrix(x2, ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
          LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x) ))
          if(num.lv == 1) LambdaLambdaT <- t(LambdaLambdaT)
        }

        eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params  + offset
        if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
        if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas, nrow = p))
        I.numlv <- diag(x = num.lv)
        SGam <- 0

        phis.mat<-matrix(phis,n,p,byrow=TRUE)

        if(num.lv > 0) {
          V <- 1/phis.mat*((2-pp)*exp((2-pp)*eta.mat)-y*(1-pp)*exp((1-pp)*eta.mat))
          G<-t(t(V %*% LambdaLambdaT) + c(I.numlv))
          if(num.lv == 1) SGam <- SGam - sum(0.5 * log(G))
          if(num.lv == 2) SGam <- SGam - sum(0.5 * log(G[,1] * G[,4] - G[,2] * G[,3]))
          if(num.lv > 2) { for(i in 1:n) {
            Ga <- matrix(G[i,], nrow = num.lv, ncol = num.lv)
            SGam <- SGam - 0.5 * log(det(Ga))}
          }
        }
        if(is.null(lw)){
          fun<-function(j){fishMod::dTweedie(y[,j], mu = exp(eta.mat[,j]), phi = phis[j], p = pp)} ## derivatives of the log density w.r.t. phi
          SGam<-SGam +sum(sapply(1:p,fun))

        } else{ # log W not evaluated:
          SGam<-SGam + sum( 1/phis.mat*(y*exp(eta.mat)^(1-pp)/(1-pp) - exp(eta.mat)^(2-pp)/(2-pp)) )
        }

        if(num.lv > 0) SGam <- SGam - 0.5 * sum(diag(lvs %*% t(lvs)))
        return(SGam)
      }


      ## Return n x num.lv^2 matrix with each row as inverse of Gamma(\theta,z_i) matrix
      ## as in eq. 9 in Huber et al.
      # Hub.Gamma.tw(y, eta.mat, lambdas.mat, phis, pp)
      Hub.Gamma.tw <- function(y, eta, lambdas.mat, phis, pp,LambdaLambdaT) {
        I.numlv <- diag(x=num.lv)
        V <- 1/matrix(phis, n, p, byrow = TRUE)*((2-pp)*exp((2-pp)*eta)-y*(1-pp)*exp((1-pp)*eta))
        s <- matrix(0, n, num.lv^2)
        for(i in 1:n) s[i,] <- c(solve(matrix(colSums(V[i,] * LambdaLambdaT), num.lv, num.lv) + I.numlv))
        return(s)
      }

      param.gradient<-function(params, lvs, y, phis, pp,lw){
        if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
        if(!row.eff) { row.params <- 0; x2 <- params }
        beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
        betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
        if(num.lv > 0) {
          lambdas.mat <- matrix(x2, ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
          LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x) ))
          if(num.lv == 1) LambdaLambdaT <- t(LambdaLambdaT)
        }

        eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params  + offset
        if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
        if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas, nrow = p))

        grad.beta0 <- numeric(p);
        grad.lambdas <- NULL; if(num.lv > 0) grad.lambdas <- matrix(0, p, num.lv)
        grad.row.params <- NULL; if(row.eff) grad.row.params <- numeric(n)
        grad.betas <- NULL; if(!is.null(X)) grad.betas <- matrix(0, p, num.X)

        if(num.lv > 0) {
          GA <- Hub.Gamma.tw(y = y, eta = eta.mat, lambdas.mat, phis,pp,LambdaLambdaT)
          Kron <- list()
          es<-diag(num.lv)
          for(s in 1:num.lv) {
            Kron[[s]]<-apply(lambdas.mat, 1, function(x) kronecker(es[s,],t(x))) + apply(lambdas.mat, 1, function(x) kronecker(t(es[s,]),x))
          }
          if(num.lv == 1) Kron[[1]] <- matrix(Kron[[1]],nrow=1)
        }


        V <- 1/matrix(phis, n, p, byrow = TRUE) * ((2-pp)^2*exp((2-pp)*eta.mat) - y*(1-pp)^2*exp((1-pp)*eta.mat))
        V2 <- 1/matrix(phis, n, p, byrow = TRUE) * ((2-pp)*exp((2-pp)*eta.mat) - y*(1-pp)*exp((1-pp)*eta.mat))## variance
        V1 <- 1/matrix(phis, n, p, byrow = TRUE) * (y*exp((1-pp)*eta.mat)-exp((2-pp)*eta.mat))
        if(row.eff) {grad.row.params <- rowSums(V1); if(num.lv>0) grad.row.params <-  grad.row.params-0.5 * rowSums(GA %*% t(LambdaLambdaT) * V);}
        grad.beta0 <- colSums(V1); if(num.lv>0) grad.beta0 <- grad.beta0 -0.5 * colSums(GA %*% t(LambdaLambdaT) * V)
        if(!is.null(X)) { grad.betas <- t(V1) %*% X; if(num.lv>0) grad.betas <- grad.betas -0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X}
        if(num.lv > 0) {
          for(s in 1:num.lv){ grad.lambdas[,s] <- -0.5 * colSums((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s], n, p, byrow = FALSE) + GA %*% Kron[[s]] * V2)) + colSums(V1 * lvs[,s]) }
        }

        if(num.lv > 0) grad.lambdas[upper.tri(grad.lambdas)] <- 0
        score.out <- c(grad.row.params, grad.beta0, grad.betas, grad.lambdas)
        return(score.out)
      }

      phi.gradient<-function(params, lvs, y, phis, pp,lw){
        if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
        if(!row.eff) { row.params <- 0; x2 <- params }
        beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
        betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
        if(num.lv > 0) {
          lambdas.mat <- matrix(x2, ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
          LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x) ))
          if(num.lv == 1) LambdaLambdaT <- t(LambdaLambdaT)
        }

        eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params  + offset
        if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
        if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas, nrow = p))
        grad.phis <- numeric(p);

        if(num.lv > 0) {
          GA <- Hub.Gamma.tw(y = y, eta = eta.mat, lambdas.mat, phis,pp,LambdaLambdaT)
          Kron <- list()
          es<-diag(num.lv)
          for(s in 1:num.lv) {
            Kron[[s]]<-apply(lambdas.mat, 1, function(x) kronecker(es[s,],t(x))) + apply(lambdas.mat, 1, function(x) kronecker(t(es[s,]),x))
          }
          if(num.lv == 1) Kron[[1]] <- matrix(Kron[[1]],nrow=1)
        }

        V2 <- 1/matrix(phis^2, n, p, byrow = TRUE) * (y*(1-pp)*exp((1-pp)*eta.mat) - (2-pp)*exp((2-pp)*eta.mat))## variance
        fun<-function(j){nd2(x0=phis[j],y=y[,j],mu=exp(eta.mat[,j]), p = pp,f=fishMod::dTweedie)} ## numeric derivatives, doesn't mgcv package but fishmod
        grad.phis <-colSums(sapply(1:p,fun))
        if(num.lv>0) grad.phis <- grad.phis -0.5 * colSums(GA %*% t(LambdaLambdaT) * V2)

        score.out <- c(grad.phis)
        return(score.out)
      }




      # Update model parameters
      system.time(update.params <- try(optim(params, fn = ll, gr = param.gradient, method = "BFGS", hessian = FALSE, control = list(trace = 0, fnscale = -1, maxit=maxit), lvs = lvs, y = y, phis = phis, pp = pp,lw=FALSE), silent = TRUE))
      if(!inherits( update.params, "try-error")) {
        new.params=update.params$par
        if(trace) cat("Model parameters updated", "\n")
      }

      # Update dispersion parameters
      system.time(update.phis <- try(optim(phis, fn = ll, gr = phi.gradient, method = "L-BFGS-B", control = list(trace = 0, fnscale = -1, maxit=maxit), lvs = lvs, y = y, params=new.params, pp = pp,lw=NULL,lower = rep(0.03,p),upper = rep(1e3,p)), silent = TRUE))
      if(!inherits( update.phis, "try-error")) {
        phis=update.phis$par
        if(trace) cat("Dispersion parameters updated", "\n")
      }


      if(!fixed.power){
        system.time(update.pow <- optim(pp, fn = ll, gr = NULL, method = "Brent",lower = 1.01,upper = 1.99, hessian = FALSE, control = list(trace = 0, fnscale = -1, maxit=maxit), lvs = lvs, y = y, params=new.params, phis=phis,lw=NULL))
        if(!inherits( update.pow, "try-error")) {
          pp=update.pow$par
          if(trace) {cat("Power parameter updated", "\n"); cat("New Power", pp, "\n")}}
      }

      if(num.lv > 0) {
        ll.i <- function(lvs.i, params, y, phis, pp, i) {
          if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
          if(!row.eff) { row.params <- rep(0,n); x2 <- params }
          beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
          betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
          lambdas.mat <- matrix(x2, ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0

          eta.i <- beta0  + lambdas.mat %*% lvs.i  + offset[i,]
          if(row.eff) eta.i <- eta.i + row.params[i]
          if(!is.null(X)) eta.i <- eta.i + matrix(betas,nrow = p) %*% X[i,]

          SGam <-sum(fishMod::dTweedie(y[i,], mu=c(exp(eta.i)), phi=phis, p=pp))
          SGam <- SGam - 0.5 * t(lvs.i) %*% lvs.i
          return(SGam)
        }

        index <- matrix(0,n,num.lv)
        for(i in 1:n) {
          update.lvs <- try(optim(lvs[i,], fn = ll.i, gr = NULL, method = "BFGS", control = list(trace = 0, fnscale = -1, maxit = maxit,reltol=1e-5), y = y, params = new.params, phis = phis, pp=pp, i = i, hessian = FALSE), silent = TRUE)
          if(!inherits(update.lvs, "try-error")){
            index[i,] <- update.lvs$par
          } else index[i,]<-lvs[i,]
        }
        if(trace) cat("Latent variables updated", "\n")

        ## Plot
        if(plot){
          plot(rbind(index), type = "n", main = "Ordination points", xlab = "LV1", ylab = "LV2");
          text(index, labels = 1:n)} ##
        lvs <- index}

      params <- new.params
      new.loglik <- ll(params, lvs, y, phis, pp,lw=NULL)
      err <- abs(new.loglik/current.loglik);
      if(trace) {cat("Iterations #", iter, "\n")
        cat("New Loglik:", new.loglik, "Current Loglik:", current.loglik, "Ratio", err,"\n")}
      current.loglik <- new.loglik;
      iter <- iter + 1
    }

    if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
    if(!row.eff) { row.params <- NULL; x2 <- params }
    beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
    betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
    if(num.lv > 0) {
      lambdas.mat <- matrix(x2, ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
      LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x) ))
      if(num.lv == 1) LambdaLambdaT <- t(LambdaLambdaT)
    }


    if(n.i==1 || out$logL < new.loglik){
      out$logL = current.loglik
      out$param=params
      if(num.lv > 0) {
        out$lvs=lvs
        out$params$theta=lambdas.mat
        rownames(out$lvs) = rownames(y);
        colnames(out$params$theta) = colnames(out$lvs) = paste("LV", 1:num.lv, sep="");
        rownames(out$params$theta) = colnames(y)
      }

      names(beta0)  = colnames(out$y); out$params$beta0 = beta0;
      if(!is.null(X)){betas = matrix(betas,ncol=ncol(X)); out$params$Xcoef=betas;
      rownames(out$params$Xcoef) <- colnames(y); colnames(out$params$Xcoef) <- colnames(X); }
      if(row.eff) {out$params$row.params = row.params; names(out$params$row.params) = rownames(y)}
      out$params$phi = phis; names(out$params$phi) <- colnames(y);
      out$Power=pp;

      out$start<-fit}

    n.i=n.i+1;
  }#

  if(info){
    if(trace) cat("Calculating standard errors for parameters...\n")
    for(i in 1:n) {
      update.lvs <- try(optim(out$lvs[i,], fn = ll.i, gr = NULL, method = "BFGS", control = list(trace = 0, fnscale = -1, maxit = 1), y = y, params = out$param, phis = out$params$phi, pp=out$Power, i = i, hessian = info), silent = TRUE)
      covM.lvs[i,,]<-solve(-update.lvs$hessian)
    }
    full.param.gradient<-function(params, lvs, y, X=NULL, pp,lw, phis=NULL){
      if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
      if(!row.eff) { row.params <- 0; x2 <- params }
      beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
      betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
      if(num.lv > 0) {
        lambdas.mat <- matrix(x2[1:(p * num.lv)], ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
        LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x) ))
        if(num.lv == 1) LambdaLambdaT <- t(LambdaLambdaT)
        x2 <- x2[-(1:(p*num.lv))]
      }
      phis <- x2[1:p]; x2 <- x2[-(1:p)]

      eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params + offset
      if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
      if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas, nrow = p))

      grad.beta0 <- numeric(p);
      grad.lambdas <- NULL; if(num.lv > 0) grad.lambdas <- matrix(0, p, num.lv)
      grad.row.params <- NULL; if(row.eff) grad.row.params <- numeric(n)
      grad.betas <- NULL; if(!is.null(X)) grad.betas <- matrix(0, p, num.X)

      if(num.lv > 0) {
        GA <- Hub.Gamma.tw(y = y, eta = eta.mat, lambdas.mat, phis,pp,LambdaLambdaT)
        Kron <- list()
        es<-diag(num.lv)
        for(s in 1:num.lv) {
          Kron[[s]]<-apply(lambdas.mat, 1, function(x) kronecker(es[s,],t(x))) + apply(lambdas.mat, 1, function(x) kronecker(t(es[s,]),x))
        }
        if(num.lv == 1) Kron[[1]] <- matrix(Kron[[1]],nrow=1)
      }


      V <- 1/matrix(phis, n, p, byrow = TRUE) * ((2-pp)^2*exp((2-pp)*eta.mat) - y*(1-pp)^2*exp((1-pp)*eta.mat))
      V2 <- 1/matrix(phis, n, p, byrow = TRUE) * ((2-pp)*exp((2-pp)*eta.mat) - y*(1-pp)*exp((1-pp)*eta.mat))## variance
      V1 <- 1/matrix(phis, n, p, byrow = TRUE) * (y*exp((1-pp)*eta.mat)-exp((2-pp)*eta.mat))
      if(row.eff) {grad.row.params <- rowSums(V1); if(num.lv>0) grad.row.params <-  grad.row.params-0.5 * rowSums(GA %*% t(LambdaLambdaT) * V);}
      grad.beta0 <- colSums(V1); if(num.lv>0) grad.beta0 <- grad.beta0 -0.5 * colSums(GA %*% t(LambdaLambdaT) * V)
      if(!is.null(X)) { grad.betas <- t(V1) %*% X; if(num.lv>0) grad.betas <- grad.betas -0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X}
      if(num.lv > 0) {
        for(s in 1:num.lv){ grad.lambdas[,s] <- -0.5 * colSums((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s], n, p, byrow = FALSE) + GA %*% Kron[[s]] * V2)) + colSums(V1 * lvs[,s]) }
      }

      if(num.lv > 0) grad.lambdas[upper.tri(grad.lambdas)] <- 0

      V3 <- 1/matrix(phis^2, n, p, byrow = TRUE) * (y*(1-pp)*exp((1-pp)*eta.mat) - (2-pp)*exp((2-pp)*eta.mat))## variance
      fun<-function(j){nd2(x0=phis[j],y=y[,j],mu=exp(eta.mat[,j]), p = pp,f=fishMod::dTweedie)} ## numeric derivatives, doesn't mgcv package but fishmod
      grad.phis <-colSums(sapply(1:p,fun))
      if(num.lv>0) grad.phis <- grad.phis -0.5 * colSums(GA %*% t(LambdaLambdaT) * V3)

      score.out <- c(grad.row.params, grad.beta0, grad.betas, grad.lambdas,grad.phis)
      return(score.out)
    }

    get.info <- try(-nd2.la(x0 = c(out$params$row.params, out$params$beta0, c(out$params$Xcoef), out$params$theta,out$params$phi), func = full.param.gradient, lvs = out$lvs, y = y, X = X, pp = out$Power,lw=FALSE), silent = TRUE) #

    if(!inherits(get.info, "try-error")) {
      find.const.lambdas <- which(is.finite(colSums(get.info)))
      get.info <- get.info[find.const.lambdas, find.const.lambdas]

      se <- try(sqrt(diag(abs(MASS::ginv(get.info)))))#, silent=TRUE)
      if(row.eff) { se.row.params <- se[1:n]; names(se.row.params) <- rownames(y); se <- se[-(1:n)] }
      se.beta0 <- se[1:p]; names(se.beta0)  = colnames(out$y); se <- se[-(1:p)];
      if(!is.null(X)){ se.Xcoef <- matrix(se[1:(p*ncol(X))],nrow=p); se<-se[-(1:(p*ncol(X)))]
      rownames(se.Xcoef) <- colnames(y); colnames(se.Xcoef) <- colnames(X); }
      if(num.lv>0) {
        se.lambdas<-matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag = TRUE)]<-se[1:(p * num.lv - sum(0:(num.lv-1)))];
        colnames(se.lambdas) = paste("LV", 1:num.lv, sep="");
        rownames(se.lambdas) = colnames(out$y)
        se <- se[-(1:(p * num.lv - sum(0:(num.lv-1))))];}
      se.phis=se;  names(se.phis) <- colnames(y);
      if(num.lv>0) out$lvs.cov=covM.lvs
      if(num.lv>0) out$sd<-list(theta=se.lambdas)
      out$sd$beta0=se.beta0
      if(!is.null(X)) out$sd$Xcoef=se.Xcoef
      if(row.eff) out$sd$row.params=se.row.params
      out$sd$phi = se.phis
    }
    if(inherits(get.info, "try-error")) { cat("Standard errors for parameters could not be calculated.\n") }
  }

  if(num.lv>1){
    eta.mat <- matrix(beta0, n, p, byrow=TRUE)
    if(row.eff) eta.mat <- eta.mat + row.params
    if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
    if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas, nrow = p))
    out$mu=exp(eta.mat)}
  out$start<-fit
  return(out)
}





