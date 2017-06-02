########################################################################################
## ZIP GLLVM, with estimation done via analytical Laplace approximation
## Original author: Jenni Niku
########################################################################################

gllvm.ZIP <- function(y, X = NULL, num.lv = 2, family = "ZIP", row.eff = TRUE, max.iter = 100, eps = 1e-4, seed = NULL,maxit = 1000, start.lvs = NULL, info = FALSE,samep0=FALSE, plot=FALSE, trace=TRUE,start.params=NULL,offset=NULL,n.init=1) {
  n <- dim(y)[1]; p <- dim(y)[2]; if(!is.null(X)) {X <- data.matrix(X); num.X <- dim(X)[2]	}
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

fit <- start.values.gllvm(y = y, X = X, T = NULL, family = "poisson", offset= offset, num.lv = num.lv, start.lvs = start.lvs, seed = seed[n.i])

if(is.null(start.params)){
  beta0 <- fit$params[,1] ## column intercepts
  betas <- NULL; if(!is.null(X)) betas <- c(fit$X.params) ## covariates coefficients
  lambdas <- NULL; if(num.lv > 0) lambdas <- c(fit$params[,2:(num.lv+1)]) ## LV coefficients
  row.params <- NULL; if(row.eff) row.params <- rnorm(n) ## row parameters
  lvs <- NULL; if(num.lv > 0) lvs <- matrix(fit$index, ncol = num.lv); ## LVs

  beta0 <- fit$params[,1]
  betas <- NULL; if(!is.null(X)) betas <- c(fit$params[,2:(num.X + 1)])
  lambdas <- NULL;
  if(num.lv > 0) {
    lambdas <- as.matrix(fit$params[,(ncol(fit$params) - num.lv + 1):ncol(fit$params)])
    if(num.lv > 0) lambdas[upper.tri(lambdas)] <- 0
    covM.lvs<- array(NA,dim=c(n,num.lv,num.lv))
  }
  row.params <- NULL; if(row.eff) row.params <- rep(0,n)
  lvs <- NULL; if(num.lv > 0) lvs <- matrix(fit$index, ncol = num.lv);
} else{
  if(dim(start.params$y)==dim(y) && is.null(X)==is.null(start.params$X) && (row.eff == start.params$row.eff)){
  beta0 <- start.params$params$beta0 ## column intercepts
  betas <- NULL; if(!is.null(X)) betas <- c(start.params$params$Xcoef) ## covariates coefficients
  lambdas <- NULL; if(num.lv > 0) lambdas <- c(start.params$params$theta) ## LV coefficients
  row.params <- NULL; if(row.eff) row.params <- start.params$params$row.params ## row parameters
  lvs <- NULL; if(num.lv > 0){ lvs <- matrix(start.params$lvs, ncol = num.lv); ## LVs
  covM.lvs<- array(NA,dim=c(n,num.lv,num.lv))}
  } else { stop("Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that attributes y, X and row.eff match to each other.");}
}

  p0 <- 0; if (family == "ZIP") p0 <- 0.5
  if (family == "ZIP" && !samep0) p0 <- colMeans(y==0)*0.98+0.01 # ZIP probability

  if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)

  new.params <- params <- c(row.params, beta0, betas, lambdas)

  current.loglik <- -1e6; iter <- 1; err <- 10;
  while((err > (1 + eps) || err < (1 - eps)) && iter <= max.iter) {
    ## LA-likelihood
    ll <- function(params, lvs, y, p0, type = family) {
      if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
      if(!row.eff) { row.params <- 0; x2 <- params }
      beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
      betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
      if(num.lv > 0) {
        lambdas.mat <- matrix(x2, ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
        LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x) ))
        if(num.lv == 1) LambdaLambdaT <- t(LambdaLambdaT)
      }

      eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params + offset
      if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
      if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas, nrow = p))
      mu<-exp(eta.mat)
      I.numlv <- diag(x = num.lv)
      SGam <- 0

      p0m=matrix(p0,n,p,byrow=TRUE)

               if(num.lv > 0) {
                 V <- exp(eta.mat)
                 A <- exp(-exp(eta.mat))
                 der1 <- (1 - p0m) * (A * mu * (mu - 1)) / (p0m + (1 - p0m) * A)
                 der2 <- (1 - p0m)^2 * (A^2 * mu^2) / (p0m + (1 - p0m) * A)^2

                 G<-t(t((V * (y > 0) - (der1 - der2) * (y == 0)) %*% LambdaLambdaT) + c(I.numlv))
                 if(num.lv == 1) SGam <- SGam - sum(0.5 * log(G))
                 if(num.lv == 2) SGam <- SGam - sum(0.5 * log(G[,1] * G[,4] - G[,2] * G[,3]))
                 if(num.lv > 2) { for(i in 1:n) {
                   Ga <- matrix(G[i,], nrow = num.lv, ncol = num.lv)
                   SGam <- SGam - 0.5 * log(det(Ga))}
                 }}

               SGam <- SGam + sum((log(p0m +  (1 - p0m) * exp(-exp(eta.mat)))) * (y == 0) + ((log(1-p0m) - mu + y * eta.mat - lfactorial(y))  * (y > 0)))       ################### ZIP likelihood; is this correct?


      if(num.lv > 0) SGam <- SGam - 0.5 * sum(diag(lvs %*% t(lvs)))
      return(SGam)
    }

    # Hub.Gamma for ZIP?  -3443.938

    ## Return n x num.lv^2 matrix with each row as inverse of Gamma(\theta,z_i) matrix
    ## as in eq. 9 in Huber et al.
    Hub.Gamma2 <- function(y, eta, lambdas.mat,p0,LambdaLambdaT) {
      I.numlv <- diag(x=num.lv)
      p0m=matrix(p0,n,p,byrow=TRUE)
      mu=exp(eta)
      A <- exp(-mu)
      der1 <- (1 - p0m) * (A * mu * (mu - 1)) / (p0m + (1 - p0m) * A)
      der2 <- (1 - p0m)^2 * (A^2 * mu^2) / (p0m + (1 - p0m) * A)^2
      V<-(mu * (y > 0) - (der1 - der2) * (y == 0))
      s <- matrix(0, n, num.lv^2)
      for(i in 1:n) s[i,] <- c(solve(matrix(colSums(V[i,] * LambdaLambdaT), num.lv, num.lv) + I.numlv))
      return(s)
    }

    ## Gradients of parameters row.params, beta0 and lambdas
    param.gradient2 <- function(params, lvs, y, p0, type = family) {
      if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
      if(!row.eff) { row.params <- 0; x2 <- params }
      beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
      betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
      if(num.lv > 0) {
        lambdas.mat <- matrix(x2, ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
        x2 <- x2
        LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x)))
        if(num.lv == 1) { LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) }
      }

      eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params + offset
      if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
      if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas, nrow = p))
      mu<-exp(eta.mat)

      grad.beta0 <- numeric(p);
      grad.lambdas <- NULL; if(num.lv > 0) grad.lambdas <- matrix(0, p, num.lv)
      grad.row.params <- NULL; if(row.eff) grad.row.params <- numeric(n)
      grad.betas <- NULL; if(!is.null(X)) grad.betas <- matrix(0, p, num.X)

      if(num.lv > 0) {
        GA <- Hub.Gamma2(y = y, eta = eta.mat, lambdas.mat,p0,LambdaLambdaT)
        Kron <- list()
        es<-diag(num.lv)
        for(s in 1:num.lv) {
          Kron[[s]]<-apply(lambdas.mat, 1, function(x) kronecker(es[s,],t(x))) + apply(lambdas.mat, 1, function(x) kronecker(t(es[s,]),x))
        }
        if(num.lv == 1) Kron[[1]] <- matrix(Kron[[1]],nrow=1)
      }

      p0m=matrix(p0,n,p,byrow=TRUE)
      A <- exp(-mu)
      V1 <- ((y-mu) * (y > 0) - (1-p0m)*(A * mu) / (p0m + (1 - p0m) * A) * (y == 0))
      der1 <- (1 - p0m) * (A * mu * (mu - 1)) / (p0m + (1 - p0m) * A)
      der2 <- (1 - p0m)^2 * (A^2 * mu^2) / (p0m + (1 - p0m) * A)^2
      V2 <- (mu * (y > 0) - (der1 - der2) * (y == 0))
      Dder1 <- (1 - p0m) *( (-A*mu^2*(mu - 1) + A*mu*(2*mu-1) ) / (p0m + (1 - p0m) * A)  + ((A*mu)^2*(mu-1)*(1-p0m))/(p0m + (1 - p0m) * A)^2 )
      Dder2 <- (1 - p0m)^2 *(  (2*(A*mu)^2*(1-mu)) / (p0m + (1 - p0m) * A)^2 + (2*(1-p0m)*(A*mu)^3)/(p0m + (1 - p0m) * A)^3 )
      V3 <-(mu * (y > 0) - (Dder1 - Dder2) * (y == 0))


      if(row.eff) {grad.row.params <- rowSums(V1); if(num.lv>0) grad.row.params <-  grad.row.params-0.5 * rowSums(GA %*% t(LambdaLambdaT) * V3);}
      grad.beta0 <- colSums(V1); if(num.lv>0) grad.beta0 <- grad.beta0 -0.5 * colSums(GA %*% t(LambdaLambdaT) * V3)
      if(!is.null(X)) { grad.betas <- t(V1) %*% X; if(num.lv>0) grad.betas <- grad.betas -0.5 * t(GA %*% t(LambdaLambdaT) * V3) %*% X}
      if(num.lv > 0) {
        for(s in 1:num.lv){ grad.lambdas[,s] <- -0.5 * colSums((GA %*% t(LambdaLambdaT) * V3 * matrix(lvs[,s], n, p, byrow = FALSE) + GA %*% Kron[[s]] * V2)) + colSums(V1 * lvs[,s]) }
      }

      if(num.lv > 0) grad.lambdas[upper.tri(grad.lambdas)] <- 0
      score.out <- c(grad.row.params, grad.beta0, grad.betas, grad.lambdas)
      return(score.out)
    }


    ## Gradients of parameters row.params, beta0 and lambdas
    p0.gradient <- function(params, lvs, y, p0, type = family) {
      if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
      if(!row.eff) { row.params <- 0; x2 <- params }
      beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
      betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
      if(num.lv > 0) {
        lambdas.mat <- matrix(x2, ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
        x2 <- x2
        LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x)))
        if(num.lv == 1) { LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) }
      }

      eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params + offset
      if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
      if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas, nrow = p))
      mu<-exp(eta.mat)

      grad.p0 <- numeric(p);

      if(num.lv > 0) {
        GA <- Hub.Gamma2(y = y, eta = eta.mat, lambdas.mat,p0,LambdaLambdaT)
        Kron <- list()
        es<-diag(num.lv)
        for(s in 1:num.lv) {
          Kron[[s]]<-apply(lambdas.mat, 1, function(x) kronecker(es[s,],t(x))) + apply(lambdas.mat, 1, function(x) kronecker(t(es[s,]),x))
        }
        if(num.lv == 1) Kron[[1]] <- matrix(Kron[[1]],nrow=1)
      }

      p0m=matrix(p0,n,p,byrow=TRUE)
      A <- exp(-mu)
      pV1 <- ((1-A) / (p0m + (1 - p0m) * A) * (y == 0)  -  1 / (1-p0m) * (y > 0))
      pD1 <- -A*mu*(mu - 1)/(p0m + (1 - p0m) * A) - ((1-p0m)*A*mu*(mu-1)*(1-A))/(p0m + (1 - p0m) * A)^2
      pD2 <- (  -(2*(1-p0m)*(A*mu)^2) / (p0m + (1 - p0m) * A)^2 - (2*(1-p0m)^2*(A*mu)^2*(1-A))/(p0m + (1 - p0m) * A)^3 )
      pV3 <-  0*(y > 0) -(pD1 - pD2) * (y == 0)

      grad.p0 <- colSums(pV1); if(num.lv>0) grad.p0 <- grad.p0 -0.5 * colSums(GA %*% t(LambdaLambdaT) * pV3)

      score.out <- c(grad.p0)
      return(score.out)
    }

    update.params <- try(optim(params, fn = ll, gr = param.gradient2, method = "BFGS", hessian = FALSE, control = list(trace = 0, fnscale = -1, maxit=maxit), lvs = lvs, y = y, p0 = p0, type = family), silent = TRUE)

    if(!inherits(update.params, "try-error")) {
      new.params <- update.params$par
      if(trace) cat("Model parameters updated", "\n")
    } else {new.params<-params}

    ## Update p0 for ZIP distribution;
    update.p0 <- try(optim(p0, fn = ll, gr = p0.gradient, method = "L-BFGS-B",lower =0.01,upper = 0.99, hessian = FALSE, control = list(trace = 0, fnscale = -1, maxit = maxit), y = y, params = new.params, lvs = lvs), silent = TRUE)
    if(!inherits(update.p0, "try-error")) {
      new.p0 <- update.p0$par
      if(trace) cat("Probabilities of zero inflation updated", "\n")
    } else {new.p0<-p0}


    ## Update dispersion params for NB distribution


    if(row.eff) { row.params <- new.params[1:n]; x2 <- new.params[-(1:n)] }
    if(!row.eff) { row.params <- NULL; x2 <- new.params }
    beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
    betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
    lambdas.mat = NULL; if(num.lv > 0)  lambdas.mat <- matrix(x2, ncol = num.lv)

    if(num.lv > 0) {
      ## Update latent variables
      ll.i <- function(lvs.i, params, y, p0, i) {
        if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
        if(!row.eff) { row.params <- rep(0,n); x2 <- params }
        beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
        betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
        lambdas.mat <- matrix(x2, ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0

        eta.i <- beta0 + row.params[i] + lambdas.mat %*% lvs.i + offset[i,]
        if(!is.null(X)) eta.i <- eta.i + matrix(betas,nrow = p) %*% X[i,]

        if(family == "ZIP") SGam <- sum((y[i,] * eta.i - exp(eta.i) - lfactorial(y[i,]) + log(1 - p0) ) * (y[i,] > 0) + (log(p0 + (1 - p0 ) * exp(-exp(eta.i)))) * (y[i,] == 0))
        SGam <- SGam - 0.5 * t(lvs.i) %*% lvs.i
        return(SGam)
      }

      index <- matrix(0,n,num.lv)
      for(i in 1:n) {
        update.lvs <- try(optim(lvs[i,], fn = ll.i, gr = NULL, method = "BFGS", hessian = FALSE, control = list(trace = 0, fnscale = -1, maxit = maxit), y = y, params = new.params, p0 = new.p0, i = i), silent = TRUE)
        if(!inherits(update.lvs, "try-error")){
          index[i,] <- update.lvs$par
        } else index[i,]<-lvs[i,]
      }
      if(trace) cat("Latent variables updated", "\n")

      ## Plot
      if(plot) {
      plot(rbind(index), type = "n", main = "Ordination points", xlab = "LV1", ylab = "LV2");
      text(index, labels = 1:n) } ## new

      lvs <- index
    }

    params <- new.params
    p0 <- new.p0
    new.loglik <- ll(params = params, lvs = lvs, y = y, p0 = p0)
    err <- abs(new.loglik/current.loglik);
    if(trace) {cat("Iterations #", iter, "\n");
    cat("New Loglik:", new.loglik, "Current Loglik:", current.loglik, "Ratio", err,"\n")}
    current.loglik <- new.loglik;
    iter <- iter + 1
  }


  if(n.i==1 || out$logL < new.loglik){

    if(num.lv > 0) {
      out$lvs <- lvs
      out$params$theta <- lambdas.mat
      rownames(out$lvs) <- rownames(out$y);
      colnames(out$params$theta) <- colnames(out$lvs) <- paste("LV", 1:num.lv, sep="");
      rownames(out$params$theta) <- colnames(out$y)
    }
    out$logL <- current.loglik
    out$param <- params
    out$params$beta0 <- beta0; names(out$params$beta0) <- colnames(out$y);
    if(!is.null(X)){betas <- matrix(betas,ncol=ncol(X)); out$params$Xcoef <- betas;
                    rownames(out$params$Xcoef) <- colnames(out$y); colnames(out$params$Xcoef) <- colnames(X); }
    if(row.eff) {out$params$row.params <- row.params; names(out$params$row.params) <- rownames(out$y)}
    out$params$p=p0; names(out$params$p) <- colnames(out$y);

    out$start<-fit}

  n.i <- n.i+1;
  }


  if(info) {
    if(trace) cat("Calculating standard errors for parameters...\n")

    if(num.lv>0){for(i in 1:n) {
      update.lvs <- try(optim(out$lvs[i,], fn = ll.i, gr = NULL, method = "BFGS", control = list(trace = 0, fnscale = -1, maxit = 1), y = y, params = out$param, p0 = out$params$p, i = i, hessian = info))
      if(info) covM.lvs[i,,]<-solve(-update.lvs$hessian)
    }}
      update.p0 <- try(optim(out$params$p, fn = ll, gr = p0.gradient, method = "BFGS", hessian = info, control = list(trace = 0, fnscale = -1, maxit = 1), y = y, params = out$param, lvs = out$lvs), silent = TRUE)
        Hessp0=update.p0$hessian


    ## Gradients of all model parameters in order: row.params, beta0, lambdas, betas
    full.param.gradient <- function(params, lvs, y, X = NULL,p0) {
      if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
      if(!row.eff) { row.params <- 0; x2 <- params }
      beta0 <- x2[1:p]; x2 <- x2[-(1:p)]

      if(num.lv>0) {
        lambdas.mat <- matrix(x2[1:(p * num.lv)], ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
        LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x)))
        if(num.lv==1) { LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) }
        x2 <- x2[-(1:(p*num.lv))]
      }
      if(!is.null(X)) { betas <- x2[1:(p*ncol(X))]; x2 <- x2[-(1:(p*ncol(X)))] }

      eta.mat <- matrix(beta0, n, p, byrow = TRUE) + row.params + offset
      if(num.lv>0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
      if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas,nrow=p))
      mu<-exp(eta.mat)

      grad.beta0 <- numeric(p);
      grad.lambdas <- NULL; if(num.lv > 0) grad.lambdas <- matrix(0, p, num.lv)
      grad.row.params <- NULL; if(row.eff) grad.row.params <- numeric(n)
      grad.betas <- NULL; if(!is.null(X)) grad.betas <- matrix(0, p, num.X)

      if(num.lv > 0) {
        GA <- Hub.Gamma2(y = y, eta = eta.mat, lambdas.mat,p0,LambdaLambdaT)
        Kron <- list()
        es<-diag(num.lv)
        for(s in 1:num.lv) {
          Kron[[s]]<-apply(lambdas.mat, 1, function(x) kronecker(es[s,],t(x))) + apply(lambdas.mat, 1, function(x) kronecker(t(es[s,]),x))
        }
        if(num.lv == 1) Kron[[1]] <- matrix(Kron[[1]],nrow=1)
      }

      p0m=matrix(p0,n,p,byrow=TRUE)
      A <- exp(-mu)
      V1 <- ((y-mu) * (y > 0) - (1-p0m)*(A * mu) / (p0m + (1 - p0m) * A) * (y == 0))
      der1 <- (1 - p0m) * (A * mu * (mu - 1)) / (p0m + (1 - p0m) * A)
      der2 <- (1 - p0m)^2 * (A^2 * mu^2) / (p0m + (1 - p0m) * A)^2
      V2 <- (mu * (y > 0) - (der1 - der2) * (y == 0))
      Dder1 <- (1 - p0m) *( (-A*mu^2*(mu - 1) + A*mu*(2*mu-1) ) / (p0m + (1 - p0m) * A)  + ((A*mu)^2*(mu-1)*(1-p0m))/(p0m + (1 - p0m) * A)^2 )
      Dder2 <- (1 - p0m)^2 *(  (2*(A*mu)^2*(1-mu)) / (p0m + (1 - p0m) * A)^2 + (2*(1-p0m)*(A*mu)^3)/(p0m + (1 - p0m) * A)^3 )
      V3 <-(mu * (y > 0) - (Dder1 - Dder2) * (y == 0))


      if(row.eff) {grad.row.params <- rowSums(V1); if(num.lv>0) grad.row.params <-  grad.row.params-0.5 * rowSums(GA %*% t(LambdaLambdaT) * V3);}
      grad.beta0 <- colSums(V1); if(num.lv>0) grad.beta0 <- grad.beta0 -0.5 * colSums(GA %*% t(LambdaLambdaT) * V3)
      if(!is.null(X)) { grad.betas <- t(V1) %*% X; if(num.lv>0) grad.betas <- grad.betas -0.5 * t(GA %*% t(LambdaLambdaT) * V3) %*% X}
      if(num.lv > 0) {
        for(s in 1:num.lv){ grad.lambdas[,s] <- -0.5 * colSums((GA %*% t(LambdaLambdaT) * V3 * matrix(lvs[,s], n, p, byrow = FALSE) + GA %*% Kron[[s]] * V2)) + colSums(V1 * lvs[,s]) }
      }
      if(num.lv>0) grad.lambdas[upper.tri(grad.lambdas)] <- 0

      score.out <- c(grad.row.params, grad.beta0, grad.lambdas,grad.betas)
      return(score.out)
    }
    get.info <- try(-nd2.la(x0 = c(out$params$row.params, out$params$beta0, c(out$params$theta), c(out$params$Xcoef)), func = full.param.gradient, lvs = out$lvs, y = y, X = X,p0=out$params$p), silent = TRUE) #

    if(!inherits(get.info, "try-error")) {
      find.const.lambdas <- which(is.finite(colSums(get.info)))
      get.info <- get.info[find.const.lambdas, find.const.lambdas]
      se <- try(sqrt(diag(abs(MASS::ginv(get.info)))))#, silent=TRUE)
      if(row.eff) { se.row.params <- se[1:n]; names(se.row.params) <- rownames(y); se <- se[-(1:n)] }
      se.beta0 <- se[1:p]; names(se.beta0)  = colnames(out$y); se <- se[-(1:p)]
      if(num.lv>0) {
        se.lambdas<-matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag = TRUE)]<-se[1:(p * num.lv - sum(0:(num.lv-1)))];
        colnames(se.lambdas) = paste("LV", 1:num.lv, sep="");
        rownames(se.lambdas) = colnames(out$y)
         se <- se[-(1:(p * num.lv - sum(0:(num.lv-1))))];}
      if(!is.null(X)){ se.Xcoef <- matrix(se[1:(p*ncol(X))],nrow=p); se<-se[-(1:(p*ncol(X)))]
      rownames(se.Xcoef) <- colnames(y); colnames(se.Xcoef) <- colnames(X);}

      if(num.lv>0) out$sd<-list(theta=se.lambdas)
      out$sd$beta0=se.beta0
      if(!is.null(X)) out$sd$Xcoef=se.Xcoef
      if(row.eff) out$sd$row.params=se.row.params
      out$sd$p = try(sqrt(diag(abs(MASS::ginv(Hessp0)))))
      names(out$sd$p) <- colnames(out$y);
    }
    if(inherits(get.info, "try-error")) { cat("Standard errors for parameters could not be calculated.\n") }
    if(num.lv>0) out$lvs.cov=covM.lvs

  }
  out$y=y

  return(out)
}


##############################################
