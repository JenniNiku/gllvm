########################################################################################
## GLLVM, with estimation done via analytical Laplace approximation
## Original author: Jenni Niku
## Authors: F.K.C. Hui, J. Niku
##########################################################################################


gllvm.LA <- function(y, X = NULL, num.lv = 2, family = "poisson", row.eff = TRUE, max.iter = 100, eps = 1e-4, seed = NULL,maxit = 100, start.lvs = NULL, offset=NULL, info = FALSE,trace=TRUE,plot=FALSE,phi.upd="phi",link="logit",n.init=1,constrOpt=FALSE,restrict=30,start.params=NULL) {
  n <- dim(y)[1]; p <- dim(y)[2]; if(!is.null(X)) {X <- data.matrix(X); num.X <- dim(X)[2]	}
  num.lv <- num.lv
  y <- as.matrix(y)
  if(!is.numeric(y)) stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
  if(is.null(rownames(y))) rownames(y) <- paste("Row",1:n,sep="")
  if(is.null(colnames(y))) colnames(y) <- paste("Col",1:p,sep="")

  ## Set initial values for model parameters (including dispersion prm) and latent variables
  if(!is.null(seed)) {set.seed(seed);}

  n.i<-1;
  out <- list(y=y,X=X, logL = -Inf)
  if(n.init>1) seed<-sample(1:10000,n.init)
  while(n.i<=n.init){
    if(n.init>1 && trace) cat("Initial run ",n.i,"\n");

  fit <- start.values.gllvm(y = y, X = X, T = NULL, family = family, offset= offset, num.lv = num.lv, start.lvs = start.lvs, seed = seed[n.i])
if(is.null(start.params)){
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
  lvs <- NULL; if(num.lv > 0){ lvs <- matrix(start.params$lvs, ncol = num.lv);
  covM.lvs<- array(NA,dim=c(n,num.lv,num.lv))}## LVs
  } else { stop("Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that attributes y, X and row.eff match to each other.");}
}
  phis <- NULL; if(family == "negative.binomial") {phis <- fit$phi; if(any(phis>100))phis[phis>100]=100; if(any(phis<0.10))phis[phis<0.10]=0.10;}# && phi.upd=="inv.phi"

  params <- c(row.params, beta0, betas, lambdas)
  new.params <- params

    if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)

  current.loglik <- -1e6; iter <- 1; err <- 10;
  while((err > (1 + eps) || err < (1 - eps)) && iter <= max.iter) {
    ## LA-likelihood
    ll <- function(params, lvs, y, phis, type = family) {
      if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
      if(!row.eff) { row.params <- 0; x2 <- params }
      beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
      betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
      if(num.lv > 0) {
        lambdas.mat <- matrix(x2, ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
        LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x) ))
        if(num.lv == 1) LambdaLambdaT <- t(LambdaLambdaT)
        lvs<-matrix(c(lvs),nrow = n,ncol = num.lv)
      }

      eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params  + offset
      if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
      if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas, nrow = p))
      mu<-exp(eta.mat)
      I.numlv <- diag(x = num.lv)
      SGam <- 0

      switch(type,
             binomial = {
               probs<-binomial(link = link)$linkinv(eta.mat);
               if(num.lv > 0) {
                 if(link=="logit"){
                   V <- binomial()$variance(probs)
                 }else {
                   V <- mu - y*(mu - mu*exp(-mu) - mu^2*exp(-mu))/((1-exp(-mu))^2)
                 }
                 G<-t(t(V %*% LambdaLambdaT) + c(I.numlv))
                 if(num.lv == 1) SGam <- SGam - sum(0.5 * log(G))
                 if(num.lv == 2) SGam <- SGam - sum(0.5 * log(G[,1] * G[,4] - G[,2] * G[,3]))
                 if(num.lv > 2){ for(i in 1:n) {
                   Ga <- matrix(G[i,], nrow = num.lv, ncol = num.lv)
                   SGam <- SGam - 0.5 * log(det(Ga))}
                 }
               }
               SGam <- SGam + sum(dbinom(y, 1, prob = probs, log = TRUE))
             },
             poisson = {
               if(num.lv > 0) {
                 V <- poisson()$variance(poisson()$linkinv(eta.mat))
                 G<-t(t(V %*% LambdaLambdaT) + c(I.numlv))
                 if(num.lv == 1) SGam <- SGam - sum(0.5 * log(G))
                 if(num.lv == 2) SGam <- SGam - sum(0.5 * log(G[,1] * G[,4] - G[,2] * G[,3]))
                 if(num.lv > 2){ for(i in 1:n) {
                   Ga <- matrix(G[i,], nrow = num.lv, ncol = num.lv)
                   SGam <- SGam - 0.5 * log(det(Ga))}
                 }
               }
               SGam <- SGam + sum(dpois(y, lambda = mu, log = TRUE))
             },
             negative.binomial = {
               if(num.lv > 0) {
                 V <- mu * (1 + y * matrix(phis,n,p,byrow=TRUE)) / (1 + matrix(phis,n,p,byrow=TRUE) * mu)^2
                 G<-t(t(V %*% LambdaLambdaT) + c(I.numlv))
                 if(num.lv == 1) SGam <- SGam - sum(0.5 * log(G))
                 if(num.lv == 2) SGam <- SGam - sum(0.5 * log(G[,1] * G[,4] - G[,2] * G[,3]))
                 if(num.lv > 2) { for(i in 1:n) {
                   Ga <- matrix(G[i,], nrow = num.lv, ncol = num.lv)
                   SGam <- SGam - 0.5 * log(det(Ga))}
                 }
               }
               SGam <- SGam + sum(dnbinom(y, mu = mu, size = matrix(1 / phis, n, p, byrow=TRUE), log = TRUE))
             })

      if(num.lv > 0) SGam <- SGam - 0.5 * sum(diag(lvs %*% t(lvs)))
      return(SGam)
    }

    ## Return n x num.lv^2 matrix with each row as inverse of Gamma(\theta,z_i) matrix
    ## as in eq. 9 in Huber et al.
    Hub.Gamma <- function(y, eta, lambdas.mat, phis,LambdaLambdaT) {
      I.numlv <- diag(x=num.lv)
      if(family == "binomial"){
        if(link=="logit"){ V <- binomial()$variance(binomial()$linkinv(eta))
        }	else { mu<-exp(eta);
        probs<-binomial(link = link)$linkinv(eta);
        V <- mu - y*( mu/(probs) - exp(-mu)*(mu/(probs))^2)#
        }
      }
      if(family == "poisson")  V <- poisson()$variance(poisson()$linkinv(eta))
      if(family == "negative.binomial")  V <- as.matrix(exp(eta) * (1 + y * matrix(phis, n, p, byrow = TRUE)) / (1 + matrix(phis, n, p, byrow = TRUE) * exp(eta))^2 )
      s <- matrix(0, n, num.lv^2)
      for(i in 1:n) s[i,] <- c(solve(matrix(c(V[i,]) %*% LambdaLambdaT, num.lv, num.lv) + I.numlv))
      return(s)
    }

    ## Gradients of parameters row.params, beta0 and lambdas
    param.gradient <- function(params, lvs, y, phis, type = family) {
      if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
      if(!row.eff) { row.params <- 0; x2 <- params }
      beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
      betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
      if(num.lv > 0) {
        lambdas.mat <- matrix(x2, ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
        LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x)))
        if(num.lv == 1) { LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) }
      }

      eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params  + offset
      if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
      if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas, nrow = p))
      mu<-exp(eta.mat)

      grad.beta0 <- numeric(p);
      grad.lambdas <- NULL; if(num.lv > 0) grad.lambdas <- matrix(0, p, num.lv)
      grad.row.params <- NULL; if(row.eff) grad.row.params <- numeric(n)
      grad.betas <- NULL; if(!is.null(X)) grad.betas <- matrix(0, p, num.X)

      if(num.lv > 0) {
        GA <- Hub.Gamma(y = y, eta = eta.mat, lambdas.mat, phis,LambdaLambdaT)
        Kron <- list()
        es<-diag(num.lv)
        for(s in 1:num.lv) { Kron[[s]]<-t((lambdas.mat) %x% t(es[s,])  ) + t(t(es[s,]) %x% (lambdas.mat)) }
        if(num.lv == 1) Kron[[1]] <- matrix(Kron[[1]],nrow=1)
      }

      switch(type,
             binomial = {
               if(link=="logit"){
                 V3 <- (mu * (1 - mu)/(1 + mu)^3) ## third derivative
                 V2 <- binomial()$variance(binomial()$linkinv(eta.mat))	## second derivative
                 V1 <- y - binomial()$linkinv(eta.mat)
               }	else {
                 probs<-binomial(link = link)$linkinv(eta.mat);
                 V3 <- mu - y*( mu/(probs) - exp(-mu)*(3*mu^2-mu^3)/(probs^2) + 2*mu^3*exp(-2*mu)/(probs^3) ) ## third derivative
                 V2 <- mu - y*( mu/(probs) - mu^2*exp(-mu)/(probs)^2)
                 V1 <- y*mu/(probs) - mu

                 }

               if(row.eff) {
                 grad.row.params <- rowSums(V1);
                 if(num.lv > 0) grad.row.params <- grad.row.params - 0.5 * rowSums(GA %*% t(LambdaLambdaT) * V3)
               }
               grad.beta0 <- colSums(V1); if(num.lv > 0) grad.beta0 <- grad.beta0 - 0.5 * colSums(GA %*% t(LambdaLambdaT) * V3);
               if(!is.null(X)) { grad.betas <- t(V1) %*% X; if(num.lv > 0) grad.betas <- grad.betas -0.5 * t(GA %*% t(LambdaLambdaT) * V3) %*% X;}
               if(num.lv > 0){
                 for(s in 1:num.lv) { grad.lambdas[,s] <- -0.5 * colSums((GA %*% t(LambdaLambdaT) * V3 * matrix(lvs[,s], n, p, byrow = FALSE) + GA %*% Kron[[s]] * V2)) + colSums((V1) * lvs[,s]) }
               }
             },

             poisson = {
               V <- mu
               if(row.eff) {grad.row.params <- rowSums(y - poisson()$linkinv(eta.mat)); if(num.lv>0) grad.row.params <- grad.row.params -0.5 * rowSums(GA %*% t(LambdaLambdaT) * V);}
               grad.beta0 <- colSums(y - poisson()$linkinv(eta.mat)); if(num.lv>0) grad.beta0 <- grad.beta0 -0.5 * colSums(GA %*% t(LambdaLambdaT) * V);
               if(!is.null(X)) { grad.betas <- t(y - poisson()$linkinv(eta.mat)) %*% X; if(num.lv>0) grad.betas<- grad.betas-0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X}
               if(num.lv > 0) {
                 for(s in 1:num.lv) { grad.lambdas[,s] <- -0.5 * colSums((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s], n, p, byrow = FALSE) + GA %*% Kron[[s]] * V)) + colSums((y - poisson()$linkinv(eta.mat)) * lvs[,s]) }
               }
             },

             negative.binomial = {
               V <- mu * (1 + y * matrix(phis, n, p, byrow = TRUE)) * (1 - matrix(phis, n, p, byrow = TRUE) * mu)/(1 + matrix(phis, n, p, byrow = TRUE) * mu)^3 ## third derivative
               V2 <- mu * (1 + y * matrix(phis, n, p, byrow = TRUE)) / (1 + matrix(phis, n, p, byrow = TRUE) * mu)^2 ## variance
               yminuseta <- (y - mu) / (1 + matrix(phis, n, p, byrow = TRUE) * mu)
               if(row.eff) {grad.row.params <- rowSums(yminuseta); if(num.lv>0) grad.row.params <-  grad.row.params-0.5 * rowSums(GA %*% t(LambdaLambdaT) * V);}
               grad.beta0 <- colSums(yminuseta); if(num.lv>0) grad.beta0 <- grad.beta0 -0.5 * colSums(GA %*% t(LambdaLambdaT) * V)
               if(!is.null(X)) { grad.betas <- t(yminuseta) %*% X; if(num.lv>0) grad.betas <- grad.betas -0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X}
               if(num.lv > 0) {
                 for(s in 1:num.lv){ grad.lambdas[,s] <- -0.5 * colSums((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s], n, p, byrow = FALSE) + GA %*% Kron[[s]] * V2)) + colSums(yminuseta * lvs[,s]) }
               }
             })

      if(num.lv > 0) grad.lambdas[upper.tri(grad.lambdas)] <- 0
      score.out <- c(grad.row.params, grad.beta0, grad.betas, grad.lambdas)
      return(score.out)
    }

    if(constrOpt) {
    A=NULL
    if(row.eff){for(i in 1:p){
      A=rbind(A,diag(n))
    }}

    B0=matrix(1,n,1)
    B1=B0
    for(i in 1:(p-1)){
      B1=as.matrix(Matrix::bdiag(B1,B0))
    }

    B2=NULL
    if(!is.null(X)){
      for(s in 1:num.X){
        B2.=X[,s]
        B0=X[,s]
        for(i in 1:(p-1)){
          B2.=as.matrix(Matrix::bdiag(B2.,B0))
        }
        B2=cbind(B2,B2.)
        }
    }

    G=NULL
    if(num.lv>0){for(s in 1:num.lv){
      G0=lvs[,s]
      Gi=lvs[,s]
      for(i in 1:(p-1)){
        Gi=as.matrix(Matrix::bdiag(Gi,G0))
      }
      G=cbind(G,Gi)
    }}

    U=cbind(A,B1,B2,G)

    update.params <- try(constrOptim(params, f = ll, grad = param.gradient, method = "BFGS", control = list(trace = 0, fnscale = -1, maxit=maxit), lvs = lvs, y = y, phis = phis, type = family,ui=rbind(U,-U),ci=rep(-restrict,NROW(U)*2)), silent = TRUE)
    } else {
      update.params <- try(optim(params, fn = ll, gr = param.gradient, method = "BFGS", control = list(trace = 0, fnscale = -1, maxit=maxit), lvs = lvs, y = y, phis = phis, type = family), silent = TRUE)
      }

    if(!inherits(update.params, "try-error")) {
      new.params <- update.params$par
      if(trace) cat("Model parameters updated", "\n")
    } else {new.params<-params}

    if(row.eff) { row.params <- new.params[1:n]; x2 <- new.params[-(1:n)] }
    if(!row.eff) { row.params <- NULL; x2 <- new.params }
    beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
    betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
    lambdas.mat <- NULL; if(num.lv > 0) lambdas.mat <- matrix(x2, ncol = num.lv)


    ## Update dispersion params for NB distribution
    if(family == "negative.binomial") {
      ## LA-likelihood
      ll.nb <- function(phis, params, lvs, y) {
        if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
        if(!row.eff) { row.params <- 0; x2 <- params }
        beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
        betas <- NULL;
        if(!is.null(X)) {
          betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))]
        }
        if(num.lv > 0) {
          lambdas.mat <- matrix(x2, ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
          LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x%*%t(x) ))
          if(num.lv == 1) {LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) }
        }

        eta.mat <- matrix(beta0, n, p, byrow = TRUE) + row.params
        if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
        if(!is.null(X)) eta.mat <- eta.mat + X %*% t(betas)
        mu<-exp(eta.mat)

        I.numlv <- diag(x = num.lv)
        SGam <- 0

        if(num.lv>0){
          V <- mu * (1 + y * matrix(phis, n, p, byrow = TRUE)) / (1 + matrix(phis, n, p, byrow = TRUE) * mu)^2
          G<-t(t(V %*% LambdaLambdaT) + c(I.numlv))
          if(num.lv == 1) SGam <- SGam - sum(0.5 * log(G))
          if(num.lv == 2) SGam <- SGam - sum(0.5 * log(G[,1] * G[,4] - G[,2] * G[,3]))
          if(num.lv > 2){ for(i in 1:n) {
            Ga <- matrix(G[i,], nrow = num.lv, ncol = num.lv)
            SGam <- SGam - 0.5 * log(det(Ga))}
          }
        }
        SGam <- SGam + sum(dnbinom(y, mu = mu, size = matrix(1 / phis, n, p, byrow = TRUE), log = TRUE))

        if(num.lv>0) SGam <- SGam - 0.5 * sum(diag(lvs %*% t(lvs)))
        return(SGam)
      }
      ## gradients for dispersion parameters \phi

      phi.gradient <- function(phis, params, lvs, y,type=family) {
        if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
        if(!row.eff) { row.params <- 0; x2 <- params }
        beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
        betas <- NULL;
        if(!is.null(X)) {
          betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))]
        }

        if(num.lv>0) {
          lambdas.mat <- matrix(x2, ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
          LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x) ))
          if(num.lv==1) { LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) }
        }

        eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params
        if(num.lv>0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
        if(!is.null(X)) eta.mat <- eta.mat + X %*% t(betas)
        mu<-exp(eta.mat)

        #
        V <- mu
        yf <- y + matrix(1 / phis,n,p,byrow=TRUE)
        V2 <- V * (y/(1+matrix(phis,n,p,byrow=TRUE)*V)^2 - 2*V*(1+y*matrix(phis,n,p,byrow=TRUE))/(1+matrix(phis,n,p,byrow=TRUE)*V)^3)
        V3 <- matrix(1 / phis^2,n,p,byrow=TRUE) * (log(1 + matrix(phis,n,p,byrow=TRUE) * V) - digamma(yf) + digamma(matrix(1 / phis,n,p,byrow=TRUE)))
        grad.phis <- colSums(V3 - V*yf/(1+matrix(phis,n,p,byrow=TRUE)*V) + y/matrix(phis,n,p,byrow=TRUE))
        if(num.lv>0) { GA <- Hub.Gamma(y = y, eta = eta.mat, lambdas.mat, phis,LambdaLambdaT); grad.phis <- grad.phis+colSums(-0.5*GA%*%t(LambdaLambdaT)*V2)}
        return(grad.phis)
      }

      ll.nb2 <- function(iphis, params, lvs, y) {
        if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
        if(!row.eff) { row.params <- 0; x2 <- params }
        beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
        betas <- NULL;
        if(!is.null(X)) {
          betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))]
        }
        if(num.lv > 0) {
          lambdas.mat <- matrix(x2, ncol = num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
          LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x%*%t(x) ))
          if(num.lv == 1) {LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) }
        }
        ff<-iphis
        phis<-1/iphis
        eta.mat <- matrix(beta0, n, p, byrow = TRUE) + row.params
        if(num.lv > 0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
        if(!is.null(X)) eta.mat <- eta.mat + X %*% t(betas)
        mu<-exp(eta.mat)

        I.numlv <- diag(x = num.lv)
        SGam <- 0

        if(num.lv>0){
          V <- mu * (1 + y * matrix(phis, n, p, byrow = TRUE)) / (1 + matrix(phis, n, p, byrow = TRUE) * mu)^2
          G<-t(t(V %*% LambdaLambdaT) + c(I.numlv))
          if(num.lv == 1) SGam <- SGam - sum(0.5 * log(G))
          if(num.lv == 2) SGam <- SGam - sum(0.5 * log(G[,1] * G[,4] - G[,2] * G[,3]))
          if(num.lv > 2){ for(i in 1:n) {
            Ga <- matrix(G[i,], nrow = num.lv, ncol = num.lv)
            SGam <- SGam - 0.5 * log(det(Ga))}
          }
        }
        SGam <- SGam + sum(dnbinom(y, mu = mu, size = matrix(1 / phis, n, p, byrow = TRUE), log = TRUE))

        if(num.lv>0) SGam <- SGam - 0.5 * sum(diag(lvs %*% t(lvs)))
        return(SGam)
      }
      phi.gradient2 <- function(iphis, params, lvs, y,type=family) {
        if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
        if(!row.eff) { row.params <- 0; x2 <- params }
        beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
        betas <- NULL;
        if(!is.null(X)) {
          betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))]
        }

        if(num.lv>0) {
          lambdas.mat <- matrix(x2, ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
          LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x) ))
          if(num.lv==1) { LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) }
        }
        ff<-iphis
        phis=1/iphis

        eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params
        if(num.lv>0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
        if(!is.null(X)) eta.mat <- eta.mat + X %*% t(betas)
        mu <- exp(eta.mat)
        grad.phis <- NULL
        #
        ff <- matrix(ff,n,p,byrow=TRUE)
        V2 <- mu *ff^(-2)* (-y/(1+mu/ff)^2 + 2*mu*(1+y/ff)/(1+mu/ff)^3)
        V3 <- (-log(1+mu/ff) + (1+y/ff)*mu/(ff+mu) -y/ff + digamma(y+ff) - digamma(ff))
        grad.phis <- colSums(V3)
        if(num.lv>0) { GA <- Hub.Gamma(y = y, eta = eta.mat, lambdas.mat, phis,LambdaLambdaT); grad.phis <- grad.phis+colSums(-0.5*GA%*%t(LambdaLambdaT)*V2)}
        return(grad.phis)
      }

     if(phi.upd=="phi"){
       update.phis <- try(optim(par = phis, fn = ll.nb, gr = phi.gradient, method = "L-BFGS-B", lower = 1e-4, upper = 1e4, control = list(trace = 0, fnscale = -1), lvs = lvs, y = y, params = new.params), silent = TRUE)
       #update.phis <- try(optim(par = phis, fn = ll.nb, gr = phi.gradient, method = "BFGS", control = list(trace = 0, fnscale = -1), lvs = lvs, y = y, params = new.params), silent = TRUE)#,reltol=1e-4,maxit=maxit
      if(!inherits( update.phis, "try-error")) {
          if(trace) cat("Dispersion parameters updated", "\n")
          phis <- (update.phis$par)
      }
      else { phis <- phis }}

      if(phi.upd=="inv.phi"){
        update.phis <- try(optim(par = 1/phis, fn = ll.nb2, gr = phi.gradient2, method = "L-BFGS-B", lower = 1e-4, upper = 1e4, control = list(trace = 0, fnscale = -1), lvs = lvs, y = y, params = new.params), silent = TRUE)
        #update.phis <- try(optim(par = 1/phis, fn = ll.nb2, gr = phi.gradient2, method = "BFGS", control = list(trace = 0, fnscale = -1), lvs = lvs, y = y, params = new.params), silent = TRUE)
      if(!inherits( update.phis, "try-error")) {
        if(trace) cat("Dispersion parameters updated", "\n")
        phis <- 1/(update.phis$par)
      }
      else { phis <- phis }}

    }



    params <- new.params


    if(num.lv > 0) {
      ## Update latent variables lvs
      ll.i <- function(lvs.i, params, y, phis, i) {
        if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
        if(!row.eff) { row.params <- rep(0,n); x2 <- params }
        beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
        betas <- NULL; if(!is.null(X)) { betas <- matrix(x2[1:(p * num.X)], p, num.X); x2 <- x2[-(1:(p * num.X))] }
        lambdas.mat <- matrix(x2, ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0

        eta.i <- beta0 + row.params[i] + lambdas.mat %*% lvs.i   + offset[i,]
        if(!is.null(X)) eta.i <- eta.i + matrix(betas,nrow = p) %*% X[i,]

        if(family == "binomial")  SGam <- sum(dbinom(y[i,], 1, prob = binomial(link=link)$linkinv(eta.i), log = TRUE))
        if(family == "poisson")  SGam <- sum(dpois(y[i,], lambda = exp(eta.i), log = TRUE))
        if(family == "negative.binomial")  SGam <- sum(dnbinom(y[i,], mu = exp(eta.i), size = 1 / phis, log = TRUE))
        SGam <- SGam - 0.5 * t(lvs.i) %*% lvs.i
        return(SGam)
      }

      if(row.eff){ rp.i=row.params }else{ rp.i=rep(0,n);}
      index <- matrix(0,n,num.lv)
      for(i in 1:n) {
        if(constrOpt){
          rp=rp.i[i]
          if(!is.null(X)) rp=rp+betas%*%X[i,];
          update.lvs <- constrOptim(lvs[i,], f = ll.i,grad=NULL,ui=rbind(lambdas.mat,-lambdas.mat),ci=-restrict+c(-(beta0+rp),(beta0+rp)), y = y, params = new.params, phis = phis, i = i , control = list(trace = 0, fnscale = -1, maxit = maxit,reltol=1e-5),hessian = FALSE)
        } else {
          update.lvs <- try(optim(lvs[i,], fn = ll.i, gr = NULL, method = "BFGS", control = list(trace = 0, fnscale = -1, maxit = maxit), y = y, params = new.params, phis = phis, i = i))#, silent = TRUE)#,reltol=1e-5
        }
        if(!inherits( update.lvs, "try-error")) {index[i,] <- update.lvs$par}else{index[i,] <- lvs[i,]}
      }
      if(trace) cat("Latent variables updated", "\n")

      ## Plot
      if(plot){
      plot(rbind(index), type = "n", main = "Ordination points", xlab = "LV1", ylab = "LV2");
      text(index, labels = 1:n) ##
      }
      lvs <- index}


    new.loglik <- ll(params = params, lvs = lvs, y = y, phis = phis)
    err <- abs(new.loglik/current.loglik);
    if(trace) {cat("Iterations #", iter, "\n")
    cat("New Loglik:", new.loglik, "Current Loglik:", current.loglik, "Ratio", err,"\n")}
    current.loglik <- new.loglik;
    iter <- iter + 1
  }


  if(n.i==1 || out$logL < new.loglik){
    out$logL <- current.loglik
    out$params <- params
  names(beta0) <- colnames(out$y); out$beta0 <- beta0;
  if(row.eff) {out$row.params <- row.params; names(out$row.params) <- rownames(out$y)}
  if(!is.null(X)){betas <- matrix(betas,ncol=ncol(X)); out$betas <- betas;
                  rownames(out$betas) <- colnames(out$y); colnames(out$betas) <- colnames(X); }
  if(num.lv > 0) {
    out$lvs <- lvs
    out$lambdas <- lambdas.mat
    rownames(out$lvs) <- rownames(out$y);
    colnames(out$lambdas) <- colnames(out$lvs) <- paste("LV", 1:num.lv, sep="");
    rownames(out$lambdas) <- colnames(out$y)
  }
  if(family == "negative.binomial") {out$phis <- phis; names(out$phis) <- colnames(out$y);}
  if(family == "binomial") out$link <- link;

  out$start<-fit}

  n.i <- n.i+1;
}

  if(info) {
    if(trace) cat("Calculating standard errors for parameters...\n")

    if(num.lv>0){for(i in 1:n) {
      update.lvs <- try(optimHess(out$lvs[i,], fn = ll.i, gr = NULL, control = list(trace = 0, fnscale = -1, maxit = 0), y = y, params = out$params, phis = out$phis, i = i))
      if(info) covM.lvs[i,,]<-solve(-update.lvs)
    }}

    ## Gradients of all model parameters in order: row.params, beta0, lambdas, betas, phis
    full.param.gradient <- function(params, lvs, y, X = NULL, type = family,phis=NULL) {
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

      if(type == "negative.binomial") { ff <- x2[1:p]; x2 <- x2[-(1:p)] ; phis <- 1/ff }

      eta.mat <- matrix(beta0, n, p, byrow = TRUE) + row.params
      if(num.lv>0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
      if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas,nrow=p))
      mu<-exp(eta.mat)

      grad.beta0 <- numeric(p);
      grad.lambdas <- NULL; if(num.lv>0) grad.lambdas <- matrix(0, p, num.lv)
      grad.row.params <- NULL; if(row.eff) grad.row.params <- numeric(n)
      grad.betas <- NULL; if(!is.null(X)) grad.betas <- matrix(0, p, ncol(X))
      grad.phis <- NULL;

      if(num.lv>0) {
        GA <- Hub.Gamma(y = y, eta = eta.mat, lambdas.mat, phis,LambdaLambdaT)
        Kron <- list()
        es<-diag(num.lv)
        for(s in 1:num.lv) { Kron[[s]]<-t((lambdas.mat) %x% t(es[s,])) + t(t(es[s,]) %x% (lambdas.mat)) }
        if(num.lv==1) Kron[[1]] <- matrix(Kron[[1]],nrow=1)
      }


      switch(type,
             binomial = {

               if(link=="logit"){
                 V3 <- (mu * (1 - mu)/(1 + mu)^3) ## third derivative
                 V2 <- binomial()$variance(binomial()$linkinv(eta.mat))	## second derivative
                 V1 <- y - binomial()$linkinv(eta.mat)
               }	else {
                 V3 <- mu - y*( mu/((1-exp(-mu))+1e-5) - exp(-mu)*(3*mu^2-mu^3)/((1-exp(-mu))^2+1e-5) + 2*mu^3*exp(-2*mu)/((1-exp(-mu))^3+1e-5) ) ## third derivative
                 V2 <- mu - y*( mu/(1-exp(-mu)+1e-5) - mu^2*exp(-mu)/(1-exp(-mu)+1e-5)^2)
                 #V2[V2>1e7]<- 1e7
                 V1 <- y*mu/(1-exp(-mu)+1e-5) - mu
               }

               if(row.eff) {
                 grad.row.params <- rowSums(V1);
                 if(num.lv > 0) grad.row.params <- grad.row.params - 0.5 * rowSums(GA %*% t(LambdaLambdaT) * V3)
               }
               grad.beta0 <- colSums(V1); if(num.lv > 0) grad.beta0 <- grad.beta0 - 0.5 * colSums(GA %*% t(LambdaLambdaT) * V3);
               if(!is.null(X)) { grad.betas <- t(V1) %*% X; if(num.lv > 0) grad.betas <- grad.betas -0.5 * t(GA %*% t(LambdaLambdaT) * V3) %*% X;}
               if(num.lv > 0){
                 for(s in 1:num.lv) { grad.lambdas[,s] <- -0.5 * colSums((GA %*% t(LambdaLambdaT) * V3 * matrix(lvs[,s], n, p, byrow = FALSE) + GA %*% Kron[[s]] * V2)) + colSums((V1) * lvs[,s]) }
               }
             },


             poisson = {
               V <- mu ## third derivative and variance
               if(row.eff){ grad.row.params <- rowSums(y - poisson()$linkinv(eta.mat)); if(num.lv>0) grad.row.params <- grad.row.params -0.5 * rowSums(GA %*% t(LambdaLambdaT) * V)}
               grad.beta0 <- colSums(y - poisson()$linkinv(eta.mat)); if(num.lv>0) grad.beta0 <- grad.beta0 -0.5 * colSums(GA %*% t(LambdaLambdaT) * V)

               if(num.lv>0) {
                 for(s in 1:num.lv) { grad.lambdas[,s] <- -0.5 * colSums((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s],n,p,byrow=FALSE) + GA %*% Kron[[s]] * V)) + colSums((y - poisson()$linkinv(eta.mat)) * lvs[,s]) }
               }
               if(!is.null(X)) { grad.betas <- t(y - poisson()$linkinv(eta.mat)) %*% X; if(num.lv>0) grad.betas <- grad.betas -0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X }
             },


             negative.binomial = {
               phimat<-matrix(phis,n,p,byrow=TRUE)
               V <- mu * (1 + y * phimat) * (1 - phimat * mu)/(1 + phimat * mu)^3 ## third derivative
               V2 <- mu * (1 + y * phimat)/(1 + phimat * mu)^2 ## variance
               yminuseta <- (y - mu) / (1 + phimat * mu)

               if(row.eff) {grad.row.params <- rowSums(yminuseta); if(num.lv>0) grad.row.params<-grad.row.params-0.5 * rowSums(GA %*% t(LambdaLambdaT) * V)}
               grad.beta0 <- colSums(yminuseta); if(num.lv>0) grad.beta0 <- grad.beta0-0.5 * colSums(GA %*% t(LambdaLambdaT) * V)
               if(num.lv>0) {
                 for(s in 1:num.lv){ grad.lambdas[,s] <- -0.5 * colSums((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s],n,p,byrow=FALSE) + GA %*% Kron[[s]] * V2)) + colSums(yminuseta * lvs[,s]) }
               }
               if(!is.null(X)) { grad.betas <-t(yminuseta) %*% X ; if(num.lv>0)grad.betas<- grad.betas-0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X }

               ff<-matrix(ff,n,p,byrow=TRUE)
               V2 <- mu *ff^(-2)* (-y/(1+mu/ff)^2 + 2*mu*(1+y/ff)/(1+mu/ff)^3)
               V3 <- (-log(1+mu/ff) + (1+y/ff)*mu/(ff+mu) -y/ff + digamma(y+ff) - digamma(ff))
               grad.phis <- colSums(V3)
               if(num.lv>0) { grad.phis <- grad.phis+colSums(-0.5*GA%*%t(LambdaLambdaT)*V2)}
             }
      )

      if(num.lv>0) grad.lambdas[upper.tri(grad.lambdas)] <- 0
      score.out <- c(grad.row.params, grad.beta0, grad.lambdas, grad.betas,grad.phis)
      return(score.out)
    }
    get.info <- try(-nd2.la(x0 = c(out$row.params, out$beta0, out$lambdas, c(out$betas),1/out$phis), func = full.param.gradient, lvs = out$lvs, y = y, X = X, type = family, phis=out$phis), silent = TRUE) #

    if(!inherits(get.info, "try-error")) {
      find.const.lambdas <- which(is.finite(colSums(get.info)))
      get.info <- get.info[find.const.lambdas, find.const.lambdas]
      se <- try(sqrt(diag(abs(MASS::ginv(get.info)))))#, silent=TRUE)
      if(row.eff) { out$se.row.params <- se[1:n]; names(out$se.row.params) <- rownames(y); se <- se[-(1:n)] }
      out$se.beta0 <- se[1:p]; names(out$se.beta0)  = colnames(out$y); se <- se[-(1:p)]
      if(num.lv>0) {
        se.lambdas<-matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag = TRUE)]<-se[1:(p * num.lv - sum(0:(num.lv-1)))];
        colnames(se.lambdas) = paste("LV", 1:num.lv, sep="");
        rownames(se.lambdas) = colnames(out$y)
        out$se.lambdas<-se.lambdas; se <- se[-(1:(p * num.lv - sum(0:(num.lv-1))))];}
      if(!is.null(X)){ out$se.Xcoef <- matrix(se[1:(p*ncol(X))],nrow=p); se<-se[-(1:(p*ncol(X)))]
      rownames(out$se.Xcoef) <- colnames(y); colnames(out$se.Xcoef) <- colnames(X);}
      if(family == "negative.binomial") { out$se.phis=se;  names(out$se.phis) <- colnames(y);
      if(phi.upd=="phi"){
        hess.phis <- try(optimHess(par = out$phis, fn = ll.nb, gr = phi.gradient, control = list(trace = 0, fnscale = -1,maxit=0), lvs = lvs, y = y, params = c(out$row.params, out$beta0, out$lambdas, c(out$betas)) ))
        if(!inherits( update.phis, "try-error")) {
        out$se.phis1 <- try(sqrt(diag(abs(MASS::ginv(-hess.phis)))))
      }
      }}
    }
    if(inherits(get.info, "try-error")) { cat("Standard errors for parameters could not be calculated.\n") }
    if(num.lv>0) out$lvs.cov=covM.lvs
  }

  return(out)
}

