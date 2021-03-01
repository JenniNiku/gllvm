start.values.gllvm.TMB <- function(y, X = NULL, TR=NULL, family, 
        offset= NULL, trial.size = 1, num.lv = 0, start.lvs = NULL, 
        seed = NULL,power=NULL,starting.val="res",formula=NULL, 
        jitter.var=0,yXT=NULL, row.eff=FALSE, TMB=TRUE, 
        link = "probit", randomX = NULL, beta0com = FALSE, zeta.struc=NULL) {
  if(!is.null(seed)) set.seed(seed)
  N<-n <- nrow(y); p <- ncol(y); y <- as.matrix(y)
  num.T <- 0; if(!is.null(TR)) num.T <- dim(TR)[2]
  num.X <- 0; if(!is.null(X)) num.X <- dim(X)[2]
  Br <- sigmaB <- sigmaij <- NULL
  mu<-NULL
  out <- list()

  sigma=1
  
  row.params <- rep(0, n);
  if(starting.val %in% c("res","random") || row.eff == "random"){
    rmeany <- rowMeans(y)
    if(family=="binomial"){
      rmeany=1e-3+0.99*rmeany
      if(row.eff %in% c("fixed",TRUE)) {
        row.params <-  binomial(link = link)$linkfun(rmeany) - binomial(link = link)$linkfun(rmeany[1])
      } else{
        row.params <-  binomial(link = link)$linkfun(rmeany) - binomial(link = link)$linkfun(mean(rmeany))
      }
    } else if(family=="gaussian"){
      rmeany=1e-3+0.99*rmeany
      if(row.eff %in% c("fixed",TRUE)) {
        row.params <-  rmeany - rmeany[1]
      } else{
        row.params <-  rmeany - mean(rmeany)
      }
    } else {
      if(row.eff %in% c("fixed",TRUE)) {
        row.params <-  row.params <- log(rmeany)-log(rmeany[1])
      } else{
        row.params <-  row.params <- log(rmeany)-log(mean(y))
      }
    }
    if(any(abs(row.params)>1.5)) row.params[abs(row.params)>1.5] <- 1.5 * sign(row.params[abs(row.params)>1.5])
    sigma=sd(row.params)
  }

  if(!is.numeric(y))
    stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")

  if(family=="ZIP") family="poisson"

  if(!(family %in% c("poisson","negative.binomial","binomial","ordinal","tweedie", "gaussian", "gamma", "exponential")))
    stop("inputed family not allowed...sorry =(")

  if(num.lv > 0) {
    unique.ind <- which(!duplicated(y))
    if(is.null(start.lvs)) {
      index <- mvtnorm::rmvnorm(N, rep(0, num.lv));
      colnames(index) <- paste("LV",1:num.lv, sep = "")
      unique.index <- as.matrix(index[unique.ind,])
    }

    if(!is.null(start.lvs)) {
      index <- as.matrix(start.lvs)
      unique.index <- as.matrix(index[unique.ind,])
    }
  }

  if(num.lv == 0) { index <- NULL }

  y <- as.matrix(y)

  if(family == "ordinal" && zeta.struc == "species") {
    max.levels <- apply(y,2,function(x) length(min(x):max(x)));
    if(any(max.levels == 1) || all(max.levels == 2)) stop("Ordinal data requires all columns to have at least has two levels. If al columns only have two levels, please use family == binomial instead. Thanks")
  }else if(family=="ordinal" && zeta.struc == "common"){
    max.levels=length(min(y):max(y))
  }

  if(is.null(rownames(y))) rownames(y) <- paste("row",1:N,sep="")
  if(is.null(colnames(y))) colnames(y) <- paste("col",1:p,sep="")

  options(warn = -1)

  if(family != "tweedie" && family!="ordinal") { ## Using logistic instead of prbit regession here for binomial, but whatever...
    if(starting.val=="res" && is.null(start.lvs) ){# && num.lv>0
      if(is.null(TR)){
        if(family!="gaussian") {
          if(!is.null(X)) fit.mva <- gllvm.TMB(y=y, X=X, family = family, num.lv=0, starting.val = "zero", sd.errors = FALSE, optimizer = "nlminb")#mvabund::manyglm(y ~ X, family = family, K = trial.size)
          if(is.null(X)) fit.mva <- gllvm.TMB(y=y, family = family, num.lv=0, starting.val = "zero", sd.errors = FALSE, optimizer = "nlminb")#mvabund::manyglm(y ~ 1, family = family, K = trial.size)
          coef <- cbind(fit.mva$params$beta0,fit.mva$params$Xcoef)
          fit.mva$phi <- fit.mva$params$phi
          resi <- NULL
          mu <- cbind(rep(1,n),fit.mva$X.design)%*%t(cbind(fit.mva$params$beta0, fit.mva$params$Xcoef))
          if(family %in% c("poisson", "negative.binomial", "gamma", "exponential")) {
            mu <- exp(mu)
          }
          if(family == "binomial") {
            mu <-  binomial(link = link)$linkinv(mu)
          }
        } else {
          if(!is.null(X)) fit.mva <- mlm(y, X = X)
          if(is.null(X)) fit.mva <- mlm(y)
          mu <- NULL
          resi <- fit.mva$residuals; resi[is.infinite(resi)] <- 0; resi[is.nan(resi)] <- 0
          coef <- t(fit.mva$coef)
          fit.mva$phi <- apply(fit.mva$residuals,2,sd)
        }
        gamma=NULL
        if(num.lv>0){
          lastart <- FAstart(mu=mu, family=family, y=y, num.lv = num.lv, phis=fit.mva$phi, resi=resi)
          gamma<-lastart$gamma
          index<-lastart$index
        }
      } else {
        n1 <- colnames(X)
        n2 <- colnames(TR)
        form1 <-NULL
        if(is.null(formula)){
          form1 <- "~ 1"
          for(m in 1:length(n2)){
            for(l in 1:length(n1)){
              ni <- paste(n1[l],n2[m],sep = "*")
              form1 <- paste(form1,ni,sep = "+")
            }
          }
          formula=form1
        }
        if(!TMB) fit.mva <- gllvm.VA(y, X = X, TR = TR, formula=formula(formula), family = family, num.lv = 0, Lambda.struc = "diagonal", trace = FALSE, plot = FALSE, sd.errors = FALSE, maxit = 1000, seed=seed,n.init=1,starting.val="zero",yXT=yXT)
        if(TMB) {
          fit.mva <- try(trait.TMB(y, X = X, TR = TR, formula=formula(formula), family = family, num.lv = 0, Lambda.struc = "diagonal", trace = FALSE, sd.errors = FALSE, maxit = 1000, seed=seed,n.init=1,starting.val="zero",yXT = yXT, row.eff = row.eff, diag.iter = 0, optimizer = "nlminb", beta0com=beta0com), silent = TRUE);
          fit.mva$method = "VA"
          if(!is.null(randomX) && !inherits(fit.mva, "try-error") && is.finite(fit.mva$logL)) {
            fit.mva <- trait.TMB(y, X = X, TR = TR, formula=formula(formula), family = family, num.lv = 0, Lambda.struc = "diagonal", trace = FALSE, sd.errors = FALSE, maxit = 1000, seed=seed,n.init=1,starting.val="zero",yXT = yXT, row.eff = row.eff, diag.iter = 0, optimizer = "nlminb", randomX = randomX, beta0com=beta0com, start.params = fit.mva);
          } else {fit.mva <- trait.TMB(y, X = X, TR = TR, formula=formula(formula), family = family, num.lv = 0, Lambda.struc = "diagonal", trace = FALSE, sd.errors = FALSE, maxit = 1000, seed=seed,n.init=1,starting.val="zero",yXT = yXT, row.eff = row.eff, diag.iter = 0, optimizer = "nlminb", randomX = randomX, beta0com=beta0com);}
          fit.mva$coef=fit.mva$params
          if(row.eff=="random") { # !!!!  
            sigma=c(max(fit.mva$params$sigma[1],sigma),fit.mva$params$sigma[-1])
            fit.mva$params$row.params <- fit.mva$params$row.params/sd(fit.mva$params$row.params)*sigma[1]
          }
          out$fitstart <- list(A=fit.mva$A, Ab=fit.mva$Ab, TMBfnpar=fit.mva$TMBfn$par) #params = fit.mva$params, 
        }
        if(!is.null(form1)){
          if(!is.null(fit.mva$coef$row.params)) row.params=fit.mva$coef$row.params
          env <- (fit.mva$coef$B)[n1]
          trait <- (fit.mva$coef$B)[n2]
          inter <- (fit.mva$coef$B)[!names(fit.mva$coef$B) %in% c(n1,n2)]
          B <- c(env,trait,inter)
        } else {
          B<-fit.mva$coef$B
          if(!is.null(fit.mva$coef$row.params)) row.params=fit.mva$coef$row.params
        }

        
        fit.mva$phi <- phi <- fit.mva$coef$phi
        ds.res <- matrix(NA, n, p)
        rownames(ds.res) <- rownames(y)
        colnames(ds.res) <- colnames(y)
        mu <- (matrix(fit.mva$X.design%*%fit.mva$coef$B,n,p)+ matrix(fit.mva$coef$beta0,n,p,byrow = TRUE))
        if(row.eff %in% c(TRUE, "random", "fixed")) {mu <- mu + row.params }
        if(!is.null(randomX)) {
          Br <- fit.mva$params$Br
          sigmaB <- fit.mva$params$sigmaB
          if(ncol(fit.mva$Xrandom)>1) sigmaij <- fit.mva$TMBfn$par[names(fit.mva$TMBfn$par)=="sigmaij"]#fit.mva$params$sigmaB[lower.tri(fit.mva$params$sigmaB)]
          mu <- mu + fit.mva$Xrandom%*%Br
        }
        if(family %in% c("poisson", "negative.binomial", "gamma", "exponential")) {
          mu <- exp(mu)
        }
        if(family == "binomial") {
          mu <-  binomial(link = link)$linkinv(mu)
        }
        
        gamma=NULL
        if(num.lv>0){
          lastart <- FAstart(mu, family=family, y=y, num.lv = num.lv, phis=fit.mva$phi)
          gamma<-lastart$gamma
          index<-lastart$index
        } 
      }
      
      if(is.null(TR)){params <- cbind(coef,gamma)
      } else { params <- cbind((fit.mva$coef$beta0),gamma)}
    } else  if(starting.val != "zero") {
      if(family!="gaussian") {
        if(is.null(TR)){
          
          if(!is.null(X) & num.lv > 0) fit.mva <- gllvm.TMB(y=y, X=cbind(X,index), family = family, num.lv=0, starting.val = "zero", sd.errors = FALSE, optimizer = "nlminb")
          if(is.null(X) & num.lv > 0) fit.mva <- gllvm.TMB(y=y, X=index, family = family, num.lv=0, starting.val = "zero", sd.errors = FALSE, optimizer = "nlminb")
          if(!is.null(X) & num.lv == 0) fit.mva <- gllvm.TMB(y=y, X=X, family = family, num.lv=0, starting.val = "zero", sd.errors = FALSE, optimizer = "nlminb")
          if(is.null(X) & num.lv == 0) fit.mva <- gllvm.TMB(y=y, X=NULL, family = family, num.lv=0, starting.val = "zero", sd.errors = FALSE, optimizer = "nlminb")

        } else {
          if(num.lv > 0) fit.mva <- gllvm.TMB(y=y, X=index, family = family, num.lv=0, starting.val = "zero", sd.errors = FALSE, optimizer = "nlminb")
          if(num.lv == 0) fit.mva <- gllvm.TMB(y=y, X=NULL, family = family, num.lv=0, starting.val = "zero", sd.errors = FALSE, optimizer = "nlminb")
          env  <-  rep(0,num.X)
          trait  <-  rep(0,num.T)
          inter <- rep(0, num.T * num.X)
          B <- c(env,trait,inter)
        }
        params <- cbind(fit.mva$params$beta0, fit.mva$params$Xcoef, fit.mva$params$theta)
        fit.mva$phi <- fit.mva$params$phi
        
      } else {
        if(is.null(TR)){
          if(!is.null(X) & num.lv > 0) fit.mva <- mlm(y, X = X, index = index)
          if(is.null(X) & num.lv > 0) fit.mva <- mlm(y, index = index)
          if(!is.null(X) & num.lv == 0) fit.mva <- mlm(y, X = X)
          if(is.null(X) & num.lv == 0) fit.mva <- mlm(y)
        } else {
          if(num.lv > 0) fit.mva <- mlm(y, index = index)
          if(num.lv == 0) fit.mva <- mlm(y)
          env  <-  rep(0,num.X)
          trait  <-  rep(0,num.T)
          inter <- rep(0, num.T * num.X)
          B <- c(env,trait,inter)
        }
        params <- t(fit.mva$coefficients)
        fit.mva$phi <- apply(fit.mva$residuals,2,sd)
      }
      #params <- t(fit.mva$coef) #!!
    } else {
      fit.mva<-list()
      fit.mva$phi <- rep(1, p)
    }
    }


  if(family == "negative.binomial") {
    phi <- fit.mva$phi  + 1e-5
  } else if(family %in% c("gaussian", "gamma")) {
    phi <- fit.mva$phi
  } else { phi <- NULL }
  
  if(family == "tweedie") {
    if(is.null(TR)){
      indeX <- cbind(index, X)
    }else{indeX <- cbind(index)}
    if(num.lv == 0) indeX <- X

    phi0<-NULL
    coefs0<-NULL
    for(j in 1:p){
      fitj <- try(glm( y[,j] ~ indeX, family=statmod::tweedie(var.power=power, link.power=0) ),silent = TRUE)
      if(is.null(X) && num.lv == 0) fitj<-try(glm( y[,j] ~ 1, family=statmod::tweedie(var.power=power, link.power=0) ),silent = TRUE)
      if(!inherits(fitj, "try-error")){
        coefs0<-rbind(coefs0,fitj$coefficients)
        phi0<-c(phi0,summary(fitj)$dispersion)
      } else {coefs0<-rbind(coefs0,rnorm(num.lv+dim(X)[2]+1)); phi0<-c(phi0,runif(1,0.2,3));}
    }

    if(!is.null(TR)) B=rep(0,num.X+num.T+num.X*num.T)
    params <- as.matrix(coefs0)
    phi <- phi0
  }

  if(family == "ordinal") {
    max.levels <- length(unique(c(y)))
    params <- matrix(NA,p,ncol(cbind(1,X))+num.lv)
    if(zeta.struc == "species"){
      zeta <- matrix(NA,p,max.levels - 1)
      zeta[,1] <- 0 ## polr parameterizes as no intercepts and all cutoffs vary freely. Change this to free intercept and first cutoff to zero
    }else{
      cw.fit <- MASS::polr(factor(y) ~ 1, method = "probit")
      zeta <- cw.fit$zeta
      zeta[1] <- 0
    }
    
    for(j in 1:p) {
      y.fac <- factor(y[,j])
      if(length(levels(y.fac)) > 2) {
        if(starting.val%in%c("zero","res") || num.lv==0){
          if(is.null(X) || !is.null(TR)) cw.fit <- MASS::polr(y.fac ~ 1, method = "probit")
          if(!is.null(X) & is.null(TR) ) cw.fit <- MASS::polr(y.fac ~ X, method = "probit")
        } else {
          if(is.null(X) || !is.null(TR)) cw.fit <- MASS::polr(y.fac ~ index, method = "probit")
          if(!is.null(X) & is.null(TR) & num.lv > 0) cw.fit <- MASS::polr(y.fac ~ X+index, method = "probit")
        }
        if(starting.val=="random"){
          params[j,] <- c(cw.fit$zeta[1],-cw.fit$coefficients)
          if(zeta.struc == "species"){
            zeta[j,2:length(cw.fit$zeta)] <- cw.fit$zeta[-1]-cw.fit$zeta[1]
          }
        }else{
          params[j,1:ncol(cbind(1,X))] <- c(cw.fit$zeta[1],-cw.fit$coefficients)
          if(zeta.struc == "species"){
            zeta[j,2:length(cw.fit$zeta)] <- cw.fit$zeta[-1]-cw.fit$zeta[1]
          }
        }
      }
      if(length(levels(y.fac)) == 2) {
        if(starting.val%in%c("zero","res") || num.lv==0){
          if(is.null(X) || !is.null(TR)) cw.fit <- glm(y.fac ~ 1, family = binomial(link = "probit"))
          if(!is.null(X) & is.null(TR) ) cw.fit <- glm(y.fac ~ X, family = binomial(link = "probit"))
        } else {
          if(is.null(X) || !is.null(TR)) cw.fit <- glm(y.fac ~ index, family = binomial(link = "probit"))
          if(!is.null(X) & is.null(TR) & num.lv > 0) cw.fit <- glm(y.fac ~ X+index, family = binomial(link = "probit"))
        }
        params[j,1:length(cw.fit$coef)] <- cw.fit$coef
      }
    }
    if(starting.val%in%c("res") && num.lv>0){
      eta.mat <- matrix(params[,1],n,p,byrow=TRUE)
      if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + (X %*% matrix(params[,2:(1+num.X)],num.X,p))
      mu <- eta.mat
      
      lastart <- FAstart(eta.mat, family=family, y=y, num.lv = num.lv, zeta = zeta, zeta.struc = zeta.struc)
      gamma<-lastart$gamma
      index<-lastart$index
      params[,(ncol(cbind(1,X))+1):ncol(params)]=gamma
    }
    env <- rep(0,num.X)
    trait <- rep(0,num.T)
    inter <- rep(0, num.T * num.X)
    B=c(env,trait,inter)
  }

  if((family!="ordinal" || (family=="ordinal" & starting.val=="res")) & starting.val!="zero"){
    if(num.lv>1 && p>2){
      gamma<-as.matrix(params[,(ncol(params) - num.lv + 1):ncol(params)])
      qr.gamma <- qr(t(gamma))
      params[,(ncol(params) - num.lv + 1):ncol(params)]<-t(qr.R(qr.gamma))
      index<-(index%*%qr.Q(qr.gamma))
    }}
  if(starting.val=="zero"){
    params=matrix(0,p,1+num.X+num.lv)
    params[,1:(ncol(params) - num.lv)] <- 0
    env <- rep(0,num.X)
    trait <- rep(0,num.T)
    inter <- rep(0, num.T * num.X)
    B=c(env,trait,inter)
    if(num.lv > 0) {
      gamma <- matrix(1,p,num.lv)
      gamma[upper.tri(gamma)]=0
      params[,(ncol(params) - num.lv + 1):ncol(params)] <- gamma
      index <- matrix(0,n,num.lv)
    }
    phi <- rep(1,p)
  }
  if(num.lv > 0) {
    index <- index+mvtnorm::rmvnorm(n, rep(0, num.lv),diag(num.lv)*jitter.var);
    try({
      gamma.new <- as.matrix(params[,(ncol(params) - num.lv + 1):ncol(params)]);
      sig <- sign(diag(gamma.new));
      params[,(ncol(params) - num.lv + 1):ncol(params)] <- t(t(gamma.new)*sig)
      index <- t(t(index)*sig)}, silent = TRUE)
  }
  
  if(num.lv == 0){
    qqs <- quantile(params)
    uppout <- qqs[4] + (qqs[4]- qqs[2])*2
    lowout <- qqs[2] - (qqs[4]- qqs[2])*2
    params[params < lowout] <- lowout
    params[params > uppout] <- uppout
  }

  out$params <- params
  out$phi <- phi
  out$mu <- mu
  if(!is.null(TR)) { out$B <- B}
  if(num.lv > 0) out$index <- index
  if(family == "ordinal") out$zeta <- zeta
  options(warn = 0)
  if(row.eff!=FALSE) {
    out$row.params=row.params
    if(row.eff=="random") out$sigma=sigma
  }
  if(!is.null(randomX)){
    out$Br <- Br
    out$sigmaB <- sigmaB
    out$sigmaij <- sigmaij
  }
  
  return(out)
}


FAstart <- function(mu, family, y, num.lv, zeta = NULL, zeta.struc = NULL, phis = NULL, 
                    jitter.var = 0, resi = NULL, row.eff = FALSE){
  n<-NROW(y); p <- NCOL(y)
  row.params <- NULL # !!!!
  
  if(is.null(resi)){
    ds.res <- matrix(NA, n, p)
    rownames(ds.res) <- rownames(y)
    colnames(ds.res) <- colnames(y)
    for (i in 1:n) {
      for (j in 1:p) {
        if (family == "poisson") {
          a <- ppois(as.vector(unlist(y[i, j])) - 1, mu[i,j])
          b <- ppois(as.vector(unlist(y[i, j])), mu[i,j])
          u <- runif(n = 1, min = a, max = b)
          ds.res[i, j] <- qnorm(u)
        }
        if (family == "negative.binomial") {
          phis <- phis + 1e-05
          a <- pnbinom(as.vector(unlist(y[i, j])) - 1, mu = mu[i, j], size = 1/phis[j])
          b <- pnbinom(as.vector(unlist(y[i, j])), mu = mu[i, j], size = 1/phis[j])
          u <- runif(n = 1, min = a, max = b)
          ds.res[i, j] <- qnorm(u)
        }
        if (family == "binomial") {
          a <- pbinom(as.vector(unlist(y[i, j])) - 1, 1, mu[i, j])
          b <- pbinom(as.vector(unlist(y[i, j])), 1, mu[i, j])
          u <- runif(n = 1, min = a, max = b)
          ds.res[i, j] <- qnorm(u)
        }
        if (family == "gaussian") {
          a <- pnorm(as.vector(unlist(y[i, j])), mu[i, j], sd = phis[j])
          b <- pnorm(as.vector(unlist(y[i, j])), mu[i, j], sd = phis[j])
          u <- runif(n = 1, min = a, max = b)
          ds.res[i, j] <- qnorm(u)
        }
        if (family == "gamma") {
          a <- pgamma(as.vector(unlist(y[i, j])), shape = phis[j], scale = mu[i, j]/phis[j])
          b <- pgamma(as.vector(unlist(y[i, j])), shape = phis[j], scale = mu[i, j]/phis[j])
          u <- runif(n = 1, min = a, max = b)
          ds.res[i, j] <- qnorm(u)
        }
        if (family == "exponential") {
          a <- pexp(as.vector(unlist(y[i, j])), rate = 1/mu[i, j])
          b <- pexp(as.vector(unlist(y[i, j])), rate = 1/mu[i, j])
          u <- runif(n = 1, min = a, max = b)
          ds.res[i, j] <- qnorm(u)
        }
      if (family == "ordinal") {
      if(zeta.struc == "species"){
        probK <- NULL
        probK[1] <- pnorm(zeta[j,1]-mu[i,j],log.p = FALSE)
        probK[max(y[,j])] <- 1 - pnorm(zeta[j,max(y[,j]) - 1] - mu[i,j])
        if(max(y[,j]) > 2) {
          j.levels <- 2:(max(y[,j])-1)
          for(k in j.levels) { probK[k] <- pnorm(zeta[j,k] - mu[i,j]) - pnorm(zeta[j,k - 1] - mu[i,j]) }
        }
        probK <- c(0,probK)
        cumsum.b <- sum(probK[1:(y[i,j]+1)])
        cumsum.a <- sum(probK[1:(y[i,j])])
        u <- runif(n = 1, min = cumsum.a, max = cumsum.b)
        if (abs(u - 1) < 1e-05)
          u <- 1
        if (abs(u - 0) < 1e-05)
          u <- 0
        ds.res[i, j] <- qnorm(u)
      }else{
        probK <- NULL
        probK[1] <- pnorm(zeta[1] - mu[i, j], log.p = FALSE)
        probK[max(y)] <- 1 - pnorm(zeta[max(y) - 1] - mu[i,j])
          levels <- 2:(max(y) - min(y))#
          for (k in levels) {
            probK[k] <- pnorm(zeta[k] - mu[i, j]) - pnorm(zeta[k - 1] - mu[i, j])
          }
        probK <- c(0, probK)
        cumsum.b <- sum(probK[1:(y[i, j] + 2 - min(y))])
        cumsum.a <- sum(probK[1:(y[i, j])])
        u <- runif(n = 1, min = cumsum.a, max = cumsum.b)
        if (abs(u - 1) < 1e-05)
          u <- 1
        if (abs(u - 0) < 1e-05)
          u <- 0
        ds.res[i, j] <- qnorm(u)
      }
  }
}
  }
    
  } else {
    ds.res <- resi
  }
  resi <- as.matrix(ds.res); resi[is.infinite(resi)] <- 0; resi[is.nan(resi)] <- 0
  
  if(p>2 && n>2){
    if(any(is.nan(resi))){stop("Method 'res' for starting values can not be used, when glms fit too poorly to the data. Try other starting value methods 'zero' or 'random' or change the model.")}
    
    if(n>p){
      fa  <-  try(factanal(resi,factors=num.lv,scores = "regression"))
      if(inherits(fa,"try-error")) stop("Factor analysis for calculating starting values failed. Maybe too many latent variables. Try smaller 'num.lv' value or change 'starting.val' to 'zero' or 'random'.")
      gamma<-matrix(fa$loadings,p,num.lv)
      index <- fa$scores
    } else if(n<p) {
      fa  <-  try(factanal(t(resi),factors=num.lv,scores = "regression"))
      if(inherits(fa,"try-error")) stop("Factor analysis for calculating starting values failed. Maybe too many latent variables. Try smaller 'num.lv' value or change 'starting.val' to 'zero' or 'random'.")
      gamma<-fa$scores
      index <- matrix(fa$loadings,n,num.lv)
    } else {
      tryfit <- TRUE; tryi <- 1
      while(tryfit && tryi<5) {
        fa  <-  try(factanal(rbind(resi,rnorm(p,0,0.01)),factors=num.lv,scores = "regression"), silent = TRUE)
        tryfit <- inherits(fa,"try-error"); tryi <- tryi + 1;
      }
      if(tryfit) {
        warning(attr(fa,"condition")$message, "\n Factor analysis for calculating starting values failed. Maybe too many latent variables. Try smaller 'num.lv' value or change 'starting.val' to 'zero' or 'random'. Using solution from Principal Component Analysis instead./n")
        pr <- princomp(resi)
        gamma<-matrix(pr$loadings[,1:num.lv],p,num.lv)
        index<-matrix(pr$scores[,1:num.lv],n,num.lv)
      }else{
        gamma<-matrix(fa$loadings,p,num.lv)
        index <- fa$scores[1:n,]
      }
    }
  } else {
    gamma <- matrix(1,p,num.lv)
    gamma[upper.tri(gamma)]=0
    index <- matrix(0,n,num.lv)
  }
  
  if(row.eff == "random") { # !!!!
    row.params <- index[,num.lv]* (sd(matrix(index[,num.lv])%*%matrix(gamma[, num.lv], nrow=1)))/sd(index[,num.lv])
    
    index <- as.matrix(index[,-num.lv])
    gamma <- as.matrix(gamma[,-num.lv])
    num.lv <- num.lv-1
  }
  
  if(num.lv>1 && p>2){
    qr.gamma <- qr(t(gamma))
    gamma.new<-t(qr.R(qr.gamma))
    sig <- sign(diag(gamma.new));
    gamma <- t(t(gamma.new)*sig)
    index<-(index%*%qr.Q(qr.gamma))
    index <- t(t(index)*sig)
  } else {
    sig <- sign(diag(gamma));
    gamma <- t(t(gamma)*sig)
    index <- t(t(index)*sig)
  }
  if(p>n) {
    sdi <- sqrt(diag(cov(index)))
    sdt <- sqrt(diag(cov(gamma)))
    indexscale <- diag(x = 0.8/sdi, nrow = length(sdi))
    index <- index%*%indexscale
    gammascale <- diag(x = 1/sdt, nrow = length(sdi))
    gamma <- gamma%*%gammascale
  } else {
    sdi <- sqrt(diag(cov(index)))
    sdt <- sqrt(diag(cov(gamma)))
    indexscale <- diag(x = 1/sdi, nrow = length(sdi))
    index <- index%*%indexscale
    gammascale <- diag(x = 0.8/sdt, nrow = length(sdi))
    gamma <- gamma%*%gammascale
  }
  index <- index + mvtnorm::rmvnorm(n, rep(0, num.lv),diag(num.lv)*jitter.var);
  return(list(index = index, gamma = gamma, row.params =row.params))
}




## Calculates information matrix based on scores of VA log-likelihood

calc.infomat <- function(theta = NULL, beta0=NULL, env = NULL, row.params = NULL, vameans = NULL, lambda=NULL, phi = NULL, zeta = NULL, num.lv, family, Lambda.struc, row.eff, y, X = NULL,TR = NULL, Xd=NULL, offset=0, B=NULL) {
  n <- nrow(y); p <- ncol(y)
  if(!is.null(X)) num.X <- ncol(X)
  if(!is.null(TR)) num.T <- ncol(TR)
  trial.size <- 1
  if(family == "ordinal") {
    max.levels <- max(y)
    if(any(max.levels == 1) || all(max.levels == 2)) stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks")
  }
  if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)
  if (!is.null(Xd))  nd <- ncol(Xd)

  ## score of VA logL wrt model parameters
  grad.all.mod <- function(x,vameans=NULL,lambda=NULL) {

    x2 <- x
    new.vameans <- new.lambda <- new.theta <- NULL
    if(num.lv > 0) {
      new.lambda <- lambda
      new.vameans <- vameans
      new.theta <- matrix(c(x2[1:(p * num.lv)]),p,num.lv); x2 <- x2[-(1:(p * num.lv))]
      new.theta[upper.tri(new.theta)] <- 0
    }
    new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
    new.env <- NULL
    if(!is.null(X)) {
      if(is.null(TR)) { new.env <- matrix(x2[1:(p * num.X)],p,num.X); x2 <- x2[-(1:(p * num.X))]
      } else { B=x2[1:nd]; x2 <- x2[-(1:nd)]}
    }
    new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
    new.phi <- NULL; if(family == "negative.binomial") { new.phi <- x2[1:p]; x2 <- x2[-(1:p)] }
    new.zeta <- NULL; if(family == "ordinal") {
      new.zeta <- matrix(NA,p,max.levels - 2)
      for(j in 1:p) { if(max(y[,j]) > 2) { new.zeta[j,1:(max(y[,j])-2)] <- x2[1:(max(y[,j])-2)]; x2 <- x2[-(1:(max(y[,j])-2))] } }
      new.zeta <- cbind(0,new.zeta);
    }

    eta.mat <- matrix(new.beta0,n,p,byrow=TRUE)  + offset
    if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
    if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% B),n,p)
    if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta)
    if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)

    grad.theta <- grad.beta0 <- grad.env <- grad.B <- grad.phi <- grad.row.params <- grad.zeta <- NULL;

    if(Lambda.struc == "unstructured" & num.lv > 0) {
      Lambda.theta <- array(NA,dim=c(p,n,num.lv))
      for(i2 in 1:n) { Lambda.theta[,i2,] <- new.theta%*%new.lambda[i2,,] }
    }

    if(family=="poisson") {
      if(num.lv > 0) eta.mat <- eta.mat + calc.quad(new.lambda,new.theta,Lambda.struc)$mat
      grad.beta0 <- colSums(y-exp(eta.mat))
      if(num.lv > 0) {
        for(l in 1:num.lv) {
          if(Lambda.struc == "unstructured") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep(t(Lambda.theta[,,l]),1,new.vameans[,l],"+") * exp(eta.mat)
          if(Lambda.struc == "diagonal") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep((new.lambda[,l]%*%t(new.theta[,l])),1,new.vameans[,l],"+")*exp(eta.mat)
          grad.theta <- c(grad.theta,colSums(sum1))
        }
      }
      if(!is.null(X) && is.null(TR)) {
        for (l in 1:num.X) {
          sum1 <- sweep(y-exp(eta.mat),1,X[,l],"*")
          grad.env <- c(grad.env,colSums(sum1))
        }
      }
      if(!is.null(TR)) {
        sum1 <- c(y-exp(eta.mat)) * Xd
        grad.B <- c(grad.B,colSums(sum1))
      }
      if(row.eff) { grad.row.params <- rowSums(y - exp(eta.mat)) }
    }

    if(family=="negative.binomial") {
      eta.mat <- matrix(new.beta0,n,p,byrow=TRUE) + offset
      if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
      if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% B),n,p)
      if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
      if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta)

      phi.mat <- matrix(new.phi,n,p,byrow=TRUE)
      mu.mat.noquad <- eta.mat
      if(num.lv > 0) eta.mat <- eta.mat - calc.quad(new.lambda,new.theta,Lambda.struc)$mat
      grad.beta0 <- colSums(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(eta.mat) + phi.mat)))
      if(num.lv > 0) {
        for(l in 1:num.lv) {
          if(Lambda.struc == "unstructured") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep(t(Lambda.theta[,,l]),1,-new.vameans[,l],"+") * (-(y + phi.mat) * (phi.mat / (exp(eta.mat) + phi.mat)))
          if(Lambda.struc == "diagonal") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep((new.lambda[,l] %*% t(new.theta[,l])),1,-new.vameans[,l],"+") * (-(y + phi.mat) * (phi.mat / (exp(eta.mat) + phi.mat)))
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

      grad.phi <- colSums(-mu.mat.noquad + 1 + log(phi.mat) - digamma(phi.mat) - log(1 + phi.mat / exp(eta.mat)) - (y + phi.mat)/(phi.mat + exp(eta.mat)) + digamma(y + phi.mat))

    }

    if(family=="binomial") {
      grad.beta0 <- colSums(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat))),na.rm=TRUE)

      if(num.lv > 0) {
        for(l in 1:num.lv) {
          if(Lambda.struc == "unstructured") sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat))/(pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
          if(Lambda.struc == "diagonal") sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,new.vameans[,l],"*") - (new.lambda[,l] %*% t(new.theta[,l]))
          grad.theta <- c(grad.theta,colSums(sum1,na.rm=TRUE)) }
      }

      if(!is.null(X) && is.null(TR)) {
        for(l in 1:num.X) {
          sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,X[,l],"*")
          grad.env <- c(grad.env,colSums(sum1))
        }
      }
      if(!is.null(TR)) {
        sum1 <- c(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10)) * Xd
        grad.B <- c(grad.B,colSums(sum1))
      }
      if(row.eff) { grad.row.params <- rowSums((dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)))),na.rm=TRUE) }
    }

    if(family=="ordinal") {
      grad.zeta <- matrix(0,p,ncol(new.zeta))
      deriv.trunnorm <- matrix(0,n,p)
      for(j in 1:p) {
        deriv.trunnorm[y[,j] == 1,j] <- -dnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j]) / pnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j])
        deriv.trunnorm[y[,j] == max(y[,j]),j] <- dnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]) / (1 - pnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]))
        j.levels <- (1:max(y[,j]))[-c(1,max(y[,j]))]
        for(k in j.levels) {
          deriv.trunnorm[y[,j] == k,j] <- (-dnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) + dnorm(new.zeta[j,k - 1] - eta.mat[y[,j] == k,j])) / (pnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) - pnorm(new.zeta[j,k-1] - eta.mat[y[,j] == k,j])) }
      }
      deriv.trunnorm[!is.finite(deriv.trunnorm)] <- 0
      grad.beta0 <- colSums(deriv.trunnorm,na.rm=TRUE)
      if(num.lv > 0) {
        for(l in 1:num.lv) {
          if(Lambda.struc == "unstructured") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
          if(Lambda.struc == "diagonal") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - (new.lambda[,l] %*% t(new.theta[,l]))
          grad.theta <- c(grad.theta,colSums(sum1,na.rm=TRUE))
        }
      }
      if(!is.null(X) && is.null(TR)) {
        for(l in 1:num.X) { sum1 <- sweep(deriv.trunnorm,1,X[,l],"*");
        grad.env <- c(grad.env,colSums(sum1,na.rm=TRUE))
        }
      }
      if(!is.null(TR)) {
        sum1 <- c(deriv.trunnorm) * Xd
        grad.B <- c(grad.B,colSums(sum1))
      }

      if(row.eff) { grad.row.params <- rowSums(deriv.trunnorm,na.rm=TRUE) }

      for(j in 1:p) {
        zeta0 <- new.zeta[j,]
        grad.zeta[j,max(y[,j]) - 1] <- grad.zeta[j,max(y[,j]) - 1] - sum(dnorm(zeta0[max(y[,j]) - 1] - eta.mat[which(y[,j] == max(y[,j])) ,j]) / (1 - pnorm(zeta0[max(y[,j]) - 1] - eta.mat[which(y[,j] == max(y[,j])),j])))
        j.levels <- (1:max(y[,j]))[-c(1,max(y[,j]))]
        for(k in j.levels) {
          grad.zeta[j,k] <- grad.zeta[j,k] + sum(dnorm(zeta0[k] - eta.mat[y[,j] == k,j]) / (pnorm(zeta0[k] - eta.mat[y[,j] == k,j]) - pnorm(zeta0[k - 1] - eta.mat[y[,j] == k,j])))
          grad.zeta[j,k - 1] <- grad.zeta[j,k - 1] - sum(dnorm(zeta0[k - 1] - eta.mat[y[,j] == k,j]) / (pnorm(zeta0[k] - eta.mat[y[,j] == k,j]) - pnorm(zeta0[k - 1] - eta.mat[y[,j] == k,j])))
        }
      }
      grad.zeta <- grad.zeta[,-1]
      if(length(unique(apply(y,2,max))) > 1) { grad.zeta <- c(grad.zeta[-which(grad.zeta == 0)]) }
    }

    if(num.lv>0){
      grad.theta=matrix(grad.theta,nrow=p,ncol=num.lv)
      grad.theta[upper.tri(grad.theta)] <- 0}
    if(row.eff) grad.row.params[1]=0;

    return(c(grad.theta, grad.beta0, grad.env, grad.B, grad.row.params, grad.phi, grad.zeta))
  }


  ## score of VA logL wrt vameans_sel.i and lambda_sel.i
  grad.all.var <- function(x,vameans=NULL,lambda=NULL,mod.x,sel.i) {
    x2 <- mod.x
    if(num.lv > 0) {
      new.vameans <- vameans; new.lambda <- lambda
      new.theta <- matrix(c(x2[1:(p * num.lv)]),p,num.lv); x2 <- x2[-(1:(p * num.lv))]
      new.theta[upper.tri(new.theta)] <- 0
    }
    new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]; new.env <- NULL
    if(!is.null(X)) {
      if(is.null(TR)) { new.env <- matrix(x2[1:(p * num.X)],p,num.X); x2 <- x2[-(1:(p * num.X))]
      } else { B=x2[1:nd]; x2 <- x2[-(1:nd)]}
    }
    new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
    new.phi <- NULL; if(family == "negative.binomial") { new.phi <- x2[1:p]; x2 <- x2[-(1:p)] }
    new.zeta <- NULL
    if(family == "ordinal") {
      new.zeta <- matrix(NA,p,max.levels-2);
      for(j in 1:p) { if(max(y[,j]) > 2) { new.zeta[j,1:(max(y[,j])-2)] <- x2[1:(max(y[,j])-2)]; x2 <- x2[-(1:(max(y[,j])-2))] } }
      new.zeta <- cbind(0,new.zeta);
    }

    ## Replace vameans and lambda with elements in x
    if(num.lv > 0) {
      new.vameans[sel.i,] <- x[1:num.lv]; va.x <- x[-(1:num.lv)]
      if(Lambda.struc == "diagonal") { new.lambda[sel.i,] <- va.x[1:(num.lv)]; va.x <- va.x[-(1:(num.lv))] }
      if(Lambda.struc == "unstructured") { ## Rebuilt lambda array from unique elements
        new.lambda[sel.i,,][lower.tri(new.lambda[sel.i,,],diag=TRUE)] <- va.x[1:(num.lv + num.lv * (num.lv - 1) / 2)]; va.x <- va.x[-(1:(num.lv + num.lv * (num.lv - 1) / 2))]
        new.lambda[sel.i,,][upper.tri(new.lambda[sel.i,,])] <- new.lambda[sel.i,,][lower.tri(new.lambda[sel.i,,])]
      }
    }

    grad.vameans <- NULL
    grad.lambda <- matrix(0,num.lv,num.lv)
    eta.mat <- matrix(new.beta0,n,p,byrow=TRUE) + new.vameans %*% t(new.theta) + offset
    if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
    if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% B),n,p)
    if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)

    if(family=="poisson") {
      eta.mat <- eta.mat + calc.quad(new.lambda,new.theta,Lambda.struc)$mat
      sum1 <- (y[sel.i,]-exp(eta.mat[sel.i,]))
      for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,sum(sum1 * new.theta[,l]) - new.vameans[sel.i,l]) }
    }

    if(family=="negative.binomial") {
      eta.mat <- eta.mat - calc.quad(new.lambda,new.theta,Lambda.struc)$mat
      sum1 <- -new.phi - (y[sel.i,] + new.phi) * (-new.phi / (exp(eta.mat[sel.i,]) + new.phi))
      for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,sum(sum1 * new.theta[,l]) - new.vameans[sel.i,l]) }
    }

    if(family=="binomial") {
      sum1 <- dnorm(eta.mat[sel.i,]) * (y[sel.i,]-trial.size * pnorm(eta.mat[sel.i,])) / (pnorm(eta.mat[sel.i,]) * (1 - pnorm(eta.mat[sel.i,])) + 1e-10)
      for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,sum(sum1 * new.theta[,l]) - new.vameans[sel.i,l]) }
    }

    if(family=="ordinal") {
      deriv.trunnorm <- matrix(NA,n,p)
      for(j in 1:p) {
        deriv.trunnorm[y[,j] == 1,j] <- -dnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j]) / pnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j])
        deriv.trunnorm[y[,j] == max(y[,j]),j] <- dnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]) / (1 - pnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]))
        if(max(y[,j])>2) {
          j.levels <- 2:(max(y[,j])-1)
          for(k in j.levels) { deriv.trunnorm[y[,j] == k,j] <- (-dnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) + dnorm(new.zeta[j,k-1] - eta.mat[y[,j] == k,j])) / (pnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) - pnorm(new.zeta[j,k - 1] - eta.mat[y[,j] == k,j])) }
        }
      }
      deriv.trunnorm[!is.finite(deriv.trunnorm)] <- 0

      for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,sum(deriv.trunnorm[sel.i,] * new.theta[,l]) - new.vameans[sel.i,l]) }
    }

    if(family %in% c("poisson","negative.binomial")) {
      if(family == "poisson") { deriv1 <- exp(eta.mat[sel.i,]) }
      if(family == "negative.binomial") {
        eta.mat <- matrix(new.beta0,n,p,byrow=TRUE) + offset + new.vameans %*% t(new.theta) - calc.quad(new.lambda,new.theta,Lambda.struc)$mat
        if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
        if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% B),n,p)
        if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
        deriv1 <- (y[sel.i,] + new.phi) * (new.phi / (exp(eta.mat[sel.i,]) + new.phi))
      }

      if(Lambda.struc == "unstructured") {
        theta2 <- sapply(1:p,function(j,theta) theta[j,]%*%t(theta[j,]), theta = new.theta)
        theta2 <- t(theta2)
        grad.lambda <- -0.5 * matrix(apply(deriv1 * theta2,2,sum),nrow=num.lv) + 0.5 * solve(new.lambda[sel.i,,]) - 0.5 * diag(nrow=num.lv)
      }

      if(Lambda.struc == "diagonal") {
        theta2 <- new.theta^2
        grad.lambda <- -0.5 * diag(apply(deriv1 * theta2,2,sum),nrow=num.lv) + 0.5 * solve(diag(x = new.lambda[sel.i,],nrow=num.lv)) - 0.5 * diag(nrow=num.lv)
      }
    }

    if(family %in% c("binomial","ordinal")) {
      if(Lambda.struc == "unstructured") {
        theta2 <- sapply(1:p,function(j,theta) theta[j,] %*% t(theta[j,]), theta = new.theta)
        theta2 <- t(theta2)
        grad.lambda <- -0.5 * matrix(apply(theta2,2,sum),nrow=num.lv) + 0.5 * solve(new.lambda[sel.i,,]) - 0.5 * diag(nrow=num.lv)
      }
      if(Lambda.struc == "diagonal") {
        theta2 <- new.theta^2
        grad.lambda <- -0.5 * diag(apply(theta2,2,sum),nrow=num.lv) + 0.5 * solve(diag(x=new.lambda[sel.i,],nrow=num.lv)) - 0.5 * diag(nrow=num.lv)
      }
    }

    grad.lambda2 <- NULL
    if(Lambda.struc == "diagonal") grad.lambda2 <- diag(x=as.matrix(grad.lambda)) ## Only extract diagonal elements
    if(Lambda.struc == "unstructured") grad.lambda2 <- grad.lambda[lower.tri(grad.lambda,diag=TRUE)]

    return(c(grad.vameans,grad.lambda2))
  }

  ## cross derivatives of VA logL - taking VA parameters vameans_sel.i and lambda_sel.i and diffing wrt to model parameters
  grad.all.cross <- function(x,vameans=NULL,lambda=NULL,mod.x,sel.i) {
    x2 <- mod.x
    if(num.lv > 0) {
      new.vameans <- vameans; new.lambda <- lambda
      new.theta <- matrix(c(x2[1:(p * num.lv)]),p,num.lv); x2 <- x2[-(1:(p * num.lv))]
      new.theta[upper.tri(new.theta)] <- 0
    }
    new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
    new.env <- NULL
    if(!is.null(X)) {
      if(is.null(TR)) { new.env <- matrix(x2[1:(p * num.X)],p,num.X); x2 <- x2[-(1:(p * num.X))]
      } else { B=x2[1:nd]; x2 <- x2[-(1:nd)]}
    }
    new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
    new.phi <- NULL; if(family == "negative.binomial") { new.phi <- x2[1:p]; x2 <- x2[-(1:p)] }
    new.zeta <- NULL
    if(family == "ordinal") {
      new.zeta <- matrix(NA,p,max.levels - 2);
      for(j in 1:p) { if(max(y[,j]) > 2) { new.zeta[j,1:(max(y[,j])-2)] <- x2[1:(max(y[,j])-2)]; x2 <- x2[-(1:(max(y[,j]) - 2))] } }
      new.zeta <- cbind(0,new.zeta) }

    ## Replace vameans and lambda with elements in x
    if(num.lv > 0) {
      new.vameans[sel.i,] <- x[1:num.lv]; va.x <- x[-(1:num.lv)]
      if(Lambda.struc == "diagonal") { new.lambda[sel.i,] <- va.x[1:(num.lv)]; va.x <- va.x[-(1:(num.lv))] }
      if(Lambda.struc == "unstructured") { ## Rebuilt lambda array from unique elements
        new.lambda[sel.i,,][lower.tri(new.lambda[sel.i,,],diag=TRUE)] <- va.x[1:(num.lv+num.lv * (num.lv - 1) / 2)]; va.x <- va.x[-(1:(num.lv + num.lv * (num.lv - 1) / 2))]
        new.lambda[sel.i,,][upper.tri(new.lambda[sel.i,,])] <- new.lambda[sel.i,,][lower.tri(new.lambda[sel.i,,])]
      }
    }

    eta.mat <- matrix(new.beta0,n,p,byrow=TRUE) + offset
    if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
    if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% B),n,p)
    if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
    if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta)
    grad.theta <- grad.beta0 <- grad.env <- grad.B <- grad.phi <- grad.row.params <- grad.zeta <- NULL

    if(row.eff) grad.row.params <- numeric(n)
    if(Lambda.struc == "unstructured" & num.lv > 0) {
      Lambda.theta <- array(NA,dim=c(p,n,num.lv))
      for(i2 in 1:n) { Lambda.theta[,i2,] <- new.theta %*% new.lambda[i2,,] }
    }

    if(family=="poisson") {
      if(num.lv > 0) eta.mat <- eta.mat + calc.quad(new.lambda,new.theta,Lambda.struc)$mat
      grad.beta0 <- colSums(y-exp(eta.mat))
      if(num.lv > 0) {
        for(l in 1:num.lv) {
          if(Lambda.struc == "unstructured") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep(t(Lambda.theta[,,l]),1,new.vameans[,l],"+") * exp(eta.mat)
          if(Lambda.struc == "diagonal") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep((new.lambda[,l]%*%t(new.theta[,l])),1,new.vameans[,l],"+")*exp(eta.mat)
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
      if(row.eff) { grad.row.params <- rowSums(y - exp(eta.mat)) }
    }


    if(family=="negative.binomial") {
      eta.mat <- matrix(new.beta0,n,p,byrow=TRUE) + offset
      if(!is.null(X) && is.null(TR)) eta.mat <- eta.mat + X %*% t(new.env)
      if(!is.null(TR)) eta.mat <- eta.mat + matrix((Xd %*% B),n,p)
      if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
      if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta)

      phi.mat <- matrix(new.phi,n,p,byrow=TRUE)
      mu.mat.noquad <- eta.mat
      if(num.lv > 0) eta.mat <- eta.mat - calc.quad(new.lambda,new.theta,Lambda.struc)$mat
      grad.beta0 <- colSums(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(eta.mat) + phi.mat)))
      if(num.lv > 0) {
        for(l in 1:num.lv) {
          if(Lambda.struc == "unstructured") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep(t(Lambda.theta[,,l]),1,-new.vameans[,l],"+") * (-(y + phi.mat) * (phi.mat / (exp(eta.mat) + phi.mat)))
          if(Lambda.struc == "diagonal") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep((new.lambda[,l] %*% t(new.theta[,l])),1,-new.vameans[,l],"+") * (-(y + phi.mat) * (phi.mat / (exp(eta.mat) + phi.mat)))
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

      grad.phi <- colSums(-mu.mat.noquad + 1 + log(phi.mat) - digamma(phi.mat) - log(1 + phi.mat / exp(eta.mat)) - (y + phi.mat)/(phi.mat + exp(eta.mat)) + digamma(y + phi.mat))
    }

    if(family=="binomial") {
      grad.beta0 <- colSums(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat))),na.rm=TRUE)

      if(num.lv > 0) {
        for(l in 1:num.lv) {
          if(Lambda.struc == "unstructured") sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat))/(pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
          if(Lambda.struc == "diagonal") sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,new.vameans[,l],"*") - (new.lambda[,l] %*% t(new.theta[,l]))
          grad.theta <- c(grad.theta,colSums(sum1,na.rm=TRUE)) }
      }

      if(!is.null(X) && is.null(TR)) {
        for(l in 1:num.X) {
          sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,X[,l],"*")
          grad.env <- c(grad.env,colSums(sum1))
        }
      }
      if(!is.null(TR)) {
        sum1 <- c(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10)) * Xd
        grad.B <- c(grad.B,colSums(sum1))
      }
      if(row.eff) { grad.row.params <- rowSums((dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)))),na.rm=TRUE) }
    }

    if(family=="ordinal") {
      grad.zeta <- matrix(0,p,ncol(new.zeta))
      deriv.trunnorm <- matrix(0,n,p)
      for(j in 1:p) {
        deriv.trunnorm[y[,j] == 1,j] <- -dnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j]) / pnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j])
        deriv.trunnorm[y[,j] == max(y[,j]),j] <- dnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]) / (1 - pnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]))
        j.levels <- (1:max(y[,j]))[-c(1,max(y[,j]))]
        for(k in j.levels) {
          deriv.trunnorm[y[,j] == k,j] <- (-dnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) + dnorm(new.zeta[j,k - 1] - eta.mat[y[,j] == k,j])) / (pnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) - pnorm(new.zeta[j,k-1] - eta.mat[y[,j] == k,j])) }
      }
      deriv.trunnorm[!is.finite(deriv.trunnorm)] <- 0
      grad.beta0 <- colSums(deriv.trunnorm,na.rm=TRUE)
      if(num.lv > 0) {
        for(l in 1:num.lv) {
          if(Lambda.struc == "unstructured") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
          if(Lambda.struc == "diagonal") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - (new.lambda[,l] %*% t(new.theta[,l]))
          grad.theta <- c(grad.theta,colSums(sum1,na.rm=TRUE))
        }
      }
      if(!is.null(X) && is.null(TR)) {
        for(l in 1:num.X) { sum1 <- sweep(deriv.trunnorm,1,X[,l],"*");
        grad.env <- c(grad.env,colSums(sum1,na.rm=TRUE))
        }
      }
      if(!is.null(TR)) {
        sum1 <- c(deriv.trunnorm) * Xd
        grad.B <- c(grad.B,colSums(sum1))
      }

      if(row.eff) { grad.row.params <- rowSums(deriv.trunnorm,na.rm=TRUE) }

      for(j in 1:p) {
        zeta0 <- new.zeta[j,]
        grad.zeta[j,max(y[,j]) - 1] <- grad.zeta[j,max(y[,j]) - 1] - sum(dnorm(zeta0[max(y[,j]) - 1] - eta.mat[which(y[,j] == max(y[,j])) ,j]) / (1 - pnorm(zeta0[max(y[,j]) - 1] - eta.mat[which(y[,j] == max(y[,j])),j])))
        j.levels <- (1:max(y[,j]))[-c(1,max(y[,j]))]
        for(k in j.levels) {
          grad.zeta[j,k] <- grad.zeta[j,k] + sum(dnorm(zeta0[k] - eta.mat[y[,j] == k,j]) / (pnorm(zeta0[k] - eta.mat[y[,j] == k,j]) - pnorm(zeta0[k - 1] - eta.mat[y[,j] == k,j])))
          grad.zeta[j,k - 1] <- grad.zeta[j,k - 1] - sum(dnorm(zeta0[k - 1] - eta.mat[y[,j] == k,j]) / (pnorm(zeta0[k] - eta.mat[y[,j] == k,j]) - pnorm(zeta0[k - 1] - eta.mat[y[,j] == k,j])))
        }
      }
      grad.zeta <- grad.zeta[,-1]
      if(length(unique(apply(y,2,max))) > 1) { grad.zeta <- c(grad.zeta[-which(grad.zeta == 0)]) }
    }

    if(num.lv>0){
      grad.theta=matrix(grad.theta,nrow=p,ncol=num.lv)
      grad.theta[upper.tri(grad.theta)] <- 0}
    if(row.eff) grad.row.params[1]=0;
    return(c(grad.theta, grad.beta0, grad.env, grad.B, grad.row.params, grad.phi, grad.zeta))
  }

  zero.cons <- which(theta==0)
  if(family == "ordinal") {
    zeta2 <- as.vector(t(zeta[,-1]))
    find.na.zeta <- which(is.na(zeta2)); if(length(find.na.zeta) > 0) zeta2 <- zeta2[-find.na.zeta]
  }
  if(family != "ordinal") { zeta2 <- NULL }

  A.mat <- -nd2(x0=c(theta,beta0,env,B,row.params,phi,zeta2), f = grad.all.mod, vameans = vameans, lambda = lambda)
  r1.off=NULL; if(row.eff) {r1.off=length(c(theta,beta0,env,B))+1}
  if(length(zero.cons) > 0) { A.mat <- A.mat[-c(zero.cons,r1.off),]; A.mat <- A.mat[,-c(zero.cons,r1.off)]
  } else if(row.eff) {A.mat <- A.mat[-r1.off,]; A.mat <- A.mat[,-r1.off] } #!!!!!

  w1=FALSE
  if(num.lv > 0) {
    D.mat.fun <- function(i3) {
      if(Lambda.struc == "unstructured") { lambda.i <- lambda[i3,,][lower.tri(lambda[i3,,],diag=TRUE)] }
      if(Lambda.struc == "diagonal") { lambda.i <- lambda[i3,]; }

      return(-nd2(f = grad.all.var, x0=c(vameans[i3,],lambda.i), vameans = vameans, lambda = lambda, mod.x = c(theta,beta0,env,B,row.params,phi,zeta2), sel.i = i3))
    }

    D.mat <- vector("list",n);
    if(p>1){
      unique.ind <- which(!duplicated(y))
      for(k in 1:length(unique.ind)) {
        subD <- D.mat.fun(i3=unique.ind[k]); match.seq <- which(apply(y,1,identical,y[unique.ind[k],]) == 1) ## Unfortunately replace() only works if you want to replace elements with numerics
        if( any(is.nan(subD)) ){ logi<-is.nan(subD); subD[logi]=diag(dim(subD)[1])[logi]; w1=TRUE}
        for(k2 in match.seq) { D.mat[[k2]] <- subD }
      }
    } else {
      unique.ind=1:n
      for(k in 1:length(unique.ind)) {
        subD <- D.mat.fun(i3=unique.ind[k]);
        if( any(is.nan(subD)) ){ logi<-is.nan(subD); subD[logi]=diag(dim(subD)[1])[logi]; w1=TRUE}
        for(k2 in 1:n) { D.mat[[k2]] <- subD }
      }
    }
    D.mat <- as.matrix(Matrix::bdiag(D.mat)); rm(subD)
    B.mat <- vector("list",n)
    for(i3 in 1:length(unique.ind)) {
      #cat("Onto row",i3,"\n")
      if(Lambda.struc == "unstructured") { lambda.i <- lambda[unique.ind[i3],,][lower.tri(lambda[unique.ind[i3],,],diag=TRUE)] }
      if(Lambda.struc == "diagonal") { lambda.i <- lambda[unique.ind[i3],]; }

      subB.mat <- -nd2(f = grad.all.cross, x0=c(vameans[unique.ind[i3],],lambda.i), vameans = vameans, lambda = lambda, mod.x = c(theta,beta0,env,B,row.params,phi,zeta2), sel.i = unique.ind[i3])
      if(length(zero.cons) > 0) { subB.mat <- subB.mat[-c(zero.cons,r1.off),]
      }else if(row.eff) {subB.mat <- subB.mat[-c(r1.off),] }
      if( any(is.nan(subB.mat)) ){ logi<-is.nan(subB.mat); subB.mat[logi]=0; w1=TRUE}

      if(p>1){
        match.seq <- which(apply(y,1,identical,y[unique.ind[i3],]) == 1)
        for(k2 in match.seq) { B.mat[[k2]] <- subB.mat }
      } else {  for(k2 in 1:n) { B.mat[[k2]] <- subB.mat }}

    }
    B.mat <- do.call(cbind,B.mat)

    cov.mat.mod <- MASS::ginv(A.mat-B.mat%*%solve(D.mat)%*%t(B.mat))
  }

  if(num.lv == 0) cov.mat.mod <- MASS::ginv(A.mat)
  if(w1) warning("Standard errors of parameters might not be reliable, due to the technical remedy.\n");
  sd.errors <- sqrt(diag(abs(cov.mat.mod)))
  theta.errors <- NULL;

  sdout=NULL
  if(num.lv > 0) {
    theta.vec <- sd.errors[1:(p*num.lv-length(zero.cons))]; theta.errors <- matrix(0,p,num.lv);
    if(length(zero.cons) > 0) theta.errors[-zero.cons] <- theta.vec
    if(length(zero.cons) == 0) theta.errors <- matrix(theta.vec,p,num.lv)
    sd.errors <- sd.errors[-(1:(p*num.lv-length(zero.cons)))]
    if(num.lv > 1) theta.errors[upper.tri(theta.errors)] <- 0
    rownames(theta.errors) <- colnames(y);
    colnames(theta.errors) <- paste("LV", 1:num.lv, sep="");
    sdout$theta <- theta.errors
  }

  beta.errors <- sd.errors[1:p]; sd.errors <- sd.errors[-(1:p)]
  names(beta.errors) <- colnames(y)
  sdout$beta0 <- beta.errors
  env.errors=NULL
  if(!is.null(X)) {
    if(is.null(TR)){
      beta.errors <- matrix(sd.errors[1:(p*ncol(X))],p,ncol(X));
      rownames(beta.errors) <- colnames(y); colnames(beta.errors) <- colnames(X);
      sd.errors <- sd.errors[-(1:(p*ncol(X)))]
      sdout$Xcoef <- beta.errors
    } else {
      coef.errors <- sd.errors[1:nd]
      names(coef.errors) <- colnames(Xd)
      sd.errors <- sd.errors[-(1:(nd))]
      sdout$B <- coef.errors
    }
  }

  row.params.errors <- NULL;
  if(row.eff) {
    row.params.errors <- c(0,sd.errors[1:(n-1)]); sd.errors <- sd.errors[-(1:(n-1))];
    names(row.params.errors) <- rownames(y)
    sdout$row.params <- row.params.errors
  }



  phi.errors <- NULL;
  if(family == "negative.binomial") {
    phi.errors <- sd.errors[1:p]; sd.errors <- sd.errors[-(1:p)];
    names(phi.errors) <- colnames(y)
    sdout$inv.phi <- phi.errors
  }

  zeta.errors <- NULL;
  if(family == "ordinal") {
    zeta.errors <- matrix(NA,p,max.levels - 2)
    for(j in 1:p) { if(max(y[,j]) > 2) { zeta.errors[j,1:(max(y[,j]) - 2)] <- sd.errors[1:(max(y[,j]) - 2)]; sd.errors <- sd.errors[-(1:(max(y[,j]) - 2))] } }
    zeta.errors <- cbind(0,zeta.errors)
    rownames(zeta.errors) <- colnames(y); colnames(zeta.errors) <- paste(1:(max(y)-1),"|",2:max(y),sep="")
    sdout$zeta <- zeta.errors
  }

  return(sdout)
}


## A function to compute highly accurate first-order derivatives
## From Fornberg and Sloan (Acta Numerica, 1994, p. 203-267; Table 1, page 213)
## Adapted by Scott Foster from code nicked off the net 2007
nd2 <- function(x0, f, m = NULL, D.accur = 2, ...) {
  D.n <- length(x0)
  if(is.null(m)) { D.f0 <- f(x0, ...); m <- length(D.f0) }
  if(D.accur == 2) { D.w <- tcrossprod(rep(1,m), c(-1/2,1/2)); D.co <- c(-1, 1) }
  else {
    D.w <- tcrossprod(rep(1,m), c(1/12, -2/3, 2/3, -1/12))
    D.co <- c(-2, -1, 1, 2)
  }
  D.n.c <- length(D.co)
  macheps <- .Machine$double.eps
  D.h <- macheps^(1/3) * abs(x0)
  D.deriv <- matrix(NA, nrow = m, ncol = D.n)

  for (ii in 1:D.n) {
    D.temp.f <- matrix(0, m, D.n.c)
    for (jj in 1:D.n.c) {
      D.xd <- x0 + D.h[ii] * D.co[jj] * (1:D.n == ii)
      D.temp.f[,jj] <- f( D.xd, ...)
    }
    D.deriv[,ii] <- rowSums(D.w * D.temp.f) / D.h[ii]
  }
  return(D.deriv)
}







calc.quad <- function(lambda,theta,Lambda.struc) {
  if(Lambda.struc == "diagonal") out <- 0.5 * (lambda) %*% t(theta^2)
  if(Lambda.struc == "unstructured") {
    if(class(lambda) == "array") { n <- dim(lambda)[1]; num.lv <- dim(lambda)[2] }
    if(class(lambda) == "matrix") { num.lv <- dim(lambda)[2]; n <- 1 }
    if(class(theta) == "matrix") { p <- dim(theta)[1] }
    if(class(theta) == "numeric") { p <- 1; theta <- matrix(theta,1) }

    out <- matrix(NA,n,p)

    if(n == 1) {
      for(j in 1:p) { out[1,j] <- 0.5 * t(theta[j,]) %*% lambda %*% theta[j,] }
    }

    if(n > 1) {
      if(n <= p) out <- t(sapply(1:n, function(x) 0.5 * rowSums(theta * (theta %*% lambda[x,,]))))
      if(n > p) {
        lambda.mat <- aperm(lambda,c(3,2,1)); dim(lambda.mat) <- c(num.lv,num.lv * n)
        f <- function(x) 0.5 * rowSums(matrix((t(lambda.mat) * theta[x,]) %*% theta[x,],ncol=num.lv,byrow=TRUE))
        out <- sapply(1:p,f)
      }
    }
  }
  return(list(mat = out, mat.sum = sum(out)))
}

lambda.convert <- function(lambda,type=1) {
  if(type == 1) { ## Convert unstructured to diagonal
    n <- dim(lambda)[1]
    lambda2 <- matrix(1,dim(lambda)[1],dim(lambda)[3])
    for(i in 1:n) { lambda2[i,] <- diag(lambda[i,,]) }
    return(lambda2)
  }
  if(type == 2) { ## Convert diagonal to unstructured
    n <- nrow(lambda)
    lambda2 <- array(0,dim=c(nrow(lambda),ncol(lambda),ncol(lambda)))
    for(i in 1:n) { diag(lambda2[i,,]) <- lambda[i,] }
    return(lambda2)
  }
}

inf.criteria <- function(fit)
{
  family=fit$family
  abund=fit$y
  num.lv=fit$num.lv
  n <- dim(abund)[1]
  p <- dim(abund)[2]
  k<-attributes(logLik.gllvm(fit))$df

  BIC <- -2*fit$logL + (k) * log(n)
  # AIC
  AIC <- -2*fit$logL + (k) * 2
  # AICc
  AICc <- AIC + 2*k*(k+1)/(n*p-k-1)
  list(BIC = BIC, AIC = AIC, AICc = AICc, k = k)
}

# Creates matrix of fourth corner terms from a vector
getFourthCorner<- function(object){
  if(is.null(object$X) || is.null(object$TR)) stop();

  n1=colnames(object$X)
  n2=colnames(object$TR)

  nams=names(object$params$B)
  fx<-cbind(apply(sapply(n1,function(x){grepl(x, nams)}),1,any), apply(sapply(n2,function(x){grepl(x, nams)}),1,any))
  fourth.index<-rowSums(fx)>1
  nams2=nams[fourth.index]
  fourth.corner=object$params$B[fourth.index]

  splitsnam <- sapply(nams2,strsplit, split = ":")
  isinvec <- function(vec, ob =""){
    ob %in% vec
  }
  n1 <- n1[(n1 %in% unlist(splitsnam))]
  n2 <- n2[(n2 %in% unlist(splitsnam))]
  i=1; j=1;
  fourth<-matrix(0,length(n1),length(n2))
  for (i in 1:length(n1)) {
    for (j in 1:length(n2)) {
      fur=(unlist(lapply(splitsnam, isinvec, ob = n1[i])) + unlist(lapply(splitsnam, isinvec, ob = n2[j]))) >1
      if(any(fur)){ fourth[i,j]=fourth.corner[fur]}
    }
  }
  
  colnames(fourth)=n2
  rownames(fourth)=n1
  return(fourth)
}


# Calculates standard errors for random effects
sdrandom<-function(obj, Vtheta, incl, ignore.u = FALSE){
  r <- obj$env$random
  par = obj$env$last.par.best
  hessian.random <- obj$env$spHess(par, random = TRUE)
  L <- obj$env$L.created.by.newton
  if (ignore.u) {
    diag.term2 <- 0
  } else {
    f <- obj$env$f
    w <- rep(0, length(par))
    reverse.sweep <- function(i) {
      w[i] <- 1
      f(par, order = 1, type = "ADGrad", rangeweight = w, doforward = 0)[r]
    }
    nonr <- setdiff(seq_along(par), r)
    tmp <- sapply(nonr, reverse.sweep)
    if (!is.matrix(tmp))
      tmp <- matrix(tmp, ncol = length(nonr))
    A <- solve(hessian.random, tmp[, incl])
    diag.term2 <- diag(rowSums((A %*% Vtheta) * A))
  }
  diag.term1 <- Matrix::chol2inv(L)
  diag.cov.random <- diag.term1 + diag.term2
  return(diag.cov.random)
}


start.values.randomX <- function(y, Xb, family, starting.val, power = NULL) {
  y <- as.matrix(y)
  Xb <- as.matrix(Xb)
  n <- NROW(y)
  p <- NCOL(y)
  tr0 <- try({
    
    if(starting.val %in% c("res", "random")){
      if(family %in% c("poisson", "negative.binomial", "binomial", "ZIP")){
        if(family == "ZIP") family <- "poisson"
        f1 <- gllvm.TMB(y=y, X=Xb, family = family, num.lv=0, starting.val = "zero")
        coefs0 <- as.matrix(scale((f1$params$Xcoef), scale = FALSE))
        Br <- coefs0/max(apply(coefs0, 2, sd))
        sigmaB <- cov(Br)
        Br <- t(Br)
      } else if(family == "tweedie"){
        coefs0 <- NULL
        for(j in 1:p){
          fitj <- try(glm( y[,j] ~ Xb, family = statmod::tweedie(var.power=power, link.power=0) ),silent = TRUE)
          if(!inherits(fitj, "try-error")){
            coefs0 <- rbind(coefs0, fitj$coefficients[-1])
          } else { coefs0 <- rbind(coefs0,rnorm(dim(Xb)[2])); }
        }
        Br <- coefs0/max(apply(coefs0, 2, sd))
        sigmaB <- cov(Br)
        Br <- t(Br)
      } else {
        Br <- matrix(0, ncol(Xb), p)
        sigmaB <- diag(ncol(Xb))
      }
    } else {
      Br <- matrix(0, ncol(Xb), p)
      sigmaB <- diag(ncol(Xb))
    }
  }, silent = TRUE)
  
  if(inherits(tr0, "try-error")){
    Br <- matrix(0, ncol(Xb), p)
    sigmaB <- diag(ncol(Xb))
  }
  
  
  return(list(Br = Br, sigmaB = sigmaB))
}



# Calculates adjusted prediction errors for random effects
sdA<-function(fit){
  n<-nrow(fit$y)
  A<- fit$Hess$cov.mat.mod  #
  B<- fit$Hess$Hess.full[fit$Hess$incl, fit$Hess$incla]
  C<- fit$Hess$Hess.full[fit$Hess$incla, fit$Hess$incl]
  D<- solve(fit$Hess$Hess.full[fit$Hess$incla, fit$Hess$incla])
  covb <- (D%*%C)%*%(A)%*%(B%*%t(D))
  se <- (diag(abs(covb)))
  CovAerr<-array(0, dim(fit$A))
  for (i in 1:ncol(fit$A)) {
    CovAerr[,i,i] <- se[1:n]; se<-se[-(1:n)]
  }
  CovAerr
}

sdB<-function(fit){
  n<-nrow(fit$y)
  p<-ncol(fit$y)
  incla<-rep(FALSE, length(fit$Hess$incl))
  incla[names(fit$TMBfn$par)%in%c("u", "Br")] <- TRUE
  fit$Hess$incla <- incla
  
  A<- fit$Hess$cov.mat.mod  #
  B<- fit$Hess$Hess.full[fit$Hess$incl, fit$Hess$incla]
  C<- fit$Hess$Hess.full[fit$Hess$incla, fit$Hess$incl]
  D<- solve(fit$Hess$Hess.full[fit$Hess$incla, fit$Hess$incla])
  covb <- (D%*%C)%*%(A)%*%(B%*%t(D))
  se <- (diag(abs(covb)))
  CovAerr<-array(0, dim(fit$A))
  for (i in 1:ncol(fit$A)) {
    CovAerr[,i,i] <- se[1:n]; se<-se[-(1:n)]
  }
  CovABerr<-array(0, dim(fit$Ab))
  for (i in 1:ncol(fit$Ab)) {
    CovABerr[,i,i] <- se[1:p]; se<-se[-(1:p)]
  }
  CovABerr
}


# draw an ellipse
ellipse<-function(center, covM, rad, col=4, lwd = 1, lty = 1){
  seg <- 51
  Qc <- chol(covM, pivot = TRUE)
  angles <- (0:seg) * 2 * pi / seg
  unit.circ <- cbind(cos(angles), sin(angles))
  order <- order(attr(Qc, "pivot"))
  ellips <- t(center + rad * t(unit.circ %*% Qc[, order]))
  lines(ellips, col = col, lwd = lwd, lty = lty)
}

gamEnvelope <- function(x, y,line.col = "red", envelope.col = c("blue","lightblue"), col = 1, envelopes = TRUE, ...){
  xSort <- sort(x, index.return = TRUE)
  gam.yx <- gam(y[xSort$ix] ~ xSort$x)
  pr.y <- predict.gam(gam.yx, se.fit = TRUE)
  n.obs <- length(xSort$ix)
  prHi <- pr.y$fit + 1.96*pr.y$se.fit
  prLow <- pr.y$fit - 1.96*pr.y$se.fit
  if(envelopes) polygon(xSort$x[c(1:n.obs,n.obs:1)], c(prHi,prLow[n.obs:1]), col = envelope.col[2], border = NA)
  lines(xSort$x, pr.y$fit, col = envelope.col[1])
  abline(h = 0, col = 1)
  points(x, y, col = col, ...)
}

mlm <- function(y, X = NULL, index = NULL){
  out <- list(coefficients=NULL, residuals = NULL)
  coefficients <- residuals <- NULL
  if(!is.null(X) || !is.null(index)){ 
    Xx<- cbind(X, index)
    for (j in 1:ncol(y)) {
      f0 <- lm(y[,j]~Xx)
      coefficients <- rbind(coefficients, f0$coef)
      residuals <- cbind(residuals, f0$residuals)
    }
  } else {
    for (j in 1:ncol(y)) {
      f0 <- lm(y[,j]~1)
      coefficients <- rbind(coefficients, f0$coef)
      residuals <- cbind(residuals, f0$residuals)
    }
  }
  out$coefficients <- t(coefficients)
  out$residuals <- residuals
  out
}
                                 
#'@export
#'@export nobs.gllvm
nobs.gllvm <- function(object){
  n <- prod(dim(object$y))
  return(n)
}

#'@export nobs
nobs <- function(object, ...)
{
  UseMethod(generic = "nobs")
}


#'@export
#'@export AICc.gllvm
AICc.gllvm <- function(object, ...){
  objectlist <- list(object, ...)
  IC<-lapply(objectlist,function(x){
    abund=x$y
    n <- dim(abund)[1]
    p <- dim(abund)[2]
    k<-attributes(logLik.gllvm(x))$df
    AICc <- -2*x$logL + (k) * 2 + 2*k*(k+1)/(n*p-k-1)
    return(AICc)
  })
  return(unlist(IC))
}

#'@export AICc
AICc <- function(object, ...)
{
  UseMethod(generic = "AICc")
}