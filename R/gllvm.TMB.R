########################################################################################
## GLLVM, with estimation done via Variational approximation using TMB-package
## Original author: Jenni Niku
########################################################################################
gllvm.TMB <- function(y, X = NULL, lv.X = NULL, formula = NULL, family = "poisson", 
      num.lv = 2, num.lv.c = 0, num.RR = 0, num.lv.cor=0, lv.formula = NULL, corWithin = TRUE, randomB = FALSE, 
      method="VA",Lambda.struc="unstructured", Ar.struc="diagonal", row.eff = FALSE, reltol = 1e-8, reltol.c = 1e-8,
      seed = NULL,maxit = 3000, max.iter=200, start.lvs = NULL, offset=NULL, sd.errors = FALSE,
      trace=FALSE,link="logit",n.init=1,n.init.max = 10, restrict=30,start.params=NULL, dr=NULL, rstruc =0, cstruc = "diag", dist =matrix(0),
      optimizer="optim",starting.val="res",Power=1.5,diag.iter=1, dependent.row = FALSE, scalmax = 10, MaternKappa = 1.5,
      Lambda.start=c(0.1,0.5), quad.start=0.01, jitter.var=0, zeta.struc = "species", quadratic = FALSE, start.struc = "LV", optim.method = "BFGS", disp.group = NULL, NN=matrix(0), setMap=NULL) { 
  # , Dthreshold=0
  # If there is no random effects/LVs set diag iter to zero:
  if(((num.lv+num.lv.c)==0) & (row.eff!="random") & (randomB==FALSE)) diag.iter <-  0
  

  if(!is.null(start.params)) starting.val <- "zero"
  ignore.u <- FALSE
  n <- nr <- nu <- dim(y)[1]
  p <- dim(y)[2]
  times <- 1
  objrFinal <- optrFinal <- NULL
  if(is.null(disp.group)) disp.group <- 1:NCOL(y)

  cstrucn = switch(cstruc, "diag" = 0, "corAR1" = 1, "corExp" = 2, "corCS" = 3, "corMatern" = 4)
  
  # Structure for row effects
  model = 0
  xr = NULL
  if(rstruc==0 & num.lv.cor==0){ # No structure
    dr <- diag(n)
  }
  
  Astruc = 0;
  scaledc = 0;
  rho.lv =NULL
  if(rstruc>0){#rstruc
    dist<-as.matrix(dist)
    if(is.null(dr)) stop("Define structure for row params if 'rstruc == ",rstruc,"'.")
    if(rstruc==1){# group specific
      nr <- dim(dr)[2]
      if((cstrucn == 2) | (cstrucn == 4)) {
        if(is.null(dist) || NROW(dist)!=nr)
          dist=matrix(1:nr)
        if(NROW(dist)!=nr)
          stop("Number of rows in 'dist' should be same as maximum number of groups when corWithin = FALSE")
      }
    }
    if(rstruc==2) { # correlated within groups
      if(is.null(dr)) stop("Define structure for row params if 'rstruc == 2'.")
      nr <- dim(dr)[2]
      times <- n/nr#dim(dr)[1]
      if((cstrucn == 2) | (cstrucn == 4)) {
        if(is.null(dist) || NROW(dist)!=times)
          dist=matrix(1:times)
        if(NROW(dist)!=times)
          stop("Number of rows in 'dist' should be same as maximum number of units within groups when corWithin = TRUE")
        
      }
    }
    if((cstrucn == 2) | (cstrucn == 4)) {
      AD1 = (apply(as.matrix(dist),2,max)-apply(as.matrix(dist),2,min))/scalmax
      scaledc = log(AD1)
      # AD1 = pmax(apply(as.matrix(dist),2,function(x) min(dist(unique(x), diag = FALSE))),1)
      # md = min(dist(as.matrix(dist)%*%diag(1/(AD1), length(AD1)), diag = FALSE))/2
      # if(md>5) AD1 = AD1*md
      # scaledc = log(AD1)
      # if(!is.null(setMap$scaledc)) {
      #   if( (length(setMap$scaledc)!= NCOL(dist))) stop("setMap$scaledc must be a numeric vector and have length that is same as the number of columns in 'dist'.")
      #   scaledc[is.na(setMap$scaledc)]=0
      # }
      # if(NCOL(dist)==1) {
      #   scaledc = scaledc*0
      #   setMap$scaledc = factor(rep(NA, length(scaledc)))
      # }
    }
    if(nr==1) Ar.struc = "diagonal"
  }
  if(num.lv.cor > 0){#rstruc
    dist<-as.matrix(dist)
    if(is.null(dr)) stop("Define structure for row params if 'rstruc == ",rstruc,"'.")
    # LVs correlated within groups
    if(is.null(dr)) stop("Define structure for row params if 'rstruc == 2'.")
    nu <- dim(dr)[2]
    times <- n/nu#dim(dr)[1]
    if((cstrucn == 2) | (cstrucn == 4)) {
      if(corWithin){
          if(is.null(dist))
            dist=matrix(1:times)
          if(NROW(dist)!=times)
            stop("Number of rows in 'dist' should be same as maximum number of units within groups when corWithin = TRUE")
        } else {
          if(is.null(dist))
            dist=matrix(1:nu)
          if(NROW(dist)!=nu)
            stop("Number of rows in 'dist' should be same as maximum number of groups when corWithin = FALSE")
        }
      AD1 = (apply(as.matrix(dist),2,max)-apply(as.matrix(dist),2,min))/scalmax
      scaledc<-log(AD1)
      # AD1<-pmax(apply(as.matrix(dist),2,function(x) min(dist(unique(x), diag = FALSE))),1)
      # md<-min(dist(as.matrix(dist)%*%diag(1/(AD1), length(AD1)), diag = FALSE))/2
      # if(md>5) AD1=AD1*md
      # scaledc<-log(AD1)
      # if(!is.null(setMap$scaledc)) {
      #   if( (length(setMap$scaledc)!= NCOL(dist))) stop("setMap$scaledc must be a numeric vector and have length that is same as the number of columns in 'dist'.")
      #   scaledc[is.na(setMap$scaledc)]=0
      # }
      # if(NCOL(dist)==1) {
      #   scaledc = scaledc*0
      #   setMap$scaledc = factor(rep(NA, length(scaledc)))
      # }
      # scaledc<-rep(log(1),ncol(dist))
    }
    rho_lvc<- matrix(rep(0, num.lv.cor))
    if(Lambda.struc == "unstructured") {Astruc=1}
    if(Lambda.struc == "bdNN") {Astruc=2}
    if(Lambda.struc %in% c("diagU","UNN","UU")) {
      if(num.lv.cor>1){
        if(Lambda.struc == "UU") {Astruc=3; }#Lambda.struc = "unstructured"}
        if(Lambda.struc == "UNN" && num.lv.cor>0) {Astruc=4; Lambda.struc = "bdNN"}
        if(Lambda.struc == "diagU" && num.lv.cor>0) {Astruc=5; Lambda.struc = "diagonal"}
      } else {
        if(Lambda.struc == "UU") {Astruc=1; }#Lambda.struc = "unstructured"}
        if(Lambda.struc == "UNN" && num.lv.cor>0) {Astruc=2; Lambda.struc = "bdNN"}
        if(Lambda.struc == "diagU" && num.lv.cor>0) {Astruc=0; Lambda.struc = "diagonal"}
      }
      }
  }

  tr <- NULL
  y <- as.matrix(y)
  if(is.null(lv.X)){
    lv.X<-matrix(0)
  }
  formula1 <- formula
  if(method=="VA" && (family =="binomial")){ link="probit"}
  jitter.var.r <- 0
  if(length(jitter.var)>1){ 
    jitter.var.r <- jitter.var[2]
    jitter.var <- jitter.var[1]
  }
  
  if (!is.numeric(y))
    stop( "y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
  if ((family %in% c("ZIP")) && (method %in% c("VA", "EVA"))) #"tweedie", 
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
    if(any(!apply(y,2,function(x)all(diff(sort(unique(x)))==1)))&zeta.struc=="species"){
      warning("Can't fit ordinal model if there are species with missing classes. Setting 'zeta.struc = `common`'")
      zeta.struc = "common"
    }
      
    if(!all(min(y)==apply(y,2,min))&zeta.struc=="species"){
      stop("For ordinal data and zeta.struc=`species` all species must have the same minimum category.Setting 'zeta.struc = `common`'.")
      zeta.struc = "common"
    }
    if(any(diff(sort(unique(c(y))))!=1)&zeta.struc=="common")
      stop("Can't fit ordinal model if there are missing response classes. Please reclassify.")
  }

  # Define design matrix for covariates
  num.X <- 0;
  Xorig <- X
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
    }
  }
  
  if (is.null(formula) && is.null(X)) {
    formula = "~ 1"
  }

  ## Set initial values for model parameters (including dispersion prm) and latent variables

  ### Seeds
  
  # If number of seeds is less than n.init, sample the seeds randomly, but using the given seed
  if((length(seed) >1) & (length(seed) < n.init)) {
    stop("Seed length doesn't match with the number of initial starts.")
  }
  if(!is.null(seed) & (length(seed) ==1) & (length(seed) < n.init)) {
    set.seed(seed)
    seed <- sample(1:10000, n.init)
  }
  # If no seed is sampled it is randomly drawn
  if(is.null(seed)&starting.val!="zero"){
    seed <- sample(1:10000, n.init)
  }


  n.i <- 1
  if(starting.val!="zero"){seed.best <- seed[n.i]}else{seed.best <- NULL}
  out <- list( y = y, X = X, logL = Inf, num.lv = num.lv, num.lv.c = num.lv.c, row.eff = row.eff, family = family, X.design = X, method = method, zeta.struc = zeta.struc)
  #if (n.init > 1)

  # n.init model fits
  while(n.i <= n.init){

    if(n.init > 1 && trace)
      cat("Initial run ", n.i, "\n")

    #### Calculate starting values
    fit <- start.values.gllvm.TMB(y = y, X = Xorig, formula = formula, lv.X = lv.X, TR = NULL, family = family, offset= offset, num.lv = num.lv, num.lv.c = num.lv.c, num.RR = num.RR, start.lvs = start.lvs, seed = seed[n.i], starting.val = starting.val, Power = Power, jitter.var = jitter.var, row.eff = row.eff, TMB=TRUE, link=link, zeta.struc = zeta.struc, disp.group = disp.group, method=method, randomB = randomB)

    ## Set initial values
        sigma <- 1
    if (is.null(start.params)) {
      beta0 <- fit$params[, 1]
      if((num.lv.c+num.RR)>0){b.lv <- fit$b.lv}else{b.lv<-matrix(0)}
      betas <- NULL
      if (!is.null(X))
        betas <- c(fit$params[, 2:(num.X + 1)])
      lambdas <- NULL

      if ((num.lv+(num.lv.c+num.RR)) > 0) {
        sigma.lv <- fit$sigma.lv
        lambdas <- as.matrix(fit$params[, (ncol(fit$params) - num.lv - (num.lv.c+num.RR) + 1):ncol(fit$params)])
        if(start.struc=="LV"&quadratic!=FALSE|quadratic=="LV"){
          lambda2 <- matrix(quad.start, ncol = num.lv + (num.lv.c+num.RR), nrow = 1)  
        }else if(start.struc=="all"&quadratic!=FALSE){
          lambda2 <- matrix(quad.start, ncol = num.lv + (num.lv.c+num.RR), nrow = p)
        }else if(quadratic==FALSE){
          lambda2 <- 0
        }
        if(randomB!=FALSE){
          sigmab_lv <- fit$sigmab_lv
        }
        
          if(num.lv>1&(num.lv.c+num.RR)==0){
            lambdas[upper.tri(lambdas)] <- 0  
          }else if(num.lv==0&(num.lv.c+num.RR)>1){
            lambdas[upper.tri(lambdas)] <- 0
          }else if(num.lv>0&num.lv.c>0){
            if((num.lv.c+num.RR)>1)lambdas[,1:(num.lv.c+num.RR)][upper.tri(lambdas[,1:(num.lv.c+num.RR)])] <- 0
            if(num.lv>1)lambdas[,((num.lv.c+num.RR)+1):ncol(lambdas)][upper.tri(lambdas[,((num.lv.c+num.RR)+1):ncol(lambdas)])] <- 0
          }
          
        if(quadratic != FALSE){
          fit$params <- cbind(fit$params, matrix(lambda2,nrow=p,ncol=num.lv+(num.lv.c+num.RR)))  
        }else{
          fit$params <- fit$params
        }
        
        if(num.lv.cor>0){ # In correlation model, sigmas are  scale parameters
          # lambdas <- lambdas%*%diag(sigma.lv, nrow = length(sigma.lv), ncol = length(sigma.lv))
          rho_lvc<- rep(0, num.lv.cor);
          if((cstrucn == 2) | (cstrucn == 4)) {
            AD1 = (apply(as.matrix(dist),2,max)-apply(as.matrix(dist),2,min))/scalmax
            scaledc<-log(AD1)
          }
        }
        if(family == "betaH"){ # Own loadings for beta distr in hurdle model
          thetaH <- t(lambdas%*%diag(sigma.lv, nrow = length(sigma.lv), ncol = length(sigma.lv)))
        }
      }

      row.params <- NULL
      
      if (row.eff != FALSE) {
        row.params <- fit$row.params
        if(rstruc==0 && row.eff=="random") row.params <- row.params[1:nr]#rstruc
        if(rstruc==1 && row.eff=="random") try(row.params <- (t(dr)%*%(row.params))/(dim(dr)[1]/dim(dr)[2]), silent = TRUE)#rstruc
        if (row.eff == "random") {
          sigma <- sd(row.params)#1;#
        }
      }#rep(0,n)
      lvs <- NULL
      if ((num.lv+num.lv.c) > 0)
        lvs <- matrix(fit$index, ncol = num.lv+num.lv.c)

    } else{
      if (all(dim(start.params$y) == dim(y)) &&
          is.null(X) == is.null(start.params$X) &&
          (row.eff == start.params$row.eff)) {
        if(class(start.params)[2]=="gllvm.quadratic" && quadratic != FALSE){
          lambda2 <- start.params$params$theta[,-c(1:(start.params$num.lv+start.params$num.lv.c+start.params$num.RR)),drop=F]
        }else if(class(start.params)[1]=="gllvm" && quadratic != FALSE){
          if(start.struc=="LV"|quadratic=="LV"){
            lambda2 <- matrix(quad.start, ncol = num.lv+num.lv.c+num.RR, nrow = 1)  
          }else if(start.struc=="all"&quadratic==TRUE){
            lambda2 <- matrix(quad.start, ncol = num.lv+num.lv.c+num.RR, nrow = p)
          }
        }else if(quadratic == FALSE){
          lambda2 <- 0
        }
        if(start.params$randomB!=FALSE && randomB !=FALSE){
          sigmab_lv <- start.params$sigmaLvXcoef
        }else if(randomB!=FALSE){
          sigmab_lv <- fit$sigmab_lv
        }
        if((start.params$num.lv.c+start.params$num.RR)==0){
          b.lv <- matrix(0)
        }else{
          b.lv <- start.params$params$LvXcoef
        }
        beta0 <- start.params$params$beta0 ## column intercepts
        betas <- NULL
        if (!is.null(X))
          if(!all((dim(X) == dim(start.params$X)))) stop( "Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that predictors X are the same in both models.")
          betas <- c(start.params$params$Xcoef) ## covariates coefficients
        lambdas <- NULL
        if ((num.lv+(num.lv.c+num.RR)) > 0){
          sigma.lv <- start.params$params$sigma.lv
          lambdas <- start.params$params$theta
          if((num.lv.c+num.RR)>1)lambdas[,1:(num.lv.c+num.RR)][upper.tri(lambdas[,1:(num.lv.c+num.RR)])] <- 0
          if(num.lv>1)lambdas[,((num.lv.c+num.RR)+1):ncol(lambdas)][upper.tri(lambdas[,((num.lv.c+num.RR)+1):ncol(lambdas)])] <- 0
          
        }

        row.params <- NULL
        if (start.params$row.eff != FALSE) {
          row.params <- start.params$params$row.params
          if(row.eff=="fixed")
            row.params[1] <- 0
          if(row.eff=="random")
            sigma <- start.params$params$sigma
        }## row parameters
        lvs <- NULL
        sigma.lv <- 0
        if ((num.lv+num.lv.c) > 0) {
          sigma.lv <- start.params$params$sigma.lv
          lvs <- matrix(start.params$lvs, ncol = num.lv+num.lv.c)
        }
        if(num.lv.cor>0){ # sigmas are scale parameters # just diagonal values, not
          if(is.numeric(start.params$params$rho.lv) & ((cstrucn == 2) | (cstrucn == 4))) {
            # if(cstrucn == 4) start.params$params$rho.lv <- start.params$params$rho.lv[,-ncol(start.params$params$rho.lv), drop=FALSE]
            scaledc = colMeans(as.matrix(start.params$params$rho.lv)); 
            if(length(scaledc) < ncol(dist) ) scaledc <- rep(scaledc, ncol(dist))[1:ncol(dist)]
          }
        }
        
      } else {
        stop( "Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that attributes y, X and row.eff match to each other.")
      }
    }

    phis <- NULL
    if (family == "negative.binomial") {
      phis <- fit$phi
      if (any(phis > 100))
        phis[phis > 100] <- 100
      if (any(phis < 0.01))
        phis[phis < 0.01] <- 0.01
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
      if(length(unique(disp.group))!=p){
      phis <- sapply(1:length(unique(disp.group)),function(x)mean(y[,which(disp.group==x)]==0))*0.98 + 0.01  
      phis <- phis[disp.group]
      }else{
        phis <- (colMeans(y == 0) * 0.98) + 0.01  
      }
      phis <- phis / (1 - phis)
    } # ZIP probability
    if (family %in% c("gaussian", "gamma", "beta", "betaH")) {
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
    if(!is.null(row.params)){ r0 <- row.params} else {r0 <- rep(0,n)}
    a <- c(beta0)
    b <- NULL; if(!is.null(X)) b <- matrix(betas, ncol(X), p,byrow = TRUE)
    lambda=0
    if((num.lv+(num.lv.c))==0)u <- matrix(0)
    if((num.lv+num.RR+num.lv.c)==0)lambda2 <- matrix(0)
    if(num.lv > 0 & (num.lv.c+num.RR) == 0) {
      # diag(lambdas) <- log(diag(lambdas)) #!!!
      lambda <- lambdas[lower.tri(lambdas,diag = F)]
      u <- lvs
    }else if(num.lv == 0 & (num.lv.c+num.RR) > 0){
      lambda <- lambdas[lower.tri(lambdas,diag = F)]
      if(num.lv.c>0)u <- lvs
    }else if(num.lv>0&(num.lv.c+num.RR)>0){
      lambda <- lambdas[,1:(num.lv.c+num.RR)][lower.tri(lambdas[,1:(num.lv.c+num.RR),drop=F],diag = F)]
      lambda <- c(lambda,lambdas[,((num.lv.c+num.RR)+1):ncol(lambdas)][lower.tri(lambdas[,((num.lv.c+num.RR)+1):ncol(lambdas),drop=F],diag = F)])
      u <- lvs
    }
    if((num.lv+num.lv.c)==0){
      sigma.lv <- 0
    }
    if(!is.null(phis)) {
      phi <- phis 
    } else { 
        phi <- rep(1, p); 
        fit$phi <- phi
    }
    
    
    q <- num.lv+(num.lv.c+num.RR)

## map.list defines parameters which are not estimated in this model

    map.list <- list()    
    if(is.list(setMap)) {
      map.list <- setMap
    }
    map.list$B <- map.list$Br <- map.list$sigmaB <- map.list$sigmaij <- map.list$Abb <- factor(NA)
    
    xb<-Br<-matrix(0); sigmaB=diag(1);sigmaij=0; lg_Ar=0; Abb=0; Ab_lv = 0;
    if(row.eff==FALSE) map.list$r0 <- factor(rep(NA,n))
    if(randomB==FALSE){
      map.list$sigmab_lv <- factor(NA)
    }
    
    if(family %in% c("poisson","binomial","ordinal","exponential")){
      map.list$lg_phi <- factor(rep(NA,p))
    } else if(family %in% c("tweedie", "negative.binomial", "gamma", "gaussian", "beta", "betaH", "ZIP")){
      map.list$lg_phi <- factor(disp.group)
    }
    if(family != "ordinal") map.list$zeta <- factor(NA)
    if((num.lv.c+num.RR)==0){
      map.list$b_lv = factor(rep(NA, length(b.lv)))
    }
    
    randoml=c(0,0,0)
    if(row.eff=="fixed"){xr <- matrix(1,1,p)} else {xr <- matrix(0,1,p)}
    if(row.eff=="random") randoml[1]=1

    if(row.eff == "random" && rstruc ==0){ nlvr<-num.lv+num.lv.c+1 } else {nlvr=num.lv+num.lv.c}
    if(row.eff=="fixed"){xr <- matrix(1,1,p)} else {xr <- matrix(0,1,p)}
    if(randomB!=FALSE){
      randoml[3]<-1
    }else{
      sigmab_lv <- 0
    }
    if(!is.null(X)){Xd <- cbind(1,X)} else {Xd <- matrix(1,n)}
    extra <- c(0,0,0)
    
    optr<-NULL
    timeo<-NULL
    se <- NULL

## Set up starting values for scale (and shape) parameters for collelated LVs
    if(num.lv.cor>0 & cstrucn>0){
      rho_lvc<- matrix(rep(0, num.lv.cor))
      if(cstruc=="corExp"){
        if(is.null(rho.lv)) {
          rho.lv=rep(0, num.lv.cor) 
        } else if(length(rho.lv)==num.lv.cor) {
          rho.lv=c(log(rho.lv))
        }
        rho_lvc<- matrix(c(rep(scaledc, each=num.lv.cor)), num.lv.cor)
      } else if(cstruc=="corMatern"){
        if(is.null(rho.lv)) {
          rho.lv=rep(log(MaternKappa), each=num.lv.cor)
        } else if(length(rho.lv)==num.lv.cor) {
          rho.lv=c(log(rho.lv))
        }
        rho_lvc<- matrix(c(rep(scaledc, each=num.lv.cor), rho.lv), num.lv.cor)
        # rho_lvc<- matrix(rho.lv,nrow = num.lv.cor)
      } 
      # else {
      #   map.list$scaledc = factor(rep(NA, length(scaledc)))
      # }
      
      if(cstrucn %in% c(2,4)){
        iv<-1:length(rho_lvc); 
        if(!is.null(setMap$rho_lvc)){
          if((length(setMap$rho_lvc)==length(rho_lvc))) 
            iv = (setMap$rho_lvc)
          map.list$rho_lvc = factor(iv)
        }
      }
      fit$rho.lv = rho_lvc
    } else {
      rho_lvc <- matrix(0)
      map.list$rho_lvc = factor(NA) 
    }
    
   ### VA method, used only if there is some random effects/LVs in the model
    
    if(((method %in% c("VA", "EVA")) && (nlvr>0 || row.eff == "random" || (randomB!=FALSE)))){

      # Variational covariances for latent variables
      if((num.lv+num.lv.c)>0){
        if(is.null(start.params) || start.params$method=="LA" || num.lv.cor>0){
          if(Lambda.struc=="diagonal" || (Lambda.struc=="bdNN") || diag.iter>0){
            Au <- log(rep(Lambda.start[1],(num.lv+num.lv.c)*n)) #1/2, 1
          } else{
            Au <- c(log(rep(Lambda.start[1],(num.lv+num.lv.c)*n)),rep(0,(num.lv+num.lv.c)*((num.lv+num.lv.c)-1)/2*n)) #1/2, 1
          }
        } else {
          Au <- NULL
          for(d in 1:(num.lv+num.lv.c)) {
            if(start.params$Lambda.struc=="unstructured" || length(dim(start.params$A))==3){
              Au <- c(Au,log(start.params$A[,d,d]))
            } else {
              Au <- c(Au,log(start.params$A[,d]))
              }
          }
          if(Lambda.struc!="diagonal" && diag.iter==0){
            Au <- c(Au,rep(0,(num.lv+num.lv.c)*((num.lv+num.lv.c)-1)/2*n))
          }
        }} else { Au <- 0}
      
      # Variational covariances for random slopes of const. ord.
      if((num.RR+num.lv.c)>0&(randomB!=FALSE)){
        if(randomB=="P"|randomB=="single"){
          ab12 <- num.RR+num.lv.c
          ab3 <- ncol(lv.X)
        }else{
          ab12 <- ncol(lv.X)
          ab3 <- num.RR+num.lv.c
        }
        if(is.null(start.params) || start.params$method=="LA" || start.params$randomB==FALSE){
          if(Lambda.struc=="diagonal" || diag.iter>0){
                Ab_lv <- log(rep(Lambda.start[1],ab12*ab3)) #1/2, 1
          } else{
            Ab_lv <- c(log(rep(Lambda.start[1],ab12*ab3)),rep(0.01,ab12*(ab12-1)/2*ab3)) #1/2, 1
          }
        } else {
          Ab_lv <- NULL
          for(d in 1:ab12) {
            if(start.params$Lambda.struc=="unstructured" || length(dim(start.params$Ab_lv))==3){
              Ab_lv <- c(Ab_lv,log(start.params$Ab_lv[,d,d]))
            } else {
              Ab_lv <- c(Ab_lv,log(start.params$Ab_lv[,d]))
            }
          }
          if(Lambda.struc!="diagonal" && diag.iter==0){
            Ab_lv <- c(Ab_lv,rep(0.01,ab12*(ab12-1)/2*ab3))
          }
        }} else { Ab_lv <- 0; map.list$Ab_lv = factor(NA)}
      
      # Variational covariances for  random rows
      if(row.eff == "random"){
        if(rstruc ==1){
          lg_Ar <- rep(log(Lambda.start[2]), nr)
        } else {
          lg_Ar <- rep(log(Lambda.start[2]), n)
        }
      
        if((rstruc == 0) && (nlvr>(num.lv+num.lv.c)) && (num.lv.cor==0) & Ar.struc!="diagonal"){
          lg_Ar<-c(lg_Ar, rep(0, (num.lv+num.lv.c)*n))
        }
        if(rstruc == 1 & (cstrucn %in% c(1,2,3,4)) & Ar.struc!="diagonal"){
          lg_Ar<-c(lg_Ar, rep(0, nr*(nr-1)/2))
        }
        if(rstruc == 2 & Ar.struc!="diagonal"){
          lg_Ar<-c(lg_Ar, rep(0, nr*times*(times-1)/2))
        }
      } else {lg_Ar <- 0}
      
      #quadratic model starting values
      if(quadratic == TRUE && start.struc == "LV"){
        start.fit <- try(gllvm.TMB(y=y, X=X, lv.X = lv.X, num.lv=num.lv, num.lv.c = num.lv.c, num.RR = num.RR, family = family, Lambda.struc = Lambda.struc, row.eff=row.eff, reltol=reltol, seed =  seed[n.i], maxit = maxit, start.lvs = start.lvs, offset = offset, n.init = 1, diag.iter=diag.iter, dependent.row=dependent.row, quadratic="LV", starting.val = starting.val, Lambda.start = Lambda.start, quad.start = quad.start, jitter.var = jitter.var, zeta.struc = zeta.struc, sd.errors = FALSE, optimizer = optimizer, optim.method = optim.method, max.iter=max.iter, start.struc="all", disp.group = disp.group, randomB = randomB),silent=T)
        if(inherits(start.fit,"try-error")&starting.val!="zero"){
          start.fit <- try(gllvm.TMB(y=y, X=X, lv.X = lv.X, num.lv=num.lv, num.lv.c = num.lv.c, num.RR = num.RR, family = family, Lambda.struc = Lambda.struc, row.eff=row.eff, reltol=reltol, seed =  seed[n.i], maxit = maxit, start.lvs = start.lvs, offset = offset, n.init = 1, diag.iter=diag.iter, dependent.row=dependent.row, quadratic="LV", starting.val = "zero", Lambda.start = Lambda.start, quad.start = quad.start, jitter.var = jitter.var, zeta.struc = zeta.struc, sd.errors = FALSE, optimizer = optimizer, optim.method = optim.method, max.iter=max.iter, start.struc="all", disp.group = disp.group, randomB = randomB),silent=T)
        }
        if(!inherits(start.fit,"try-error")&starting.val!="zero"){
          if(is.null(start.fit$lvs)){
            start.fit <- try(gllvm.TMB(y=y, X=X, lv.X = lv.X, num.lv=num.lv, num.lv.c = num.lv.c, num.RR = num.RR, family = family, Lambda.struc = Lambda.struc, row.eff=row.eff, reltol=reltol, seed =  seed[n.i], maxit = maxit, start.lvs = start.lvs, offset = offset, n.init = 1, diag.iter=diag.iter, dependent.row=dependent.row, quadratic="LV", starting.val = "zero", Lambda.start = Lambda.start, quad.start = quad.start, jitter.var = jitter.var, zeta.struc = zeta.struc, sd.errors = FALSE, optimizer = optimizer, optim.method = optim.method, max.iter=max.iter, start.struc="all", disp.group = disp.group, randomB = randomB),silent=T)
          }
        }
        if(!inherits(start.fit,"try-error")){
          if(!is.null(start.fit$lvs)){
            u <- start.fit$lvs
            fit$index <- u
          }
        }
        start.struc="all"
      }
    ### family settings
      extra[1] <- 0
      if(family == "poisson") { familyn <- 0}
      if(family == "negative.binomial") { familyn <- 1}
      if(family == "binomial") { 
        familyn <- 2
        if(link=="probit") extra[1]=1
      }
      if(family == "gaussian") {familyn=3}
      if(family == "gamma") {familyn=4}
      if(family == "tweedie"){ familyn=5; extra[1]=Power}
      if(family == "ordinal") {familyn=7}
      if(family == "exponential") {familyn=8}
      if(family == "beta"){ 
        familyn=9
        if(link=="probit") extra[1]=1
      }
      # if(family == "betaH"){ # Not used for VA atm
      #   familyn = 10
      #   bH <- rbind(a,b)
      #   if(num.lv>0) {
      #     mapLH<-factor(1:length(thetaH))
      #     mapLH[lower.tri(thetaH)] <- NA
      #     map.list$thetaH <- factor(mapLH)
      #   } else {
      #     thetaH<- matrix(0); 
      #     map.list$thetaH = factor(NA)
      #   }
      # } else {
        thetaH<- matrix(0)
        map.list$thetaH = factor(NA)
        bH <- matrix(0)
        map.list$bH = factor(NA)
      # }

      # Set up parameter.list, data.list and map.list
      # latent vars
      if((num.lv+num.lv.c)>0){
        u<-cbind(u)
      } else {
        u<-matrix(0)
        if(num.RR==0)lambda = 0
        if(num.RR==0)map.list$lambda = factor(NA)
        if(num.RR==0&quadratic==F)map.list$lambda2 = factor(NA)
        map.list$u = factor(NA) 
        map.list$Au = factor(NA) 
      }
        
      if(num.lv.cor>0){
        if(corWithin) {
          if(diag.iter>0){
            if(Astruc>=3){
              Au <- c(Au[1:(n)])
              AQ<-diag(rep(log(Lambda.start[1]),num.lv.cor),num.lv.cor)
              Au<-c(Au,AQ[lower.tri(AQ, diag = TRUE)])
            }
          } else {
          if(Lambda.struc == "unstructured" && Astruc==1) {
            Au <- c(Au[1:(n*num.lv.cor)], rep(0,sum(lower.tri(matrix(0,n,n)))*num.lv.cor) )
          } else if(Lambda.struc == "bdNN" && Astruc==2){
            Au <- c(Au[1:(n*num.lv.cor)], rep(0,nrow(NN)*num.lv.cor*nu) )
          } else if(Astruc==3) {
            Au <- c(Au[1:(n)], rep(0,sum(lower.tri(matrix(0,n,n)))) )
            AQ<-diag(rep(log(Lambda.start[1]),num.lv.cor),num.lv.cor)
            Au<-c(Au,AQ[lower.tri(AQ, diag = TRUE)])
          } else if(Astruc==4) {
            Au <- c(Au[1:(n)], rep(0,nrow(NN)*nu) )
            AQ<-diag(rep(log(Lambda.start[1]),num.lv.cor),num.lv.cor)
            Au<-c(Au,AQ[lower.tri(AQ, diag = TRUE)])
          } else if(Astruc==5) {
            Au <- c(Au[1:(n)])
            AQ<-diag(rep(log(Lambda.start[1]),num.lv.cor),num.lv.cor)
            Au<-c(Au,AQ[lower.tri(AQ, diag = TRUE)])
          }}
        } else {
          # u <- as.matrix(u[1:nu,])
          if(nrow(lvs) != nu){
            u=as.matrix((t(dr)%*%lvs/colSums(dr))[1:nu,, drop=FALSE])
          } else {
            u<-lvs
          }
          if(diag.iter>0){
            if(Astruc<3){
              Au <- c(Au[1:(nu*num.lv.cor)])
            } else {
              Au <- c(Au[1:(nu)])
              AQ<-diag(rep(log(Lambda.start[1]),num.lv.cor),num.lv.cor)
              Au<-c(Au,AQ[lower.tri(AQ, diag = TRUE)])
            }
          } else {
            if(Lambda.struc == "unstructured" && Astruc==1 & cstrucn==0 & diag.iter==0){
              Au <- c(Au[1:(nu*num.lv.cor)], rep(0, nu*num.lv.cor*(num.lv.cor-1)/2))
            } else  if(Astruc==2){
              Au <- c(Au[1:(nu*num.lv.cor)], rep(0,nrow(NN)*num.lv.cor) )
            } else  if(Astruc==3){
              Au <- c(Au[1:(nu)], rep(0,sum(lower.tri(matrix(0,nu,nu)))) )
              AQ<-diag(rep(log(Lambda.start[1]),num.lv.cor),num.lv.cor)
              Au<-c(Au,AQ[lower.tri(AQ, diag = TRUE)])
            } else  if(Astruc==4){
              Au <- c(Au[1:(nu)], rep(0,nrow(NN)) )
              AQ<-diag(rep(log(Lambda.start[1]),num.lv.cor),num.lv.cor)
              Au<-c(Au,AQ[lower.tri(AQ, diag = TRUE)])
            } else  if(Astruc==5){
              Au <- c(Au[1:(nu)] )
              AQ<-diag(rep(log(Lambda.start[1]),num.lv.cor),num.lv.cor)
              Au<-c(Au,AQ[lower.tri(AQ, diag = TRUE)])
            } else  if(Astruc==0){
              Au <- c(Au[1:(nu*num.lv.cor)])
            }
          }
        }
        # u <- u*0.5
        # lambda <- lambda *0.5
        if(is.null(start.params)) sigma.lv <- (sigma.lv*0.5) #
        # sigma.lv <- log(sigma.lv*0.5) #
        # sigma.lv <- log(sigma.lv) #
        Au = Au + 1e-3
      }

      if(num.RR==0 && num.lv.c==0) map.list$b_lv = factor(NA)
      
    ## Row effect settings
      if(row.eff=="random"){
        if(dependent.row&quadratic==F|dependent.row&starting.val=="zero") 
          sigma<-c(log(sigma), rep(0, num.lv+num.lv.c))
        if((rstruc %in% 1:2)) {
          if(cstrucn %in% c(1,3)) {
            sigma = c(log(sigma[1]),0)
          } else if(cstrucn %in% c(2)){
            sigma = c(log(sigma[1]),scaledc)
          } else if(cstrucn %in% c(4)){
            sigma = c(log(sigma[1]),scaledc,0)
          } else {
            sigma = c(log(sigma[1]))
          }
        }
        if(cstruc=="corMatern"){
          sigma<- c(sigma,log(MaternKappa))
        }
      } else {
        sigma=0
        map.list$log_sigma <- factor(NA)
        map.list$lg_Ar <- factor(NA)
        # if(row.eff != "fixed") map.list$r0 <- factor(rep(NA, length(r0)))
      }
      
      if(quadratic == FALSE){
        # if(num.RR==0&quadratic==F) map.list$lambda2 = factor(NA)
        map.list$lambda2 = factor(NA)
      }

      
  ## generate starting values quadratic coefficients in some cases
      if(starting.val!="zero" && quadratic != FALSE && (num.lv+num.lv.c+num.RR)>0){
      data.list = list(y = y, x = Xd, x_lv = lv.X, xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, num_lv_c = num.lv.c, num_RR = num.RR, num_corlv=num.lv.cor, quadratic = 1, family=familyn,extra=extra,method=switch(method, VA=0, EVA=2),model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist, Astruc=Astruc, NN = NN)
      
      # if(row.eff=="random"){
      #   if(dependent.row) sigma<-c(log(sigma), rep(0, num.lv))
      #     #parameter.list = list(r0 = matrix(r0), b = rbind(a,b), B = matrix(0), Br=Br,lambda = lambda, u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=sigma,Au=Au, lg_Ar=lg_Ar,Abb=0, zeta=zeta)
      # } else {
      #   sigma = 0 
      #   map.list$log_sigma = factor(NA)
      #   #parameter.list = list(r0 = matrix(r0), b = rbind(a,b), B = matrix(0), Br=Br,lambda = lambda, u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=0,Au=Au, lg_Ar=lg_Ar,Abb=0, zeta=zeta)
      # }
        map.list2 <- map.list 
        map.list2$sigmaLV = factor(rep(NA,length(sigma.lv)))
        map.list2$r0 = factor(rep(NA, length(r0)))
        map.list2$b_lv = factor(rep(NA, length(b.lv)))
        map.list2$Ab_lv = factor(rep(NA, length(Ab_lv)))
        map.list2$sigmab_lv = factor(rep(NA, length(sigmab_lv)))
        map.list2$b = factor(rep(NA, length(rbind(a, b))))
        map.list2$B = factor(rep(NA, 1))
        map.list2$Br = factor(rep(NA,length(Br)))
        #map.list2$lambda = factor(rep(NA, length(lambda)))
        map.list2$u = factor(rep(NA, length(u)))
        map.list2$lg_phi = factor(rep(NA, p))
        map.list2$log_sigma = factor(rep(NA, length(sigma)))
        map.list2$sigmaB = factor(rep(NA,length(sigmaB)))
        map.list2$sigmaij = factor(rep(NA,length(sigmaij)))
        map.list2$Au = factor(rep(NA, length(Au)))
        map.list2$zeta = factor(rep(NA, length(zeta)))
        map.list2$r0 = factor(rep(NA, length(r0)))
        map.list2$lg_Ar = factor(rep(NA, length(lg_Ar)))


        parameter.list = list(r0 = matrix(r0), b = rbind(a,b), bH=bH, b_lv = b.lv, sigmab_lv = sigmab_lv, Ab_lv = Ab_lv, B = matrix(0), Br=Br,lambda = lambda, lambda2 = t(lambda2), thetaH = thetaH, sigmaLV = (sigma.lv), u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=sigma, rho_lvc=rho_lvc, Au=Au, lg_Ar =lg_Ar, Abb=0, zeta=zeta) #, scaledc=scaledc

      objr <- TMB::MakeADFun(
        data = data.list, silent=TRUE,
        parameters = parameter.list, map = map.list2,
        inner.control=list(maxit = maxit), #mgcmax = 1e+200,
        DLL = "gllvm")##GLLVM
      
      if(optimizer=="nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit)),silent = TRUE))
      }
      if(optimizer=="optim" || !(optimizer %in%c("optim","nlminb") )) {
        if(optimizer == "optim" && optim.method != "BFGS")
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE))
        else
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
      }
    
      if(!inherits(optr,"try-error")){
        try({
        lamba <- optr$par[names(optr$par)=="lambda"]

        if(start.struc=="LV"|quadratic=="LV"){
          lambda2 <- matrix(optr$par[names(optr$par)=="lambda2"], byrow = T, ncol = num.lv+(num.lv.c+num.RR), nrow = 1)
        }else if(quadratic==TRUE){
          lambda2 <- matrix(optr$par[names(optr$par)=="lambda2"], byrow = T, ncol = num.lv+(num.lv.c+num.RR), nrow = p)
        }
        },silent=T)
      }
      }
    

    ## Set up data and parameters
      
      data.list = list(y = y, x = Xd, x_lv = lv.X , xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, num_lv_c = num.lv.c, num_RR = num.RR, num_corlv=num.lv.cor, quadratic = ifelse(quadratic!=FALSE,1,0), family=familyn,extra=extra,method=switch(method, VA=0, EVA=2),model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist,Astruc=Astruc, NN = NN)
      
      parameter.list = list(r0 = matrix(r0), b = rbind(a,b), bH=bH, b_lv = b.lv, sigmab_lv = sigmab_lv, Ab_lv = Ab_lv, B = matrix(0), Br=Br,lambda = lambda, lambda2 = t(lambda2),thetaH = thetaH, sigmaLV = (sigma.lv), u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=sigma, rho_lvc=rho_lvc, Au=Au, lg_Ar=lg_Ar,Abb=0, zeta=zeta) #, scaledc=scaledc

#### Call makeADFun
      objr <- TMB::MakeADFun(
        data = data.list, silent=TRUE,
        parameters = parameter.list, map = map.list,
        inner.control=list(maxit = maxit), #mgcmax = 1e+200,
        DLL = "gllvm")##GLLVM

#### Fit model 
      if((num.lv.c+num.RR)<=1|randomB!=FALSE){
      if(optimizer=="nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit)),silent = TRUE))
      }
      if(optimizer=="optim") {
        if(optim.method != "BFGS")# Due the memory issues, "BFGS" should not be used for Tweedie
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE))
        else
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
      }
      }else{
        if(optimizer == "alabama"){
          if(!optim.method%in%c("L-BFGS-B","nlminb")){
            control.optim <- list(maxit=maxit, reltol = reltol.c)
          }else if(optim.method == "L-BFGS-B"){
            control.optim <- list(maxit=maxit, factr = 1/reltol.c)
          }else if(optim.method == "nlminb"){
            control.optim <-  list(rel.tol=reltol.c,iter.max=max.iter,eval.max=maxit)
          }
          suppressWarnings(timeo <- system.time(optr <- try(auglag(objr$par, objr$fn, objr$gr, heq = eval_eq_c, heq.jac = eval_eq_j, control.optim=control.optim, control.outer = list(eps = reltol.c, itmax=maxit, trace = FALSE, kkt2.check = FALSE, method = optim.method), obj = objr),silent = TRUE)))
        }else{
          local_opts <- list( "algorithm" = optim.method,
                              "xtol_rel" = reltol,
                              "maxeval" = maxit,
                              "tol_constraints_eq" = rep(reltol.c,(num.lv.c+num.RR)*(num.lv.c+num.RR-1)/2))
          
          opts <- list( "algorithm" = optimizer,
                        "xtol_rel" = reltol,
                        "maxeval" = maxit,
                        "tol_constraints_eq" = rep(reltol.c,(num.lv.c+num.RR)*(num.lv.c+num.RR-1)/2),
                        "local_opts" = local_opts)
          timeo <- system.time(optr <- try(nloptr(x0 = objr$par, eval_f=eval_f, eval_g_eq=eval_g_eq, opts=opts, obj = objr),silent = TRUE))
          if(!inherits(optr,"try-error")){
            optr$convergence <- as.integer(optr$status<0&optr$status!=5)
            #need to return objr$env$last.par.best, because when nloptr hits maxeval it doesn't return the last set of estimates
            optr$par <- objr$env$last.par.best;names(objr$env$last.par.best) = names(optr$par) = names(objr$par);   
            if(optr$status<0){
              optr[1] <- optr$message
              class(optr) <- "try-error"
            }
          }
        }
        
      }
      if(inherits(optr,"try-error")) warning(optr[1]);

      
 ### Now diag.iter, improves the model fit sometimes
      if((diag.iter>0) && !(Lambda.struc %in% c("diagonal", "diagU")) && ((nlvr+randoml[3]*num.RR)>1) && !inherits(optr,"try-error")){
        objr1 <- objr
        optr1 <- optr
        param1 <- optr$par
        nam <- names(param1)
        if(length(param1[nam=="r0"])>0){ r1 <- matrix(param1[nam=="r0"])} else {r1 <- matrix(r0)}
        b1 <- matrix(param1[nam=="b"],num.X+1,p)
        if(randomB!=FALSE){
          sigmab_lv1 <- param1[nam=="sigmab_lv"]
          Ab_lv1<- c(pmax(param1[nam=="Ab_lv"],rep(log(1e-6), ab12*ab3)), rep(0.01,ab12*(ab12-1)/2*ab3))
          
        }else{
            sigmab_lv1<-0
            Ab_lv1 <- 0
            }
        if((num.lv.c+num.RR)>0){b.lv1 <- matrix(param1[nam=="b_lv"],ncol(lv.X),(num.lv.c+num.RR))}else{b.lv1<-matrix(0)}
        lambda1 <- param1[nam=="lambda"]
        if (quadratic=="LV" | quadratic == T && start.struc == "LV"){
          lambda2 <- matrix(param1[nam == "lambda2"], byrow = T, ncol = num.lv+(num.lv.c+num.RR), nrow = 1)#In this scenario we have estimated two quadratic coefficients before
        }else if(quadratic == T){
          lambda2 <- matrix(param1[nam == "lambda2"], byrow = T, ncol = num.lv+(num.lv.c+num.RR), nrow = p)
        }

        sigma.lv1 <- param1[nam=="sigmaLV"]
        if((num.lv+num.lv.c)>0){u1 <- matrix(param1[nam=="u"],nrow(u),num.lv+num.lv.c)}else{u1<-u}
        if(family %in% c("poisson","binomial","ordinal","exponential")){ lg_phi1 <- log(phi)} else {lg_phi1 <- param1[nam=="lg_phi"][disp.group]} #cat(range(exp(param1[nam=="lg_phi"])),"\n")
        sigmaB1 <- param1[nam=="sigmaB"]
        sigmaij1 <- param1[nam=="sigmaij"]
        if(row.eff == "random"){
          lg_Ar<- param1[nam=="lg_Ar"]
          log_sigma1 <- param1[nam=="log_sigma"]
          } else {log_sigma1 = 0}

        if(num.lv.cor>0){
          Au1<- c(param1[nam=="Au"])
          if(corWithin) {
            if(Lambda.struc == "unstructured" && Astruc==1) {
              Au1 <- c(pmax(Au1[1:(n*num.lv.cor)],log(1e-2)), rep(1e-3,sum(lower.tri(matrix(0,n,n)))*num.lv.cor) )
            } else if(Lambda.struc == "bdNN" && Astruc==2){
              Au1 <- c(pmax(Au1[1:(n*num.lv.cor)],log(1e-2)), rep(1e-3,nrow(NN)*num.lv.cor*nu) )
            } else if(Astruc==3) {
              Au1 <- c(log(exp(Au1[1:(n)])+1e-2), rep(1e-3,sum(lower.tri(matrix(0,n,n)))), Au1[-(1:n)])
            } else if(Astruc==4) {
              Au1 <- c(log(exp(Au1[1:(n)])+1e-2), rep(1e-3,nrow(NN)*nu), Au1[-(1:n)])
            }
          } else {
            if(Lambda.struc == "unstructured" && Astruc==1 & cstrucn==0){
              Au1 <- c(pmax(Au1[1:(nu*num.lv.cor)],log(1e-2)), rep(1e-3, nu*num.lv.cor*(num.lv.cor-1)/2))
            } else  if(Astruc==2){
              Au1 <- c(pmax(Au1[1:(nu*num.lv.cor)],log(1e-2)), rep(1e-3,nrow(NN)*num.lv.cor) )
            } else  if(Astruc==3){
              Au1 <- c(log(exp(Au1[1:(nu)])+1e-2), rep(1e-3,sum(lower.tri(matrix(0,nu,nu)))), Au1[-(1:nu)])
            } else  if(Astruc==4){
              Au1 <- c(log(exp(Au1[1:(nu)])+1e-2), rep(1e-3,nrow(NN)), Au1[-(1:nu)])
            }
          }
          if(cstrucn>0){
            if(cstrucn %in% c(2,4)){ #cstruc=="corExp" || cstruc=="corMatern"
              if(num.lv.cor>0){
                rho_lvc <- matrix((param1[nam=="rho_lvc"])[map.list$rho_lvc],nrow(rho_lvc),ncol(rho_lvc)); rho_lvc[is.na(rho_lvc)]=0 
              } #rho_lvc[-1]<- param1[nam=="rho_lvc"]
            } else {
              rho_lvc[1:length(rho_lvc)]<- param1[nam=="rho_lvc"]
            }
          }
        } else if((num.lv+num.lv.c)>0) {
          Au1<- c(pmax(param1[nam=="Au"],rep(log(1e-6), (num.lv+num.lv.c)*nrow(u1))), rep(0,(num.lv+num.lv.c)*((num.lv+num.lv.c)-1)/2*nrow(u1)))
        } else {Au1<-Au}

        if(family == "ordinal"){ zeta <- param1[nam=="zeta"] } else { zeta <- 0 }
        #Because then there is no next iteration
        data.list = list(y = y, x = Xd,  x_lv = lv.X, xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, num_lv_c = num.lv.c, num_RR = num.RR, num_corlv=num.lv.cor, quadratic = ifelse(quadratic!=FALSE,1,0), family=familyn,extra=extra,method=switch(method, VA=0, EVA=2),model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist, Astruc=Astruc, NN = NN)

        parameter.list = list(r0=r1, b = b1, bH=bH, b_lv = b.lv1, sigmab_lv = sigmab_lv1, Ab_lv = Ab_lv1, B=matrix(0), Br=Br,lambda = lambda1, lambda2 = t(lambda2),thetaH = thetaH, sigmaLV = sigma.lv1, u = u1,lg_phi=lg_phi1,sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=log_sigma1, rho_lvc=rho_lvc, Au=Au1, lg_Ar=lg_Ar, Abb=0, zeta=zeta) #, scaledc=scaledc
        
        objr <- TMB::MakeADFun(
          data = data.list, silent=TRUE,
          parameters = parameter.list, map = map.list,
          inner.control=list(maxit = maxit), #mgcmax = 1e+200,
          DLL = "gllvm")

        if((num.lv.c+num.RR)<=1|randomB!=FALSE){
        if(optimizer=="nlminb") {
          timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit)),silent = TRUE))
        }
        if(optimizer=="optim") {
          if(optim.method != "BFGS")
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE))
          else
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
        }
        }else{
          if(optimizer == "alabama"){
              if(!optim.method%in%c("L-BFGS-B","nlminb")){
                control.optim <- list(maxit=maxit, reltol = reltol.c)
              }else if(optim.method == "L-BFGS-B"){
                control.optim <- list(maxit=maxit, factr = 1/reltol.c)
              }else if(optim.method == "nlminb"){
                control.optim <-  list(rel.tol=reltol.c,iter.max=max.iter,eval.max=maxit)
              }
              suppressWarnings(timeo <- system.time(optr <- try(auglag(objr$par, objr$fn, objr$gr, heq = eval_eq_c, heq.jac = eval_eq_j, control.optim=control.optim, control.outer = list(eps = reltol.c, itmax=maxit, trace = FALSE, kkt2.check = FALSE, method = optim.method), obj = objr),silent = TRUE)))
              }else{
            local_opts <- list( "algorithm" = optim.method,
                                "xtol_rel" = reltol,
                                "maxeval" = maxit,
                                "tol_constraints_eq" = rep(reltol.c,(num.lv.c+num.RR)*(num.lv.c+num.RR-1)/2))
            
            opts <- list( "algorithm" = optimizer,
                          "xtol_rel" = reltol,
                          "maxeval" = maxit,
                          "tol_constraints_eq" = rep(reltol.c,(num.lv.c+num.RR)*(num.lv.c+num.RR-1)/2),
                          "local_opts" = local_opts)
            timeo <- system.time(optr <- try(nloptr(x0 = objr$par, eval_f=eval_f, eval_g_eq=eval_g_eq, opts=opts, obj = objr),silent = TRUE))
            if(!inherits(optr,"try-error")){
              optr$convergence <- as.integer(optr$status<0&optr$status!=5)
              #need to return objr$env$last.par.best, because when nloptr hits maxeval it doesn't return the last set of estimates
              optr$par <- objr$env$last.par.best; names(objr$env$last.par.best) = names(optr$par) = names(objr$par);   
              if(optr$status<0){
                optr[1] <- optr$message
                class(optr) <- "try-error"
              }
            }
          }
          
        }
        if(optimizer%in%c("nlminb","NLOPT_LD_AUGLAG","NLTOPT_LD_SLSQP")){
          if(inherits(optr, "try-error") || is.nan(optr$objective) || is.na(optr$objective)|| is.infinite(optr$objective) || optr$objective < 0){optr=optr1; objr=objr1; Lambda.struc="diagonal"}
        }else if(optimizer%in%c("optim","alabama")){
          if(inherits(optr, "try-error") || is.nan(optr$value) || is.na(optr$value)|| is.infinite(optr$value) || optr$value < 0){optr=optr1; objr=objr1; Lambda.struc="diagonal"}
        }
        if(inherits(optr,"try-error")) warning(optr[1]);
      }

      
#### Extract estimated values
      
      param<-objr$env$last.par.best
      if(family %in% c("negative.binomial", "tweedie", "gaussian", "gamma", "beta", "betaH")) {
        phis <- exp(param[names(param)=="lg_phi"])[disp.group]
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
        out$zeta.struc = zeta.struc
      }
      if((num.lv.c+num.RR)>0){
        bi.lv <- names(param)=="b_lv"
        if(randomB!=FALSE)sib <- names(param)=="sigmab_lv"
      }
      bi <- names(param)=="b"
      li <- names(param)=="lambda"
      si <- names(param) == "sigmaLV"
      li2 <- names(param)=="lambda2"
      ui <- names(param)=="u"

      if(num.lv.cor > 0){ # Correlated latent variables
        if(corWithin){
          lvs<-(matrix(param[ui],n,num.lv.cor))
        } else {
          lvs = matrix(param[ui],nu,num.lv.cor)
          rownames(lvs) =colnames(dr)
          # lvs = dr%*%matrix(param[ui],nu,num.lv.cor)
        }
        sigma.lv <- abs(param[si])
        theta <- matrix(0,p,num.lv.cor)
        if(num.lv.cor>1){
          diag(theta) <- 1 #sigma.lv 
        } else if(num.lv.cor==1) {
          theta[1,1] <- 1 #sigma.lv[1]
        }
        
        if(p>1) {
          theta[lower.tri(theta[,1:num.lv.cor,drop=F],diag=FALSE)] <- param[li];
        } else {
          theta <- as.matrix(1)
        }
        rho_lvc = param[names(param)=="rho_lvc"]
        if((cstrucn %in% c(1,3))) rho.lv<- param[names(param)=="rho_lvc"] / sqrt(1.0 + param[names(param)=="rho_lvc"]^2);
        if((cstrucn %in% c(2,4))) {
          rho.lv<- exp(param[names(param)=="rho_lvc"]);
          # scaledc<- exp(param[names(param)=="scaledc"]);
        }
      } else if((num.lv+num.lv.c+num.RR) > 0){
        if((num.lv+num.lv.c)>0)lvs<-(matrix(param[ui],n,num.lv+num.lv.c))
        theta <- matrix(0,p,num.lv+num.lv.c+num.RR)  
        if((num.lv.c+num.RR)>1){diag(theta[,1:(num.lv.c+num.RR)])<-1}else if((num.lv.c+num.RR)==1){theta[1,1]<-1}
        if(num.lv>1){diag(theta[,((num.lv.c+num.RR)+1):((num.lv.c+num.RR)+num.lv)])<-1}else if(num.lv==1){theta[1,((num.lv.c+num.RR)+1):((num.lv.c+num.RR)+num.lv)]<-1}
        if(nlvr>0)sigma.lv <- abs(param[si])
        if(num.lv>0&(num.lv.c+num.RR)==0){
        
        if(p>1) {
          theta[lower.tri(theta[,1:num.lv,drop=F],diag=FALSE)] <- param[li];
          if(quadratic!=FALSE){
            theta<-cbind(theta,matrix(-abs(param[li2]),ncol=num.lv,nrow=p,byrow=T))
          }
        } else {
          if(quadratic==FALSE){
            theta <- as.matrix(1)
          }else{
            theta <- c(as.matrix(1),-abs(param[li2]))}  
          }
        }else if(num.lv==0&(num.lv.c+num.RR)>0){
          if(p>1) {
            theta[lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE)] <- param[li];
            if(quadratic!=FALSE){
              theta<-cbind(theta,matrix(-abs(param[li2]),ncol=(num.lv.c+num.RR),nrow=p,byrow=T))
            }
          } else {
            if(quadratic==FALSE){
              theta <- as.matrix(1)
            }else{
              theta <- c(as.matrix(1),-abs(param[li2]))}  
          }
        }else if(num.lv>0&(num.lv.c+num.RR)>0){
          if(p>1) {
            theta[,1:(num.lv.c+num.RR)][lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE)] <- param[li][1:sum(lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE))];
            theta[,((num.lv.c+num.RR)+1):ncol(theta)][lower.tri(theta[,((num.lv.c+num.RR)+1):ncol(theta),drop=F],diag=FALSE)] <- param[li][(sum(lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE))+1):length(param[li])];
            if(quadratic!=FALSE){
              theta<-cbind(theta,matrix(-abs(param[li2]),ncol=num.lv+(num.lv.c+num.RR),nrow=p,byrow=T))
            }
          } else {
            if(quadratic==FALSE){
              theta <- as.matrix(1)
            }else{
              theta <- c(as.matrix(1),-abs(param[li2]))}  
          }
        }
        #diag(theta) <- exp(diag(theta)) # !!!
      }
      
      if(row.eff!=FALSE) {
        ri <- names(param)=="r0"
        row.params <- param[ri]#c(0,param[ri])
        if(row.eff=="random"){
          sigma<-exp(param[names(param)=="log_sigma"])[1]
          if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,3))) rho<- param[names(param)=="log_sigma"][2] / sqrt(1.0 + param[names(param)=="log_sigma"][2]^2);
          if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(2,4))) {
            rho<- exp(param[names(param)=="log_sigma"][-1]);
            # scaledc<- exp(param[names(param)=="scaledc"]);
          }
          if(num.lv>0 && dependent.row && rstruc==0) sigma <- c(sigma,(param[names(param)=="log_sigma"])[-1])
        }
      }
      betaM <- matrix(param[bi],p,num.X+1,byrow=TRUE)
      beta0 <- betaM[,1]
      if(!is.null(X)) betas <- betaM[,-1]
      if((num.lv.c+num.RR)>0)b.lv <- matrix(param[bi.lv],ncol(lv.X),(num.lv.c+num.RR))
      # if(family %in% "betaH"){
      #   bHi <- names(param)=="bH"
      #   betaH <- matrix(param[bHi],p,num.X+1,byrow=TRUE)
      # }
      if((randomB!=FALSE)&(num.lv.c+num.RR)>0)sigmab_lv <- exp(param[sib])
      new.loglik <- objr$env$value.best[1]

    }
    
    
## Laplace method / nlvr==0
    if(method=="LA" || (nlvr==0 && (method %in% c("VA", "EVA")) && row.eff!="random" && (randomB==FALSE))){
      if(!is.null(X)){Xd=cbind(1,X)} else {Xd=matrix(1,n)}
### Family settings
  
      extra[1]=0
      if(family == "poisson") {familyn=0}
      if(family == "negative.binomial") {familyn=1}
      if(family == "binomial") {
        familyn=2;
        if(link=="probit") extra[1]=1
      }
      if(family == "gaussian") {familyn=3}
      if(family == "gamma") {familyn=4}
      if(family == "tweedie"){ familyn=5; extra[1]=Power}
      if(family == "ZIP"){ familyn=6;}
      if(family == "ordinal"){ familyn=7}
      if(family == "exponential"){ familyn=8}
      if(family == "beta"){ 
        familyn=9
        if(link=="probit") extra[1]=1
      }
      if(family == "betaH"){
        familyn = 10
        if(link=="probit") extra[1]=1
        
        bH <- rbind(a,b)
        if(num.lv>0) {
          mapLH<-factor(1:length(thetaH))
          mapLH[lower.tri(thetaH)] <- NA
          map.list$thetaH <- factor(mapLH)
        } else {
          thetaH<- matrix(0); 
          map.list$thetaH = factor(NA)
        }
      } else {
        thetaH<- matrix(0)
        map.list$thetaH = factor(NA)
        bH <- matrix(0)
        map.list$bH = factor(NA)
      }
      
      ## generate starting values quadratic coefficients in some cases
      if(starting.val!="zero" && quadratic == TRUE && num.RR>0&(num.lv+num.lv.c)==0 && start.struc=="LV"){
        data.list = list(y = y, x = Xd, x_lv = lv.X, xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, num_lv_c = num.lv.c, num_RR = num.RR, num_corlv=num.lv.cor, quadratic = 1, family=familyn,extra=extra,method=switch(method, VA=0, EVA=2),model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist, Astruc=Astruc, NN = NN)
        
        map.list2 <- map.list 
        map.list2$log_sigma = factor(NA)
        map.list2$Au <- map.list2$Abb <- map.list2$Ab_lv <- factor(NA)# map.list$lambda2 <- 
        map.list2$b_lv <- factor(rep(NA,length(b.lv)))

        parameter.list = list(r0 = matrix(r0), b = rbind(a,b), bH=bH, b_lv = b.lv, sigmab_lv = 0, Ab_lv = Ab_lv, B = matrix(0), Br=Br,lambda = lambda, lambda2 = t(lambda2), thetaH = thetaH, sigmaLV = (sigma.lv), u = u,lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=sigma,rho_lvc=rho_lvc, Au=0, lg_Ar =0, Abb=0, zeta=zeta) #, scaledc=scaledc
        
        objr <- TMB::MakeADFun(
          data = data.list, silent=TRUE,
          parameters = parameter.list, map = map.list2,
          inner.control=list(maxit = maxit), #mgcmax = 1e+200,
          DLL = "gllvm")##GLLVM

        if(optimizer=="nlminb") {
          timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit)),silent = TRUE))
        }
        if(optimizer=="optim" | !(optimizer %in% c("optim","nlminb"))) {
          if( optimizer == "optim" && optim.method != "BFGS")
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE))
          else
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
        }

        if(!inherits(optr,"try-error")){
          # lambda <- optr$par[names(optr$par)=="lambda"]
          try({lambda2 <- matrix(optr$par[names(optr$par)=="lambda2"],ncol=num.RR,nrow=p,byrow=T)},silent=T)
          # b.lv <- matrix(objr$par[names(objr$par)=="b_lv"],ncol=num.RR)
          # fit$params[,2:(1+num.RR)][lower.tri(fit$params[,2:(1+num.RR)],diag=F)] <- lambda
          fit$params[,(ncol(fit$params)-num.RR+1):ncol(fit$params)] <- lambda2
          # fit$b.lv <- b.lv
        }
      }
      data.list = list(y = y, x = Xd, x_lv = lv.X, xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, num_lv_c = num.lv.c, num_RR = num.RR, num_corlv=num.lv.cor, quadratic = ifelse(quadratic!=FALSE,1,0), family=familyn,extra=extra,method=1,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist, Astruc=Astruc, NN = NN)
      
      if(family == "ordinal"){
        data.list$method = 0
      }
      
      randomp <- "u"
      map.list$Au <- map.list$lg_Ar <- map.list$Abb
      map.list$Ab_lv = factor(NA)
      if(quadratic==FALSE) map.list$lambda2 <- factor(NA)
      
      randomp <- NULL
      # latent vars
      if((num.lv+num.lv.c)>0){
        u<-cbind(u)
        randomp <- c(randomp,"u")
      } else {
        u = matrix(0)
        if(num.RR==0)lambda = 0
        if(num.RR==0)map.list$lambda = factor(NA)
        if(num.RR==0&quadratic==F)map.list$lambda2 = factor(NA)
        map.list$u = factor(NA) 
      }
      
      if(num.lv.cor>0 & (!corWithin)){
        u <- as.matrix(u[1:nu,])
      }
      
      # Row parameter settings
      if(row.eff=="random"){
        randoml[1] <- 1
        randomp <- c(randomp,"r0")
        
        if(dependent.row && (rstruc == 0)) 
          sigma<-c(log(sigma[1]), rep(0, num.lv+num.lv.c))
        if((rstruc %in% 1:2) & (cstrucn %in% c(1,2,3,4))) {
          sigma = c(log(sigma[1]),0)
        }
        if(cstruc=="corMatern"){
          sigma<- c(sigma,log(MaternKappa))
        }
      } else {
        sigma=0
        map.list$log_sigma <- factor(NA)
        # if(row.eff != "random") map.list$r0 <- factor(rep(NA, length(r0)))
      }
      if(randomB!=FALSE){
        randoml[3] <- 1
        randomp <- c(randomp, "b_lv")
      }else{
        sigmab_lv <- 0
      }
      
#### Set up data and parameters
      
      if(family == "ordinal"){
        data.list$method = 0
      }
      parameter.list = list(r0=matrix(r0), b = rbind(a,b), bH=bH, b_lv = b.lv, sigmab_lv = sigmab_lv, Ab_lv = Ab_lv, B=matrix(0), Br=Br,lambda = lambda, lambda2 = t(lambda2),thetaH = thetaH, sigmaLV = (sigma.lv), u = u, lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=c(sigma), rho_lvc=rho_lvc, Au=0, lg_Ar=0, Abb=0, zeta=zeta) #, scaledc=scaledc
      
      #### Call makeADFun
      objr <- TMB::MakeADFun(
        data = data.list, silent=!trace,
        parameters = parameter.list, map = map.list,
        inner.control=list(mgcmax = 1e+200,maxit = maxit,tol10=0.01),
        random = randomp, DLL = "gllvm")
     #### Fit model 
      
      # Not used for now
      # if(family=="ZIP" && FALSE) {
      #   m <- length(objr$par)
      #   low <- rep(-restrict,m); upp=rep(restrict,m);
      #   low[names(objr$par)=="lg_phi"]=0.0; upp[names(objr$par)=="lg_phi"]=1#0.99
      #   timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit),lower = low,upper = upp),silent = TRUE))
      # }
      if((num.lv.c+num.RR)<=1|randomB!=FALSE){
        if(optimizer=="nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit)),silent = TRUE))
      }
      if(optimizer=="optim") {
        if(optim.method != "BFGS")
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE))
        else
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
      }
      }else{
        if(optimizer == "alabama"){
          if(!optim.method%in%c("L-BFGS-B","nlminb")){
            control.optim <- list(maxit=maxit, reltol = reltol.c)
          }else if(optim.method == "L-BFGS-B"){
            control.optim <- list(maxit=maxit, factr = 1/reltol.c)
          }else if(optim.method == "nlminb"){
            control.optim <-  list(rel.tol=reltol.c,iter.max=max.iter,eval.max=maxit)
          }
          suppressWarnings(timeo <- system.time(optr <- try(auglag(objr$par, objr$fn, objr$gr, heq = eval_eq_c, heq.jac = eval_eq_j, control.optim=control.optim, control.outer = list(eps = reltol.c, itmax=maxit, trace = FALSE, kkt2.check = FALSE, method = optim.method), obj = objr),silent = TRUE)))
        }else{
          local_opts <- list( "algorithm" = optim.method,
                              "xtol_rel" = reltol,
                              "maxeval" = maxit,
                              "tol_constraints_eq" = rep(reltol.c,(num.lv.c+num.RR)*(num.lv.c+num.RR-1)/2))
          
          opts <- list( "algorithm" = optimizer,
                        "xtol_rel" = reltol,
                        "maxeval" = maxit,
                        "tol_constraints_eq" = rep(reltol.c,(num.lv.c+num.RR)*(num.lv.c+num.RR-1)/2),
                        "local_opts" = local_opts)
          timeo <- system.time(optr <- try(nloptr(x0 = objr$par, eval_f=eval_f, eval_g_eq=eval_g_eq, opts=opts, obj = objr),silent = TRUE))
          if(!inherits(optr,"try-error")){
            optr$convergence <- as.integer(optr$status<0&optr$status!=5)
            #need to return objr$env$last.par.best, because when nloptr hits maxeval it doesn't return the last set of estimates
            optr$par <- objr$env$last.par.best[!objr$env$lrandom()]; names(objr$env$last.par.best) <- rep(objr$env$parNameOrder,unlist(lapply(objr$env$parameters,length))); names(optr$par) = names(objr$par);   
            if(optr$status<0){
              optr[1] <- optr$message
              class(optr) <- "try-error"
            }
          }
        }
        
      }
      
      if(inherits(optr,"try-error")) warning(optr[1]);
      
      
      if(quadratic == TRUE && starting.val=="zero" && start.struc=="LV" & num.RR>0){
        if(family == "ordinal"){
          data.list$method = 0
        }
        
        lambda <- objr$env$last.par.best[names(objr$env$last.par.best)=="lambda"]
        lambda2 <- matrix(objr$env$last.par.best[names(objr$env$last.par.best)=="lambda2"],ncol=num.RR,nrow=p,byrow=T)
        b.lv <- matrix(objr$env$last.par.best[names(objr$env$last.par.best)=="b_lv"],ncol=num.RR,nrow=ncol(lv.X))
        if(randomB!=FALSE){
          sigmab_lv <- objr$env$last.par.best[names(objr$env$last.par.best)=="sigmab_lv"]
        }else{
          sigmab_lv <- 0
        }
        b <- matrix(objr$env$last.par.best[names(objr$env$last.par.best)=="b"],num.X+1,p)

        if(!(family %in% c("poisson","binomial","ordinal","exponential"))) phi <- exp(objr$env$last.par.best[names(objr$env$last.par.best)=="lg_phi"])[disp.group]
        parameter.list = list(r0=matrix(r0), b = b, bH=bH, b_lv = b.lv, sigmab_lv = sigmab_lv, Ab_lv = Ab_lv, B=matrix(0), Br=Br,lambda = lambda, lambda2 = t(lambda2),thetaH = thetaH, sigmaLV = (sigma.lv), u = u, lg_phi=log(phi),sigmaB=log(diag(sigmaB)),sigmaij=sigmaij,log_sigma=c(sigma), rho_lvc=rho_lvc, Au=0, lg_Ar=0, Abb=0, zeta=zeta) #, scaledc=scaledc

        #### Call makeADFun
        objr <- TMB::MakeADFun(
          data = data.list, silent=!trace,
          parameters = parameter.list, map = map.list,
          inner.control=list(mgcmax = 1e+200,maxit = maxit,tol10=0.01),
          random = randomp, DLL = "gllvm")

        if((num.lv.c+num.RR)<=1|randomB!=FALSE){
        if(optimizer=="nlminb") {
          timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit)),silent = TRUE))
        }
        if(optimizer=="optim") {
          if(optim.method != "BFGS")
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE))
          else
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
        }
        }else{
          if(optimizer == "alabama"){
            if(!optim.method%in%c("L-BFGS-B","nlminb")){
              control.optim <- list(maxit=maxit, reltol = reltol.c)
            }else if(optim.method == "L-BFGS-B"){
              control.optim <- list(maxit=maxit, factr = 1/reltol.c)
            }else if(optim.method == "nlminb"){
              control.optim <-  list(rel.tol=reltol.c,iter.max=max.iter,eval.max=maxit)
            }
            suppressWarnings(timeo <- system.time(optr <- try(auglag(objr$par, objr$fn, objr$gr, heq = eval_eq_c, heq.jac = eval_eq_j, control.optim=control.optim, control.outer = list(eps = reltol.c, itmax=maxit, trace = FALSE, kkt2.check = FALSE, method = optim.method), obj = objr),silent = TRUE)))
          }else{
            local_opts <- list( "algorithm" = optim.method,
                                "xtol_rel" = reltol,
                                "maxeval" = maxit,
                                "tol_constraints_eq" = rep(reltol.c,(num.lv.c+num.RR)*(num.lv.c+num.RR-1)/2))
            
            opts <- list( "algorithm" = optimizer,
                          "xtol_rel" = reltol,
                          "maxeval" = maxit,
                          "tol_constraints_eq" = rep(reltol.c,(num.lv.c+num.RR)*(num.lv.c+num.RR-1)/2),
                          "local_opts" = local_opts)
            timeo <- system.time(optr <- try(nloptr(x0 = objr$par, eval_f=eval_f, eval_g_eq=eval_g_eq, opts=opts, obj = objr),silent = TRUE))
            if(!inherits(optr,"try-error")){
              optr$convergence <- as.integer(optr$status<0&optr$status!=5)
              #need to return objr$env$last.par.best, because when nloptr hits maxeval it doesn't return the last set of estimates
              optr$par <- objr$env$last.par.best; names(objr$env$last.par.best) <- rep(objr$env$parNameOrder,unlist(lapply(objr$env$parameters,length))); names(optr$par) = names(objr$par);   
              if(optr$status<0){
                optr[1] <- optr$message
                class(optr) <- "try-error"
              }
            }
          }
          
        }
        
        if(inherits(optr,"try-error")) warning(optr[1]);
        
      }

      
#### Extract estimated values
      param <- objr$env$last.par.best
      if((num.lv.c+num.RR)>0){
        bi.lv <- names(param)=="b_lv"
        if(randomB!=FALSE)sib <- names(param)=="sigmab_lv"
      }
      bi <- names(param)=="b"
      li <- names(param)=="lambda"
      li2 <- names(param)=="lambda2"
      si <- names(param) == "sigmaLV"
      ui <- names(param)=="u"
      
      if(num.lv.cor > 0){
        if((num.lv.cor)>0 & corWithin){
          lvs<-(matrix(param[ui],n,num.lv.cor))
        } else{ 
          lvs = matrix(param[ui],nu,num.lv.cor)
          rownames(lvs) =colnames(dr)
          # lvs = dr%*%matrix(param[ui],nu,num.lv.cor)
        }
        
        theta <- matrix(0,p,num.lv.cor)
        sigma.lv <- abs(param[si])
        
        if(num.lv.cor>1){
          diag(theta) <- 1 #sigma.lv 
        } else if(num.lv.cor==1) {
          theta[1,1] <- 1 #sigma.lv[1]
        }

        if(p>1) {
          theta[lower.tri(theta[,1:num.lv.cor,drop=F],diag=FALSE)] <- param[li];
        } else {
          theta <- as.matrix(1)
        }
        rho_lvc = as.matrix(param[names(param)=="rho_lvc"])
        if((cstrucn %in% c(1,3))) rho.lv<- param[names(param)=="rho_lvc"] / sqrt(1.0 + param[names(param)=="rho_lvc"]^2);
        if((cstrucn %in% c(2,4))) {
          rho.lv<- exp(param[names(param)=="rho_lvc"]);
          # scaledc<- exp(param[names(param)=="scaledc"]);
        }
        
      } else if((num.lv+num.lv.c+num.RR) > 0){
        if((num.lv.c+num.lv)>0)lvs<-(matrix(param[ui],n,num.lv+num.lv.c))
        theta <- matrix(0,p,num.lv+(num.lv.c+num.RR))
        if((num.lv.c+num.RR)>1){diag(theta[,1:(num.lv.c+num.RR)])<-1}else if((num.lv.c+num.RR)==1){theta[1,1]<-1}
        if(num.lv>1){diag(theta[,((num.lv.c+num.RR)+1):((num.lv.c+num.RR)+num.lv)])<-1}else if(num.lv==1){theta[1,((num.lv.c+num.RR)+1):((num.lv.c+num.RR)+num.lv)]<-1}
        if((num.lv+num.lv.c)>0) sigma.lv <- abs(param[si])
        if(num.lv>0&(num.lv.c+num.RR)==0){
          if(p>1) {
            theta[,1:num.lv][lower.tri(theta[,1:num.lv,drop=F],diag=FALSE)] <- param[li];
            if(quadratic!=FALSE){
              theta<-cbind(theta,matrix(-abs(param[li2]),ncol=num.lv,nrow=p,byrow=T))
            }
          } else {
            if(quadratic==FALSE){
              theta <- as.matrix(1)
            }else{
              theta <- c(as.matrix(1),-abs(param[li2]))}  
          }
        }else if(num.lv==0&(num.lv.c+num.RR)>0){
          if(p>1) {
            theta[,1:(num.lv.c+num.RR)][lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE)] <- param[li];
            if(quadratic!=FALSE){
              theta<-cbind(theta,matrix(-abs(param[li2]),ncol=(num.lv.c+num.RR),nrow=p,byrow=T))
            }
          } else {
            if(quadratic==FALSE){
              theta <- as.matrix(1)
            }else{
              theta <- c(as.matrix(1),-abs(param[li2]))}  
          }
        }else if(num.lv>0&(num.lv.c+num.RR)>0){
          if(p>1) {
            theta[,1:(num.lv.c+num.RR)][lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE)] <- param[li][1:sum(lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE))];
            theta[,((num.lv.c+num.RR)+1):ncol(theta)][lower.tri(theta[,((num.lv.c+num.RR)+1):ncol(theta),drop=F],diag=FALSE)] <- param[li][(sum(lower.tri(theta[,1:(num.lv.c+num.RR),drop=F],diag=FALSE))+1):length(param[li])];
            if(quadratic!=FALSE){
              theta<-cbind(theta,matrix(-abs(param[li2]),ncol=num.lv+(num.lv.c+num.RR),nrow=p,byrow=T))
            }
          } else {
            if(quadratic==FALSE){
              theta <- as.matrix(1)
            }else{
              theta <- c(as.matrix(1),-abs(param[li2]))}  
          }
        }
        #diag(theta) <- exp(diag(theta)) # !!!

      }
      
      if(row.eff!=FALSE) {
        ri <- names(param)=="r0"
        row.params=param[ri]
        if(row.eff=="random"){
          sigma<-exp(param[names(param)=="log_sigma"])[1]
          if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,3))) rho<- param[names(param)=="log_sigma"][2] / sqrt(1.0 + param[names(param)=="log_sigma"][2]^2);
          if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(2,4))) {
            rho<- exp(param[names(param)=="log_sigma"][-1]);
            # scaledc<- exp(param[names(param)=="scaledc"]);
            }
          if((num.lv+num.lv.c)>0 && dependent.row && rstruc==0) sigma <- c(sigma, (param[names(param)=="log_sigma"])[-1])
        }
      }

      betaM <- matrix(param[bi],p,num.X+1,byrow=TRUE)
      beta0 <- betaM[,1]
      if(!is.null(X)) betas=betaM[,-1]
      if((num.lv.c+num.RR)>0){
        b.lv <- matrix(param[bi.lv],ncol(lv.X),(num.lv.c+num.RR))
        if(randomB!=FALSE)sigmab_lv <- exp(param[sib])
      }
      
      new.loglik <- objr$env$value.best[1]
      
      if(family %in% c("negative.binomial", "tweedie", "ZIP", "gaussian", "gamma", "beta", "betaH")) {
        phis <- exp(param[names(param)=="lg_phi"])[disp.group]
        if(family=="ZIP") {
          lp0 <- param[names(param)=="lg_phi"][disp.group]; out$lp0 <- lp0
          phis <- exp(lp0)/(1+exp(lp0));
        }
      }
      if(family %in% "betaH"){
        bHi <- names(param)=="bH"
        betaH <- matrix(param[bHi],p,num.X+1,byrow=TRUE)
        if(num.lv>0) {
          thetaH[!is.na(map.list$thetaH)] <- param[names(param)=="thetaH"]
        }
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
    }


#### Check if model fit succeeded/improved on this iteration n.i
    out$start <- fit
    
    # Gradient check with n.i >2 so we don't get poorly converged models - relatively relaxed tolerance
    if(n.i>1){
      if(!is.null(objrFinal)){
        gr1 <- objrFinal$gr()
        gr1 <- gr1/length(gr1)
        norm.gr1 <- norm(as.matrix(gr1))
      }else{
        gr1 <- NaN
        norm.gr1 <- NaN
      }
     
      gr2 <- objr$gr()
      gr2 <- gr2/length(gr2)
      norm.gr2 <- norm(as.matrix(gr2))
      n.i.i <- n.i.i +1
      grad.test1 <- all.equal(norm.gr1, norm.gr2, tolerance = 1, scale = 1)#check if gradients are similar when accepting on log-likelihood
      grad.test2 <- all.equal(norm.gr1, norm.gr2, tolerance = .1, scale = 1)#check if gradient are (sufficiently) different from each other, when accepting on gradient. Slightly more strict for norm(gr2)<norm(gr1)
      }else{
      n.i.i <- 0
    }
    if(n.i.i>n.init.max){
      n.init <- n.i
      warning("n.init.max reached after ", n.i, " iterations.")
    }

    if((n.i==1 || ((is.nan(norm.gr1) && !is.nan(norm.gr2)) || !is.nan(norm.gr2) && ((isTRUE(grad.test1) && out$logL > (new.loglik)) || (!isTRUE(grad.test2) && norm.gr2<norm.gr1))))  && is.finite(new.loglik) && !inherits(optr, "try-error")){
      objrFinal<-objr1 <- objr; optrFinal<-optr1<-optr;n.i.i<-0;
      out$start <- fit
      out$logL <- new.loglik
      if((num.lv+(num.lv.c+num.RR)) > 0) {
        if((num.lv+num.lv.c)>0)out$lvs <- lvs
        out$params$theta <- theta
        if((num.lv+num.lv.c)>0) out$params$sigma.lv  <- sigma.lv
        if((num.lv.c+num.RR)>0){
          out$params$LvXcoef <- b.lv
          colnames(out$params$LvXcoef) <- paste("CLV",1:(num.lv.c+num.RR), sep="")
          row.names(out$params$LvXcoef) <- colnames(lv.X)
          if(randomB!=FALSE){
            out$params$sigmaLvXcoef <- sigmab_lv
            if(randomB=="LV")names(out$params$sigmaLvXcoef) <- paste("CLV",1:(num.lv.c+num.RR), sep="")
            if(randomB=="P")names(out$params$sigmaLvXcoef) <- colnames(lv.X)
            # if(randomB=="all")names(out$params$sigmaLvXcoef) <- paste(paste("CLV",1:(num.lv.c+num.RR),sep=""),rep(colnames(lv.X),each=num.RR+num.lv.c),sep=".")
            if(randomB=="single")names(out$params$sigmaLvXcoef) <- NULL
            
          }
        }
        
        if((num.lv+num.lv.c)>0 & !is.null(out$lvs)) if((nrow(out$lvs)==nrow(out$y))) rownames(out$lvs) <- rownames(out$y);
        if(num.lv>0&(num.lv.c+num.RR)==0) {
          if(quadratic==FALSE){
            colnames(out$params$theta)<- paste("LV", 1:num.lv, sep="")
            colnames(out$lvs) <- paste("LV", 1:num.lv, sep="")};
          if(quadratic!=FALSE){
            colnames(out$lvs) <- paste("LV", 1:num.lv, sep="");
            colnames(out$params$theta)<- c(paste("LV", 1:num.lv, sep=""),paste("LV", 1:num.lv, "^2",sep=""));
          }
        rownames(out$params$theta) <- colnames(out$y)
        }else if((num.lv.c+num.RR)>0&num.lv==0) {
          if(quadratic==FALSE){
            if(num.lv.c>0){
              colnames(out$lvs) <- paste("CLV", 1:num.lv.c, sep="")
            }
          colnames(out$params$theta) <- paste("CLV", 1:(num.lv.c+num.RR), sep="")
          }
          if(quadratic!=FALSE){
            if(num.lv.c>0)colnames(out$lvs) <- paste("CLV", 1:num.lv.c, sep="");
            colnames(out$params$theta)<- c(paste("CLV", 1:(num.lv.c+num.RR), sep=""),paste("CLV", 1:(num.lv.c+num.RR), "^2",sep=""));
          }
          rownames(out$params$theta) <- colnames(out$y)
        }else if(num.lv>=1&(num.lv.c+num.RR)>=1){
          if(quadratic==FALSE){
            colnames(out$params$theta)<- c(paste("CLV", 1:(num.lv.c+num.RR), sep=""),paste("LV", 1:num.lv, sep=""))
            if((num.lv+num.lv.c)>0){
              if(num.lv>0&num.lv.c>0){
                colnames(out$lvs)<- c(paste("CLV", 1:num.lv.c, sep=""),paste("LV", 1:num.lv, sep=""))
              }else if(num.lv>0&num.lv.c==0){
                colnames(out$lvs)<- paste("LV", 1:num.lv, sep="")
              }
              
            }
            };
          if(quadratic!=FALSE){
            if(num.lv.c>0&num.lv==0){colnames(out$lvs) <-  paste("CLV", 1:num.lv.c, sep="")
            }else if(num.lv>0&num.lv.c==0){colnames(out$lvs) <- paste("LV", 1:num.lv, sep="")
            }else if(num.lv>0&num.lv.c>0){colnames(out$lvs) <- c(paste("CLV", 1:num.lv.c, sep=""),paste("LV", 1:num.lv, sep=""))}
            
            colnames(out$params$theta)<- c(paste("CLV", 1:(num.lv.c+num.RR), sep=""),paste("LV", 1:num.lv, sep=""),paste("CLV", 1:(num.lv.c+num.RR), "^2",sep=""),paste("LV", 1:num.lv, "^2",sep=""));
          }
          rownames(out$params$theta) <- colnames(out$y)
        }
        if((num.lv+num.lv.c)>0){
          if(num.lv>0&num.lv.c>0){
            names(out$params$sigma.lv)<- c(paste("CLV", 1:num.lv.c, sep=""),paste("LV", 1:num.lv, sep=""))
          }else if(num.lv>0&num.lv.c==0){
            names(out$params$sigma.lv)<- paste("LV", 1:num.lv, sep="")
          }else if(num.lv.c>0&num.lv==0){
            names(out$params$sigma.lv)<- paste("LV", 1:num.lv.c, sep="")
          }
        }
      }
      names(beta0) <- colnames(out$y); out$params$beta0 <- beta0;
      if(!is.null(X)){
        betas <- matrix(betas,ncol=ncol(X)); out$params$Xcoef <- betas;
        rownames(out$params$Xcoef) <- colnames(out$y); colnames(out$params$Xcoef) <- colnames(X); 
      }

      if(family %in% "betaH"){
        out$params$betaH <- betaH;
        rownames(out$params$betaH) <- colnames(out$y); 
        colnames(out$params$betaH) <- paste("c", 1:ncol(betaH))
        # colnames(out$params$betaH)[1] <- "Intercept"; 
        if(!is.null(X)){ colnames(out$params$betaH) <- c("Intercept",colnames(X)); }
        if(num.lv>0) {
          out$params$thetaH <- thetaH
        }      
      }
      if(family =="negative.binomial") {
        out$params$inv.phi <- phis;
        out$params$phi <- 1/phis;
          names(out$params$phi) <- colnames(y);
        if(!is.null(names(disp.group))){
          try(names(out$params$phi) <- names(disp.group),silent=T)
        }
        names(out$params$inv.phi) <-  names(out$params$phi)
      }
      if(family %in% c("gaussian", "tweedie", "gamma","beta", "betaH")) {
        out$params$phi <- phis;
        names(out$params$phi) <- colnames(y);
        if(!is.null(names(disp.group))){
          try(names(out$params$phi) <- names(disp.group),silent=T)
        }
      }
      if(family =="ZIP") {
        out$params$phi <- phis;
        names(out$params$phi) <- colnames(y);
        
        if(!is.null(names(disp.group))){
          try(names(out$params$phi) <- names(disp.group),silent=T)
        }
        
      }
      if(row.eff!=FALSE) {
        if(row.eff=="random"){ 
          out$params$sigma=sigma; 
          names(out$params$sigma)="sigma"
          if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,2,3,4))){ 
            out$params$rho <- rho
            # if(cstrucn %in% c(2,4)){ out$params$scaledc=scaledc}
          }
          if((num.lv+num.lv.c)>1 && dependent.row) names(out$params$sigma) <- paste("sigma",c("",1:(num.lv+num.lv.c)), sep = "")
        }
        out$params$row.params <- row.params; 
        if(length(row.params) == n) names(out$params$row.params) <- rownames(out$y)
        if((length(row.params) == ncol(dr)) && (rstruc==1)) try(names(out$params$row.params) <- colnames(dr), silent = TRUE)
      }
      if(num.lv.cor>0 & cstrucn>0){
        out$params$rho.lv <- rho.lv; 
        if(cstrucn %in% c(2,4)){ 
          if(length(out$params$rho.lv)>0) 
            names(out$params$rho.lv) <- paste("rho.lv",1:length(out$params$rho.lv), sep = "") #[!is.na(map.list$rho_lvc)]
        } else {
          names(out$params$rho.lv) <- paste("rho.lv",1:num.lv.cor, sep = "") 
        }
      }

      if(family %in% c("binomial", "beta")) out$link <- link;
      if(family == "tweedie") out$Power <- Power;
      if(family=="ordinal"){
        out$params$zeta <- zetas
      }
      out$row.eff <- row.eff
      out$time <- timeo
      pars <- optr$par

## Collect VA covariances
      if((method %in% c("VA", "EVA"))){
        param <- objr$env$last.par.best
        
        if(num.lv.cor>0 && !corWithin){
          Au <- param[names(param)=="Au"]
          AQ <- NULL
          
          if(cstrucn==0){
            A <- array(0, dim=c(nu, num.lv.cor, num.lv.cor))
            for (d in 1:(num.lv.cor)){
              for(i in 1:nu){
                A[i,d,d] <- exp(Au[(d-1)*nu+i]);
              }
            }
            if(Astruc>0 & (length(Au)>((num.lv.cor)*nu))){ # var cov Unstructured
              k=0;
              for (d in  1:num.lv.cor){
                r=d+1
                while (r <= num.lv.cor){
                  for(i in 1:nu){
                    A[i,r,d]=Au[nu*num.lv.cor+k*nu+i];
                  }
                  k=k+1; r=r+1
                }}
            }
            
            for(i in 1:nu){
              A[i,,] <- A[i,,]%*%t(A[i,,])
            }
          } else {
            # A <- array(0, dim=c(nu, nu, num.lv.cor))
            if(Astruc<3){
              nMax<- num.lv.cor
            } else {
              nMax<- 1
            }
            A <- array(0, dim=c(nu, nu, nMax))
            Alvm <- objr$report()$Alvm
            
            
            for (d in 1:nMax) {
              # for (i in 1:nu){
              #   A[i,i,d]=exp(Au[(d-1)*nu+i]);
              # }
              # if(Astruc>0 & length(Au)>(nu*nMax)){
              #   k=0;
              #   if(Astruc %in% c(1,3)){ # var cov struct unstructured
              #     for (i in  1:nu){
              #       r=i+1
              #       while (r <= nu){
              #         A[r,i,d]=Au[nu*nMax+k*nMax+d];
              #         k=k+1; r=r+1;
              #       }
              #     }
              #   } else if(Astruc %in% c(2,4)) { # var cov struct NN
              #     arank = nrow(NN);
              #     for (r in (1:arank)){
              #       A[NN[r,1],NN[r,2],d]=Au[nu*nMax+k*nMax+d];
              #       k=k+1;
              #     }
              #   }
              # 
              # }
              if(Astruc %in% c(3,4)){
                A[,,d] <- Alvm #%*%t(Alvm)
              } else {
                A[,,d] <- Alvm[,,d]%*%t(Alvm[,,d])
              }
            }
            if(Astruc %in% c(3,4)){
              AQ <- matrix(0,num.lv.cor,num.lv.cor)
              # for (d in 1:num.lv.cor){
              #   AQ[d,d]=exp(Au[(nu+k)]);
              #   k=k+1;
              #   r=d+1
              #   while (r <= num.lv.cor){
              #     AQ[r,d]=Au[(nu+k)];
              #     k=k+1; r=r+1
              #   }
              # }
              AQ <- objr$report()$AQ
              # AQ<-AQ%*%t(AQ)
            }
            
            # for(d in 1:nMax){ #num.lv.cor
            #   A[,,d] <- A[,,d]%*%t(A[,,d])
            # }
          }
          out$A <- A
          out$AQ <- AQ
          
        } else if(num.lv.cor>0 && corWithin){
          Au <- param[names(param)=="Au"]
          # A <- array(0, dim=c(times*nu,times*nu,num.lv.cor))
          if(Astruc<3){ 
            nMax<- num.lv.cor
          } else {
            nMax<- 1
          }
          A <- array(0, dim=c(times*nu, times*nu, nMax))
          Alvm <- objr$report()$Alvm
          
          AQ <- NULL
          
          for (q in 1:nMax) {
            # for (i in  1:nu){
            #   for (r in (1:times)){
            #     A[(i-1)*times+r,(i-1)*times+r,q]=exp(Au[(q-1)*n+(i-1)*times+r]);
            #   }
            # }
            # if(Astruc>0){#var cov
            #   k=0;
            #   if(Astruc %in% c(1,3)){ # var cov struct unstructured
            #     for(i in 1:nu){
            #       for (d in 1:times){
            #         r=d+1
            #         while (r<=(times)){
            #           A[(i-1)*times+r,(i-1)*times+d,q]=Au[nu*times*nMax+k*nMax+q];
            #           k=k+1; r=r+1
            #         }
            #       }
            #     }
            #   } else if(Astruc %in% c(2,4)) { # var cov struct NN
            #     arank = nrow(NN);
            #     for(i in 1:nu){
            #       for (r in (1:arank)){
            #         A[(i-1)*times+NN[r,1],(i-1)*times+NN[r,2],q]=Au[nu*times*nMax+k*nMax+q];
            #         k=k+1;
            #       }
            #     }
            #   }
            # }
            # A[,,q] <- A[,,q]%*%t(A[,,q])
            if(Astruc %in% c(3,4)){
              A[,,q] <- Alvm%*%t(Alvm)
            } else {
              A[,,q] <- Alvm[,,q]%*%t(Alvm[,,q])
            }
          }
          
          if(Astruc %in% c(3,4)){
            AQ <- matrix(0,num.lv.cor,num.lv.cor)
            # for (d in 1:num.lv.cor){
            #   AQ[d,d]=exp(Au[(nu*times+k+1)]);
            #   k=k+1;
            #   r=d+1
            #   while (r <= num.lv.cor){
            #     AQ[r,d]=Au[(nu*times+k+1)];
            #     k=k+1; r=r+1
            #   }
            # }
            AQ <- objr$report()$AQ
            AQ<-AQ%*%t(AQ)
          }
          
          
          out$AQ <- AQ
          out$A <- A
          
        } else if(nlvr>0){
          param <- objr$env$last.par.best
          A <- array(0, dim=c(n, nlvr, nlvr))
          if(nlvr>(num.lv+num.lv.c)){
            lg_Ar <- param[names(param)=="lg_Ar"]
            for(i in 1:n){
              A[i,1,1]=exp(lg_Ar[i]);
            }
            if(length(lg_Ar)>n){
              for (r in 2:nlvr){
                for(i in 1:n){
                  A[i,r,1]=lg_Ar[((r-1)*n+i)];
                }}
            }
          }
          
          if((num.lv+num.lv.c)>0){
            Au <- param[names(param)=="Au"]
            for (d in 1:(num.lv+num.lv.c)){
              for(i in 1:n){
                A[i,(nlvr-(num.lv+num.lv.c))+ d,(nlvr-(num.lv+num.lv.c))+ d] <- exp(Au[(d-1)*n+i]);
              }
            }
            if(length(Au) > (num.lv+num.lv.c)*n){
              k <- 0;
              for (c1 in 1:(num.lv+num.lv.c)){
                r <- c1 + 1;
                while (r <= (num.lv+num.lv.c)){
                  for(i in 1:n){
                    A[i,(nlvr-(num.lv+num.lv.c))+ r,(nlvr-(num.lv+num.lv.c))+ c1] <- Au[(num.lv+num.lv.c)*n+k*n+i];
                    # A[i,c1,r] <- A[i,r,c1];
                  }
                  k <- k+1; r <- r+1;

                }
              }
            }
            for(i in 1:n){
              A[i,,] <- A[i,,]%*%t(A[i,,])
            }
            out$A <- A
          } else {
            out$Ar <- A
          }
        }
        
        # For random slopes constr. ord.
        if((num.RR+num.lv.c)>0&(randomB!=FALSE)){
          param <- objr$env$last.par.best
          AB_lv <- array(0, dim=c(ab3, ab12, ab12))

          Ab_lv <- param[names(param)=="Ab_lv"]
            for (d in 1:ab12){
              for(i in 1:ab3){
                AB_lv[i,d, d] <- exp(Ab_lv[(d-1)*ab3+i]);
              }
            }
            if(length(Ab_lv) > ab12*ab3){
              k <- 0;
              for (c1 in 1:ab12){
                r <- c1 + 1;
                while (r <= ab12){
                  for(i in 1:ab3){
                    AB_lv[i,r,c1] <- Ab_lv[ab12*ab3+k*ab3+i];
                    # A[i,c1,r] <- A[i,r,c1];
                  }
                  k <- k+1; r <- r+1;
                  
                }
              }
            }
            for(i in 1:ab3){
              AB_lv[i,,] <- AB_lv[i,,]%*%t(AB_lv[i,,])
            }
            out$Ab.lv <- AB_lv
          } 
        
        
        if((num.lv+num.lv.c) == nlvr && row.eff=="random"){
          lg_Ar <- param[names(param)=="lg_Ar"]
          Ar <- exp((lg_Ar)[1:length(out$params$row.params)])
          out$Ar <- Ar
          if(rstruc == 1 && cstrucn>0){
            Arm <- matrix(0,nr,nr)
            diag(Arm)<-Ar
            if(length(lg_Ar)>nr){
              k=0;
              for(d in 1:nr){
                r <- d + 1;
                while (r <= nr){
                  Arm[r,d] = lg_Ar[nr+k];
                  k=k+1; r=r+1;
                }}
            }
            Arm <- Arm %*% t(Arm)
            out$Ar <- diag(Arm)
          }
          if(rstruc == 2){
            
            Arm <- array(0, dim = c(times,times,nr));
            for(i in 1:nr){
              for(d in 1:times){
                Arm[d,d,i]=Ar[(i-1)*times+d];
              }
            }
            if(length(lg_Ar)>(nr*times)){
              k=0;
              for(d in 1:times){
                r <- d + 1;
                while (r <= times){
                  for(i in 1:nr){
                    Arm[r,d,i]=lg_Ar[nr*times+k*nr+i];
                  }
                  k=k+1; r=r+1;
                }}
            }
            for (i in 1:nr) {
              Arm[,,i] <- Arm[,,i] %*% t(Arm[,,i])
            }
            out$Ar <- c(apply(Arm,3,diag))
          }
        }
        
      }
      seed.best <- seed[n.i]
    }
    
    n.i <- n.i+1;
  
    }
  
  #Store the seed that gave the best results, so that we may reproduce results, even if a seed was not explicitly provided
  out$seed <- seed.best
  
  if(is.null(formula1)){ out$formula <- formula} else {out$formula <- formula1}
  
  
  # DW, 7/5/19: adding TMBfn to output:
  out$TMBfn <- objrFinal
  out$TMBfn$par <- optrFinal$par #ensure params in this fn take final values
  out$convergence <- optrFinal$convergence == 0
  out$logL <- -out$logL
  
  # if((method %in% c("VA", "EVA"))){ # These have been moved to gllvm.cpp
  #   if(num.lv > 0) out$logL = out$logL + n*0.5*num.lv
  #   if(row.eff == "random") out$logL = out$logL + n*0.5
  #   if(family=="gaussian") {
  #     out$logL <- out$logL - n*p*log(pi)/2
  #   }
  # }
  
#### Try to calculate sd errors
  
  tr<-try({
    if(sd.errors && !is.infinite(out$logL)) {
      if(trace) cat("Calculating standard errors for parameters...\n")
      out$TMB <- TRUE
      # out <- c(out, se.gllvm(out))
      
      if(family=="ZIP") {
        p0i <- names(pars)=="lg_phi"
        p0 <- pars[p0i][disp.group]
        p0 <- p0+runif(length(unique(disp.group)),0,0.001)[disp.group]
        pars[p0i] <- p0
      }
      if((method %in% c("VA", "EVA"))){
        sdr <- objrFinal$he(pars)
      }
      if(method == "LA"){
        sdr <- optimHess(pars, objrFinal$fn, objrFinal$gr, control = list(reltol=reltol,maxit=maxit))#maxit=maxit
      }
      # makes a small correction to the partial derivatives of LvXcoef if fixed-effect
      # because of constrained objective function
      # assumes L(x) = f(x) + lambda*c(x) for constraint function c(x)
      # though this is not (exactly) how we are fitting the model.
      if((num.RR+num.lv.c)>1 && randomB==FALSE){
        b_lvHE <- sdr[names(pars)=="b_lv",names(pars)=="b_lv"]
        Lmult <- lambda(pars,objrFinal) #estimates  lagranian multiplier
        sdr[names(pars)=="b_lv",names(pars)=="b_lv"] = b_lvHE + b_lvHEcorrect(Lmult,K = ncol(lv.X), d = num.lv.c+num.RR)
      }
      
      m <- dim(sdr)[1]; incl <- rep(TRUE,m); incld <- rep(FALSE,m); inclr <- rep(FALSE,m)
      # Not used for this model
      incl[names(objrFinal$par)=="B"] <- FALSE
      incl[names(objrFinal$par)%in%c("Br","sigmaB","sigmaij")] <- FALSE

      # Variational params not included for incl
      incl[names(objrFinal$par)=="Ab_lv"] <- FALSE;
      incl[names(objrFinal$par)=="Abb"] <- FALSE;
      incl[names(objrFinal$par)=="lg_Ar"] <- FALSE;
      incl[names(objrFinal$par)=="Au"] <- FALSE;
      incl[names(objrFinal$par)=="u"] <- FALSE;
      
      #slopes for reduced rank predictors
      if((num.lv.c+num.RR)==0){
        incl[names(objrFinal$par)=="sigmab_lv"] <- FALSE
        incl[names(objrFinal$par)=="b_lv"] <- FALSE
      }else if((num.lv.c+num.RR)>0&randomB==FALSE){
        incl[names(objrFinal$par)=="b_lv"] <- TRUE
        incl[names(objrFinal$par)=="sigmab_lv"] <- FALSE
      }else if((num.lv.c+num.RR>0)&randomB!=FALSE){
        incl[names(objrFinal$par)=="b_lv"] <- FALSE
        
        inclr[names(objrFinal$par)=="b_lv"] <- TRUE
        incld[names(objrFinal$par)=="b_lv"] <- TRUE
        
        incld[names(objrFinal$par)=="Ab_lv"] <- TRUE
        
        incl[names(objrFinal$par)=="sigmab_lv"] <- TRUE
      }
      
      #loadings for quadratic models
      if(quadratic == FALSE){incl[names(objrFinal$par)=="lambda2"]<-FALSE}
      
      if(familyn!=7) incl[names(objrFinal$par)=="zeta"] <- FALSE
      if(familyn==0 || familyn==2 || familyn==7 || familyn==8) incl[names(objrFinal$par)=="lg_phi"] <- FALSE

      if((num.lv+num.lv.c)>0){
        inclr[names(objrFinal$par)=="u"] <- TRUE;
        incld[names(objrFinal$par)=="u"] <- TRUE;
        incld[names(objrFinal$par)=="Au"] <- TRUE;
      } else {
        #only set lambda to FALSE if all dimensions are zero
        if(num.RR==0)incl[names(objrFinal$par)=="lambda"] <- FALSE;
        #set sigmaLV to zero if no latent variables
        incl[names(objrFinal$par)=="sigmaLV"] <- FALSE;
        
      }
      
      if(row.eff=="random") {
        incld[names(objrFinal$par) == "lg_Ar"] <- TRUE
        incld[names(objrFinal$par) == "r0"] <- TRUE
        inclr[names(objrFinal$par) == "r0"] <- TRUE;
        incl[names(objrFinal$par) == "r0"] <- FALSE; 
      } else {
        incl[names(objrFinal$par)=="log_sigma"] <- FALSE
        if(row.eff==FALSE) { incl[names(objrFinal$par)=="r0"] <- FALSE }
        if(row.eff=="fixed"){ incl[1] <- FALSE }
      }
      
      if(method=="LA" || ((num.lv+num.lv.c)==0 && (method %in% c("VA", "EVA")) && row.eff!="random" && randomB == FALSE)){
        covM <- try(MASS::ginv(sdr[incl,incl]))
        se <- try(sqrt(diag(abs(covM))))
        trpred<-try({
        if(((num.lv+num.lv.c+num.RR) > 0 || row.eff == "random")&& method == "LA") {
          sd.random <- sdrandom(objrFinal, covM, incl,ignore.u = ignore.u)
          prediction.errors <- list()
          
          if(row.eff=="random"){
            prediction.errors$row.params <- sd.random$row
          }
          if((num.lv+num.lv.c+num.RR)>0){
            # cov.lvs <- array(0, dim = c(n, num.lv+num.lv.c, num.lv+num.lv.c))
            # for (i in 1:n) {
            #   cov.lvs[i,,] <- as.matrix(sd.random[(0:(num.lv+num.lv.c-1)*n+i),(0:(num.lv+num.lv.c-1)*n+i)])
              # cov.lvs[i,,] <- as.matrix(sd.random[(0:(nlvr-1)*n+i),(0:(nlvr-1)*n+i)])
            #}
            # if(row.eff=="random"){
            #   prediction.errors$row.params <- cov.lvs[,1,1]
            #   if(num.lv > 0) cov.lvs <- array(cov.lvs[,-1,-1], dim = c(n, num.lv, num.lv))
            # }
            prediction.errors$lvs <- sd.random$A
            if(randomB!=FALSE){
              prediction.errors$Ab.lv <- sd.random$Ab_lv
            }
          }
          
          out$prediction.errors <- prediction.errors
        }
        }, silent = TRUE)
        if(inherits(trpred, "try-error")) { cat("Prediction errors for latent variables could not be calculated.\n") }
        out$Hess <- list(Hess.full=sdr, incla = NULL, incl=incl, incld=NULL, cov.mat.mod=covM)
        
      } else {

        A.mat <- sdr[incl, incl] # a x a
        D.mat <- sdr[incld, incld] # d x d
        B.mat <- sdr[incl, incld] # a x d
        
        try({
          cov.mat.mod<- MASS::ginv(A.mat-B.mat%*%solve(D.mat)%*%t(B.mat))
          se <- sqrt(diag(abs(cov.mat.mod)))
          
          incla<-rep(FALSE, length(incl))
          incla[names(objrFinal$par)=="u"] <- TRUE
          
          out$Hess <- list(Hess.full=sdr, incla = incla, incl=incl, incld=incld, cov.mat.mod=cov.mat.mod)
        }, silent = TRUE)

      }

      if(row.eff == "fixed") {
        se.row.params <- c(0,se[1:(n-1)]); 
        names(se.row.params) <- rownames(out$y); se <- se[-(1:(n-1))] 
        }
      sebetaM <- matrix(se[1:((num.X+1)*p)],p,num.X+1,byrow=TRUE);  
      se <- se[-(1:((num.X+1)*p))]
      
      if(family %in% "betaH"){
        out$sd$betaH <- matrix(se[1:((num.X+1)*p)],p,num.X+1,byrow=TRUE);  se <- se[-(1:((num.X+1)*p))]
        rownames(out$sd$betaH) <- rownames(out$params$betaH); 
        colnames(out$sd$betaH) <- colnames(out$params$betaH)
      }
      

      if((num.lv.c+num.RR)>0 & (randomB==FALSE)){
        se.LvXcoef <- matrix(se[1:((num.lv.c+num.RR)*ncol(lv.X))],ncol=(num.lv.c+num.RR),nrow=ncol(lv.X))
        se <- se[-c(1:((num.lv.c+num.RR)*ncol(lv.X)))]
        colnames(se.LvXcoef) <- paste("CLV",1:(num.lv.c+num.RR),sep="")
        row.names(se.LvXcoef) <- colnames(lv.X)
        out$sd$LvXcoef <- se.LvXcoef
      }
      
      if((num.lv.c+num.lv)>0){
        se.sigma.lv <- se[1:(num.lv+num.lv.c)]; ##*out$params$sigma.lv; 
        se<-se[-c(1:(num.lv+num.lv.c))]
      }
      
      if(num.lv > 0&(num.lv.c+num.RR)==0) {
        se.lambdas <- matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag=FALSE)] <- se[1:(p * num.lv - sum(0:num.lv))];
        colnames(se.lambdas) <- paste("LV", 1:num.lv, sep="");
        rownames(se.lambdas) <- colnames(out$y)
        out$sd$theta <- se.lambdas; 
        if((p * num.lv - sum(0:num.lv))>0) se <- se[-(1:(p * num.lv - sum(0:num.lv)))];
        if(quadratic==TRUE){
          se.lambdas2 <- matrix(se[1:(p * num.lv)], p, num.lv, byrow = T)  
          colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
          se <- se[-(1:(num.lv*p))]
          out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
        }else if(quadratic=="LV"){
          se.lambdas2 <- matrix(se[1:num.lv], p, num.lv, byrow = T)
          colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
          se <- se[-(1:num.lv)]
          out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
        }
      }else if(num.lv==0&(num.lv.c+num.RR)>0){
        se.lambdas <- matrix(0,p,(num.lv.c+num.RR)); se.lambdas[lower.tri(se.lambdas, diag=FALSE)] <- se[1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR)))];

        colnames(se.lambdas) <- paste("CLV", 1:(num.lv.c+num.RR), sep="");
        rownames(se.lambdas) <- colnames(out$y)
        out$sd$theta <- se.lambdas; se <- se[-(1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR))))];
        
        if(quadratic==TRUE){
          se.lambdas2 <- matrix(se[1:(p * (num.lv.c+num.RR))], p, (num.lv.c+num.RR), byrow = T)  
          colnames(se.lambdas2) <- paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = "")
          se <- se[-(1:((num.lv.c+num.RR)*p))]
          out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
        }else if(quadratic=="LV"){
          se.lambdas2 <- matrix(se[1:(num.lv.c+num.RR)], p, (num.lv.c+num.RR), byrow = T)
          colnames(se.lambdas2) <- paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = "")
          se <- se[-(1:(num.lv.c+num.RR))]
          out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
        }
        
      }else if(num.lv>0&(num.lv.c+num.RR)>0){
        se.lambdas <- matrix(0,p,num.lv+(num.lv.c+num.RR));
        se.lambdas[,1:(num.lv.c+num.RR)][lower.tri(se.lambdas[,1:(num.lv.c+num.RR),drop=F], diag=FALSE)] <- se[1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR)))];
        se <- se[-c(1:(p * (num.lv.c+num.RR) - sum(0:(num.lv.c+num.RR))))];
        se.lambdas[,((num.lv.c+num.RR)+1):ncol(se.lambdas)][lower.tri(se.lambdas[,((num.lv.c+num.RR)+1):ncol(se.lambdas),drop=F], diag=FALSE)] <- se[1:(p * num.lv - sum(0:num.lv))];
        if((p * num.lv - sum(0:num.lv))>0) se <- se[-c(1:(p * num.lv - sum(0:num.lv)))]
        colnames(se.lambdas) <- c(paste("CLV", 1:(num.lv.c+num.RR), sep=""),paste("LV", 1:num.lv, sep=""));
        rownames(se.lambdas) <- colnames(out$y)
        out$sd$theta <- se.lambdas;
        
        if(quadratic==TRUE){
          se.lambdas2 <- matrix(se[1:(p * ((num.lv.c+num.RR)+num.lv))], p, (num.lv.c+num.RR)+num.lv, byrow = T)  
          colnames(se.lambdas2) <- c(paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = ""),paste("LV", 1:num.lv, "^2", sep = ""))
          se <- se[-(1:((num.lv+(num.lv.c+num.RR))*p))]
          out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
        }else if(quadratic=="LV"){
          se.lambdas2 <- matrix(se[1:((num.lv.c+num.RR)+num.lv)], p, (num.lv.c+num.RR)+num.lv, byrow = T)
          colnames(se.lambdas2) <- c(paste("CLV", 1:(num.lv.c+num.RR), "^2", sep = ""),paste("LV", 1:num.lv, "^2", sep = ""))
          se <- se[-(1:((num.lv.c+num.RR)+num.lv))]
          out$sd$theta <- cbind(out$sd$theta,se.lambdas2)
        }
        
      }
      
      if(num.lv>0 & family == "betaH"){
        se.thetaH <- matrix(0,num.lv,p); 
        se.thetaH[upper.tri(se.thetaH, diag=TRUE)] <- se[1:(sum(upper.tri(se.thetaH, diag=TRUE)))];
        se <- se[-( 1:sum(upper.tri(se.thetaH, diag=TRUE)) )];
        out$sd$se.thetaH <- se.thetaH
      }
      
      if((num.lv+num.lv.c)>0){
        out$sd$sigma.lv  <- se.sigma.lv
        names(out$sd$sigma.lv) <- colnames(out$params$theta[,1:(num.lv+(num.lv.c))])
      }

      out$sd$beta0 <- sebetaM[,1]; names(out$sd$beta0) <- colnames(out$y);
      if(!is.null(X)){
        out$sd$Xcoef <- matrix(sebetaM[,-1],nrow = nrow(sebetaM));
        rownames(out$sd$Xcoef) <- colnames(y); colnames(out$sd$Xcoef) <- colnames(X);
      }
      if(row.eff=="fixed") {out$sd$row.params <- se.row.params}

      if(family %in% c("negative.binomial")) {
        se.lphis <- se[1:length(unique(disp.group))][disp.group];  out$sd$inv.phi <- se.lphis*out$params$inv.phi;
        out$sd$phi <- se.lphis*out$params$phi;
        names(out$sd$phi) <- colnames(y);
        
        if(!is.null(names(disp.group))){
          try(names(out$sd$phi) <- names(disp.group),silent=T)
        }
        names(out$sd$inv.phi) <-  names(out$sd$phi)
        se <- se[-(1:length(unique(disp.group)))]
      }
      if(family %in% c("tweedie", "gaussian", "gamma", "beta", "betaH")) {
        se.lphis <- se[1:length(unique(disp.group))][disp.group];
        out$sd$phi <- se.lphis*out$params$phi;
        names(out$sd$phi) <- colnames(y);
         if(!is.null(names(disp.group))){
          try(names(out$sd$phi) <- names(disp.group),silent=T)
        }
        
        se <- se[-(1:length(unique(disp.group)))]
      }
      if(family %in% c("ZIP")) {
        se.phis <- se[1:length(unique(disp.group))][disp.group];
        out$sd$phi <- se.phis*exp(lp0)/(1+exp(lp0))^2;#
        names(out$sd$phi) <- colnames(y);
        
        if(!is.null(names(disp.group))){
          try(names(out$sd$phi) <- names(disp.group),silent=T)
        }
        se <- se[-(1:length(unique(disp.group)))]
      }
      if((randomB!=FALSE)&(num.lv.c+num.RR)>0){
        se.lsigmab.lv <-  se[1:length(out$paramssigmaLvXcoef)];
        se <-  se[-c(1:length(out$params$sigmaLvXcoef))]
        out$sd$sigmaLvXcoef <- se.lsigmab.lv*out$params$sigmaLvXcoef
        if(randomB=="LV")names(out$sd$sigmaLvXcoef) <- paste("CLV",1:(num.lv.c+num.RR), sep="")
        if(randomB=="P")names(out$sd$sigmaLvXcoef) <- colnames(lv.X)
        # if(randomB=="all")names(out$sd$sigmaLvXcoef) <- paste(paste("CLV",1:(num.lv.c+num.RR),sep=""),rep(colnames(lv.X),each=num.RR+num.lv.c),sep=".")
        if(randomB=="single")names(out$sd$sigmaLvXcoef) <- NULL
      }
      
      if(row.eff=="random") {
        out$sd$sigma <- se[1:length(out$params$sigma)]*c(out$params$sigma[1],rep(1,length(out$params$sigma)-1));
        se=se[-(1:length(out$sd$sigma))]
        names(out$sd$sigma) <- "sigma"
        if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,3))) {out$sd$rho <- se[1]*(1-out$params$rho^2)^1.5; if(length(out$params$rho)) se = se[-(1:length(out$params$rho))]}#*sqrt(1-out$params$rho^2)
        if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(2,4))) {out$sd$rho <- se[1]*out$params$rho; if(length(out$params$rho)) se = se[-(1:length(out$params$rho))]}
      }
      if(num.lv.cor>0 & cstrucn>0){ 
        if((cstrucn %in% c(1,3))) {out$sd$rho.lv <- se[(1:length(out$params$rho.lv))]*(1-out$params$rho.lv^2)^1.5; se = se[-(1:length(out$params$rho.lv))]}
        if((cstrucn %in% c(2,4))) {out$sd$rho.lv <- se[(1:length(out$params$rho.lv))]*out$params$rho.lv; se = se[-(1:length(out$params$rho.lv))]}
        names(out$sd$rho.lv) <- names(out$params$rho.lv)
        # if((cstrucn %in% c(2,4))) {out$sd$scaledc <- se[(1:length(out$params$scaledc))]*out$params$scaledc; se = se[-(1:length(out$sd$scaledc))]}
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


  return(out)
}

