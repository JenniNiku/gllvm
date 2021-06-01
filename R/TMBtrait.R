########################################################################################
## GLLVM fourth corner model, with estimation done via Laplace and Variational approximation using TMB-package
## Original author: Jenni Niku
##########################################################################################

trait.TMB <- function(
      y, X = NULL,TR=NULL,formula=NULL, num.lv = 2, family = "poisson",
      Lambda.struc = "unstructured", Ab.struct = "unstructured", row.eff = FALSE, reltol = 1e-6, seed = NULL,
      maxit = 3000, max.iter=200, start.lvs = NULL, offset=NULL, sd.errors = FALSE,trace=FALSE,
      link="logit",n.init=1,start.params=NULL,start0=FALSE,optimizer="optim", dr=NULL, rstruc =0, cstruc = "diag", dist =0,
      starting.val="res",method="VA",randomX=NULL,Power=1.5,diag.iter=1, Ab.diag.iter = 0, dependent.row = FALSE,
      Lambda.start=c(0.2, 0.5), jitter.var=0, yXT = NULL, scale.X = FALSE, randomX.start = "zero", beta0com = FALSE
      ,zeta.struc = "species", quad.start=0.01, start.struc="LV",quadratic=FALSE, optim.method = "BFGS") {
  if(is.null(X) && !is.null(TR)) stop("Unable to fit a model that includes only trait covariates")
  if(!is.null(start.params)) starting.val <- "zero"
  
  objrFinal <- optrFinal <- NULL
  cstrucn = switch(cstruc, "diag" = 0, "corAR1" = 1, "corExp" = 2, "corCS" = 3)
  
    term <- NULL
  n <- nr <- dim(y)[1]; p <- dim(y)[2];
  times = 1
  
  # Structure for row effects
  model = 1
  xr = NULL
  if(rstruc==0){ # No structure
    dr <- diag(n)
  }
  if(rstruc>0){#rstruc
    if(is.null(dr)) stop("Define structure for row params if 'rstruc == ",rstruc,"'.")
    if(rstruc==1){# group specific
      nr <- dim(dr)[2]
      if((cstrucn == 2)) {
        if(is.null(dist) || length(dist)!=nr)
          dist=1:nr
      }
    }
    if(rstruc==2) { # correlated within groups
      nr <- dim(dr)[2]
      times <- n/nr#dim(dr)[1]
      if((cstrucn == 2)) {
        if(is.null(dist) || length(dist)!=times)
          dist=1:times
      }
    }
  }
  
  
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
        dum <- model.matrix( ~ X[,i]-1)
        dum <- as.matrix(dum[, !(colnames(dum) %in% c("(Intercept)"))])
        # colnames(dum) <- paste(colnames(X)[i], levels(X[,i])[ - 1], sep = "")
        colnames(dum) <- paste(colnames(X)[i], levels(X[,i]), sep = "")
        X.new <- cbind(X.new, dum)
      }
    }
    X.new <- data.frame(X.new);
  }
  
  num.T <- 0
  T.new <- NULL
  if(!is.null(TR)) {
    num.T <- dim(TR)[2]
    T.new <- matrix(0, p, 0)
    if(num.T > 0){
    for (i in 1 : num.T) {
      #if(!is.factor(TR[,i])  && length(unique(TR[,i])) > 2) { #!!!
      if(is.numeric(TR[,i])  && length(unique(TR[,i])) > 2) {
        TR[,i] <- scale(TR[,i])
        T.new <- cbind(T.new,scale(TR[,i], scale = scale.X, center = scale.X)); colnames(T.new)[dim(T.new)[2]] <- colnames(TR)[i]
      } else {
        if(!is.factor(TR[,i])) TR[,i] <- factor(TR[,i]) #!!!
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

  # Define design matrix for covariates
  if(!is.null(X) || !is.null(TR)){
    yX <- cbind(cbind(X,id = 1:nrow(y))[rep(1:nrow(X), times=ncol(y)),],  species = rep(1:ncol(y), each= nrow(y)), y = c(as.matrix(y))) #reshape(data.frame(cbind(y, X)), direction = "long", varying = colnames(y), v.names = "y")
    TR2 <- data.frame(species = 1:p, TR)
    if(is.null(yXT)){
      yXT <- merge(yX, TR2, by = "species")
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
    if(length(nxd)>1) for(i in 2:length(nxd)) formulab <- paste(formulab,nxd[i],sep = "+")
    formula1 <- formulab
  }


  if(!(family %in% c("poisson","negative.binomial","binomial","tweedie","ZIP", "gaussian", "ordinal", "gamma", "exponential")))
    stop("Selected family not permitted...sorry!")
  if(!(Lambda.struc %in% c("unstructured","diagonal")))
    stop("Lambda matrix (covariance of variational distribution for latent variable) not permitted...sorry!")
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

  out <-  list(y = y, X = X1, TR = TR1, num.lv = num.lv, row.eff = row.eff, logL = Inf, family = family, offset=offset,randomX=randomX,X.design=Xd,terms=term, method = method)
  if(is.null(formula) && is.null(X) && is.null(TR)){formula ="~ 1"}
  
  n.i <- 1;
  if(n.init > 1) seed <- sample(1:10000, n.init)

  # n.init model fits
  while(n.i <= n.init){

    randomXb <- NULL
    
    # Design for random slopes
    if(!is.null(randomX)){
      #
      if(num.lv>0 && randomX.start == "res" && starting.val == "res") {randomXb <- randomX}
      #
      xb <- as.matrix(model.matrix(randomX, data = data.frame(X)))
      rnam <- colnames(xb)[!(colnames(xb) %in% c("(Intercept)"))]
      xb <- as.matrix(xb[, rnam]); #as.matrix(X.new[, rnam])
      if(NCOL(xb) == 1) colnames(xb) <- rnam
      bstart <- start.values.randomX(y, xb, family, starting.val = randomX.start, power = Power, link = link)
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
    
    num.X <- dim(X)[2]
    num.T <- dim(TR)[2]
    phi<-phis <- NULL
    sigma <- 1
    phi <- phis <- NULL;
    
    if(n.init > 1 && trace) cat("initial run ",n.i,"\n");

    #### Calculate starting values

    res <- start.values.gllvm.TMB(y = y, X = X1, TR = TR1, family = family, offset=offset, trial.size = trial.size, num.lv = num.lv, start.lvs = start.lvs, seed = seed[n.i],starting.val=starting.val,power=Power,formula = formula, jitter.var=jitter.var, #!!!
                                  yXT=yXT, row.eff = row.eff, TMB=TRUE, link=link, randomX=randomXb, beta0com = beta0com0, zeta.struc = zeta.struc)

    ## Set initial values
    
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
        if(rstruc<2 && row.eff=="random") row.params <- row.params[1:nr] #rstruc
        if (row.eff == "random") {
          sigma <- sd(row.params);
        }
      }
      vameans <- theta <- lambda <- NULL

      if(num.lv > 0) {
        if(!is.null(randomXb) && family != "ordinal"){
          Br <- res$Br
          sigmaB <- (res$sigmaB)
          if(length(sigmaB)>1) sigmaij <- rep(0,length(res$sigmaij))
          if(randomX.start == "res" && !is.null(res$fitstart)) { ##!!!
            res$sigmaij <- sigmaij <- res$fitstart$TMBfnpar[names(res$fitstart$TMBfnpar) == "sigmaij"] 
          }
        }
        if(start.struc=="LV"&quadratic!=FALSE){
          lambda2 <- matrix(quad.start, ncol = num.lv, nrow = 1)  
        }else if(start.struc=="all"&quadratic!=FALSE){
          lambda2 <- matrix(quad.start, ncol = num.lv, nrow = p)
        }else if(quadratic==FALSE){
          lambda2 <- 0
        }
        if(quadratic != FALSE){
          res$params <- cbind(res$params, matrix(lambda2,nrow=p,ncol=num.lv))  
        }else{
          res$params <- res$params
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
            lastart <- FAstart(res$mu, family=family, y=y, num.lv = num.lv, phis = res$phi, jitter.var = jitter.var[1], zeta.struc=zeta.struc, zeta = res$zeta)
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
          if(class(start.params)[2]=="gllvm.quadratic" && quadratic != FALSE){
            lambda2 <- start.params$params$theta[,-c(1:start.params$num.lv),drop=F]
          }else if(class(start.params)[1]=="gllvm" && quadratic != FALSE){
            if(start.struc=="LV"|quadratic=="LV"){
              lambda2 <- matrix(quad.start, ncol = num.lv, nrow = 1)  
            }else if(start.struc=="all"&quadratic=="all"){
              lambda2 <- matrix(quad.start, ncol = num.lv, nrow = p)
            }
          }
        }
        if(family == "negative.binomial" && start.params$family == "negative.binomial" && !is.null(start.params$params$phi)) {res$phi<-start.params$params$phi}
      } else { stop("Model which is set as starting parameters isn't the suitable you are trying to fit. Check that attributes y, X, TR and row.eff match to each other.");}
    }
    if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)

    if(family == "negative.binomial") {
      phis <- res$phi
      if (any(phis > 10))
        phis[phis > 50] <- 50
      if (any(phis < 0.02))
        phis[phis < 0.02] <- 0.02
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
    # if (family %in% c("gaussian", "gamma")) {
    #   phis <- res$phi
    # }
    if(family=="ordinal"){
      K = max(y00)-min(y00)
      if(zeta.struc=="species"){
        zeta <- c(t(res$zeta[,-1]))
        zeta <- zeta[!is.na(zeta)]
      }else{
        zeta <- res$zeta[-1]
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
    if(row.eff == "random" && rstruc ==0){ nlvr<-num.lv+1 } else {nlvr=num.lv}
    if(row.eff=="fixed"){xr <- matrix(1,1,p)} else {xr <- matrix(0,1,p)}

    # set starting values for variational distribution covariances
    # Variational covariances for latent variables
    if(num.lv > 0){
      if(Lambda.struc=="diagonal" || diag.iter>0){
        Au <- log(rep(Lambda.start[1],num.lv*n)) #
      } else{
        Au <- c(log(rep(Lambda.start[1],num.lv*n)),rep(0,num.lv*(num.lv-1)/2*n))
      }
    } else { Au <- 0}

    # Variational covariances for  random rows
    if(row.eff == "random"){
      if(rstruc ==1){
        lg_Ar <- rep(log(Lambda.start[2]), nr)
      } else {
        lg_Ar <- rep(log(Lambda.start[2]), n)
      }
      
      if(rstruc == 0 && nlvr>num.lv && (quadratic == FALSE)){
        lg_Ar<-c(lg_Ar, rep(0, num.lv*n))
      }
      if(rstruc == 1 & (cstrucn %in% c(1,2,3))){
        lg_Ar<-c(lg_Ar, rep(0, nr*(nr-1)/2))
      }
      if(rstruc == 2){
        lg_Ar<-c(lg_Ar, rep(0, nr*times*(times-1)/2))
      }
    } else {lg_Ar <- 0}
    
    # Variational covariances for  random slopes of envs
    if(!is.null(randomX)){
      if(length(Lambda.start)>2) { 
        a.var <- Lambda.start[3];
      } else {a.var <- 0.5;}
      if(randomX.start == "res" && !is.null(res$fitstart$Ab)){ # !!!! && !is.null(res$fitstart$Ab)
        if(Ab.struct == "diagonal" || Ab.diag.iter>0){
          Abb <- c(log(c(apply(res$fitstart$Ab,1, diag))))
        } else {
          Abb <- c(log(c(apply(res$fitstart$Ab,1, diag))), rep(0, ncol(xb) * (ncol(xb) - 1) / 2 * p))
        }
        res$Br <- Br
        res$Ab <- c(apply(res$fitstart$Ab,1, diag))
      } else{ #!!!
        if(Ab.struct == "diagonal" || Ab.diag.iter>0){
          Abb <- c(log(rep(a.var, ncol(xb) * p)))
        } else {
          Abb <- c(log(rep(a.var, ncol(xb) * p)), rep(0, ncol(xb) * (ncol(xb) - 1) / 2 * p))
        }
      } #!!!
    } else { Abb <- 0 }

    optr<-NULL
    timeo<-NULL
    se <- NULL

## map.list defines parameters which are not estimated in this model
    
    map.list <- list()    
    #    if(row.eff==FALSE) map.list$r0 <- factor(rep(NA,n))
    if(family %in% c("poisson","binomial","ordinal","exponential")) map.list$lg_phi <- factor(rep(NA,p))
    if(family != "ordinal") map.list$zeta <- factor(NA)
    
  ### family settings
    
    extra <- c(0,1)
    if(family == "poisson") { familyn=0}
    if(family == "negative.binomial") { familyn=1}
    if(family == "binomial") {
      familyn <- 2;
      if(link=="probit") extra[1] <- 1
    }
    if(family == "gaussian") {familyn=3}
    if(family == "gamma") {familyn=4}
    if(family == "tweedie"){ familyn <- 5; extra[1] <- Power}
    if(family == "ZIP"){ familyn <- 6;}
    if(family == "ordinal") {familyn=7}
    if(family == "exponential") {familyn=8}
    if(family == "beta"){
      familyn=9
      if(link=="probit") extra[1] <- 1
    }
    
  # Common intercept
    if(beta0com){
      extra[2] <- 0
      Xd<-cbind(1,Xd)
      a <- a*0
      B<-c(mean(a),B)
    }

    
    # For Laplace method, specify random paramteters to randomp
    randomp= NULL #c("u","r0,"Br")
    randoml=c(0,0)
    
# Specify parameter.list, data.list and map.list
    # latent vars
    if(num.lv>0){
      u<-cbind(u)
      randomp <- c(randomp,"u")
    } else {
      u<-matrix(0)
      theta = 0; 
      lambda2 <- 0
      map.list$lambda = factor(NA)
      map.list$lambda2 = factor(NA)
      map.list$u = factor(NA) 
      map.list$Au = factor(NA) 
    }
    
    ## Row effect settings
    if(row.eff=="random"){
      randoml[1] <- 1
      randomp <- c(randomp,"r0")
      
      if(dependent.row && (rstruc == 0)) 
        sigma<-c(sigma[1], rep(0, num.lv))
      if((rstruc %in% 1:2) & (cstrucn %in% c(1,2,3))) {
        sigma = c(log(sigma[1]),0)
        }
    } else {
      sigma=0
      map.list$log_sigma <- factor(NA)
      map.list$lg_Ar <- factor(NA)
      # if(row.eff != "fixed") map.list$r0 <- factor(NA, length(r0))
    }
    
    # Random slopes
    if(!is.null(randomX)){
      randoml[2]=1
      randomp <- c(randomp,"Br")
      res$Br <- Br
      res$sigmaB <- sigmaB
    } else {
      map.list$Br = factor(NA)
      map.list$sigmaB = factor(NA)
      map.list$sigmaij = factor(NA)
      map.list$Abb = factor(NA)
    }
    if(quadratic==FALSE){
      map.list$lambda2 <- factor(NA)
    }
    
    
## Something for quadratic model, what?
    if(starting.val!="zero" && start.struc != "LV" && quadratic == TRUE && num.lv>0 && method == "VA"){
      map.list2 <- map.list
      map.list2$r0 = factor(rep(NA, length(r0)))
      map.list2$b = factor(rep(NA, length(rbind(a))))
      map.list2$B = factor(rep(NA, length(B)))
      map.list2$Br = factor(rep(NA,length(Br)))
      map.list2$lambda = factor(rep(NA, length(theta)))
      map.list2$u = factor(rep(NA, length(u)))
      map.list2$lg_phi = factor(rep(NA, length(phi)))
      map.list2$sigmaB = factor(rep(NA,length(sigmaB)))
      map.list2$sigmaij = factor(rep(NA,length(sigmaij)))
      map.list2$log_sigma = factor(rep(NA, length(sigma)))
      map.list2$Au = factor(rep(NA, length(Au)))
      map.list2$lg_Ar = factor(rep(NA, length(lg_Ar)))
      map.list2$Abb = factor(rep(NA, length(Abb)))
      map.list2$zeta = factor(rep(NA, length(zeta)))
      parameter.list = list(r0=matrix(r0), b = rbind(a), B=matrix(B), Br=Br, lambda = theta, lambda2 = t(lambda2), u = u, lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij, log_sigma=c(sigma), Au=Au, lg_Ar=lg_Ar, Abb=Abb, zeta=zeta)
      objr <- TMB::MakeADFun(
        data = list(y = y, x = Xd, xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, family=familyn, extra=extra, quadratic = 1, method=0, model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist), silent=!trace,
        parameters = parameter.list, map = map.list2,
        inner.control=list(mgcmax = 1e+200,maxit = maxit),
        DLL = "gllvm")
      
      if(optimizer=="nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit)),silent = TRUE))
      }
      if(optimizer=="optim") {
        if(optim.method != "BFGS")
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE))
        else
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
      }
      lambda2 <- matrix(optr$par, byrow = T, ncol = num.lv, nrow = p)
      
      if(inherits(optr,"try-error")) warning(optr[1]);
    }
      
      
#### Call makeADFun
    
    if(method == "VA" && (num.lv>0 || row.eff=="random" || !is.null(randomX))){
      parameter.list = list(r0=matrix(r0), b = rbind(a), B=matrix(B), Br=Br, lambda = theta, lambda2 = t(lambda2), u = u, lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij, log_sigma=c(sigma), Au=Au, lg_Ar=lg_Ar, Abb=Abb, zeta=zeta)
      objr <- TMB::MakeADFun(
        data = list(y = y, x = Xd, xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, quadratic = ifelse(quadratic!=FALSE,1,0), family=familyn, extra=extra, method=0, model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist), silent=!trace,
        parameters = parameter.list, map = map.list,
        inner.control=list(mgcmax = 1e+200,maxit = maxit),
        DLL = "gllvm")
    } else {
      Au=0; Abb=0; lg_Ar=0;
      map.list$Au <- map.list$Abb <- map.list$lg_Ar <- factor(NA)
      parameter.list = list(r0=matrix(r0), b = rbind(a), B=matrix(B), Br=Br, lambda = theta, lambda2 = t(lambda2), u = u, lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij, log_sigma=c(sigma), Au=Au, lg_Ar=lg_Ar, Abb=Abb, zeta=zeta)
      data.list <- list(y = y, x = Xd,xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, quadratic = 0, family=familyn,extra=extra,method=1,model=1,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist)
      if(family == "ordinal"){
        data.list$method = 0
      }
      objr <- TMB::MakeADFun(
        data = data.list, silent=!trace,
        parameters = parameter.list, map = map.list,
        inner.control=list(mgcmax = 1e+200,maxit = maxit,tol10=0.01),
        random = randomp, DLL = "gllvm")
    }

#### Fit model 
    
    if(optimizer=="nlminb") {
      timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit)),silent = TRUE))
    }
    if(optimizer=="optim") {
      if(optim.method != "BFGS") # Due the memory issues, "BFGS" should not be used for Tweedie
        timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE))
      else
        timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
    }
    if(inherits(optr,"try-error")) warning(optr[1]);

 ### Now diag.iter, improves the model fit sometimes
    
    if(diag.iter>0 && Lambda.struc=="unstructured" && method =="VA" && (num.lv>1 || !is.null(randomX)) && !inherits(optr,"try-error")){
      objr1 <- objr
      optr1 <- optr
      param1 <- optr$par
      nam <- names(param1)
      r1 <- matrix(param1[nam=="r0"])
      b1 <- rbind(param1[nam=="b"])
      B1 <- matrix(param1[nam=="B"])
      if(!is.null(randomX)) {
        Br1 <- matrix(param1[nam=="Br"], ncol(xb), p) #!!!
        sigmaB1 <- param1[nam=="sigmaB"]
        sigmaij1 <- param1[nam=="sigmaij"]*0
        Abb <- param1[nam=="Abb"]
        if(Ab.diag.iter>0 && Ab.struct == "unstructured")
          Abb <- c(Abb, rep(0,ncol(xb)*(ncol(xb)-1)/2*p))
      } else {
        Br1 <- Br
        sigmaB1 <- sigmaB
        sigmaij1 <- sigmaij
      }
      if(num.lv>0) {
        lambda1 <- param1[nam=="lambda"]; 
        u1 <- matrix(param1[nam=="u"],n,num.lv) 
        Au<- c(pmax(param1[nam=="Au"],rep(log(1e-6), num.lv*n)), rep(0,num.lv*(num.lv-1)/2*n))
        
        if (quadratic=="LV" | quadratic == T && start.struc == "LV"){
          lambda2 <- matrix(param1[nam == "lambda2"], byrow = T, ncol = num.lv, nrow = 1)#In this scenario we have estimated two quadratic coefficients before
        }else if(quadratic == T){
          lambda2 <- matrix(param1[nam == "lambda2"], byrow = T, ncol = num.lv, nrow = p)
        }
        
        
      } else {u1 <- u}
      if(num.lv==0) {lambda1 <- 0; }
      if(family %in% c("poisson","binomial","ordinal","exponential")){ lg_phi1 <- log(phi)} else {lg_phi1 <- param1[nam=="lg_phi"]}
      if(row.eff == "random"){
        lg_sigma1 <- param1[nam=="log_sigma"]
        lg_Ar<- param1[nam=="lg_Ar"]
      } else {lg_sigma1 = 0}

      if(family == "ordinal"){ zeta <- param1[nam=="zeta"] } else { zeta <- 0 }
      

      parameter.list = list(r0=r1, b = b1, B=B1, Br=Br1, lambda = lambda1, lambda2 = t(lambda2), u = u1, lg_phi=lg_phi1, sigmaB=sigmaB1, sigmaij=sigmaij1, log_sigma=lg_sigma1, Au=Au, lg_Ar=lg_Ar, Abb=Abb, zeta=zeta)

      data.list = list(y = y, x = Xd, xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, quadratic = ifelse(quadratic!=FALSE&num.lv>0,1,0), family=familyn, extra=extra, method=0, model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist)
      objr <- TMB::MakeADFun(
        data = data.list, silent=!trace,
        parameters = parameter.list, map = map.list,
        inner.control=list(mgcmax = 1e+200,maxit = 1000),
        DLL = "gllvm")

      if(optimizer=="nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit)),silent = TRUE))
      }
      if(optimizer=="optim") {
        if(optim.method != "BFGS")
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE))
        else
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
      }
      if(inherits(optr, "try-error")){optr <- optr1; objr <- objr1; Lambda.struc <- "diagonal"}

    }

    if(!inherits(optr,"try-error") && start.struc=="LV" && quadratic == TRUE && method == "VA"){
      objr1 <- objr
      optr1 <- optr
      param1 <- optr$par
      nam <- names(param1)
      r1 <- matrix(param1[nam=="r0"])
      b1 <- rbind(param1[nam=="b"])
      B1 <- matrix(param1[nam=="B"])
      if(!is.null(randomX)) {
        Br1 <- matrix(param1[nam=="Br"], ncol(xb), p) #!!!
        sigmaB1 <- param1[nam=="sigmaB"]
        sigmaij1 <- param1[nam=="sigmaij"]*0
        Abb <- param1[nam=="Abb"]
        if(Ab.diag.iter>0 && Ab.struct == "unstructured")
          Abb <- c(Abb, rep(0,ncol(xb)*(ncol(xb)-1)/2*p))
      } else {
        Br1 <- Br
        sigmaB1 <- sigmaB
        sigmaij1 <- sigmaij
      }
        lambda1 <- param1[nam=="lambda"]; 
        u1 <- matrix(param1[nam=="u"],n,num.lv) 
        Au<- param1[nam=="Au"]
        lambda2 <- abs(matrix(param1[nam == "lambda2"], byrow = T, ncol = num.lv, nrow = p))

      if(family %in% c("poisson","binomial","ordinal","exponential")){ lg_phi1 <- log(phi)} else {lg_phi1 <- param1[nam=="lg_phi"]}
      if(row.eff == "random"){
        lg_sigma1 <- param1[nam=="log_sigma"]
        lg_Ar<- param1[nam=="lg_Ar"]
      } else {lg_sigma1 = 0}

      if(family == "ordinal"){ zeta <- param1[nam=="zeta"] } else { zeta <- 0 }
      
      parameter.list = list(r0=r1, b = b1, B=B1, Br=Br1, lambda = lambda1, lambda2 = t(lambda2), u = u1, lg_phi=lg_phi1, sigmaB=sigmaB1, sigmaij=sigmaij1, log_sigma=lg_sigma1, Au=Au, lg_Ar=lg_Ar, Abb=Abb, zeta=zeta)

      data.list = list(y = y, x = Xd, xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, quadratic = 1, family=familyn, extra=extra, method=0, model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist)
      objr <- TMB::MakeADFun(
        data = data.list, silent=!trace,
        parameters = parameter.list, map = map.list,
        inner.control=list(mgcmax = 1e+200,maxit = 1000),
        DLL = "gllvm")
      
      if(optimizer=="nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit)),silent = TRUE))
      }
      if(optimizer=="optim") {
        if(optim.method != "BFGS")
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method, control = list(maxit=maxit),hessian = FALSE),silent = TRUE))
        else
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS", control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE))
      }
      
      #quick check to see if something actually happened
      flag <- 1
      if(all(round(lambda2,0)==round(matrix(abs(optr$par[names(optr$par)=="lambda2"]),byrow=T,ncol=num.lv,nrow=p),0))){
        flag <- 0
        warning("Full quadratic model did not properly converge or all quadratic coefficients are close to zero. Try changing 'start.struc' in 'control.start'. /n")
      }
      if(inherits(optr, "try-error") || flag == 0){optr <- optr1; objr <- objr1; quadratic <- "LV";}
      
    }

        
#### Extract estimated values
    
    param <- objr$env$last.par.best
    if(family %in% c("negative.binomial", "tweedie", "gaussian", "gamma", "beta")) {
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
    li2 <- names(param)=="lambda2"
    ui<-names(param)=="u"
    if(num.lv > 0){
      lvs <- (matrix(param[ui],n,num.lv))
      theta <- matrix(0,p,num.lv)
      if(p>1) {
        theta[lower.tri(theta,diag=TRUE)] <- param[li];
        if(quadratic!=FALSE){
          theta<-cbind(theta,matrix(-abs(param[li2]),ncol=num.lv,nrow=p,byrow=T))
        }
      } else {theta <- c(param[li],-abs(param[li2]))}
      # diag(theta) <- exp(diag(theta))#!!!
    }
    if(row.eff!=FALSE) {
      ri <- names(param)=="r0"
      row.params=param[ri] 
      if(row.eff=="random"){
        sigma<-exp(param[names(param)=="log_sigma"])[1]
        if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,3))) rho<- param[names(param)=="log_sigma"][2] / sqrt(1.0 + param[names(param)=="log_sigma"][2]^2);
        if((rstruc ==2 | (rstruc == 1)) & (cstrucn == 2)) rho<- exp(param[names(param)=="log_sigma"][2]);
        if(num.lv>0 && dependent.row && rstruc==0) sigma <- c(sigma, (param[names(param)=="log_sigma"])[-1])
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

    
#### Check if model fit succeeded/improved on this iteration n.i
    
    if( (n.i==1 || out$logL > (new.loglik)) && is.finite(new.loglik) && !inherits(optr, "try-error")){ #
      objrFinal<-objr1 <- objr; optrFinal<-optr1 <- optr;
      out$logL <- new.loglik
      if(num.lv > 0) {
        out$lvs <- lvs
        out$params$theta <- theta
        rownames(out$lvs) <- rownames(out$y);
        rownames(out$params$theta) <- colnames(out$y)
        if(quadratic==FALSE)colnames(out$params$theta) <- colnames(out$lvs) <- paste("LV", 1:num.lv, sep="");
        if(quadratic!=FALSE){
          colnames(out$lvs) <- paste("LV", 1:num.lv, sep="");
          colnames(out$params$theta)<- c(paste("LV", 1:num.lv, sep=""),paste("LV", 1:num.lv, "^2",sep=""));
        }
      }
      if(!beta0com) names(beta0) <- colnames(out$y); 
      if(beta0com) names(beta0) <- "Community intercept";
      out$params$beta0 <- beta0;
      out$params$B <- B; names(out$params$B)=colnames(Xd)

      if(row.eff!=FALSE) {
        if(row.eff=="random"){  
          out$params$sigma <- sigma; 
          names(out$params$sigma) <- "sigma"
          if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,2,3))) out$params$rho <- rho
          if(num.lv>0 && dependent.row) names(out$params$sigma) <- paste("sigma",c("",1:num.lv), sep = "")
          }
        out$params$row.params <- row.params; 
        if(length(row.params) == n) names(out$params$row.params) <- rownames(out$y)
        if((length(row.params) == ncol(dr)) && (rstruc==1)) try(names(out$params$row.params) <- colnames(dr), silent = TRUE)
      }
      if(family %in% c("negative.binomial")) {
        out$params$phi <- 1/phis; names(out$params$phi) <- colnames(out$y);
        out$params$inv.phi <- phis; names(out$params$inv.phi) <- colnames(out$y);
      }
      if(family %in% c("gaussian","tweedie","gamma","beta")) {
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
        out$corr <- sigmaB_ #!!!!
        rownames(out$params$Br) <- rownames(out$params$sigmaB) <- colnames(out$params$sigmaB) <- colnames(xb)
      }
      if(family %in% c("binomial", "beta")) out$link <- link;
      out$row.eff <- row.eff
      out$time <- timeo
      out$start <- res
      out$Power <- Power
      pars <- optr$par

## Collect VA covariances
      if(method=="VA"){
        param <- objr$env$last.par.best
        
        if(nlvr>0){
          param <- objr$env$last.par.best
          A <- array(0, dim=c(n, nlvr, nlvr))
          if(nlvr>num.lv){
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
          
          if(num.lv>0){
            Au <- param[names(param)=="Au"]
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
        
        if(num.lv == nlvr && row.eff=="random"){
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

      if(method == "VA" && !is.null(randomX)){
        Abb <- param[names(param) == "Abb"]
        xdr <- ncol(xb)
        Ab <- array(0,dim=c(p,xdr,xdr))
        for (d in 1:xdr){
          for(j in 1:p){
            Ab[j,d,d] <- exp(Abb[(d-1)*p + j]);
          }
        }
        if(length(Abb)>xdr*p){
          k <- 0;
          for (c1 in 1:xdr){
            r <- c1+1;
            while (r <= xdr){
              for(j in 1:p){
                Ab[j,r,c1] <- Abb[xdr*p+k*p+j];
                # Ab[j,c1,r] <- Ab[j,r,c1];
              }
              k <- k+1; r <- r+1;
            }
          }
        }
        for(j in 1:p){
          Ab[j,,] <- Ab[j,,]%*%t(Ab[j,,])
        }
        out$Ab <- Ab
      }
    }

    n.i <- n.i+1;
  }

  
  
  if(is.null(formula1)){ out$formula <- formula} else {out$formula <- formula1}
  
  out$Xrandom <- xb
  out$D <- Xd
  out$TMBfn <- objrFinal
  out$TMBfn$par <- optrFinal$par #ensure params in this fn take final values
  out$convergence <- optrFinal$convergence == 0
  out$quadratic <- quadratic
  out$logL <- -out$logL
  out$zeta.struc <- zeta.struc
  out$beta0com <- beta0com
  
  # if(method == "VA"){ # These have been moved to gllvm.cpp
  #   if(num.lv > 0) out$logL = out$logL + n*0.5*num.lv
  #   if(row.eff == "random") out$logL = out$logL + n*0.5
  #   if(!is.null(randomX)) out$logL = out$logL + p*0.5*ncol(xb)
  #   if(family=="gaussian") {
  #     out$logL <- out$logL - n*p*log(pi)/2
  #   }
  # }
  
#### Try to calculate sd errors
  
    tr<-try({
      if(sd.errors && is.finite(out$logL)) {
        if(trace) cat("Calculating standard errors for parameters...\n")
        out$TMB <- TRUE
        # out <- c(out, se.gllvm(out))

        if(method == "VA"){
          sdr <- objrFinal$he(optrFinal$par)
        }
        if(method == "LA"){
          sdr <- optimHess(optrFinal$par, objrFinal$fn, objrFinal$gr,control = list(reltol=reltol,maxit=maxit))
        }

        m <- dim(sdr)[1]; incl <- rep(TRUE,m); incld <- rep(FALSE,m)
        
        # Variational params not included for incl
        incl[names(objrFinal$par)=="Abb"] <- FALSE;
        incl[names(objrFinal$par)=="lg_Ar"] <- FALSE;
        incl[names(objrFinal$par)=="Au"] <- FALSE;
        incl[names(objrFinal$par)=="u"] <- FALSE; 
        
        if(quadratic == FALSE){incl[names(objrFinal$par)=="lambda2"]<-FALSE}
        if(beta0com){incl[names(objrFinal$par)=="b"] <- FALSE}
        if(familyn!=7) incl[names(objrFinal$par)=="zeta"] <- FALSE
        if(familyn==0 || familyn==2 || familyn==7 || familyn==8) incl[names(objrFinal$par)=="lg_phi"] <- FALSE

        if(num.lv>0) {
          incld[names(objrFinal$par)=="Au"] <- TRUE
          incld[names(objrFinal$par)=="u"] <- TRUE
        } else {
          incl[names(objrFinal$par)=="lambda"] <- FALSE;
          incl[names(objrFinal$par)=="lambda2"] <- FALSE;
        }
          
        if(row.eff=="random") {
          incld[names(objrFinal$par)=="lg_Ar"] <- TRUE
          incld[names(objrFinal$par)=="r0"] <- TRUE
          incl[names(objrFinal$par)=="r0"] <- FALSE; 
        } else {
          incl[names(objrFinal$par)=="log_sigma"] <- FALSE
          if(row.eff==FALSE) incl[names(objrFinal$par)=="r0"] <- FALSE
          if(row.eff=="fixed") incl[1] <- FALSE
        }
        
        if(is.null(randomX)) {
          incl[names(objrFinal$par)%in%c("Br","sigmaB","sigmaij")] <- FALSE
        } else {
          incl[names(objrFinal$par)=="Abb"] <- FALSE; incld[names(objrFinal$par)=="Abb"] <- TRUE
          incl[names(objrFinal$par)=="Br"] <- FALSE; incld[names(objrFinal$par)=="Br"] <- TRUE
          if(NCOL(xb)==1) incl[names(objrFinal$par) == "sigmaij"] <- FALSE
        }


        if(method=="LA" || (num.lv==0 && (row.eff!="random" && is.null(randomX)))){

          covM <- try(MASS::ginv(sdr[incl,incl]))
          se <- try(sqrt(diag(abs(covM))))
          if(num.lv > 0 || row.eff == "random" || !is.null(randomX)) {
            sd.random <- sdrandom(objrFinal, covM, incl)
            prediction.errors <- list()
            
            if(row.eff=="random"){
              prediction.errors$row.params <- diag(as.matrix(sd.random))[1:length(out$params$row.params)];
              sd.random <- sd.random[-(1:length(out$params$row.params)),-(1:length(out$params$row.params))]
            }
            if(!is.null(randomX)){
              prediction.errors$Br  <- matrix(diag(as.matrix(sd.random))[1:(ncol(xb)*ncol(y))], ncol(xb), ncol(y));
              sd.random <- sd.random[-(1:(ncol(xb)*ncol(y))),-(1:(ncol(xb)*ncol(y)))]
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
          A.mat <- sdr[incl,incl] # a x a
          D.mat <- sdr[incld,incld] # d x d
          B.mat <- sdr[incl,incld] # a x d
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
        if(family %in% c("gaussian","tweedie","gamma","beta")) {
          se.lphis <- se[1:p];
          out$sd$phi <- se.lphis*out$params$phi;
          names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
        }
        if(family %in% c("ZIP")) {
          se.phis <- se[1:p];
          out$sd$phi <- se.phis*exp(lp0)/(1+exp(lp0))^2;#
          names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
        }
        if(!is.null(randomX)){
          nr <- ncol(xb)
          out$sd$sigmaB <- se[1:ncol(xb)]*c(sqrt(diag(out$params$sigmaB)));
          names(out$sd$sigmaB) <- c(paste("sd",colnames(xb),sep = "."))
          se <- se[-(1:ncol(xb))]
          if(nr>1){
            out$sd$corrpar <- se[1:(nr*(nr-1)/2)]
            se <- se[-(1:(nr*(nr-1)/2))]
          }
        }
        if(row.eff=="random") {
          out$sd$sigma <- se[1:length(out$params$sigma)]*c(out$params$sigma[1],rep(1,length(out$params$sigma)-1));
          names(out$sd$sigma) <- "sigma";
          se=se[-(1:(length(out$params$sigma)))]
          if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,3))) {out$sd$rho <- se[1]*(1-out$params$rho^2)^1.5; se = se[-1]}
          if((rstruc ==2 | (rstruc == 1)) & (cstrucn ==2)) {out$sd$rho <- se[1]*out$params$rho; se = se[-1]}
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

      }
      
      })
  if(inherits(tr, "try-error")) { cat("Standard errors for parameters could not be calculated.\n") }


  return(out)
}

