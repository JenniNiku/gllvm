########################################################################################
## GLLVM fourth corner model, with estimation done via Laplace and Variational approximation using TMB-package
## Original author: Jenni Niku
##########################################################################################
trait.TMB <- function(
      y, X = NULL,TR=NULL,formula=NULL, num.lv = 2, family = "poisson", num.lv.cor=0, corWithin = FALSE,
      Lambda.struc = "unstructured", Ab.struct = "unstructured", Ar.struc="diagonal", row.eff = FALSE, reltol = 1e-6, seed = NULL,
      maxit = 3000, max.iter=200, start.lvs = NULL, offset=NULL, sd.errors = FALSE,trace=FALSE,
      link="logit",n.init=1,n.init.max = 10, start.params=NULL,start0=FALSE,optimizer="optim", dr=NULL, rstruc =0, cstruc = "diag", dist = matrix(0),  scalmax=10, MaternKappa = 1.5,
      starting.val="res",method="VA",randomX=NULL,Power=1.5,diag.iter=1, Ab.diag.iter = 0, dependent.row = FALSE,
      Lambda.start=c(0.2, 0.5), jitter.var=0, yXT = NULL, scale.X = FALSE, randomX.start = "zero", beta0com = FALSE, rangeP = NULL,
      zeta.struc = "species", quad.start=0.01, start.struc="LV",quadratic=FALSE, optim.method = "BFGS", disp.group = NULL, NN=matrix(0), setMap = NULL, Ntrials = 1) {
  if(is.null(X) && !is.null(TR)) stop("Unable to fit a model that includes only trait covariates")
  if(!is.null(start.params)) starting.val <- "zero"
  
  if(!(family %in% c("poisson","negative.binomial","binomial","tweedie","ZIP", "ZINB", "gaussian", "ordinal", "gamma", "exponential", "beta", "betaH", "orderedBeta")))
    stop("Selected family not permitted...sorry!")
  if(!(Lambda.struc %in% c("unstructured","diagonal","bdNN","UNN")))
    stop("Lambda matrix (covariance of variational distribution for latent variable) not permitted...sorry!")
  
  objrFinal <- optrFinal <- NULL
  cstrucn = switch(cstruc, "diag" = 0, "corAR1" = 1, "corExp" = 2, "corCS" = 3, "corMatern" = 4)
  
  term <- NULL
  n <- nr <- nu <- dim(y)[1]; p <- dim(y)[2];
  times = 1
  if(is.null(disp.group)) disp.group <- 1:NCOL(y)
  if(family=="binomial" && length(Ntrials) != 1 && length(Ntrials) != p){
    stop("Supplied Ntrials is of the wrong length, should be of length 1 or the number of columns in y.")
  }else if(family=="binomial" && length(Ntrials) == 1){
    Ntrials <- rep(Ntrials, p)
  }
  
  # Structure for row effects
  model = 1
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
      if(is.null(rangeP)) {
        rangeP = AD1 = (apply(as.matrix(dist),2,max)-apply(as.matrix(dist),2,min))/scalmax
      } else {
        AD1 = rep(rangeP, ncol(dist))[1:ncol(dist)]
      }
      scaledc<-log(AD1)
      # AD1<-pmax(apply(as.matrix(dist),2,function(x) min(dist(unique(x), diag = FALSE))),1)
      # md<-min(dist(as.matrix(dist)%*%diag(1/(AD1), length(AD1)), diag = FALSE))/2
      # if(md>5) AD1=AD1*md
      # if(!is.null(setMap$scaledc)) {
      #   if( (length(setMap$scaledc)!= NCOL(dist))) stop("setMap$scaledc must be a numeric vector and have length that is same as the number of columns in 'dist'.")
      #   scaledc[is.na(setMap$scaledc)]=0
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
      if(is.null(rangeP)) {
        rangeP = AD1 = (apply(as.matrix(dist),2,max)-apply(as.matrix(dist),2,min))/scalmax
      } else {
        AD1 = rangeP
      }
      scaledc<-log(AD1)
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
  
  
  y <- as.data.frame(y)
  formula1 <- formula
  beta0com0 = beta0com
  if(method=="VA" && (family =="binomial")){ link <- "probit"}
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

  if(num.lv == 1) Lambda.struc <- "diagonal" ## Prevents it going to "unstructured" loops and causing chaos
  trial.size <- 1
  
  y <- as.matrix(y)
  if(!is.numeric(y)) stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
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
  if(is.null(rownames(y))) rownames(y) <- paste("Row",1:n,sep="")
  if(is.null(colnames(y))) colnames(y) <- paste("Col",1:p,sep="")
  if(!is.null(X)) { if(is.null(colnames(X))) colnames(X) <- paste("x",1:ncol(X),sep="") }
  


  if(family == "orderedBeta") {
    if (!(method %in% c("VA", "EVA"))) #"tweedie", 
      stop("family=\"", family, "\" : family not implemented with LA method, change the method to 'VA'")
    
    if((sum(y==1) + sum(y==0))==0){
      stop("No zeros or ones in the data, so use 'family = `beta` '")
    }
    if(!all(colSums(y==1)>0) & !all(colSums(y==0)>0)){
      warning("All species do not have zeros and ones. Setting 'zeta.struc = `common`'")
      zeta.struc = "common"
    }
  }
  
  out <-  list(y = y, X = X1, TR = TR1, num.lv = num.lv, row.eff = row.eff, logL = Inf, family = family, offset=offset,randomX=randomX,X.design=Xd,terms=term, method = method, Ntrials = Ntrials)

  if(is.null(formula) && is.null(X) && is.null(TR)){formula ="~ 1"}
  
  n.i <- 1;
  
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
  # if(n.init > 1) seed <- sample(1:10000, n.init)
  
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
      bstart <- start.values.randomX(y, X, family, formula=randomX, starting.val = randomX.start, Power = Power, link = link)
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
    ZINBphi <- ZINBphis <- NULL
    sigma <- 1
    
    if(n.init > 1 && trace) cat("initial run ",n.i,"\n");
    
    #### Calculate starting values
    res <- start.values.gllvm.TMB(y = y, X = X1, TR = TR1, family = family, offset=offset, trial.size = trial.size, num.lv = num.lv, start.lvs = start.lvs, seed = seed[n.i],starting.val=starting.val,Power=Power,formula = formula, jitter.var=jitter.var, #!!!
                                  yXT=yXT, row.eff = row.eff, TMB=TRUE, link=link, randomX=randomXb, beta0com = beta0com0, zeta.struc = zeta.struc, disp.group = disp.group, method=method, Ntrials = Ntrials)
    if(is.null(res$Power) && family == "tweedie")res$Power=1.1
    if(family=="tweedie"){
      Power = res$Power
      ePower = log((Power-1)/(1-(Power-1)))
      if(ePower==0)ePower=ePower-0.01
    }else{
      ePower = 0
    }
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
        if(rstruc==0 && row.eff=="random") row.params <- row.params[1:nr]#rstruc
        if(rstruc==1 && row.eff=="random") try(row.params <- (t(dr)%*%(row.params))/(dim(dr)[1]/dim(dr)[2]), silent = TRUE)#rstruc
        # if(rstruc<2 && row.eff=="random") row.params <- row.params[1:nr] #rstruc
        if (row.eff == "random") {
          sigma <- sd(row.params);
        }
      }
      vameans <- theta <- lambda <- NULL
      
      if(num.lv > 0) {
        sigma.lv <- res$sigma.lv
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
          for(i in 1:n) { lambda[i,,] <- diag(rep(1,num.lv),num.lv) }
        }
        if(Lambda.struc == "diagonal") {
          lambda <- matrix(1,n,num.lv)
        }
        zero.cons <- which(theta == 0)
        if(num.lv.cor>0){ # In correlation model, 
          rho_lvc<- rep(0, num.lv.cor);
          if((cstrucn == 2) | (cstrucn == 4)) {
            if(is.null(rangeP)) {
              rangeP = AD1 = (apply(as.matrix(dist),2,max)-apply(as.matrix(dist),2,min))/scalmax
            } else {
              AD1 = rangeP
            }
            scaledc<-log(AD1)
          }
        }
        # if(family == "betaH"){ # Own loadings for beta distr in hurdle model
        #   thetaH <- t(theta%*%diag(sigma.lv, nrow = length(sigma.lv), ncol = length(sigma.lv)))
        # }
        
        if(n.init > 1 && !is.null(res$mu) && starting.val == "res" && family != "tweedie") {
          if(family %in% c("ZIP","ZINB")) {
            lastart <- FAstart(res$mu, family="poisson", y=y, num.lv = num.lv, jitter.var = jitter.var[1], disp.group=disp.group)
          } else {
            lastart <- FAstart(res$mu, family=family, y=y, num.lv = num.lv, phis = res$phi, jitter.var = jitter.var[1], zeta.struc=zeta.struc, zeta = res$zeta, disp.group=disp.group, link = link)
          }
          theta <- lastart$gamma#/lastart$gamma
          vameans<-lastart$index#/max(lastart$index)
        }
      }else{
        sigma.lv <- matrix(0)
      }
      
    } else{
      if(all(dim(start.params$y)==dim(y)) && is.null(X)==is.null(start.params$X) && is.null(T)==is.null(start.params$TR) && row.eff == start.params$row.eff){
        beta0 <- start.params$params$beta0
        # common env params or different env response for each spp
        B <- NULL
        if(!is.null(TR) && !is.null(X)) {
          B <- start.params$params$B;
        }
        b.lv <- matrix(0)
        
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
        sigma.lv <- 0
        if(num.lv > 0) {
          sigma.lv <- start.params$params$sigma.lv
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
        if(num.lv.cor>0){ # sigmas are scale parameters # just diagonal values, not
          if(is.numeric(start.params$params$rho.lv) & ((cstrucn == 2) | (cstrucn == 4))) {
            # if(cstrucn == 4) start.params$params$rho.lv <- start.params$params$rho.lv[,-ncol(start.params$params$rho.lv), drop=FALSE]
            scaledc = colMeans(as.matrix(start.params$params$rho.lv)); 
            if(length(scaledc) < ncol(dist) ) scaledc <- rep(scaledc, ncol(dist))[1:ncol(dist)]
          }
        }
        if(family == "negative.binomial" && start.params$family == "negative.binomial" && !is.null(start.params$params$phi)) {res$phi<-start.params$params$phi}
      } else { stop("Model which is set as starting parameters isn't the suitable you are trying to fit. Check that attributes y, X, TR and row.eff match to each other.");}
    }
    if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)
    
    
### Starting values for dispersion/shape parameters
    
    if(family == "negative.binomial") {
      phis <- res$phi
      if (any(phis > 10))
        phis[phis > 50] <- 50
      if (any(phis < 0.02))
        phis[phis < 0.02] <- 0.02
      res$phi <- phis
      phis <- 1/phis
      
    }
    if (family == "ZIP" && starting.val=="res") {
      phis <- res$phi
      phis <- phis / (1 - phis)
    }
    if (family == "ZINB" && starting.val=="res") {
      phis <- res$phi
      phis <- phis / (1 - phis)
      
      ZINBphis <- res$ZINB.phi
      if (any(ZINBphis > 100))
        ZINBphis[ZINBphis > 100] <- 100
      if (any(ZINBphis < 0.01))
        ZINBphis[ZINBphis < 0.01] <- 0.01
      res$ZINB.phi <- ZINBphis
      ZINBphis <- 1/ZINBphis
    }
    
    if(family == "tweedie") {
      phis <- res$phi; 
      if(any(phis>10)) phis[phis>10]=10; 
      if(any(phis<0.10))phis[phis<0.10]=0.10; 
      phis= (phis)
    }
    
    if (family %in% c("ZIP","ZINB") && is.null(phis)) {
      if(length(unique(disp.group))!=p){
        phis <- (sapply(1:length(unique(disp.group)),function(x)mean(y[,which(disp.group==x)]==0))*0.98 + 0.01)[disp.group]
      }else{
        phis <- (colMeans(y == 0) * 0.98) + 0.01  
      }
      phis <- phis / (1 - phis)
      } # ZIP probability
    
    if (family %in% c("gaussian", "gamma", "beta", "betaH", "orderedBeta")) {
      phis <- res$phi
      if (family %in% c("betaH", "orderedBeta") & is.null(res$phi)) {
        phis <- rep(5,p)
      }
    }
    
### Starting values for cut-off parameters
    
    if(family=="ordinal"){
      K = max(y00)-min(y00)
      if(zeta.struc=="species"){
        zeta <- c(t(res$zeta[,-1]))
        zeta <- zeta[!is.na(zeta)]
      }else{
        zeta <- res$zeta[-1]
      }
      
    } else if(family=="orderedBeta") {
      zeta <- rep(0,p)
      # if(any(y==1)) 
      zeta <- c(zeta,rep(3,p))
    } else {
      zeta = 0
    }
    
### Jittering for row effs/random coefs
    if(jitter.var.r>0){ 
      if(row.eff == "random") row.params <- row.params + rnorm(n, 0, sd = sqrt(jitter.var.r));
      if(!is.null(randomX)) Br <- Br + t(MASS::mvrnorm(p, rep(0, nrow(Br)),diag(nrow(Br))*jitter.var.r));
    }
    
    q <- num.lv
    
    a <- c(beta0)
    if(num.lv > 0) {
      # diag(theta) <- log(diag(theta)) # !!!
      theta <- theta[lower.tri(theta, diag = F)]
      u <- vameans
    }

    if(!is.null(phis)) {
      phi=(phis)
    } else {
      phi <- rep(1,p)+runif(p,0,0.001) 
      if (family %in% c("betaH", "orderedBeta")) {
        phi <- rep(5,p)
      }
      res$phi <- phi
    }
    
    if(!is.null(ZINBphis)) {
      ZINBphi <- ZINBphis 
    } else { 
      ZINBphi <- rep(1, p)+runif(p,0,0.001)  
      if(family=="ZINB") res$ZINBphi <- ZINBphi
    }
    
    if(!is.null(row.params)){ r0 <- row.params} else {r0 <- rep(0, n)}
    if(row.eff == "random" && rstruc ==0){ nlvr<-num.lv+1 } else {nlvr=num.lv}
    if(row.eff=="fixed"){xr <- matrix(1,1,p)} else {xr <- matrix(0,1,p)}
    
    
    optr<-NULL
    timeo<-NULL
    se <- NULL
    
## map.list defines parameters which are not estimated in this model
    
    map.list <- list()    
    if(is.list(setMap)) map.list <- setMap
    # thetaH = matrix(0)
    # map.list$thetaH = factor(NA)
    # map.list$bH <- factor(NA) # not used
    map.list$b_lv <- factor(NA) # not used
    map.list$sigmab_lv = factor(NA)
    map.list$Ab_lv = factor(NA)
    if(family %in% c("poisson","binomial","ordinal","exponential")) {
      map.list$lg_phi <- factor(rep(NA,p))
    } else if(family %in% c("tweedie", "negative.binomial", "gamma", "gaussian", "beta", "betaH", "orderedBeta", "ZIP", "ZINB")){
      map.list$lg_phi <- factor(disp.group)
      if(family=="tweedie" && !is.null(Power))map.list$ePower = factor(NA)
      if(family=="ZINB")map.list$lg_phiZINB <- factor(disp.group)
    }
    
    if(!(family %in% c("ordinal", "orderedBeta"))) map.list$zeta <- factor(NA)
    if((family %in% c("orderedBeta"))){
      if(zeta.struc=="species"){
        zetamap = c(1:length(zeta))
        if(!all(colSums(y==0)>0))
          zetamap[1:p] <- 1
        if(!all(colSums(y==1)>0))
          zetamap[-(1:p)] <- max(zetamap[1:p])+1
        map.list$zeta = factor( zetamap)
        
      }else{
        zetamap <- c(rep(1,p))
        # if(any(y==1))
        zetamap <- c(zetamap,rep(max(zetamap)+1,p))
        map.list$zeta <- factor( c(zetamap) )
      }
    }
    if(family != "tweedie"){map.list$ePower = factor(NA)}
    if(family!="ZINB")map.list$lg_phiZINB <- factor(rep(NA,p))
    if(row.eff==FALSE) map.list$r0 <- factor(rep(NA,n))

    
    extra <- c(0,1,0)
    
    # Common intercept
    if(beta0com){
      extra[2] <- 0
      Xd<-cbind(1,Xd)
      a <- a*0
      B<-c(mean(a),B)
      map.list$b<-factor(rep(NA,length(a)))
    }
    
    
    ## Set up starting values for scale (and shape) parameters for correlated LVs
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
        } else if(cstruc=="corMatern"){
          # Fix matern smoothness by default
          maprho = matrix(iv, nrow(rho_lvc), ncol(rho_lvc))
          maprho[, ncol(maprho)] = NA
          map.list$rho_lvc = factor(c(maprho))
        }
      }
      res$rho.lv = rho_lvc
    } else {
      rho_lvc <- matrix(0)
      map.list$rho_lvc = factor(NA) 
    }
    
    
### set starting values for variational distribution covariances
    # Variational covariances for latent variables
    if(num.lv > 0){
      if(is.null(start.params) || start.params$method=="LA" || num.lv.cor>0){
        if(Lambda.struc=="diagonal" || (Lambda.struc=="bdNN") || (Lambda.struc=="LR") || diag.iter>0){
          Au <- log(rep(Lambda.start[1],num.lv*n))
        } else{
          Au <- c(log(rep(Lambda.start[1],num.lv*n)),rep(0,num.lv*(num.lv-1)/2*n)) #1/2, 1
        }
      } else{
        Au <- NULL
        for(d in 1:num.lv) {
          if(start.params$Lambda.struc=="unstructured" || length(dim(start.params$A))==3){
            Au <- c(Au,log(start.params$A[,d,d]))
          } else {
            Au <- c(Au,log(start.params$A[,d]))
          }
        }
        if(Lambda.struc!="diagonal" && diag.iter==0){
          Au <- c(Au,rep(0,num.lv*(num.lv-1)/2*n))
        }
      }
    } else { Au <- 0}
    
    # Variational covariances for structured/correlated LVs
    if((num.lv.cor>0) & (method %in% c("VA", "EVA"))){
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
        if(diag.iter>0){
          if(Astruc<3){
            Au <- c(Au[1:(nu*num.lv.cor)])
          } else {
            Au <- c(Au[1:(nu)])
            AQ<-diag(rep(log(Lambda.start[1]),num.lv.cor),num.lv.cor)
            Au<-c(Au,AQ[lower.tri(AQ, diag = TRUE)])
          }
        } else {
          if(Lambda.struc == "unstructured" && Astruc==1 & cstrucn==0){
            Au <- c(Au[1:(nu*num.lv.cor)], rep(0, nu*num.lv.cor*(num.lv.cor-1)/2))
          } else  if(Astruc==1){
            Au <- c(Au[1:(nu*num.lv.cor)], rep(0, num.lv.cor*nu*(nu-1)/2) )
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
      # if(corWithin) {
      #   if(Lambda.struc == "unstructured" && Astruc==1) {
      #     Au <- c(Au[1:(n*num.lv.cor)], rep(0,sum(lower.tri(matrix(0,n,n))[,1:2])*num.lv.cor) )
      #   } else if(Lambda.struc == "bdNN" && Astruc==2){
      #     Au <- c(Au[1:(n*num.lv.cor)], rep(0,nrow(NN)*num.lv.cor*nu) )
      #     # Au <- c(Au[1:(n*num.lv.cor)], rep(0,length(NN)*num.lv.cor) )
      #   }
      # } else {
      #   u <- as.matrix(u[1:nu,])
      #   Au <- Au[1:(nu*num.lv.cor)]
      #   if(Lambda.struc == "unstructured" && Astruc==1 & cstrucn==0 & diag.iter==0){
      #     Au <- c(Au[1:(nu*num.lv.cor)], rep(0, nu*num.lv.cor*(num.lv.cor-1)/2))
      #   } else {
      #     Au <- Au[1:(nu*num.lv.cor)]
      #   }
      # }
    }
      
    # Variational covariances for  random rows
    if(row.eff == "random"){
      if(rstruc ==1){
        lg_Ar <- rep(log(Lambda.start[2]), nr)
      } else {
        lg_Ar <- rep(log(Lambda.start[2]), n)
      }
      
      if(rstruc == 0 && nlvr>num.lv && (num.lv.cor==0) & (Ar.struc!="diagonal")){
        lg_Ar<-c(lg_Ar, rep(0, num.lv*n))
      }
      if(rstruc == 1 & (cstrucn %in% c(1,2,3,4)) & Ar.struc!="diagonal"){
        lg_Ar<-c(lg_Ar, rep(0, nr*(nr-1)/2))
      }
      if(rstruc == 2 & Ar.struc!="diagonal"){
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
    
      
  ### Specify parameter.list, data.list and map.list
    

    # For Laplace method, specify random parameters to randomp
    randomp= NULL #c("u","r0,"Br")
    randoml=c(0,0,0)
    
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
      map.list$sigmaLV = factor(NA)
    }
    if(num.lv.cor>0){
      if(!corWithin) {
        if(nrow(u) != nu){
          u=as.matrix((t(dr)%*%u/colSums(dr))[1:nu,, drop=FALSE])
        }
      }
    }

    ## Row effect settings
    if(row.eff=="random"){
      randoml[1] <- 1
      randomp <- c(randomp,"r0")
      
      if(dependent.row && (rstruc == 0)) 
        sigma<-c(sigma[1], rep(0, num.lv))
      if((rstruc %in% 1:2)) {
        if(cstrucn %in% c(1,3)) {
          sigma = c(log(sigma[1]),0)
        } else if(cstrucn %in% c(2)){
          sigma = c(log(sigma[1]),scaledc)
          if(is.null(setMap$log_sigma)) map.list$log_sigma = factor( c(1, rep(2,length(sigma)-1) ) )
        } else if(cstrucn %in% c(4)){
          sigma = c(log(sigma[1]),scaledc)
          # Fix matern smoothness by default
          if(is.null(setMap$log_sigma)) map.list$log_sigma = factor( c(1, rep(2,length(sigma)-1), NA) )
          sigma = c(sigma,log(MaternKappa))
        } else {
          sigma = c(log(sigma[1]))
        }
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
    
    
    
    
    ### family settings
    
    if(family == "poisson") { familyn=0}
    if(family == "negative.binomial") { familyn=1}
    if(family == "binomial") {
      familyn <- 2;
      if(link=="probit") extra[1] <- 1
    }
    if(family == "gaussian") {familyn=3}
    if(family == "gamma") {familyn=4}
    if(family == "tweedie"){ familyn <- 5}
    if(family == "ZIP"){ familyn <- 6;}
    if(family == "ordinal") {familyn=7}
    if(family == "exponential") {familyn=8}
    if(family == "beta"){
      familyn=9
      if(link=="probit") extra[1] <- 1
    }
    if(family == "betaH"){ 
      familyn = 10
      if(link=="probit") extra[1]=1
      # bH <- rbind(a,b)
      # extra[2] <- 0
      # Xd<-cbind(1,Xd)
      # bH<-matrix(B)
      # if(num.lv>0) {
      #   mapLH<-factor(1:length(thetaH))
      #   mapLH[lower.tri(thetaH)] <- NA
      #   map.list$thetaH <- factor(mapLH)
      # } else {
      #   thetaH<- matrix(0);
      #   map.list$thetaH = factor(NA)
      # }
    }
    if(family == "ZINB"){ familyn <- 11;}
    if(family == "orderedBeta") {familyn=12}
    
    
    
    ## To improve starting values for quadratic model
    if(starting.val!="zero" && start.struc != "LV" && quadratic == TRUE && num.lv>0 && method == "VA"){
      map.list2 <- map.list
      map.list2$r0 = factor(rep(NA, length(r0)))
      map.list2$b = factor(rep(NA, length(rbind(a))))
      map.list2$B = factor(rep(NA, length(B)))
      map.list2$Br = factor(rep(NA,length(Br)))
      map.list2$lambda = factor(rep(NA, length(theta)))
      map.list2$sigmaLV = factor(rep(NA, length(theta)))
      map.list2$u = factor(rep(NA, length(u)))
      map.list2$lg_phi = factor(rep(NA, length(phi)))
      map.list2$lg_phiZINB = factor(rep(NA, length(ZINBphi)))
      map.list2$sigmaB = factor(rep(NA,length(sigmaB)))
      map.list2$sigmaij = factor(rep(NA,length(sigmaij)))
      map.list2$log_sigma = factor(rep(NA, length(sigma)))
      map.list2$Au = factor(rep(NA, length(Au)))
      map.list2$lg_Ar = factor(rep(NA, length(lg_Ar)))
      map.list2$Abb = factor(rep(NA, length(Abb)))
      map.list2$zeta = factor(rep(NA, length(zeta)))

      parameter.list = list(r0=matrix(r0), b = rbind(a), b_lv = matrix(0), sigmab_lv = 0, Ab_lv = 0, B=matrix(B), Br=Br, lambda = theta, lambda2 = t(lambda2), sigmaLV = (sigma.lv), u = u, lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij, log_sigma=c(sigma), rho_lvc=rho_lvc, Au=Au, lg_Ar=lg_Ar, Abb=Abb, zeta=zeta, ePower = ePower, lg_phiZINB = log(ZINBphi)) #, scaledc=scaledc, bH=bH, thetaH = thetaH

      objr <- TMB::MakeADFun(
        data = list(y = y, x = Xd, x_lv = matrix(0), xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, num_RR = 0, num_lv_c = 0, num_corlv=num.lv.cor, family=familyn, extra=extra, quadratic = 1, method=switch(method, VA=0, EVA=2), model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist, Astruc=Astruc, NN = NN, Ntrials = Ntrials), silent=!trace,
        parameters = parameter.list, map = map.list2,
        inner.control=list(mgcmax = 1e+200),
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
    
    if( (method %in% c("VA", "EVA")) && (num.lv>0 || row.eff=="random" || !is.null(randomX) || (family =="orderedBeta")) ){
      parameter.list = list(r0=matrix(r0), b = rbind(a),  b_lv = matrix(0), sigmab_lv = 0, Ab_lv = 0, B=matrix(B), Br=Br, lambda = theta, lambda2 = t(lambda2), sigmaLV = (sigma.lv), u = u, lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij, log_sigma=c(sigma), rho_lvc=rho_lvc, Au=Au, lg_Ar=lg_Ar, Abb=Abb, zeta=zeta, ePower = ePower, lg_phiZINB = log(ZINBphi)) #, scaledc=scaledc, bH=bH, thetaH = thetaH

      objr <- TMB::MakeADFun(
        data = list(y = y, x = Xd, x_lv = matrix(0), xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, num_RR = 0, num_lv_c = 0, num_corlv=num.lv.cor, quadratic = ifelse(quadratic!=FALSE,1,0), family=familyn, extra=extra, method=switch(method, VA=0, EVA=2), model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist, Astruc=Astruc, NN = NN, Ntrials = Ntrials), silent=!trace,
        parameters = parameter.list, map = map.list,
        inner.control=list(mgcmax = 1e+200),
        DLL = "gllvm")
      
    } else {
      Au=0; Abb=0; lg_Ar=0;
      map.list$Au <- map.list$Abb <- map.list$lg_Ar <- factor(NA)
      
      parameter.list = list(r0=matrix(r0), b = rbind(a), b_lv = matrix(0), sigmab_lv = 0, Ab_lv = 0, B=matrix(B), Br=Br, lambda = theta, lambda2 = t(lambda2), sigmaLV = (sigma.lv), u = u, lg_phi=log(phi), sigmaB=log(sqrt(diag(sigmaB))), sigmaij=sigmaij, log_sigma=c(sigma), rho_lvc=rho_lvc, Au=Au, lg_Ar=lg_Ar, Abb=Abb, zeta=zeta, ePower = ePower, lg_phiZINB = log(ZINBphi)) #, scaledc=scaledc, thetaH = thetaH, bH=bH
      data.list <- list(y = y, x = Xd, x_lv = matrix(0), xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, num_RR = 0, num_lv_c = 0, num_corlv=num.lv.cor, quadratic = 0, family=familyn,extra=extra,method=1,model=1,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist, Astruc=Astruc, NN = NN, Ntrials = Ntrials)

      if(family == "ordinal"){
        data.list$method = 0
      }
      
      objr <- TMB::MakeADFun(
        data = data.list, silent=!trace,
        parameters = parameter.list, map = map.list,
        inner.control=list(mgcmax = 1e+200,tol10=0.01),
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
    
    if(diag.iter>0 && !(Lambda.struc %in% c("diagonal", "diagU")) && (method %in% c("VA", "EVA")) && (num.lv>1 || !is.null(randomX)) && !inherits(optr,"try-error")){
      objr1 <- objr
      optr1 <- optr
      param1 <- optr$par
      nam <- names(param1)
      if(length(param1[nam=="r0"])>0){ r1 <- matrix(param1[nam=="r0"])} else {r1 <- matrix(r0)}
      if(length(param1[nam=="b"])>0){ b1 <- rbind(param1[nam=="b"])} else {b1 <- rbind(rep(0,p))}
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
        u1 <- matrix(param1[nam=="u"], nrow(u), num.lv)
        Au<- c(pmax(param1[nam=="Au"],rep(log(1e-6), num.lv*nrow(u1))), rep(0,num.lv*(num.lv-1)/2*nrow(u1)))
        
        if (quadratic=="LV" | quadratic == T && start.struc == "LV"){
          lambda2 <- matrix(param1[nam == "lambda2"], byrow = T, ncol = num.lv, nrow = 1)#In this scenario we have estimated two quadratic coefficients before
        }else if(quadratic == T){
          lambda2 <- matrix(param1[nam == "lambda2"], byrow = T, ncol = num.lv, nrow = p)
        }
        
      } else {u1 <- u}
      if((num.lv)>0){sigma.lv1 <- param1[nam=="sigmaLV"]}else{sigma.lv1<-0}

      
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
          } else  if(Astruc==1){
            Au1 <- c(pmax(Au1[1:(nu*num.lv.cor)],log(1e-2)), rep(1e-3, num.lv.cor*nu*(nu-1)/2) )
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
              rho_lvc <- matrix((param1[nam=="rho_lvc"])[map.list$rho_lvc],nrow(rho_lvc),ncol(rho_lvc)); 
              rho_lvc[is.na(rho_lvc)]=0 
            } #rho_lvc[-1]<- param1[nam=="rho_lvc"]
          } else {
            rho_lvc[1:length(rho_lvc)]<- param1[nam=="rho_lvc"]
          }
        }
      } else if((num.lv)>0) {
        Au1<- c(pmax(param1[nam=="Au"],rep(log(1e-6), (num.lv)*nrow(u1))), rep(0,(num.lv)*((num.lv)-1)/2*nrow(u1)))
      } else {Au1<-Au}
      
      if(num.lv==0) {lambda1 <- 0; }
      if(family %in% c("poisson","binomial","ordinal","exponential", "betaH", "orderedBeta")){ lg_phi1 <- log(phi)} else {lg_phi1 <- param1[nam=="lg_phi"][disp.group]} #cat(range(exp(param1[nam=="lg_phi"])),"\n")
      if(family=="ZINB"){lg_phiZINB1 <- param1[nam=="lg_phiZINB"][disp.group]}else{lg_phiZINB1<-log(ZINBphi)}
      if(family=="tweedie" && is.null(Power))ePower = param1[nam == "ePower"]
      if(row.eff == "random"){
        log_sigma1 <- log(exp(param1[nam=="log_sigma"])+1e-3)
        if(!is.null(map.list$log_sigma)) log_sigma1 = log_sigma1[map.list$log_sigma]
        lg_Ar<- log(exp(param1[nam=="lg_Ar"])+1e-3)
      } else {log_sigma1 = 0}

      if(family %in% c("ordinal")){
        zeta <- param1[nam=="zeta"] 
      } else if(family %in% c("orderedBeta")){
        zeta <- c(rep(0,p),rep(param1[nam=="zeta"] ,p)[1:p])
      } else {
        zeta <- 0 
      }
      parameter.list = list(r0=r1, b = b1, b_lv = matrix(0), sigmab_lv = 0, Ab_lv = 0, B=B1, Br=Br1, lambda = lambda1, lambda2 = t(lambda2), sigmaLV = sigma.lv1, u = u1, lg_phi=lg_phi1, sigmaB=sigmaB1, sigmaij=sigmaij1, log_sigma=log_sigma1, rho_lvc=rho_lvc, Au=Au1, lg_Ar=lg_Ar, Abb=Abb, zeta=zeta, ePower = ePower, lg_phiZINB = lg_phiZINB1) #, scaledc=scaledc, thetaH = thetaH, bH=bH

      data.list = list(y = y, x = Xd, x_lv = matrix(0), xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, num_RR = 0, num_lv_c = 0, num_corlv=num.lv.cor, quadratic = ifelse(quadratic!=FALSE&num.lv>0,1,0), family=familyn, extra=extra, method=switch(method, VA=0, EVA=2), model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist, Astruc=Astruc, NN = NN, Ntrials = Ntrials)

      objr <- TMB::MakeADFun(
        data = data.list, silent=!trace,
        parameters = parameter.list, map = map.list,
        inner.control=list(mgcmax = 1e+200),
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
      if(length(param1[nam=="r0"])>0){ r1 <- matrix(param1[nam=="r0"])} else {r1 <- matrix(r0)}
      b1 <- rbind(param1[nam=="b"])
      B1 <- matrix(param1[nam=="B"])
      sigma.lv1 <- param1[nam=="sigmaLV"]
      
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
      
      if(family %in% c("poisson","binomial","ordinal","exponential")){ lg_phi1 <- log(phi)} else {lg_phi1 <- param1[nam=="lg_phi"][disp.group]}
      if(family=="ZINB"){lg_phiZINB1 <- param1[nam=="lg_ZINBphi"][disp.group]}else{lg_phiZINB1<-log(ZINBphi)}
      if(row.eff == "random"){
        log_sigma1 <- param1[nam=="log_sigma"]
        if(!is.null(map.list$log_sigma)) log_sigma1 = log_sigma1[map.list$log_sigma]
        lg_Ar<- param1[nam=="lg_Ar"]
      } else {log_sigma1 = 0}
      
      if(family == "ordinal"){ zeta <- param1[nam=="zeta"] } else { zeta <- 0 }
      
        parameter.list = list(r0=r1, b = b1, b_lv = matrix(0), sigmab_lv = 0, Ab_lv = 0, B=B1, Br=Br1, lambda = lambda1, lambda2 = t(lambda2), sigmaLV = sigma.lv1, u = u1, lg_phi=lg_phi1, sigmaB=sigmaB1, sigmaij=sigmaij1, log_sigma=log_sigma1, rho_lvc=rho_lvc, Au=Au, lg_Ar=lg_Ar, Abb=Abb, zeta=zeta, ePower = ePower, lg_phiZINB = lg_phiZINB1) #, scaledc=scaledc, thetaH = thetaH, bH=bH

      data.list = list(y = y, x = Xd, x_lv = matrix(0), xr=xr, xb=xb, dr0 = dr, offset=offset, num_lv = num.lv, num_RR = 0, num_lv_c = 0, quadratic = 1, family=familyn, extra=extra, method=switch(method, VA=0, EVA=2), model=1, random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), rstruc = rstruc, times = times, cstruc=cstrucn, dc=dist, Astruc=Astruc, NN = NN, Ntrials = Ntrials)

      objr <- TMB::MakeADFun(
        data = data.list, silent=!trace,
        parameters = parameter.list, map = map.list,
        inner.control=list(mgcmax = 1e+200),
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
    if(family %in% c("negative.binomial", "tweedie", "gaussian", "gamma", "beta", "betaH", "orderedBeta")) {
      phis=exp(param[names(param)=="lg_phi"])[disp.group]
      if(family=="tweedie" && is.null(Power)){
        Power = exp(param[names(param)=="ePower"])/(1+exp(param[names(param)=="ePower"]))+1
        names(Power) = "Power"
      }
    }
    if(family %in% c("ZIP","ZINB")) {
      if(family == "ZINB")ZINBphis <- exp(param[names(param)=="lg_phiZINB"])[disp.group]
      lp0 <- param[names(param)=="lg_phi"][disp.group]; out$lp0=lp0
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
    if(family == "orderedBeta") {
      zetas <- matrix((param[names(param)=="zeta"])[map.list$zeta],p,2)
    }
    
    bi<-names(param)=="b"
    Bi<-names(param)=="B"
    li<-names(param)=="lambda"
    si <- names(param)=="sigmaLV"
    li2 <- names(param)=="lambda2"
    ui<-names(param)=="u"
    
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
        diag(theta)<- 1 #sigma.lv 
      } else if(num.lv.cor==1) {
        theta[1,1]<- 1 #sigma.lv[1]
      }
      
      if(p>1) {
        theta[lower.tri(theta[,1:num.lv.cor,drop=F],diag=FALSE)] <- param[li];
      } else {
        theta <- as.matrix(1)
      }
      rho_lvc = param[names(param)=="rho_lvc"]
      if((cstrucn %in% c(1,3))) rho.lv<- param[names(param)=="rho_lvc"] / sqrt(1.0 + param[names(param)=="rho_lvc"]^2);
      if((cstrucn %in% c(2,4))){ 
        rho.lv<- exp(param[names(param)=="rho_lvc"]);
        # scaledc<- exp(param[names(param)=="scaledc"]);
      }
    } else if(num.lv > 0){
      sigma.lv <- abs(param[si])
      lvs <- (matrix(param[ui],n,num.lv))
      theta <- matrix(0,p,num.lv)
      diag(theta)<-1
      
      if(p>1) {
        theta[lower.tri(theta,diag=F)] <- param[li];
        if(quadratic!=FALSE){
          theta<-cbind(theta,matrix(-abs(param[li2]),ncol=num.lv,nrow=p,byrow=T))
        }
      } else {theta <- c(as.matrix(1),-abs(param[li2]))}
      # diag(theta) <- exp(diag(theta))#!!!
    }
    
    if(row.eff!=FALSE) {
      ri = names(param)=="r0"
      row.params = param[ri] 
      if(row.eff=="random"){
        sigma = exp(param[names(param)=="log_sigma"])[1]
        if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,3))) rho = param[names(param)=="log_sigma"][2] / sqrt(1.0 + param[names(param)=="log_sigma"][2]^2);
        if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(2,4))) {
          rho = exp(param[names(param)=="log_sigma"][-1]);
        }
        if(num.lv>0 && dependent.row && rstruc==0) sigma = c(sigma, (param[names(param)=="log_sigma"])[-1])
      }
    }
    if(!is.null(randomX)){
      Bri <- names(param)=="Br"
      Br <- matrix(param[Bri],ncol(xb),p)
      Sri <- names(param)=="sigmaB"
      L <- diag(ncol(xb))
      if(ncol(xb)>1){
        sigmaB <- diag(exp(param[Sri]), length(param[Sri]))
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
    # if(family %in% "betaH"){
    #   bHi <- names(param)=="bH"
    #   betaH <- (param[bHi])
    #   if(num.lv>0) {
    #     thetaH[!is.na(map.list$thetaH)] <- param[names(param)=="thetaH"]
    #   }
    # }
    
    cn<-colnames(Xd)
    if(beta0com){
      beta0=B[1]
      B = B[-1]
      cn<-colnames(Xd)
      Xd<-as.matrix(Xd[,-1])
      colnames(Xd)<-cn[-1]
    }
    new.loglik<-objr$env$value.best[1]
    
    
    #### Check if model fit succeeded/improved on this iteration n.i
    # Gradient check with n.i >2 so we don't get poorly converged models - relatively relaxed tolerance
    if(n.i>1){
      if(!is.null(objrFinal)){
        gr1 <- objrFinal$gr()
        gr1 <- as.matrix(gr1/length(gr1))
        norm.gr1 <- norm(gr1)
      }else{
        gr1 <- NaN
        norm.gr1 <- NaN
      }
      
      gr2 <- objr$gr()
      gr2 <- as.matrix(gr2/length(gr2))
      norm.gr2 <- norm(gr2)
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
      objrFinal<-objr1 <- objr; optrFinal<-optr1 <- optr;n.i.i<-0;
      out$logL <- new.loglik
      if(num.lv > 0) {
        out$lvs <- lvs
        out$params$theta <- theta
        if(num.lv>0)out$params$sigma.lv  <- sigma.lv
        if(nrow(out$lvs)==nrow(out$y)) rownames(out$lvs) <- rownames(out$y);
        rownames(out$params$theta) <- colnames(out$y)
        if(quadratic==FALSE)colnames(out$params$theta) <- colnames(out$lvs) <- paste("LV", 1:num.lv, sep="");
        if(quadratic!=FALSE){
          colnames(out$lvs) <- paste("LV", 1:num.lv, sep="");
          colnames(out$params$theta)<- c(paste("LV", 1:num.lv, sep=""),paste("LV", 1:num.lv, "^2",sep=""));
        }
        names(out$params$sigma.lv) <-  paste("LV", 1:num.lv, sep="");
      }
      if(!beta0com) names(beta0) <- colnames(out$y); 
      if(beta0com) names(beta0) <- "Community intercept";
      out$params$beta0 <- beta0;
      out$params$B <- B; names(out$params$B)=colnames(Xd)
      
      # row params
      if(row.eff!=FALSE) {
        if(row.eff=="random"){  
          out$params$sigma <- sigma; 
          names(out$params$sigma) <- "sigma"
          if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,2,3,4))){
            out$params$rho <- rho
          }
          if(num.lv>0 && dependent.row) names(out$params$sigma) <- paste("sigma",c("",1:num.lv), sep = "")
        }
        out$params$row.params <- row.params; 
        if(length(row.params) == n) names(out$params$row.params) <- rownames(out$y)
        if((length(row.params) == ncol(dr)) && (rstruc==1)) try(names(out$params$row.params) <- colnames(dr), silent = TRUE)
      }
      
      # LV correlation matrix parameters
      if(num.lv.cor>0 & cstrucn>0){
        out$params$rho.lv <- rho.lv; 
        if(cstrucn %in% c(2,4)){ 
          names(out$params$rho.lv) <- paste("rho.lv",1:length(out$params$rho.lv), sep = "") #[!is.na(map.list$sigma_lvc)]
        } else {
          names(out$params$rho.lv) <- paste("rho.lv",1:num.lv.cor, sep = "") 
        }
      }
      
      # if(family %in% "betaH"){
      #   out$params$betaH <- betaH;
      #   names(out$params$betaH)=cn; #colnames(Xd)
      #   if(num.lv>0) {
      #     out$params$thetaH <- thetaH
      #   }      
      # }
      
      # Dispersion parameters
      if(family =="negative.binomial") {
        out$params$inv.phi <- phis;
        out$params$phi <- 1/phis;
        names(out$params$phi) <- colnames(y);
        
        if(!is.null(names(disp.group))){
          try(names(out$params$phi) <- names(disp.group),silent=T)
        }
        names(out$params$inv.phi) <-  names(out$params$phi)
      }

      if(family %in% c("gaussian", "tweedie", "gamma","beta", "betaH", "orderedBeta")) {
        out$params$phi <- phis;
        names(out$params$phi) <- colnames(y);
        
        if(!is.null(names(disp.group))){
          try(names(out$params$phi) <- names(disp.group),silent=T)
        }
      }
      if(family %in% c("ZIP","ZINB")) {
        out$params$phi <- phis;
        names(out$params$phi) <- colnames(y);
        
        if(!is.null(names(disp.group))){
          try(names(out$params$phi) <- names(disp.group),silent=T)
        }
        if(family =="ZINB") {
          out$params$ZINB.inv.phi <- ZINBphis;
          out$params$ZINB.phi <- 1/ZINBphis;
          names(out$params$ZINB.phi) <- colnames(y);
          if(!is.null(names(disp.group))){
            try(names(out$params$ZINB.phi) <- names(disp.group),silent=T)
          }
          names(out$params$ZINB.inv.phi) <-  names(out$params$ZINB.phi)
        }
      }
      if (family %in% c("ordinal", "orderedBeta")) {
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
      if(family == "tweedie")out$Power = Power
      
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
            
            if(Astruc==0){
              Au <- exp(param[names(param)=="Au"])^2
              for (d in 1:(num.lv.cor)){
                A[,,d] <- diag(Au[(d-1)*nu+1:nu],nu,nu)
              }
            } else {
              # Alvm <- array(objr$report()$Alvm, dim=c(nu, nu, nMax))
              for (d in 1:nMax) {
                if(Astruc %in% c(3,4)){
                  A[,,d] <- objr$report()$Alvm[[d]]
                  # A[,,d] <- Alvm #%*%t(Alvm)
                } 
                # else {
                #   A[,,d] <- Alvm[,,d]%*%t(Alvm[,,d])
                # }
              }
            }
            if(Astruc %in% c(3,4)){
              AQ <- matrix(0,num.lv.cor,num.lv.cor)
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
          if(Astruc<3){ 
            nMax<- num.lv.cor
          } else {
            nMax<- 1
          }
          A <- array(0, dim=c(times*nu, times*nu, nMax))
          Alvm <- objr$report()$Alvm
          
          AQ <- NULL
          
          for (q in 1:nMax) {
            if(Astruc %in% c(3,4)){
              A[,,q] <- Alvm%*%t(Alvm)
            } else {
              A[,,q] <- Alvm[,,q]%*%t(Alvm[,,q])
            }
          }
          
          if(Astruc %in% c(3,4)){
            AQ <- matrix(0,num.lv.cor,num.lv.cor)
            AQ <- objr$report()$AQ
            AQ<-AQ%*%t(AQ)
          }
          
          
          out$AQ <- AQ
          out$A <- A
          
        } else if(nlvr>0){
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
            for (d in 1:num.lv){
              for(i in 1:n){
                A[i,(nlvr-num.lv)+ d,(nlvr-num.lv)+ d] <- exp(Au[(d-1)*n+i]);
              }
            }
            if(length(Au) > num.lv*n){
              k <- 0;
              for (c1 in 1:num.lv){
                r <- c1 + 1;
                while (r <= num.lv){
                  for(i in 1:n){
                    A[i,(nlvr-num.lv)+ r,(nlvr-num.lv)+ c1] <- Au[num.lv*n+k*n+i];
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
      
      if((method %in% c("VA", "EVA")) && !is.null(randomX)){
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
    seed.best <- seed[n.i]
    n.i <- n.i+1;
  }
  
  #Store the seed that gave the best results, so that we may reproduce results, even if a seed was not explicitly provided
  out$seed <- seed.best
  
  
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


  return(out)
}

