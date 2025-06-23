########################################################################################
## GLLVM, with estimation done via Variational approximation using TMB-package
## Original author: Jenni Niku
########################################################################################
gllvm.TMB <- function(y, X = NULL, lv.X = NULL, xr = matrix(0), formula = NULL, family = "poisson", 
                      num.lv = 2, num.lv.c = 0, num.RR = 0, num.lv.cor=0, lv.formula = NULL, corWithinLV = FALSE, randomB = FALSE, 
                      method = "VA",Lambda.struc = "unstructured", Ar.struc = "diagonal", sp.Ar.struc = "diagonal",  sp.Ar.struc.rank = NULL, Ab.diag.iter = 1, row.eff = FALSE, col.eff = FALSE, colMat = matrix(0), nn.colMat = NULL, colMat.approx = "NNGP", colMat.rho.struct = "single", randomX.start = "res", reltol = 1e-8, reltol.c = 1e-8,
                      maxit = 3000, max.iter = 200, start.lvs = NULL, offset = NULL,
                      trace = FALSE, link = "logit", n.init = 1, n.init.max = 10, restrict = 30, start.params = NULL, RElist = NULL, dr = matrix(0), dLV=NULL, cstruc = "diag", cstruclv = "diag", dist = list(matrix(0)), distLV = matrix(0),
                      optimizer = "optim", starting.val = "res", Power = 1.5, diag.iter = 1, dependent.row = FALSE, scalmax = 10, MaternKappa = 1.5, rangeP = NULL, zetacutoff = NULL,
                      Lambda.start = c(0.1,0.5), quad.start = 0.01, jitter.var = 0, jitter.var.br = 0, zeta.struc = "species", quadratic = FALSE, start.struc = "LV", optim.method = "BFGS", disp.group = NULL, NN = matrix(0), setMap = NULL, Ntrials = 1, beta0com = FALSE, csBlv = matrix(0)) { 

  # , Dthreshold=0
  # If there is no random effects/LVs set diag iter to zero:
  # if(!is.null(dr) && ncol(dr) != length(Ar.struc) && length(Ar.struc==1)){
  #   Ar.struc <- rep(Ar.struc, ncol(dr))
  # }else if(length(Ar.struc) != length(Ar.struc))stop("'Ar.struc' should be of the same length as the number of row effects.")
  n <- nr <- nu <- dim(y)[1]
  p <- dim(y)[2]
  
  if(((num.lv+num.lv.c)==0) & ((nrow(dr)!=n) || Ar.struc == "diagonal") & (randomB==FALSE)) diag.iter <-  0
  if(col.eff != "random" || sp.Ar.struc%in%c("diagonal","MNdiagonal","diagonalCL1")) Ab.diag.iter <- 0
  
  if(is.null(colMat) && !(sp.Ar.struc %in% c("diagonal","blockdiagonal")))sp.Ar.struc <- "blockdiagonal"
  
  # if(!is.null(start.params)) starting.val <- "zero"
  ignore.u <- FALSE
  
  times <- 1
  if(is.null(disp.group)) disp.group <- 1:NCOL(y)
  if(family=="binomial" && length(Ntrials) != 1 && length(Ntrials) != p){
    stop("Supplied Ntrials is of the wrong length, should be of length 1 or the number of columns in y.")
  } else if(family=="binomial" && length(Ntrials) == 1){
    Ntrials <- rep(Ntrials, p)
  }
  
  cstrucn = 0
  for (i in 1:length(cstruc)) {
    cstrucn[i] = switch(cstruc[i], "ustruc" = 0, "diag" = 0, "corAR1" = 1, "corExp" = 2, "corCS" = 3, "corMatern" = 4)
  }
  cstruclvn = switch(cstruclv, "ustruc" = 0, "diag" = 0, "corAR1" = 1, "corExp" = 2, "corCS" = 3, "corMatern" = 4)
  
  # Structure for row effects
  model = 0
  # if(rstruc==0){ # No structure
  #   dr <- diag(n)
  # }
  if(num.lv.cor==0 || is.null(dLV)){ # No structure
    dLV <- as(matrix(0), "TsparseMatrix")
  }
  if(col.eff != "random"){
    # nsp <- 0
    spdr <- colMat <- cs <- matrix(0)
    RElist <- colMat.old <- NULL
    Abstruc <- Abranks <- 0
    blocks = list(matrix(0))
    nncolMat <- matrix(0)
  }else{
    spdr <- Matrix::t(RElist$Zt)
    cs <- RElist$cs
    Xt <- Matrix::t(RElist$Xt)
    colMat.old <- colMat
    if(!is.null(colMat) && is.list(colMat)){
      if(colMat.approx == "NNGP" && (length(colMat)!=2 && !is.null(nn.colMat) || !"dist"%in%names(colMat))){
        stop("if nn.colMat<p 'colMat' must be a list of length 2: one Phylogenetic covariance matrix, and one (dist-named) distance matrix.")
      }else if(length(colMat)==1 && is.null(nn.colMat)){
         colMat <- colMat[[1]]
      }else{
      colMat.dist <- colMat$dist  
      colMat <- colMat[[which(names(colMat)!="dist")]]
      
      if(is.null(nn.colMat)){
        nn.colMat <- round(.3*p)
      }
      }
    }else if(is.matrix(colMat) && colMat.approx != "band"){
      nn.colMat <- p
    }
    
    spdr <- as.matrix(spdr)
    # nsp <- table(factor(colnames(spdr),levels=unique(colnames(spdr))))
    
    if(!is.null(colMat)  && all(dim(colMat)!=1)){
      if(!all(colnames(colMat) == colnames(y)))stop("Please make sure that the column names for 'y' and 'colMat' are the same.")
      colMat <- colMat[colnames(y), colnames(y)]
      if(exists("colMat.dist"))colMat.dist <- colMat.dist[colnames(y), colnames(y)]
      
      if(sp.Ar.struc%in%c("diagonal","blockdiagonal")){
        sp.Ar.struc.rank = 0
      }else if(is.null(sp.Ar.struc.rank) && sp.Ar.struc != "unstructured"){
        sp.Ar.struc.rank <- p
      }else if(is.null(sp.Ar.struc.rank) && sp.Ar.struc == "unstructured"){
        sp.Ar.struc.rank <- ncol(spdr)*p
      }
      if(ncol(colMat)!=nrow(colMat)){
        stop("Matrix for column effects must be square.")
      }
      if(ncol(colMat)!=p){
        stop("Matrix for column effects is of incorrect size.")
      }
      if(!isSymmetric(colMat, tol = 1e-12)){
        stop("Matrix for column effects is not symmetric.")
      }

      colMat <- try(cov2cor(colMat),silent = TRUE)
      # find blockstructure in colMat
      blocks = list()
      B = 1
      E = B
      if(nn.colMat == p)nncolMat <- matrix(0)
      if(nn.colMat < p)nncolMat <- NULL
      if(colMat.approx == "NNGP"){
        while(B<=p){
          while(E<p && (any(colMat[(E+1):p,B:E]!=0)|any(colMat[B:E,(E+1):p]!=0))){
            # expand block
            E = E+1;
          }
          # save block
          # here we work with blocks of the inverse of the correlation matrix
          if(nn.colMat==p)blocks[[length(blocks)+1]] = solve(colMat[B:E,B:E,drop=FALSE])
          if(nn.colMat<p){
            # here we work with blocks of the correlation matrix
            blocks[[length(blocks)+1]] = colMat[B:E,B:E,drop=FALSE]
            nncolMat <- cbind(nncolMat, sapply(1:ncol(colMat.dist[B:E,B:E,drop=FALSE]),function(i)head(order(colMat.dist[B:E,B:E,drop=FALSE][i,])[order(colMat.dist[B:E,B:E,drop=FALSE][i,])<i],min(i, nn.colMat))[1:p]))
          }
          E = E+1;
          B = E;
        }
      }else{
        while(B<=p){
          while(E<p && (any(colMat[(E+1):p,B:E]!=0)|any(colMat[B:E,(E+1):p]!=0))){
            # expand block
            E = E+1;
          }
          # save block
          if(nn.colMat<p && colMat.approx == "band"){
            blocks[[length(blocks)+1]] = colMat[B:E,B:E,drop=FALSE]
            nncolMat <- cbind(nncolMat, sapply(1:ncol(colMat[B:E,B:E,drop=FALSE]),function(i)(((i-1):(i-nn.colMat))[(i-1):(i-nn.colMat)>0])[1:p]))
          }
          E = E+1;
          B = E;
        }
      }
      nncolMat[is.na(nncolMat)] <- 0 ## using zeros to represent an empty cell
      blocksp <- unlist(lapply(blocks, ncol))
      # store total species and nr of species per block in first column, 0 and log determinants of each block in second column
      blocks = append(list(cbind(c(p,blocksp),c(0,unlist(lapply(blocks,function(x)-determinant(x)$modulus))))), blocks)
      if(sp.Ar.struc!="unstructured"){
        Abranks <- ifelse(sp.Ar.struc.rank>blocksp,blocksp,sp.Ar.struc.rank)
      }else{
        Abranks <- ifelse(sp.Ar.struc.rank>(ncol(spdr)*blocksp),ncol(spdr)*blocksp,sp.Ar.struc.rank)
      }
    }else{
      nncolMat <- matrix(0)
      colMat <- matrix(0)
      blocks <- list(matrix(0))
      Abstruc <- Abranks <- 0
    }
    if(is.null(cs)) cs <- matrix(0)
  }
  
  # number of random effects in each row effect
  # factor is used with table to keep the order in tact
  Astruc = 0;
  scaledc = 0;
  rho.lv =NULL  
  if(nrow(dr)==n){
    nr <- table(factor(colnames(dr),levels=unique(colnames(dr))))
    # distance matrix checks
    if(any(cstruc%in%c("corExp","corMatern"))){
      if(length(dist)!=sum(cstruc%in%c("corExp","corMatern"))){
        stop("Number of provided distance matrices should equal the number of spatially structured row effects.")
      }else{
        if(!all(unlist(lapply(dist, nrow))==nr[cstruc%in%c("corExp","corMatern")])){
          stop("Number of rows in 'dist' matrices should be same as number of units in the corresponding spatial row effect.")
        }
      }
    }
    if(any(cstruc%in%c("corExp","corMatern"))) {
      if(is.null(rangeP)) {
        rangeP = AD1 = unlist(mapply("/", lapply(mapply('-', lapply(dist,function(x)apply(x,2,max)), lapply(dist,function(x)apply(x,2,min)), SIMPLIFY = FALSE), mean), scalmax, SIMPLIFY = FALSE))
      } else {
        if(length(rangeP) >1 && length(rangeP) != sum(cstruc%in%c("corExp","corMatern"))){
          stop("The length of rangeP should be equal to the number of correlated structured row effects, or of length one.")
        }else if(length(rangeP)==1){
          rangeP = AD1 <- rep(rangeP,sum(cstruc%in%c("corExp","corMatern")))
        }else if(length(rangeP) == sum(cstruc%in%c("corExp","corMatern"))){
          AD1 = rangeP
        }
      }
      scaledc = lapply(AD1, log)
      # AD1 = pmax(apply(as.matrix(dist),2,function(x) min(dist(unique(x), diag = FALSE))),1)
      # md = min(dist(as.matrix(dist)%*%diag(1/(AD1), length(AD1)), diag = FALSE))/2
      # if(md>5) AD1 = AD1*md
      # scaledc = log(AD1)
      # if(!is.null(setMap$scaledc)) {
      #   if( (length(setMap$scaledc)!= NCOL(dist))) stop("setMap$scaledc must be a numeric vector and have length that is same as the number of columns in 'dist'.")
      #   scaledc[is.na(setMap$scaledc)]=0
      # }
    }
    # Ar.struc <- ifelse(nr==1, "diagonal", Ar.struc)
  }else{
    dr <- as(matrix(0), "TsparseMatrix")  
    # dimnames(dr) <- list(rep("site", n), rep("site", n))
    nr <- n
    # names(nr) = "site"
  }
  
  if(num.lv.cor > 0){#rstruc
    distLV<-as.matrix(distLV)
    if(is.null(dLV)) stop("Define structure for LVs'.")
    # LVs correlated within groups
    if(is.null(dLV)) stop("Define structure for LVs.")
    nu <- dim(dLV)[2]
    times <- n/nu#dim(dLV)[1]
    if((cstruclvn == 2) | (cstruclvn == 4)) {
      if(corWithinLV){
        if(is.null(distLV))
          distLV=matrix(1:times)
        if(NROW(distLV)!=times)
          stop("Number of rows in 'distLV' should be same as maximum number of units within groups when corWithinLV = TRUE")
      } else {
        if(is.null(distLV))
          distLV=matrix(1:nu)
        if(NROW(distLV)!=nu)
          stop("Number of rows in 'distLV' should be same as maximum number of groups when corWithinLV = FALSE")
      }
      if(is.null(rangeP)) {
        rangeP = AD1 = (apply(as.matrix(distLV),2,max)-apply(as.matrix(distLV),2,min))/scalmax
      } else {
        AD1 = rep(rangeP, ncol(distLV))[1:ncol(distLV)]
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
  
  tr <- NULL
  y <- as.matrix(y)
  if(is.null(lv.X)){
    lv.X <- matrix(0)
  }
  formula1 <- formula
  # if(method=="VA" && (family =="binomial")){ link="probit"}
  jitter.var.r <- 0
  if(length(jitter.var)>1){ 
    jitter.var.r <- jitter.var[2]
    jitter.var <- jitter.var[1]
  }
  
  if (!is.numeric(y))
    stop( "y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1.")
  # if ((family %in% c("ZIP")) && (method %in% c("VA", "EVA"))) #"tweedie", 
  #   stop("family=\"", family, "\" : family not implemented with VA method, change the method to 'LA'")
  if (is.null(rownames(y)))
    rownames(y) <- paste("Row", 1:n, sep = "")
  if (is.null(colnames(y)))
    colnames(y) <- paste("Col", 1:p, sep = "")
  if(family == "ordinal") {
    y00 <- y
    if(min(y) == 0){ y = y+1}
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
  
  out <- list( y = y, X = X, logL = Inf, num.lv = num.lv, num.lv.c = num.lv.c, col.eff = col.eff, colMat = colMat, family = family, X.design = X, method = method, zeta.struc = zeta.struc, Ntrials = Ntrials)
  
    #### Calculate starting values
    if((num.lv.c+num.lv+num.RR)==0 && !is.null(RElist) || randomX.start=="zero") RElist <- NULL # calculating starting values for REs and LVs
    fit <- start_values_gllvm_TMB(y = y, xr = xr, dr = dr, cstruc = cstruc, X = Xorig, formula = formula, lv.X = lv.X, TR = NULL, family = family, offset= offset, num.lv = num.lv, num.lv.c = num.lv.c, num.RR = num.RR, start.lvs = start.lvs, starting.val = starting.val, Power = Power, jitter.var = jitter.var, TMB=TRUE, link=link, zeta.struc = zeta.struc, disp.group = disp.group, method=method, randomB = randomB, Ntrials = Ntrials, Ab.struct = sp.Ar.struc, Ab.struct.rank = sp.Ar.struc.rank, colMat = colMat.old, nn.colMat = nn.colMat, RElist = RElist, beta0com = beta0com)

    if(is.null(fit$Power) && family == "tweedie")fit$Power=1.1
    if(family=="tweedie"){
      ePower = log((fit$Power-1)/(1-(fit$Power-1)))
      if(ePower==0)ePower=ePower-0.01
    }else{
      ePower = 0
    }
    
    ## Set initial values
    sigma <- 1;Br <- matrix(0);sigmaB <- 0;rho.sp.start = 0.5;
    
    if(!is.null(start.params)){
      if (all(dim(start.params$y) == dim(y))) {
        if(class(start.params)[2]=="gllvm.quadratic" && !isFALSE(quadratic)){
          if(start.params$num.lv>0)fit$params[,num.X+tail(1:ncol(fit$params), num.lv)] <- fit$params$theta[, tail(1:ncol(start.params$params$theta), start.params$num.lv)]
          if((start.params$num.lv.c+start.params$num.RR>0)){
           fit$params[,num.X+(2+num.lv+num.lv.c+num.RR):(1+num.lv+2*num.lv.c+2*num.RR)] <- start.params$params$theta[,head(tail(1:ncol(start.params$params$theta), num.lv+num.lv.c+num.RR), num.lv.c+num.RR)]
          }
        }
        if((num.lv.c+num.RR)>0 && (start.params$num.lv.c+start.params$num.RR>0) && ((num.lv.c+num.RR) == (start.params$num.lv.c+start.params$num.RR)) && !isFALSE(start.params$randomB) && !isFALSE(randomB) && randomB!="iid" && all.equal(ncol(start.params$lv.X.design), lv.X)){
          fit$sigmab_lv <- start.params$params$sigmaLvXcoef
          if(ncol(csBlv)==2)sigmab_lv <- c(sigmab_lv, start.params$params$corsLvXcoef[csBlv])
        }
        
        if((start.params$num.lv.c+start.params$num.RR)>0 && (num.lv.c+num.RR>0) && (num.lv.c+num.RR)==(start.params$num.lv.c+start.params$num.RR)){
          b.lv <- start.params$params$LvXcoef
        }
        
        fit$params[,1] <- start.params$params$beta0 ## column intercepts
        if (!is.null(X) && !is.null(start.params$X) && ncol(X)==ncol(start.params$X))
          fit$params[,2:(num.X + 1)] <- c(start.params$params$Xcoef) ## covariates coefficients
        
        if ((num.lv+(num.lv.c)) > 0 && (start.params$num.lv.c+start.params$num.lv)>0 && (start.params$num.lv.c+start.params$num.lv)==(num.lv+num.lv.c)){
          fit$sigma.lv <- sigma.lv <- start.params$params$sigma.lv
          fit$params[,(num.X+2):(num.lv+num.lv.c+num.RR+num.X+1)] <- start.params$params$theta
          if((num.lv.c+num.RR)>1)fit$params[,1+num.X+1:(num.lv.c+num.RR)][upper.tri(fit$params[,1+num.X+1:(num.lv.c+num.RR)])] <- 0
          if(num.lv>1)fit$params[,1+num.X+((num.lv.c+num.RR)+1)+(0:(num.lv-1))][upper.tri(fit$params[,1+num.X+((num.lv.c+num.RR)+1)+(0:(num.lv-1))])] <- 0
          
        }
        
        if (!isFALSE(row.eff) && !isFALSE(start.params$row.eff) && all.equal(row.eff, start.params$row.eff)) {
          fit$row.params.fixed <- start.params$params$row.params.fixed
          fit$row.params.random <- start.params$params$row.params.random
          if(nrow(dr)==n)
            fit$sigma <- start.params$params$sigma
        }## row parameters
        if ((num.lv+num.lv.c) > 0 && (start.params$num.lv+start.params$num.lv.c) >0 && (num.lv+num.lv.c) == (start.params$num.lv+start.params$num.lv.c)) {
          fit$sigma.lv <- start.params$params$sigma.lv
          fit$index <- matrix(start.params$lvs, ncol = num.lv+num.lv.c)
        }
        if(num.lv.cor>0){ # sigmas are scale parameters # just diagonal values, not
          if(is.numeric(start.params$params$rho.lv) & ((cstruclvn == 2) | (cstruclvn == 4))) {
            # if(cstruclvn == 4) start.params$params$rho.lv <- start.params$params$rho.lv[,-ncol(start.params$params$rho.lv), drop=FALSE]
            scaledc = colMeans(as.matrix(start.params$params$rho.lv)); 
            if(length(scaledc) < ncol(distLV) ) scaledc <- rep(scaledc, ncol(distLV))[1:ncol(distLV)]
          }
        }
        if(col.eff == "random" && start.params$col.eff$col.eff == "random" && ncol(spdr) == ncol(start.params$col.eff$spdr)){
          fit$Br <- start.params$params$Br
          fit$sigmaB <- fit$params$sigmaB#log(diag(start.params$params$sigmaB))
          # if(!is.null(cs) && (ncol(cs) == 2)){
          #   if(any(sigmaB[cs]==0)){
          #     sigmaB[cs] <- sigmaB[cs]+1e-5
          #   }
          #   sigmaB <- c(sigmaB, sigmaB[cs])
          # }
          if(!is.null(start.params$params$rho.sp))rho.sp.start <- start.params$params$rho.sp
        }
      } else {
        # findproblem = c("y","X", "row.eff", "col.eff")[!c(dim(start.params$y)==dim(y), is.null(X)==is.null(start.params$X), isTRUE(all.equal(row.eff, start.params$row.eff)), (col.eff == start.params$col.eff$col.eff)  )]
        
        # stop( "Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that attributes y, X and row.eff match to each other. Problem in: ", findproblem)
        stop("Starting model needs to be fitted to the same data.")
      }
    }
    
      beta0 <- fit$params[, 1]
      if((num.lv.c+num.RR)>0){b.lv <- fit$b.lv}else{b.lv<-matrix(0)}
      betas <- NULL
      if (!is.null(X))
        betas <- c(fit$params[, 2:(num.X + 1)])
      lambdas <- NULL
      
      if ((num.lv+(num.lv.c+num.RR)) > 0) {
        sigma.lv <- (fit$sigma.lv)
        lambdas <- as.matrix(fit$params[, (ncol(fit$params) - num.lv - (num.lv.c+num.RR) + 1):ncol(fit$params)])
        if(start.struc=="LV"&quadratic!=FALSE|quadratic=="LV"){
          lambda2 <- matrix(quad.start, ncol = num.lv + (num.lv.c+num.RR), nrow = 1)  
        }else if(start.struc=="all"&quadratic!=FALSE){
          lambda2 <- matrix(quad.start, ncol = num.lv + (num.lv.c+num.RR), nrow = p)
        }else if(quadratic==FALSE){
          lambda2 <- 0
        }
        if(randomB!=FALSE & randomB!="iid"){
          sigmab_lv <- fit$sigmab_lv
          if(ncol(csBlv)==2 && starting.val == "res"){
            sigmab_lv <- c(sigmab_lv, fit$sigmab_cors[csBlv])
          }else if(ncol(csBlv)==2){
            sigmab_lv <- c(sigmab_lv, rep(1e-3, nrow(csBlv)))
          }
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
          if((cstruclvn == 2) | (cstruclvn == 4)) {
            if(is.null(rangeP)) {
              rangeP = AD1 = (apply(as.matrix(distLV),2,max)-apply(as.matrix(distLV),2,min))/scalmax
            } else {
              AD1 = rep(rangeP, ncol(distLV))[1:ncol(distLV)]
            }
            scaledc<-log(AD1)
          }
        }
        # if(family == "betaH"){ # Own loadings for beta distr in hurdle model
        #   thetaH <- t(lambdas%*%diag(sigma.lv, nrow = length(sigma.lv), ncol = length(sigma.lv)))
        # }
      }
      
      if(col.eff == "random"){
        if(!is.null(RElist) && starting.val == "res" && randomX.start=="res" && (num.lv.c+num.RR+num.lv)>0){ # getting some improved starting values
          if(nrow(Xt)==n){
          B <- rep(0, ncol(spdr))
          B[colnames(spdr)%in%colnames(Xt)] <- fit$fitstart$B
          }
          sigmaB <- log(sqrt(diag(fit$fitstart$sigmaB)))
          if(ncol(cs)==2){
            sigmaB <- c(sigmaB, fit$fitstart$TMBfnpar[names(fit$fitstart$TMBfnpar) == "sigmaB"][(ncol(spdr)+1):(ncol(spdr)+nrow(cs))])
          }
          Br <- fit$fitstart$Br
        }else{
        bstart <- start_values_randomX(y, as.matrix(spdr), family, formula=formula(paste0("~",paste0(make.unique(colnames(spdr)),collapse="+"))), starting.val = randomX.start, Power = Power, link = link)
        B <- bstart$B
        Br <- bstart$Br
        if(jitter.var.br>0)Br <- Br + matrix(rnorm(prod(dim(B)), sd=sqrt(jitter.var.br)), nrow(Br), ncol(Br))
        if(Matrix::rankMatrix(bstart$sigmaB)<ncol(bstart$sigmaB)){
          bstart$sigmaB <- Matrix::nearPD(bstart$sigmaB)$mat
        }
        bstart$sigmaB <- t(chol(bstart$sigmaB))
        sigmaB <- log(diag(bstart$sigmaB))
          if(ncol(cs)==2){
            sigmaB <- c(sigmaB,bstart$sigmaB[cs])
          }
        }
        
      fit$Br <- Br
      if(nrow(Xt)==n)fit$B <- B
      # colMat signal strength
      if(any(colMat[row(colMat)!=col(colMat)]!=0))sigmaB <- c(sigmaB, rep(log(-log(rho.sp.start)),ifelse(colMat.rho.struct == "single", 1, ncol(spdr))))
      fit$sigmaB <- sigmaB
      }
      if(col.eff != "random"){
        sigmaB <- 0;Br <- matrix(0);B<-matrix(0)
      }else if(nrow(Xt)!=n && col.eff == "random"){
        B <- matrix(0)
      }
      row.params <- row.params.fixed <- row.params.random <- NULL
      
      # if ((nrow(xr)==n) || (nrow(dr)==n)) {
      #   row.params <- fit$row.params
      #   if(nrow(xr)==n){
      #     row.lm <- lm(row.params~0+xr)
      #     row.params.fixed <- coef(row.lm)
      #     row.params <- residuals(row.lm)
      #   }
      #   if (nrow(dr)==n) {
      #     row.params<-row.params
      #     try(row.params.random <- (Matrix::t(dr)%*%(row.params))/(dim(dr)[1]/dim(dr)[2]), silent = TRUE)
      #     sigma <- aggregate(as.matrix(row.params), by = list(row.names(row.params)), FUN = sd)[,2]
      #   }
      # }
      if ((nrow(xr)==n) || (nrow(dr)==n)) {
        row.params.random <- fit$row.params.random
        row.params.fixed <- fit$row.params.fixed
        sigma <- fit$sigma
      }
      
      lvs <- NULL
      if ((num.lv+num.lv.c) > 0)
        lvs <- matrix(fit$index, ncol = num.lv+num.lv.c)

    phis <- NULL
    ZINBphis <- NULL
    if (family == "negative.binomial") {
      phis <- fit$phi
      if (any(phis > 100))
        phis[phis > 100] <- 100
      if (any(phis < 0.01))
        phis[phis < 0.01] <- 0.01
      fit$phi <- phis
      phis <- 1/phis
    }
    if (family == "ZIP" && starting.val=="res") {
      phis <- fit$phi
      phis <- phis / (1 - phis)
    }
    if (family == "ZINB" && starting.val=="res") {
      phis <- fit$phi
      phis <- phis / (1 - phis)
      
      ZINBphis <- fit$ZINB.phi
      if (any(ZINBphis > 100))
        ZINBphis[ZINBphis > 100] <- 100
      if (any(ZINBphis < 0.01))
        ZINBphis[ZINBphis < 0.01] <- 0.01
      fit$ZINB.phi <- ZINBphis
      ZINBphis <- 1/ZINBphis
    }
    if (family == "tweedie") {
      phis <- fit$phi
      if (any(phis > 10))
        phis[phis > 10] <- 10
      if (any(phis < 0.10))
        phis[phis < 0.10] <- 0.10
      phis = (phis)
    }
    if (family %in%c("ZIP","ZINB") && is.null(phis)) {
      if(length(unique(disp.group))!=p){
        phis <- sapply(1:length(unique(disp.group)),function(x)mean(y[,which(disp.group==x)]==0))*0.98 + 0.01  
        phis <- phis[disp.group]
      }else{
        phis <- (colMeans(y == 0) * 0.98) + 0.01  
      }
      phis <- phis / (1 - phis)
    } # ZIP probability
    if (family %in% c("gaussian", "gamma", "beta", "betaH", "orderedBeta")) {
      phis <- fit$phi
      if (family %in% c("betaH", "orderedBeta")) {
        phis <- rep(5,p)
      }
    }
    if(family=="ordinal"){
      K = max(y00)-min(y00)
      if(zeta.struc=="species"){
        zeta <- c(t(fit$zeta[,-1]))
        zeta <- zeta[!is.na(zeta)]
      }else{
        zeta <- fit$zeta[-1]
      }
      
    } else if(family=="orderedBeta") {
      if(is.null(zetacutoff)){
        zeta <- rep(0,p)
        zeta <- c(zeta,rep(log(3),p))
      } else {
        zetacutoff<- matrix(zetacutoff, ncol=2)
        zeta <- rep(zetacutoff[,1],p)[1:p]
        zeta <- c(zeta,rep(log(zetacutoff[,2]),p)[1:p])
      }
    } else {
      zeta = 0
    }
    
    if (is.null(offset))
      offset <- matrix(0, nrow = n, ncol = p)
    
    current.loglik <- -1e6; iter <- 1; err <- 10;
    if(!is.null(row.params.fixed)){ r0f <- row.params.fixed} else {r0f <- rep(0,ncol(xr))}
    if(nrow(dr)!=n)r0r <- 0
    if(nrow(xr)!=n)r0f <- 0
    if(!is.null(row.params.random)){
      r0r <- row.params.random
    }
    if(beta0com) a <- rep(mean(beta0), p)
    if(!beta0com) a <- c(beta0)
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
      phi <- rep(1, p)+runif(p,0,0.001); 
      if (family %in% c("betaH", "orderedBeta")) {
        phi <- rep(5,p)
      }
      fit$phi <- phi
    }
    if(!is.null(ZINBphis)) {
      ZINBphi <- ZINBphis 
    } else { 
      ZINBphi <- rep(1, p)+runif(p,0,0.001) 
      if(family=="ZINB")fit$ZINBphi <- ZINBphi
    }
    
    
    q <- num.lv+(num.lv.c+num.RR)
    
    ## map.list defines parameters which are not estimated in this model
    
    map.list <- list()    
    if(is.list(setMap)) {
      map.list <- setMap
    }
    map.list$sigmaij <- factor(NA)

    xb<-matrix(0);sigmaij=0; lg_Ar=0; Abb=0; Ab_lv = 0;
    if(nrow(xr)!=n){
      map.list$r0f <- factor(NA)
    }
    if(nrow(dr)!=n){
      map.list$r0r <- factor(NA)
    }
    if(col.eff==FALSE) {map.list$Br <- factor(NA);map.list$sigmaB <- factor(NA); map.list$Abb <- factor(NA);map.list$B <- factor(NA)}
    if(randomB==FALSE){
      map.list$sigmab_lv <- factor(NA)
    }else if(randomB=="iid"){
      map.list$sigmab_lv <- factor(NA)
      sigmab_lv <- 0
    }
    if(family %in% c("poisson","binomial","ordinal","exponential")){
      map.list$lg_phi <- factor(rep(NA,p))
    } else if(family %in% c("tweedie", "negative.binomial", "gamma", "gaussian", "beta", "betaH", "orderedBeta", "ZIP","ZINB")){
      map.list$lg_phi <- factor(disp.group)
      if(family=="tweedie" && !is.null(Power))map.list$ePower = factor(NA)
      if(family=="ZINB" & is.null(map.list$lg_phiZINB)) map.list$lg_phiZINB <- factor(disp.group)
    }
    
    if(!(family %in% c("ordinal", "orderedBeta"))) map.list$zeta <- factor(NA)
    if((family %in% c("orderedBeta"))){
      if(zeta.struc=="species"){
        zetamap = c(1:length(zeta))
        if(!all(colSums(y==0, na.rm = TRUE)>0))
          zetamap[1:p] <- 1
        if(!all(colSums(y==1, na.rm = TRUE)>0))
          zetamap[-(1:p)] <- max(zetamap[1:p])+1
        map.list$zeta = factor( zetamap)
        if("zeta" %in% names(setMap)) map.list$zeta = factor(setMap$zeta)
        # map.list$zeta <- factor(c(rep(NA,p),1:p))
      }else{
        zetamap <- c(rep(1,p))
        # if(any(y==1))
        zetamap <- c(zetamap,rep(max(zetamap)+1,p))
        map.list$zeta <- factor( c(zetamap) )
        if("zeta" %in% names(setMap)) map.list$zeta = factor(setMap$zeta)
      }
    }
    if(family != "tweedie"){map.list$ePower = factor(NA)}
    if(family!="ZINB")map.list$lg_phiZINB <- factor(rep(NA,p))
    if((num.lv.c+num.RR)==0){
      map.list$b_lv = factor(rep(NA, length(b.lv)))
    }
    if((num.lv+num.lv.c)==0)map.list$sigmaLV = factor(NA)
    
    randoml=c(0,0,0, 0)
    if(nrow(dr)==n) randoml[1]=1
    if(col.eff == "random") randoml[4]  <- 1
    nlvr=num.lv+num.lv.c
    if(randomB!=FALSE){
      randoml[3]<-1
    }
    if(isFALSE(randomB)){
      sigmab_lv <- 0
    }else if(randomB=="iid"){
      sigmab_lv <- 0
    }
    if(!is.null(X)){Xd <- cbind(1,X)} else {Xd <- matrix(1,n)}
  
    # map species common effects for REs
    if(col.eff == "random" && nrow(Xt)==n){
      map.list$B <- 1:ncol(spdr)
      if(any(!colnames(Xt)%in%colnames(spdr))){
        stop("There was a problem with the model. Did you use ordered constrasts in the random effect perhaps?")
      }
      map.list$B[!colnames(spdr)%in%colnames(Xt)] <- NA
      map.list$B <- factor(map.list$B)
      B <- as.matrix(B)
      B[!colnames(spdr)%in%colnames(Xt)] <- 0
    }
    if(beta0com){
      if(is.null(map.list[["b"]])){
      map.list$b <- 1:(p*ncol(Xd))
      map.list$b[1:p] <- 1
      map.list$b <- factor(map.list$b)
      map.list$b <- map.list$b[order(rep(1:p,num.X+1))] # back to correct ordering
      }else{
        map.list$b[1:p]<-1
        map.list$b <- factor(map.list$b)
        map.list$b <- map.list$b[order(rep(1:p,num.X+1))] # back to correct ordering
      }
    }
    
    b <- NULL; if(!is.null(X)) b <- matrix(betas, ncol(X), p,byrow = TRUE)
    extra <- c(0,0,0)
    
    optr <- timeo <- NULL
    
    ## Set up starting values for scale (and shape) parameters for correlated LVs
    if(num.lv.cor>0 & cstruclvn>0){
      rho_lvc<- matrix(rep(0, num.lv.cor))
      if(cstruclvn==2){ #"corExp"
        if(is.null(rho.lv)) {
          rho.lv=rep(0, num.lv.cor) 
        } else if(length(rho.lv)==num.lv.cor) {
          rho.lv=c(log(rho.lv))
        }
        rho_lvc<- matrix(c(rep(scaledc, each=num.lv.cor)), num.lv.cor)
      } else if(cstruclvn==4){#"corMatern"
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
      
      if(cstruclvn %in% c(2,4)){
        iv<-rep(1:nrow(rho_lvc), ncol(rho_lvc)); 
        if(!is.null(setMap$rho_lvc)){
          if((length(setMap$rho_lvc)==length(rho_lvc))) 
            iv = (setMap$rho_lvc)
          map.list$rho_lvc = factor(iv)
        } else if(cstruclvn==2){ #cstruc=="corExp"
          maprho = matrix(iv, nrow(rho_lvc), ncol(rho_lvc))
          map.list$rho_lvc = factor(c(maprho))
        } else if(cstruclvn==4){ #cstruc=="corMatern"
          # Fix matern smoothness by default
          maprho = matrix(iv, nrow(rho_lvc), ncol(rho_lvc))
          maprho[, ncol(maprho)] = NA
          map.list$rho_lvc = factor(c(maprho))
        }
      }
      fit$rho.lv = rho_lvc
    } else {
      rho_lvc <- matrix(0)
      map.list$rho_lvc = factor(NA) 
    }
    
    ### VA method, used only if there is some random effects/LVs in the model
    
    if(((method %in% c("VA", "EVA")) && (nlvr>0 || (nrow(dr)==n) || !isFALSE(randomB) || col.eff == "random")) ){
      
      # Variational covariances for latent variables
      if((num.lv+num.lv.c)>0){
        if(is.null(start.params) || start.params$method=="LA" || num.lv.cor>0){
          if(Lambda.struc=="diagonal" || (Lambda.struc=="bdNN") || (Lambda.struc=="LR") || diag.iter>0){
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
        }
      } else { Au <- 0}
      
      # Variational covariances for structured/correlated LVs
      if(num.lv.cor>0){
        if(corWithinLV) {
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
              Au <- c(Au[1:(n*num.lv.cor)], rep(0,nrow(NN)*num.lv.cor*n) )
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
            } 
            #   else if(Astruc==6){
            #   Au <- c(Au[1:(n*num.lv.cor)], rep(0.001,NN[1]*num.lv.cor*n) )
            # }
          }
        } else {
          # u <- as.matrix(u[1:nu,])
          if(diag.iter>0){
            if(Astruc<3){
              Au <- c(Au[1:(nu*num.lv.cor)])
            } else{ # if(Astruc<6)
              Au <- c(Au[1:(nu)])
              AQ<-diag(rep(log(Lambda.start[1]),num.lv.cor),num.lv.cor)
              Au<-c(Au,AQ[lower.tri(AQ, diag = TRUE)])
            } 
            # else {
            #   Au <- c(Au[1:(nu*num.lv.cor)])
            # }
          } else {
            if(Lambda.struc == "unstructured" && Astruc==1 & cstruclvn==0){
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
            # else if(Astruc==6){
            #   Au <- c(Au[1:(nu*num.lv.cor)], rep(0.001,NN[1]*num.lv.cor*nu) )
            # }
          }
        }
        if(is.null(start.params)) sigma.lv <- (sigma.lv*0.5) #
        Au = Au + 1e-3
      }
      
      # Variational covariances for random slopes of const. ord.
      if((num.RR+num.lv.c)>0&!isFALSE(randomB)){
        if(randomB=="LV"){
          ab12 <- num.RR+num.lv.c
          ab3 <- ncol(lv.X)
        }else{
          ab12 <- ncol(lv.X)
          ab3 <- num.RR+num.lv.c
        }
        if(is.null(start.params) || start.params$method=="LA" || isFALSE(start.params$randomB)){
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
      
      # Variational covariances for species effects
      if(col.eff == "random"){
        if(sp.Ar.struc == "diagonal" || sp.Ar.struc== "blockdiagonal"){
          Abstruc <- 0
          spAr <- rep(log(Lambda.start[2]), sum(p*ncol(spdr)))
          if(!is.null(fit$fitstart$Ab)){
          spAr <- log(sqrt(unlist(lapply(fit$fitstart$Ab,diag))))
          fit <- fit[names(fit)!="fitstart"]
          fit$Ab <- exp(spAr)
          }
          if(sp.Ar.struc == "blockdiagonal" && Ab.diag.iter == 0){
            spAr<-c(spAr, rep(1e-3, p*ncol(spdr)*(ncol(spdr)-1)/2))
          }
        }else if(sp.Ar.struc == "MNdiagonal" || sp.Ar.struc == "MNunstructured"  || (sp.Ar.struc=="diagonalCL2" && Ab.diag.iter == 1)  || (sp.Ar.struc=="CL1" && Ab.diag.iter == 1) || (sp.Ar.struc=="CL2" && Ab.diag.iter == 1)){
          Abstruc <- 1
          #matrix normal VA matrix
          spAr <- rep(log(Lambda.start[2]), ncol(spdr)+p-1)
          if(!is.null(fit$fitstart$Ab)){
            spAr <- log(sqrt(unlist(lapply(fit$fitstart$Ab,diag))[c(1:ncol(spdr),(ncol(spdr)+2):(ncol(spdr)+p))]))
            fit <- fit[names(fit)!="fitstart"]
            fit$Ab <- exp(spAr)
          }
          if(sp.Ar.struc == "MNunstructured" && Ab.diag.iter == 0){
              spAr<-c(spAr, c(rep(1e-2, ncol(spdr)*(ncol(spdr)-1)/2)))
          }
          spAr <- c(spAr, rep(1e-3, sum(blocksp*Abranks-Abranks*(Abranks+1)/2)))
        }else if(sp.Ar.struc %in%c("CL1","diagonalCL1")){
          Abstruc <- 3
          spAr <- rep(log(Lambda.start[2]), p*ncol(spdr)+p-length(blocksp))# variances
          if(!is.null(fit$fitstart$Ab)){
            spAr <- log(sqrt(unlist(lapply(fit$fitstart$Ab,diag))))
            fit <- fit[names(fit)!="fitstart"]
            fit$Ab <- exp(spAr)
          }
          if(sp.Ar.struc=="CL1" && Ab.diag.iter == 0)spAr <- c(spAr, rep(1e-3, p*ncol(spdr)*(ncol(spdr)-1)/2)) # rest blockdiagonal
          spAr <- c(spAr, rep(1e-3, sum(blocksp*Abranks-Abranks*(Abranks+1)/2))) # rest p*p
        }else if(sp.Ar.struc %in%c("CL2")){
          Abstruc <- 4
          spAr <- rep(log(Lambda.start[2]), ncol(spdr)-1)# variances
          if(!is.null(fit$fitstart$Ab)){
            spAr2 <- log(sqrt(unlist(lapply(fit$fitstart$Ab,diag))))
            # need to reorder these. CL2 returns in covariate-wise for all species
            # and in c++ it is per block
            for(i in 1:length(blocksp)){
              spAr <- c(spAr, spAr2[rep(c(1,blocksp[-length(blocksp)]+1)[i]:cumsum(blocksp)[i], ncol(spdr))+rep(rep(p,ncol(spdr))*(0:(ncol(spdr)-1)),each=blocksp[i])])
            }
            fit <- fit[names(fit)!="fitstart"]
            fit$Ab <- exp(spAr)
          }else{
            spAr <- c(spAr, rep(log(Lambda.start[2]), sum(p*ncol(spdr))))
          }
          if(sp.Ar.struc=="CL2" && Ab.diag.iter == 0)spAr <- c(spAr, rep(1e-3, ncol(spdr)*(ncol(spdr)-1)/2)) # rest blockdiagonal
          spAr <- c(spAr, rep(1e-3, sum(blocksp*Abranks*ncol(spdr)-ncol(spdr)*Abranks*(Abranks+1)/2))) # rest p*p
        }else if(sp.Ar.struc == "diagonalCL2" ||  (sp.Ar.struc=="unstructured" && Ab.diag.iter==1)){
          Abstruc <- 2
          spAr <- rep(log(Lambda.start[2]), p*ncol(spdr))
          if(!is.null(fit$fitstart$Ab)){
            spAr <- log(sqrt(unlist(lapply(fit$fitstart$Ab,diag))))
            fit <- fit[names(fit)!="fitstart"]
            fit$Ab <- exp(spAr)
          }
          spAr <- c(spAr, rep(1e-3, ncol(spdr)*sum(blocksp*Abranks-Abranks*(Abranks+1)/2)))
      }else if(sp.Ar.struc == "unstructured"){
            Abstruc <- 5
            spAr <- rep(log(Lambda.start[2]), p*ncol(spdr))
            if(!is.null(fit$fitstart$Ab)){
              spAr <- log(sqrt(unlist(lapply(fit$fitstart$Ab,diag))))
              fit <- fit[names(fit)!="fitstart"]
              fit$Ab <- exp(spAr)
            }
            spAr <- c(spAr,rep(1e-3, sum(ncol(spdr)*blocksp*Abranks-Abranks*(Abranks+1)/2)))
      }
          } else {spAr <- 0}
      # Variational covariances for  random rows
      if(nrow(dr)==n){
        lg_Ar <- rep(log(Lambda.start[2]), sum(nr))
        
        if(Ar.struc!="diagonal" && diag.iter == 0){
          lg_Ar<-c(lg_Ar, rep(1e-3, sum(nr*(nr-1)/2)))
        }
      } else {lg_Ar <- 0}

      #quadratic model starting values
      if(quadratic == TRUE && start.struc == "LV"){
        start.fit <- try(gllvm.TMB(y=y, X=X, lv.X = lv.X, num.lv=num.lv, num.lv.c = num.lv.c, num.RR = num.RR, family = family, Lambda.struc = Lambda.struc, reltol=reltol, maxit = maxit, start.lvs = start.lvs, offset = offset, n.init = 1, diag.iter=diag.iter, dependent.row=dependent.row, quadratic="LV", starting.val = starting.val, Lambda.start = Lambda.start, quad.start = quad.start, jitter.var = jitter.var, zeta.struc = zeta.struc, optimizer = optimizer, optim.method = optim.method, max.iter=max.iter, start.struc="all", disp.group = disp.group, randomB = randomB, Ntrials = Ntrials),silent=T)
        if(inherits(start.fit,"try-error")&starting.val!="zero"){
          start.fit <- try(gllvm.TMB(y=y, X=X, lv.X = lv.X, num.lv=num.lv, num.lv.c = num.lv.c, num.RR = num.RR, family = family, Lambda.struc = Lambda.struc, reltol=reltol, maxit = maxit, start.lvs = start.lvs, offset = offset, n.init = 1, diag.iter=diag.iter, dependent.row=dependent.row, quadratic="LV", starting.val = "zero", Lambda.start = Lambda.start, quad.start = quad.start, jitter.var = jitter.var, zeta.struc = zeta.struc, optimizer = optimizer, optim.method = optim.method, max.iter=max.iter, start.struc="all", disp.group = disp.group, randomB = randomB, Ntrials = Ntrials),silent=T)
        }
        if(!inherits(start.fit,"try-error")&starting.val!="zero"){
          if(is.null(start.fit$lvs)){
            start.fit <- try(gllvm.TMB(y=y, X=X, lv.X = lv.X, num.lv=num.lv, num.lv.c = num.lv.c, num.RR = num.RR, family = family, Lambda.struc = Lambda.struc, reltol=reltol, maxit = maxit, start.lvs = start.lvs, offset = offset, n.init = 1, diag.iter=diag.iter, dependent.row=dependent.row, quadratic="LV", starting.val = "zero", Lambda.start = Lambda.start, quad.start = quad.start, jitter.var = jitter.var, zeta.struc = zeta.struc, optimizer = optimizer, optim.method = optim.method, max.iter=max.iter, start.struc="all", disp.group = disp.group, randomB = randomB, Ntrials = Ntrials),silent=T)
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
      
      
      ### Set up parameter.list, data.list and map.list
      
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
        if(!corWithinLV) {
          if(nrow(u) != nu){
            u=as.matrix((Matrix::t(dLV)%*%u/Matrix::colSums(dLV))[1:nu,, drop=FALSE])
          }
        }
      }
      if(num.RR==0 && num.lv.c==0) map.list$b_lv = factor(NA)
      
      ## Row effect settings

      if(nrow(dr)==n){
        if(any(cstrucn==4)){
          map.list$log_sigma <- if(cstrucn[1] == 0){1}else if(cstrucn[1]==4){c(1:2,NA)}else{1:2}
          if(length(cstrucn)>1){
            for(i in 2:length(cstrucn)){
              map.list$log_sigma <- c(map.list$log_sigma, if(!cstrucn[i]%in%c(0,4)){
                c(max(map.list$log_sigma, na.rm = TRUE)+1, max(map.list$log_sigma, na.rm = TRUE)+2)
              }else if(cstrucn[i] %in% c(0)){
                max(map.list$log_sigma, na.rm = TRUE)+1
              }else if(cstrucn[i]==4){
                c(max(map.list$log_sigma, na.rm = TRUE)+1, max(map.list$log_sigma, na.rm = TRUE)+2,NA)
              })
            }
          }
          map.list$log_sigma <- factor(map.list$log_sigma)
        }
        sigmanew <- NULL
        iter = 1 # keep track of # spatial structures
        for(re in cstrucn){
          if(re %in% c(1,3)) {
            sigmanew = c(sigmanew, log(sigma[1]),0)
          } else if(re %in% c(2)){
            sigmanew = c(sigmanew, log(sigma[1]),scaledc[[iter]])
            iter <- iter + 1
          } else if(re %in% c(4)){
            sigmanew = c(sigmanew, log(sigma[1]),scaledc[[iter]])
            iter <- iter + 1
            # Fix matern smoothness by default
            sigmanew = c(sigmanew, log(MaternKappa))
          } else {
            sigmanew = c(sigmanew, log(sigma[1]))
          }
        }
        sigma <- sigmanew
      } else {
        sigma=0
        map.list$log_sigma <- factor(NA)
        map.list$lg_Ar <- factor(NA)
      }
      
      if(quadratic == FALSE){
        # if(num.RR==0&quadratic==F) map.list$lambda2 = factor(NA)
        map.list$lambda2 = factor(NA)
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
      if(family == "tweedie"){ familyn=5}
      if(family == "ZIP"){familyn=6}
      if(family == "ordinal") {
        familyn=7
        if(link=="probit")extra[1]=1
        }
      if(family == "exponential") {familyn=8}
      if(family == "beta"){ 
        familyn=9
        if(link=="probit") extra[1]=1
      }
      if(family == "ZINB"){familyn=11}
      if(family == "orderedBeta") {familyn=12}
      
      if(family == "betaH"){ # EVA
        familyn = 10
        if(link=="probit") extra[1]=1
      }
      
      ## generate starting values quadratic coefficients in some cases
      if(starting.val!="zero" && quadratic != FALSE && (num.lv+num.lv.c+num.RR)>0){
        data.list = list(y = y, x = Xd, x_lv = lv.X, xr = xr, dr0 = dr, dLV = dLV, colMatBlocksI = blocks, Abranks = Abranks, Abstruc = Abstruc, xb = spdr, cs = cs, offset=offset, nr = nr, num_lv = num.lv, num_lv_c = num.lv.c, num_RR = num.RR, num_corlv=num.lv.cor, quadratic = 1, randomB = as.integer(randomB=="LV"), family=familyn,extra=extra,method=switch(method, VA=0, EVA=2),model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), times = times, cstruc=cstrucn, cstruclv = cstruclvn, cstruclv = cstruclvn, dc=dist, dc_lv = distLV, Astruc=Astruc, NN = NN, Ntrials = Ntrials, nncolMat = nncolMat, csb_lv = csBlv)
        
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
        map.list2$r0r = factor(rep(NA, length(r0r)))
        map.list2$r0f = factor(rep(NA, length(r0f)))
        map.list2$b_lv = factor(rep(NA, length(b.lv)))
        map.list2$Ab_lv = factor(rep(NA, length(Ab_lv)))
        map.list2$sigmab_lv = factor(rep(NA, length(sigmab_lv)))
        map.list2$b = factor(rep(NA, length(rbind(a, b))))
        map.list2$B = factor(rep(NA, length(B)))
        map.list2$Br = factor(rep(NA,length(Br)))
        #map.list2$lambda = factor(rep(NA, length(lambda)))
        map.list2$u = factor(rep(NA, length(u)))
        map.list2$lg_phi = factor(rep(NA, p))
        map.list2$lg_phiZINB = factor(rep(NA, p))
        map.list2$log_sigma = factor(rep(NA, length(sigma)))
        map.list2$sigmaB = factor(rep(NA, length(sigmaB)))
        map.list2$sigmaij = factor(rep(NA,length(sigmaij)))
        map.list2$Au = factor(rep(NA, length(Au)))
        map.list2$zeta = factor(rep(NA, length(zeta)))
        map.list2$lg_Ar = factor(rep(NA, length(lg_Ar)))
        
        parameter.list = list(r0r = matrix(r0r), r0f = matrix(r0f), b = rbind(a,b), b_lv = b.lv, sigmab_lv = sigmab_lv, Ab_lv = Ab_lv, B = B, Br=Br,lambda = lambda, lambda2 = t(lambda2), sigmaLV = (sigma.lv), u = u,lg_phi=log(phi),sigmaij=sigmaij,log_sigma=sigma, sigmaB = sigmaB, rho_lvc=rho_lvc, Au=Au, lg_Ar =lg_Ar, Abb = spAr, zeta=zeta, ePower = ePower, lg_phiZINB = log(ZINBphi)) #, scaledc=scaledc, thetaH = thetaH, bH=bH
        
        objr <- TMB::MakeADFun(
          data = data.list, silent=TRUE,
          parameters = parameter.list, map = map.list2,
          DLL = "gllvm")##GLLVM
        
        
        if(optimizer=="nlminb") {
          timeo <- system.time(optr <- try(suppressWarnings(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit))),silent = TRUE), gcFirst = FALSE)
        }
        if(optimizer=="optim" || !(optimizer %in%c("optim","nlminb") )) {
          if(optimizer == "optim" && optim.method != "BFGS")
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
          else
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
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
      
      
      ### Set up data and parameters
      
      data.list <- list(y = y, x = Xd, x_lv = lv.X , xr = xr, dr0 = dr, dLV = dLV, colMatBlocksI = blocks, Abranks = Abranks, Abstruc = Abstruc, xb = spdr, cs =  cs, offset=offset, nr = nr, num_lv = num.lv, num_lv_c = num.lv.c, num_RR = num.RR, num_corlv=num.lv.cor, quadratic = ifelse(quadratic!=FALSE,1,0), randomB = as.integer(randomB=="LV"), family=familyn,extra=extra,method=switch(method, VA=0, EVA=2),model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), times = times, cstruc=cstrucn, cstruclv = cstruclvn, dc=dist, dc_lv = distLV, Astruc=Astruc, NN = NN, Ntrials = Ntrials, nncolMat = nncolMat, csb_lv = csBlv)

      parameter.list <- list(r0r = matrix(r0r), r0f = matrix(r0f), b = rbind(a,b), sigmaB = sigmaB, Abb = spAr, b_lv = b.lv, sigmab_lv = sigmab_lv, Ab_lv = Ab_lv, B = B, Br=Br,lambda = lambda, lambda2 = t(lambda2), sigmaLV = (sigma.lv), u = u,lg_phi=log(phi),sigmaij=sigmaij,log_sigma=sigma, rho_lvc=rho_lvc, Au=Au, lg_Ar=lg_Ar, zeta=zeta, ePower = ePower, lg_phiZINB = log(ZINBphi)) #, scaledc=scaledc,thetaH = thetaH, bH=bH

      #### Call makeADFun
      objr <- TMB::MakeADFun(
        data = data.list, silent=TRUE,
        parameters = parameter.list, map = map.list,
        DLL = "gllvm")##GLLVM
 
      #### Fit model 
      if((num.lv.c+num.RR)<=1|!isFALSE(randomB)){
        if(optimizer=="nlminb") {
          timeo <- system.time(optr <- try(suppressWarnings(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit))),silent = TRUE), gcFirst = FALSE)
        }
        if(optimizer=="optim") {
          if(optim.method != "BFGS")# Due the memory issues, "BFGS" should not be used for Tweedie
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
          else
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
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
          suppressWarnings(timeo <- system.time(optr <- try(auglag(objr$par, objr$fn, objr$gr, heq = eval_eq_c, heq.jac = eval_eq_j, control.optim=control.optim, control.outer = list(eps = reltol.c, itmax=maxit, trace = FALSE, kkt2.check = FALSE, method = optim.method), obj = objr),silent = TRUE), gcFirst = FALSE))
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
          timeo <- system.time(optr <- try(nloptr(x0 = objr$par, eval_f=eval_f, eval_g_eq=eval_g_eq, opts=opts, obj = objr),silent = TRUE), gcFirst = FALSE)
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
      if((diag.iter>0) && (!(Lambda.struc %in% c("diagonal", "diagU")) && (((nlvr+randoml[3]*num.RR)>1) | (num.lv.cor>0)) && !inherits(optr,"try-error") | ((nrow(dr)==n) & Ar.struc=="unstructured")) | ((Ab.diag.iter>0) && (col.eff=="random" && sp.Ar.struc%in%c("blockdiagonal","MNunstructured","unstructured","diagonalCL2","CL1","CL2")))){
        objr1 <- objr
        optr1 <- optr
        param1 <- optr$par
        nam <- names(param1)
        if(length(param1[nam=="r0r"])>0){ r0r1 <- matrix(param1[nam=="r0r"])} else {r0r1 <- matrix(0)}
        if(length(param1[nam=="r0f"])>0){ r0f1 <- matrix(param1[nam=="r0f"])} else {r0f1 <- matrix(0)}
        if(nrow(dr)==n){
          log_sigma1 <- ifelse(param1[nam=="log_sigma"]==0,1e-3,param1[nam=="log_sigma"])
          if(!is.null(map.list$log_sigma)) {
            # We need to maintain the fixed parameter for Matern smoothness
            # Which is omitted in the optimiser
            log_sigma <- sigma
            log_sigma <- log_sigma1[map.list$log_sigma[!is.na(map.list$log_sigma)]]
            log_sigma1 <- log_sigma
          }
          lg_Ar<- log(exp(param1[nam=="lg_Ar"][1:sum(nr)])+1e-3)
          if(Ar.struc=="unstructured"){
            lg_Ar <- c(lg_Ar, rep(1e-3, sum(nr*(nr-1)/2)))
          }
        } else {log_sigma1 = 0}
      
        if(beta0com){
          b1 <- matrix(paramb1[bi][model2$TMBfn$env$map$b], num.X+1, p)
        }else{
          b1 <- matrix(param1[nam=="b"],num.X+1,p)  
        }
        
        if(!isFALSE(randomB)){
          if(randomB!="iid"){
            sigmab_lv1 <- param1[nam=="sigmab_lv"]
          }else{
            sigmab_lv1 <- 0
          }
          Ab_lv1<- c(pmax(param1[nam=="Ab_lv"],rep(log(1e-6), ab12*ab3)), rep(0.01,ab12*(ab12-1)/2*ab3))
        }else{
          sigmab_lv1<-0
          Ab_lv1 <- 0
        }
        if(col.eff == 'random'){
          if(nrow(Xt)==n){
          B1 <- matrix(0, ncol = 1, nrow = ncol (spdr))
          B1[!is.na(map.list$B)] <- param1[names(param1)=="B"]
          }else{
            B1 <- matrix(0)
          }
          Br1 <- matrix(param1[nam=="Br"],nrow=ncol(spdr))
          sigmaB1 <- ifelse(round(param1[nam=="sigmaB"],8)==0,1e-3,param1[nam=="sigmaB"])
          # if(!is.null(colMat))param1[nam=="sigmaB"][length(param1[nam=="sigmaB"])] <- ifelse(tail(param1[nam=="sigmaB"], ifelse(colMat.rho.struct == "single",1,ncol(spdr)))>0.5,0.5,tail(param1[nam=="sigmaB"], ifelse(colMat.rho.struct == "single",1,ncol(spdr))))
          spAr1<- log(exp(param1[nam=="Abb"])+1e-3)
          if(Ab.diag.iter>0){
          if(sp.Ar.struc=="blockdiagonal"){# "diagonal" was previous iteration
           spAr1 <- c(spAr1, rep(1e-3, p*ncol(spdr)*(ncol(spdr)-1)/2))
          }else if(sp.Ar.struc == "MNunstructured"){# "MNdiagonal" was previous iteration
              spAr1 <- c(spAr1[1:(ncol(spdr)+p-1)], c(rep(1e-2, ncol(spdr)*(ncol(spdr)-1)/2), rep(1e-3, sum(blocksp*Abranks-Abranks*(Abranks-1)/2-Abranks))))
            }else if(sp.Ar.struc == "diagonalCL2"){# "MNdiagonal" was previous iteration
              Abstruc <- 2
              spAr1 <- c(rep(log(Lambda.start[2]), ncol(spdr)*p),rep(1e-3, ncol(spdr)*sum(blocksp*Abranks-Abranks*(Abranks-1)/2-Abranks)))
            }else if(sp.Ar.struc == "CL1"){# "MNdiagonal" was previous iteration
              Abstruc <- 3
              spAr1 <- c(rep(log(Lambda.start[2]), p*ncol(spdr)+p-length(blocksp)), rep(1e-3, p*ncol(spdr)*(ncol(spdr)-1)/2),rep(1e-3, sum(blocksp*Abranks-Abranks*(Abranks-1)/2-Abranks))) # rest blockdiagonal
            }else if(sp.Ar.struc == "CL2"){# "MNdiagonal" was previous iteration
              Abstruc <- 4
              spAr1 <- c(rep(log(Lambda.start[2]), ncol(spdr)-1+p*ncol(spdr)), rep(1e-3, ncol(spdr)*(ncol(spdr)-1)/2), rep(1e-3, sum(blocksp*Abranks*ncol(spdr)-ncol(spdr)*Abranks-ncol(spdr)*Abranks*(Abranks-1)/2)))
            }else if(sp.Ar.struc == "unstructured"){# "diagonalCL2" was previous iteration
              Abstruc <- 5
              spAr1 <- c(spAr1[1:(ncol(spdr)*p)],rep(1e-3, sum(ncol(spdr)*blocksp*Abranks-Abranks*(Abranks-1)/2-Abranks)))
            }
          }
        }else{
          B1 <- matrix(0)
          Br1 <- matrix(0)
          sigmaB1 <- 0
          spAr1 <- 0
        }
        
        if((num.lv.c+num.RR)>0){b.lv1 <- matrix(param1[nam=="b_lv"],ncol(lv.X),(num.lv.c+num.RR))}else{b.lv1<-matrix(0)}
        if((num.lv+num.lv.c+num.RR+num.lv.cor)>0){lambda1 <- param1[nam=="lambda"]}else{lambda1<-lambda}
        if (quadratic=="LV" | isTRUE(quadratic) && start.struc == "LV"){
          lambda2 <- matrix(param1[nam == "lambda2"], byrow = TRUE, ncol = num.lv+(num.lv.c+num.RR), nrow = 1)#In this scenario we have estimated two quadratic coefficients before
        }else if(isTRUE(quadratic)){
          lambda2 <- matrix(param1[nam == "lambda2"], byrow = TRUE, ncol = num.lv+(num.lv.c+num.RR), nrow = p)
        }
        
        if((num.lv+num.lv.c)>0){sigma.lv1 <- param1[nam=="sigmaLV"]}else{sigma.lv1<-0}
        if((num.lv+num.lv.c)>0){u1 <- matrix(param1[nam=="u"],nrow(u),num.lv+num.lv.c)}else{u1<-u}
        if(family %in% c("poisson","binomial","ordinal","exponential", "betaH", "orderedBeta")){ lg_phi1 <- log(phi)} else {lg_phi1 <- param1[nam=="lg_phi"][disp.group]} #cat(range(exp(param1[nam=="lg_phi"])),"\n")
        if(family=="ZINB"){lg_phiZINB1 <- param1[nam=="lg_phiZINB"][map.list$lg_phiZINB]}else{lg_phiZINB1<-log(ZINBphi)}
        if(family=="tweedie" && is.null(Power)) ePower = param1[nam == "ePower"]
        sigmaij1 <- param1[nam=="sigmaij"]

        if(num.lv.cor>0){
          Au1<- c(param1[nam=="Au"])
          if(corWithinLV) {
            if(Lambda.struc == "unstructured" && Astruc==1) {
              Au1 <- c(pmax(Au1[1:(n*num.lv.cor)],log(1e-2)), rep(1e-3,sum(lower.tri(matrix(0,n,n)))*num.lv.cor) )
            } else if(Lambda.struc == "bdNN" && Astruc==2){
              Au1 <- c(pmax(Au1[1:(n*num.lv.cor)],log(1e-2)), rep(1e-3,nrow(NN)*num.lv.cor*n) )
            } else if(Astruc==3) {
              Au1 <- c(log(exp(Au1[1:(n)])+1e-2), rep(1e-3,sum(lower.tri(matrix(0,n,n)))), Au1[-(1:n)])
            } else if(Astruc==4) {
              Au1 <- c(log(exp(Au1[1:(n)])+1e-2), rep(1e-3,nrow(NN)*nu), Au1[-(1:n)])
            } 
            # else if(Astruc==6){
            #   Au1 <- c(pmax(Au1[1:(n*num.lv.cor)],log(1e-2)), rep(1e-3,NN[1]*num.lv.cor*n) )
            # }
          } else {
            if(Lambda.struc == "unstructured" && Astruc==1 & cstruclvn==0){
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
            # else if(Astruc==6){
            #   Au1 <- c(pmax(Au1[1:(nu*num.lv.cor)],log(1e-2)), rep(1e-3,NN[1]*num.lv.cor*nu) )
            # }
          }
          if(cstruclvn>0){
            if(cstruclvn %in% c(2,4)){ #cstruc=="corExp" || cstruc=="corMatern"
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
        
        if(family %in% c("ordinal")){
          zeta <- param1[nam=="zeta"] 
        } else if(family %in% c("orderedBeta")){
          zeta <- matrix((param1[names(param1)=="zeta"])[map.list$zeta],p,2)
          if(any(is.na(map.list$zeta))) zeta[is.na(map.list$zeta)] = attr(objr1$env$parameters$zeta, "shape")[is.na(map.list$zeta)]
        } else {
          zeta <- 0 
        }
        
        #Because then there is no next iteration
        data.list = list(y = y, x = Xd,  x_lv = lv.X, xr = xr, dr0 = dr, dLV = dLV, colMatBlocksI = blocks, Abranks = Abranks, Abstruc = Abstruc, xb = spdr, cs = cs, offset=offset, nr = nr, num_lv = num.lv, num_lv_c = num.lv.c, num_RR = num.RR, num_corlv=num.lv.cor, quadratic = ifelse(quadratic!=FALSE,1,0), randomB = as.integer(randomB=="LV"), family=familyn,extra=extra,method=switch(method, VA=0, EVA=2),model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), times = times, cstruc=cstrucn, cstruclv = cstruclvn, dc=dist, dc_lv = distLV, Astruc=Astruc, NN = NN, Ntrials = Ntrials, nncolMat = nncolMat, csb_lv = csBlv)
        
        parameter.list <- list(r0r = r0r1, r0f = r0f1, b = b1, b_lv = b.lv1, sigmaB = sigmaB1, Abb = spAr1, sigmab_lv = sigmab_lv1, Ab_lv = Ab_lv1, B = B1, Br=Br1,lambda = lambda1, lambda2 = t(lambda2), sigmaLV = sigma.lv1, u = u1,lg_phi=lg_phi1,sigmaij=sigmaij,log_sigma=log_sigma1, rho_lvc=rho_lvc, Au=Au1, lg_Ar=lg_Ar, zeta=zeta, ePower = ePower, lg_phiZINB = lg_phiZINB1) #, scaledc=scaledc,thetaH = thetaH, bH=bH
        
        objr <- TMB::MakeADFun(
          data = data.list, silent=TRUE,
          parameters = parameter.list, map = map.list,
          DLL = "gllvm")
        
        if((num.lv.c+num.RR)<=1|!isFALSE(randomB)){
          if(optimizer=="nlminb") {
            timeo <- system.time(optr <- try(suppressWarnings(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit))),silent = TRUE), gcFirst = FALSE)
          }
          if(optimizer=="optim") {
            if(optim.method != "BFGS")
              timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
            else
              timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
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
            suppressWarnings(timeo <- system.time(optr <- try(auglag(objr$par, objr$fn, objr$gr, heq = eval_eq_c, heq.jac = eval_eq_j, control.optim=control.optim, control.outer = list(eps = reltol.c, itmax=maxit, trace = FALSE, kkt2.check = FALSE, method = optim.method), obj = objr),silent = TRUE), gcFirst = FALSE))
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
            timeo <- system.time(optr <- try(nloptr(x0 = objr$par, eval_f=eval_f, eval_g_eq=eval_g_eq, opts=opts, obj = objr),silent = TRUE), gcFirst = FALSE)
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
      
      param <- objr$env$last.par.best
      if(family %in% c("negative.binomial", "tweedie", "gaussian", "gamma", "beta", "betaH", "orderedBeta","ZIP","ZINB")) {
        phis <- exp(param[names(param)=="lg_phi"])[disp.group]
        if(family == "ZINB")ZINBphis <- exp(param[names(param)=="lg_phiZINB"])[map.list$lg_phiZINB]
        if(family %in% c("ZIP","ZINB")) {
          lp0 <- param[names(param)=="lg_phi"][disp.group]; out$lp0 <- lp0
          phis <- exp(lp0)/(1+exp(lp0));
        }
        if(family=="tweedie" && is.null(Power)){
          Power = exp(param[names(param)=="ePower"])/(1+exp(param[names(param)=="ePower"]))+1
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
        out$zeta.struc = zeta.struc
      }
      if(family == "orderedBeta"){
        zetas <- matrix((param[names(param)=="zeta"])[map.list$zeta],p,2)
        if(any(is.na(map.list$zeta))) zetas[is.na(map.list$zeta)] = attr(objr$env$parameters$zeta, "shape")[is.na(map.list$zeta)]
        zetas[,2] = exp(zetas[,2])
        colnames(zetas) = c("cutoff0","cutoff1")
      }
      
      if((num.lv.c+num.RR)>0){
        bi.lv <- names(param)=="b_lv"
        if(!isFALSE(randomB))sib <- names(param)=="sigmab_lv"
      }
      bi <- names(param)=="b"
      li <- names(param)=="lambda"
      si <- names(param) == "sigmaLV"
      li2 <- names(param)=="lambda2"
      ui <- names(param)=="u"
      
      if(num.lv.cor > 0){ # Correlated latent variables
        if(corWithinLV){
          lvs<-(matrix(param[ui],n,num.lv.cor))
        } else {
          lvs = matrix(param[ui],nu,num.lv.cor)
          rownames(lvs) =colnames(dLV)
          # lvs = dLV%*%matrix(param[ui],nu,num.lv.cor)
        }
        sigma.lv <- (param[si])
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
        if((cstruclvn %in% c(1,3))) rho.lv<- param[names(param)=="rho_lvc"] / sqrt(1.0 + param[names(param)=="rho_lvc"]^2);
        if((cstruclvn %in% c(2,4))) {
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
      
      if((nrow(dr) ==n) || (nrow(xr) == n)) {
        rir = names(param)=="r0r"
        rif = names(param)=="r0f"
        if(nrow(dr)==n){
          row.params.random <- param[rir]
          sigma = param[names(param)=="log_sigma"]
          # if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,3))) rho = param[names(param)=="log_sigma"][2] / sqrt(1.0 + param[names(param)=="log_sigma"][2]^2);
          # if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(2,4))) {
          #   rho = exp(param[names(param)=="log_sigma"][-1]);
          #   # scaledc<- exp(param[names(param)=="scaledc"]);
          # }
          # if(num.lv>0 && dependent.row && rstruc==0) sigma = c(sigma,(param[names(param)=="log_sigma"])[-1])
        }
        if(nrow(xr)==n){
          row.params.fixed <- param[rif]
        }
      }
      
      if(col.eff=="random"){
        sigma.sp = exp(param[names(param)=="sigmaB"])[1:ncol(spdr)]
        covsigma.sp  = param[names(param)=="sigmaB"][-c(1:ncol(spdr))]
        if(any(colMat[row(colMat)!=col(colMat)]!=0)){
          rho.sp = exp(-exp(tail(param[names(param)=="sigmaB"], ifelse(colMat.rho.struct=="single",1,ncol(spdr)))))
          if(nrow(nncolMat)<p)rho.sp = pmax(rho.sp, 1e-12)
          covsigma.sp  = head(covsigma.sp, -ifelse(colMat.rho.struct=="single",1,ncol(spdr)))
        }
        Bri = names(param)=="Br"
        Br = matrix(param[Bri], nrow = ncol(spdr))#c(0,param[ri])
        if(nrow(Xt)==n)B = param[names(param)=="B"]
      }
      
      paramb <- param[bi]
      if(beta0com)paramb <- param[bi][model2$TMBfn$env$map$b]
      betaM <- matrix(paramb,p,num.X+1,byrow=TRUE)

      beta0 <- betaM[,1]
      if(!is.null(X)) betas <- betaM[,-1]
      if((num.lv.c+num.RR)>0)b.lv <- matrix(param[bi.lv],ncol(lv.X),num.lv.c+num.RR)
      # if(family %in% "betaH"){
      #   bHi <- names(param)=="bH"
      #   betaH <- matrix(param[bHi],p,num.X+1,byrow=TRUE)
      # }
      if(!isFALSE(randomB)&(num.lv.c+num.RR)>0&randomB!="iid")sigmab_lv <- param[sib]
      new.loglik <- objr$env$value.best[1]
      
    } else if(method=="LA" || (nlvr==0 && (method %in% c("VA", "EVA")) && (nrow(dr)!=n) && isFALSE(randomB) && isFALSE(col.eff))){
      ## Laplace method / nlvr==0
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
      if(family == "tweedie"){ familyn=5}
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
        
        # bH <- rbind(a,b)
        # if(num.lv>0) {
        #   mapLH<-factor(1:length(thetaH))
        #   mapLH[lower.tri(thetaH)] <- NA
        #   map.list$thetaH <- factor(mapLH)
        # } else {
        #   thetaH<- matrix(0); 
        #   map.list$thetaH = factor(NA)
        # }
      } 
      # else {
      #   thetaH<- matrix(0)
      #   map.list$thetaH = factor(NA)
      #   bH <- matrix(0)
      #   map.list$bH = factor(NA)
      # }
      if(family == "ZINB"){familyn=11}
      if(family == "orderedBeta") {familyn=12}
      
      
      ## generate starting values quadratic coefficients in some cases
      if(starting.val!="zero" && quadratic == TRUE && num.RR>0&(num.lv+num.lv.c)==0 && start.struc=="LV"){
        data.list = list(y = y, x = Xd, x_lv = lv.X, xr = xr, dr0 = dr, dLV = dLV, colMatBlocksI = blocks, Abranks = Abranks, Abstruc = 0, xb = spdr, cs = cs, offset=offset, nr = nr, num_lv = num.lv, num_lv_c = num.lv.c, num_RR = num.RR, num_corlv=num.lv.cor, quadratic = 1, randomB = as.integer(randomB=="LV"), family=familyn,extra=extra,method=switch(method, VA=0, EVA=2),model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), times = times, cstruc=cstrucn, cstruclv = cstruclvn, dc=dist, dc_lv = distLV, Astruc=Astruc, NN = NN, Ntrials = Ntrials, nncolMat = nncolMat, csb_lv = csBlv)
        
        map.list2 <- map.list 
        map.list2$log_sigma = factor(NA)
        map.list2$Au <- map.list2$Abb <- map.list2$Ab_lv <- factor(NA)# map.list$lambda2 <- 
        map.list2$b_lv <- factor(rep(NA,length(b.lv)))
        map.list2$B <- factor(rep(NA, length(B)))
        
        parameter.list = list(r0r = matrix(r0r), r0f = matrix(r0f), b = rbind(a,b), sigmaB = sigmaB, Abb = spAr, b_lv = b.lv, sigmab_lv = 0, Ab_lv = Ab_lv, B = B, Br=Br,lambda = lambda, lambda2 = t(lambda2), sigmaLV = (sigma.lv), u = u,lg_phi=log(phi),sigmaij=sigmaij,log_sigma=sigma,rho_lvc=rho_lvc, Au=0, lg_Ar =0, zeta=zeta, ePower = ePower, lg_phiZINB = log(ZINBphi)) #, scaledc=scaledc, thetaH = thetaH, bH=bH
        
        objr <- TMB::MakeADFun(
          data = data.list, silent=TRUE,
          parameters = parameter.list, map = map.list2,
          DLL = "gllvm")##GLLVM
        
        if(optimizer=="nlminb") {
          timeo <- system.time(optr <- try(suppressWarnings(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit))),silent = TRUE), gcFirst = FALSE)
        }
        if(optimizer=="optim" | !(optimizer %in% c("optim","nlminb"))) {
          if( optimizer == "optim" && optim.method != "BFGS")
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
          else
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
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
      data.list = list(y = y, x = Xd, x_lv = lv.X, xr = xr, dr0 = dr, dLV = dLV, colMatBlocksI = blocks, Abranks = Abranks, Abstruc = 0, xb = spdr, cs = cs, offset=offset, nr = nr, num_lv = num.lv, num_lv_c = num.lv.c, num_RR = num.RR, num_corlv=num.lv.cor, quadratic = ifelse(quadratic!=FALSE,1,0), randomB = as.integer(randomB=="LV"), family=familyn,extra=extra,method=1,model=0,random=randoml, zetastruc = ifelse(zeta.struc=="species",1,0), times = times, cstruc=cstrucn, cstruclv = cstruclvn, dc=dist, dc_lv = distLV, Astruc=Astruc, NN = NN, Ntrials = Ntrials, nncolMat = nncolMat, csb_lv = csBlv)
      
      if(family %in% c("ordinal", "orderedBeta")){
        data.list$method = 0
      }
      
      randomp <- "u"
      map.list$Au <- map.list$lg_Ar <- map.list$Abb <- factor(NA)
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
      
      if(num.lv.cor>0 & (!corWithinLV)){
        u <- as.matrix(u[1:nu,])
      }
      
      # Row parameter settings
      if(nrow(dr)==n){
        randoml[1] <- 1
        randomp <- c(randomp,"r0r")
        
        sigmanew <- NULL
        iter <- 1
        if(any(cstrucn==4)){
          map.list$log_sigma <- if(cstrucn[1]==0){1}else if(cstrucn[1]==4){c(1:2, NA)}else{1:2}
          if(length(cstrucn)>1){
         for(i in 2:length(cstrucn)){
           map.list$log_sigma <- c(map.list$log_sigma, if(!cstrucn[i]%in%c(0,4)){
            c(max(map.list$log_sigma, na.rm = TRUE)+1, max(map.list$log_sigma, na.rm = TRUE)+2)
             }else if(cstrucn[i]==0){
               max(map.list$log_sigma, na.rm = TRUE)+1
               }else if(cstrucn[i]==4){
                 c(max(map.list$log_sigma, na.rm = TRUE)+1, max(map.list$log_sigma, na.rm = TRUE)+2,NA)
                 })
         }
          }
          map.list$log_sigma <- factor(map.list$log_sigma)
        }
        for(re in cstrucn){
          if(re %in% c(1,3)) {
            sigmanew = c(sigmanew, log(sigma[1]),0)
          } else if(re %in% c(2)){
            sigmanew = c(sigmanew, log(sigma[1]),scaledc[[iter]])
            iter <- iter + 1
          } else if(re %in% c(4)){
            sigmanew = c(sigmanew, log(sigma[1]),scaledc[[iter]])
            iter <- iter + 1
            # Fix matern smoothness by default
            sigmanew = c(sigmanew, log(MaternKappa))
          } else {
            sigmanew = c(sigmanew, log(sigma[1]))
          }
        }
        sigma <- sigmanew
      } else {
        sigma=0
        map.list$log_sigma <- factor(NA)
        # if(row.eff != "fixed") map.list$r0 <- factor(rep(NA, length(r0)))
      }
      if(col.eff == "random"){
        randoml[4] <- 1
        randomp <- c(randomp,"Br")
      }
      
      if(!isFALSE(randomB)){
        randoml[3] <- 1
        randomp <- c(randomp, "b_lv")
      }
      if(isFALSE(randomB)){
        sigmab_lv <- 0
      }else if(randomB=="iid"){
        sigmab_lv <- 0
      }
      
      #### Set up data and parameters
      
      if(family == "ordinal"){ # || family == "orderedBeta"
        data.list$method = 0
      }
      
      parameter.list = list(r0r= matrix(r0r), r0f = matrix(r0f), b = rbind(a,b), sigmaB = sigmaB, b_lv = b.lv, sigmab_lv = sigmab_lv, Ab_lv = Ab_lv, B = B, Br=Br,lambda = lambda, lambda2 = t(lambda2), sigmaLV = (sigma.lv), u = u, lg_phi=log(phi),sigmaij=sigmaij,log_sigma=c(sigma), rho_lvc=rho_lvc, Au=0, lg_Ar=0, Abb=0, zeta=zeta, ePower = ePower, lg_phiZINB = log(ZINBphi)) #, scaledc=scaledc,thetaH = thetaH, bH=bH
      
      #### Call makeADFun
      objr <- TMB::MakeADFun(
        data = data.list, silent=!trace,
        parameters = parameter.list, map = map.list,
        inner.control=list(mgcmax = 1e+200,tol10=0.01),
        random = randomp, DLL = "gllvm")
      #### Fit model 
      
      # Not used for now
      # if(family=="ZIP" && FALSE) {
      #   m <- length(objr$par)
      #   low <- rep(-restrict,m); upp=rep(restrict,m);
      #   low[names(objr$par)=="lg_phi"]=0.0; upp[names(objr$par)=="lg_phi"]=1#0.99
      #   timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit),lower = low,upper = upp),silent = TRUE))
      # }
      if((num.lv.c+num.RR)<=1|!isFALSE(randomB)){
        if(optimizer=="nlminb") {
          timeo <- system.time(optr <- try(suppressWarnings(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit))),silent = TRUE), gcFirst = FALSE)
        }
        if(optimizer=="optim") {
          if(optim.method != "BFGS")
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
          else
            timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
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
          suppressWarnings(timeo <- system.time(optr <- try(auglag(objr$par, objr$fn, objr$gr, heq = eval_eq_c, heq.jac = eval_eq_j, control.optim=control.optim, control.outer = list(eps = reltol.c, itmax=maxit, trace = FALSE, kkt2.check = FALSE, method = optim.method), obj = objr),silent = TRUE), gcFirst = FALSE))
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
          timeo <- system.time(optr <- try(nloptr(x0 = objr$par, eval_f=eval_f, eval_g_eq=eval_g_eq, opts=opts, obj = objr),silent = TRUE), gcFirst = FALSE)
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
        if(!isFALSE(randomB) & randomB!="iid"){
          sigmab_lv <- objr$env$last.par.best[names(objr$env$last.par.best)=="sigmab_lv"]
        }else if(isFALSE(randomB)){
          sigmab_lv <- 0
        }else if(randomB=="iid"){
          sigmab_lv <- 0
        }
        b <- matrix(objr$env$last.par.best[names(objr$env$last.par.best)=="b"],num.X+1,p)

        if(!(family %in% c("poisson","binomial","ordinal","exponential"))) phi <- exp(objr$env$last.par.best[names(objr$env$last.par.best)=="lg_phi"])[disp.group]
        if(!(family %in% c("ZINB"))) ZINBphi <- exp(objr$env$last.par.best[names(objr$env$last.par.best)=="lg_phiZINB"])[map.list$lg_phiZINB]
        parameter.list = list(r0r = matrix(r0r), r0f = matrix(r0f), b = b, sigmaB = sigmaB, b_lv = b.lv, sigmab_lv = sigmab_lv, Ab_lv = Ab_lv, B = B, Br=Br,lambda = lambda, lambda2 = t(lambda2), sigmaLV = (sigma.lv), u = u, lg_phi=log(phi),sigmaij=sigmaij,log_sigma=c(sigma), rho_lvc=rho_lvc, Au=0, lg_Ar=0, Abb=0, zeta=zeta, ePower = ePower, lg_phiZINB = log(ZINBphi)) #, scaledc=scaledc,thetaH = thetaH, bH=bH
        
        #### Call makeADFun
        objr <- TMB::MakeADFun(
          data = data.list, silent=!trace,
          parameters = parameter.list, map = map.list,
          inner.control=list(mgcmax = 1e+200,tol10=0.01),
          random = randomp, DLL = "gllvm")
        
        if((num.lv.c+num.RR)<=1|!isFALSE(randomB)){
          if(optimizer=="nlminb") {
            timeo <- system.time(optr <- try(suppressWarnings(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol,iter.max=max.iter,eval.max=maxit))),silent = TRUE), gcFirst = FALSE)
          }
          if(optimizer=="optim") {
            if(optim.method != "BFGS")
              timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = optim.method,control = list(maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
            else
              timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit),hessian = FALSE),silent = TRUE), gcFirst = FALSE)
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
            suppressWarnings(timeo <- system.time(optr <- try(auglag(objr$par, objr$fn, objr$gr, heq = eval_eq_c, heq.jac = eval_eq_j, control.optim=control.optim, control.outer = list(eps = reltol.c, itmax=maxit, trace = FALSE, kkt2.check = FALSE, method = optim.method), obj = objr),silent = TRUE), gcFirst = FALSE))
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
            timeo <- system.time(optr <- try(nloptr(x0 = objr$par, eval_f=eval_f, eval_g_eq=eval_g_eq, opts=opts, obj = objr),silent = TRUE), gcFirst = FALSE)
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
        if(!isFALSE(randomB) & randomB!="iid")sib <- names(param)=="sigmab_lv"
      }
      bi <- names(param)=="b"
      li <- names(param)=="lambda"
      li2 <- names(param)=="lambda2"
      si <- names(param) == "sigmaLV"
      ui <- names(param)=="u"
      
      if(num.lv.cor > 0){
        if((num.lv.cor)>0 & corWithinLV){
          lvs<-(matrix(param[ui],n,num.lv.cor))
        } else{ 
          lvs = matrix(param[ui],nu,num.lv.cor)
          rownames(lvs) =colnames(dLV)
          # lvs = dLV%*%matrix(param[ui],nu,num.lv.cor)
        }
        
        theta <- matrix(0,p,num.lv.cor)
        sigma.lv <- abs(param[si])
        # sigma.lv <- exp(param[si])
        
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
        if((cstruclvn %in% c(1,3))) rho.lv<- param[names(param)=="rho_lvc"] / sqrt(1.0 + param[names(param)=="rho_lvc"]^2);
        if((cstruclvn %in% c(2,4))) {
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
      
      if((nrow(dr)==n) || (nrow(xr)==n)) {
        rir = names(param)=="r0r"
        rif = names(param)=="r0f"
        if(nrow(dr)==n){
          row.params.random <- param[rir]
          sigma = param[names(param)=="log_sigma"]
        }
        if(nrow(xr)==n){
          row.params.fixed <- param[rif]
        }
      }
      if(col.eff=="random"){
        sigma.sp = exp(param[names(param)=="sigmaB"])[1:ncol(spdr)]
        covsigma.sp  = param[names(param)=="sigmaB"][-c(1:ncol(spdr))]
        if(any(colMat[row(colMat)!=col(colMat)]!=0)){
          rho.sp = exp(-exp(tail(param[names(param)=="sigmaB"], ifelse(colMat.rho.struct=="single",1,ncol(spdr)))))
          if(nrow(nncolMat)<p)rho.sp = pmax(rho.sp, 1e-12)
          covsigma.sp  = head(covsigma.sp, -ifelse(colMat.rho.struct=="single",1,ncol(spdr)))
        }
        Bri = names(param)=="Br"
        Br = matrix(param[Bri], nrow = ncol(spdr))#c(0,param[ri])
        if(nrow(Xt)==n)B <- param[names(param)=="B"]
      }
      
      paramb <- param[bi]
      if(beta0com)paramb <- param[bi][model2$TMBfn$env$map$b]
      betaM <- matrix(paramb,p,num.X+1,byrow=TRUE)

      beta0 <- betaM[,1]
      if(!is.null(X)) betas=betaM[,-1]
      if((num.lv.c+num.RR)>0){
        b.lv <- matrix(param[bi.lv],ncol(lv.X),(num.lv.c+num.RR))
        if(!isFALSE(randomB) & randomB!="iid")sigmab_lv <- param[sib]
      }
      
      new.loglik <- objr$env$value.best[1]
      
      if(family %in% c("negative.binomial", "tweedie", "ZIP", "ZINB", "gaussian", "gamma", "beta", "betaH", "orderedBeta")) {
        phis <- exp(param[names(param)=="lg_phi"])[disp.group]
        if(family == "ZINB")ZINBphis <- exp(param[names(param)=="lg_phiZINB"])[map.list$lg_phiZINB]
        if(family %in% c("ZIP","ZINB")) {
          lp0 <- param[names(param)=="lg_phi"][disp.group]; out$lp0 <- lp0
          phis <- exp(lp0)/(1+exp(lp0));
        }
        if(family=="tweedie" && is.null(Power)){
          Power = exp(param[names(param)=="ePower"])/(1+exp(param[names(param)=="ePower"]))+1
        }
      }
      # if(family %in% "betaH"){
      #   bHi <- names(param)=="bH"
      #   betaH <- matrix(param[bHi],p,num.X+1,byrow=TRUE)
      #   if(num.lv>0) {
      #     thetaH[!is.na(map.list$thetaH)] <- param[names(param)=="thetaH"]
      #   }
      # }
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
      if(family == "orderedBeta"){
        zetas <- matrix((param[names(param)=="zeta"])[map.list$zeta],p,2)
        if(any(is.na(map.list$zeta))) zetas[is.na(map.list$zeta)] = attr(objr$env$parameters$zeta, "shape")[is.na(map.list$zeta)]
        zetas[,2] = exp(zetas[,2])
        colnames(zetas) = c("cutoff0","cutoff1")
      }
      
    }
    
    
    out$start <- fit
    

    if(!inherits(optr, "try-error")){
      objrFinal<-objr1 <- objr; optrFinal<-optr1<-optr;
      out$logL <- objrFinal$env$value.best[1]
      if((num.lv+(num.lv.c+num.RR)) > 0) {
        if((num.lv+num.lv.c)>0)out$lvs <- lvs
        out$params$theta <- theta
        if((num.lv+num.lv.c)>0) out$params$sigma.lv  <- sigma.lv
        if((num.lv.c+num.RR)>0){
          out$params$LvXcoef <- b.lv
          colnames(out$params$LvXcoef) <- paste("CLV",1:(num.lv.c+num.RR), sep="")
          row.names(out$params$LvXcoef) <- colnames(lv.X)
          if(!isFALSE(randomB)){
          if(randomB=="LV"|randomB=="single"){
            out$params$sigmaLvXcoef <- exp(head(sigmab_lv, num.lv.c+num.RR))
          }else if(randomB=="P"){
            out$params$sigmaLvXcoef <- exp(head(sigmab_lv, ncol(lv.X)+num.lv.c+num.RR-1))
            }else if(randomB=="iid"){out$params$sigmaLvXcoef <- 1}
            if(randomB=="LV")names(out$params$sigmaLvXcoef) <- paste("CLV",1:(num.lv.c+num.RR), sep="")
            if(randomB=="P")names(out$params$sigmaLvXcoef) <- c(colnames(lv.X), paste("CLV",2:(num.lv.c+num.RR), sep=""))
            # if(randomB=="all")names(out$params$sigmaLvXcoef) <- paste(paste("CLV",1:(num.lv.c+num.RR),sep=""),rep(colnames(lv.X),each=num.RR+num.lv.c),sep=".")
            if(randomB=="single")names(out$params$sigmaLvXcoef) <- NULL
            if(ncol(csBlv)==2){
              covsigmaB <- tail(sigmab_lv, -ifelse(randomB=="P", ncol(lv.X), num.lv.c+num.RR))
              sigmaBij <- rep(0,(ncol(lv.X)^2-ncol(lv.X))/2)
              for(i in 1:nrow(csBlv)){
                sigmaBij[(csBlv[i,1] - 1) * (csBlv[i,1] - 2) / 2 + csBlv[i,2]] = covsigmaB[i]
              }
              bL <-  constructL(sigmaBij)
              out$params$corsLvXcoef <- bL%*%t(bL)
            }
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
      
      # if(family %in% "betaH"){
      #   out$params$betaH <- betaH;
      #   rownames(out$params$betaH) <- colnames(out$y); 
      #   colnames(out$params$betaH) <- paste("c", 1:ncol(betaH))
      #   # colnames(out$params$betaH)[1] <- "Intercept"; 
      #   if(!is.null(X)){ colnames(out$params$betaH) <- c("Intercept",colnames(X)); }
      #   if(num.lv>0) {
      #     out$params$thetaH <- thetaH
      #   }      
      # }
      
      if(family =="negative.binomial") {
        out$params$inv.phi <- phis;
        out$params$phi <- 1/phis;
        names(out$params$phi) <- colnames(y);
        if(!is.null(names(disp.group))){
          try(names(out$params$phi) <- names(disp.group),silent=T)
        }
        names(out$params$inv.phi) <-  names(out$params$phi)
      }
      
      if(family =="ZINB") {
        out$params$ZINB.inv.phi <- ZINBphis;
        out$params$ZINB.phi <- 1/ZINBphis;
        names(out$params$ZINB.phi) <- colnames(y);
        if(!is.null(names(disp.group))){
          try(names(out$params$ZINB.phi) <- names(map.list$lg_phiZINB),silent=T)
        }
        names(out$params$ZINB.inv.phi) <-  names(out$params$ZINB.phi)
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
        
      }
      if((nrow(dr)==n) || (nrow(xr == n))) {
        if(nrow(dr)==n){ 
          out$dr=dr
          iter = 1 # keep track of index
          for(re in 1:length(cstrucn)){
            if(cstrucn[re] %in% c(1,3)) {
              sigma[iter] <- exp(sigma[iter])
              names(sigma)[iter] = names(nr)[re]
              names(sigma)[iter+1] = paste0(names(nr)[re],"rho")
              sigma[iter+1] <- sigma[iter+1] / sqrt(1.0 + sigma[iter+1]^2);
              iter <- iter +2
            } else if(cstrucn[re] %in% c(2)){
              sigma[iter:(iter+1)] <- exp(sigma[iter:(iter+1)])
              names(sigma)[iter] = "Scale"
              names(sigma)[iter+1] = names(nr)[re]
              iter <- iter + 2
            } else if(cstrucn[re] %in% c(4)){
              # sigma[iter:(iter+2)] <- exp(sigma[iter:(iter+2)]) # maternKappa fixed
              sigma[iter:(iter+1)] <- exp(sigma[iter:(iter+1)])
              names(sigma)[iter] = "Scale"
              names(sigma)[iter+1] = names(nr)[re]
              iter <- iter + 2
              # Matern smoothness
              # names(sigma)[iter+1] = "Matern kappa"
              # iter <- iter +1
            } else {
              sigma[iter] <- exp(sigma[iter])
              names(sigma)[iter] = names(nr)[re]
              iter <- iter +1
            }
          }
          out$params$sigma=sigma; 
          out$params$row.params.random <- row.params.random; 
          try(names(out$params$row.params.random) <- colnames(dr), silent = TRUE)
          # if((rstruc ==2 | (rstruc == 1)) & (cstrucn %in% c(1,2,3,4))){ 
          #   out$params$rho <- rho
          #   names(out$params$rho)="rho"
          #   # if(cstrucn %in% c(2,4)){ out$params$scaledc=scaledc}
          # }
          # if((num.lv+num.lv.c)>1 && dependent.row) names(out$params$sigma) <- paste("sigma",c("",1:(num.lv+num.lv.c)), sep = "")
        }
        if(nrow(xr)==n){
          out$params$row.params.fixed <- row.params.fixed
          try(names(out$params$row.params.fixed) <- colnames(xr), silent = TRUE)
        }
        
      }
      if(col.eff == "random"){
        row.names(Br) <- colnames(spdr)
        if(!is.null(colnames(y))) colnames(Br) <- colnames(y)
        out$params$Br <- Br
        if(nrow(Xt)==n){
        out$params$B <- B
        names(out$params$B) <- colnames(Xt)
        }
        out$params$sigmaB <- diag(sigma.sp[1:ncol(spdr)], ncol(spdr))
        if(any(colMat[row(colMat)!=col(colMat)]!=0)){
          if(colMat.rho.struct == "term"){
            names(rho.sp) <- colnames(spdr)
          }else{
            names(rho.sp) <- NULL
          }
          out$params$rho.sp <- rho.sp
        }
        if(ncol(cs)==2){
          sigmaSPij <- rep(0,(ncol(spdr)^2-ncol(spdr))/2)
            for(i in 1:nrow(cs)){
              sigmaSPij[(cs[i,1] - 1) * (cs[i,1] - 2) / 2 + cs[i,2]] = covsigma.sp[i]
            }
          SprL <- out$params$sigmaB%*%constructL(sigmaSPij)
          out$params$sigmaB <- SprL%*%t(SprL)
        }
        colnames(out$params$sigmaB) <- colnames(spdr)
        out$spdr <- spdr
      }
      
      if(num.lv.cor>0 & cstruclvn>0){
        out$params$rho.lv <- rho.lv; 
        if(cstruclvn %in% c(2,4)){ 
          if(length(out$params$rho.lv)>0) 
            names(out$params$rho.lv) <- paste("rho.lv",1:length(out$params$rho.lv), sep = "") #[!is.na(map.list$rho_lvc)]
        } else if(!is.null(rho.lv)){
          names(out$params$rho.lv) <- paste("rho.lv",1:num.lv.cor, sep = "") 
        }
      }
      
      if(family %in% c("binomial", "beta", "ordinal")) out$link <- link;
      if(family == "tweedie") out$Power <- Power;
      if(family %in% c("ordinal", "orderedBeta")){
        out$params$zeta <- zetas
      }
      
      out$time <- timeo
      pars <- optr$par
      
      ## colMatect VA covariances
      if((method %in% c("VA", "EVA"))){
        param <- objr$env$last.par.best
        
        if(num.lv.cor>0 && !corWithinLV){
          Au <- param[names(param)=="Au"]
          AQ <- NULL
          
          if(cstruclvn==0){
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
            
            if(Astruc<3) {
              Au <- param[names(param)=="Au"]
              # Au <- exp(param[names(param)=="Au"])^2
              for (d in 1:(num.lv.cor)){
                A[,,d] <- diag(exp(Au[(d-1)*nu+1:nu]),nu,nu)
                
                k=0;
                if((Astruc==1) & (length(Au) > nu*num.lv.cor) ){ # unstructured variational covariance
                  for (i in 1:(nu-1)){
                    for (r in (i+1):nu){
                      A[r,i,d]=Au[nu*num.lv.cor+k*num.lv.cor+d];
                      k=k+1;
                    }
                  }
                } else if((Astruc==2) & (length(Au) > nu*num.lv.cor)) { # bdNN variational covariance
                  arank = nrow(NN);
                  for (r in  1:arank){
                    A[NN[r,1],NN[r,2],d]=Au[nu*num.lv.cor+k*num.lv.cor+d];
                    k=k+1;
                  }
                }
                A[,,d]=A[,,d]%*%t(A[,,d])
              }
              
            } else {
              # Alvm <- array(objr$report()$Alvm, dim=c(nu, nu, nMax))
              for (d in 1:nMax) {
                if(Astruc %in% c(3,4)){
                  A[,,d] <- objr$report()$Alvm
                  # A[,,d] <- Alvm %*%t(Alvm)
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
          
        } else if(num.lv.cor>0 && corWithinLV){
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
            # Uncomm
            # if(Astruc %in% c(3,4)){
            #   A[,,q] <- Alvm%*%t(Alvm)
            # } else {
            #   A[,,q] <- Alvm[,,q]%*%t(Alvm[,,q])
            # }
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
          }
        }
        
        # For random slopes constr. ord.
        if((num.RR+num.lv.c)>0&!isFALSE(randomB)){
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
        
        if(nrow(dr)==n){
          lg_Ar <- param[names(param)=="lg_Ar"]
          Ar <- vector("list", length(nr))
          Ar.sds <- exp((lg_Ar)[1:sum(nr)])
          lg_Ar <- lg_Ar[-c(1:sum(nr))]
          for(re in 1:length(nr)){
            Ar[[re]] <- diag(Ar.sds[1:nr[re]])
          }
          if(Ar.struc == "unstructured"){
            if(length(lg_Ar)>0){
              k=1;
              for(re in 1:length(nr)){
                for(d in 1:(nr[re]-1)){
                  for(r in (d+1):nr[re]){
                    if(cstruc[re]>0)Ar[[re]][r,d] = lg_Ar[k];
                    k=k+1;
                  }}
              }
            }
          }
          for(re in 1:length(nr)){
            Ar[[re]] <- Ar[[re]]%*%t(Ar[[re]])
          }
          out$Ar <- Ar
        }
        
        if(col.eff=="random"){
          spAr <- param[names(param)=="Abb"]
          if(sp.Ar.struc%in%c("blockdiagonal","diagonal")){
            spArs <- vector("list", p) #p*ncol(spdr)
            Ar.sds <- exp((spAr)[1:(p*ncol(spdr))])
            spAr <- spAr[-c(1:(p*ncol(spdr)))]
            k=1;
            
            for(j in 1:p){
              spArs[[j]] <- diag(Ar.sds[1:ncol(spdr)], ncol(spdr))
              Ar.sds <- Ar.sds[-c(1:ncol(spdr))]
              if(sp.Ar.struc == "blockdiagonal" && ncol(spdr) > 1){
                if(length(spAr)>0){
                  for(d in 1:(ncol(spdr)-1)){
                    for(r in (d+1):ncol(spdr)){
                      spArs[[j]][r,d] = spAr[k];
                      k=k+1;
                    }}
                }
              }
              spArs[[j]] <- spArs[[j]]%*%t(spArs[[j]])
            }
          }else if(sp.Ar.struc %in% c("MNdiagonal","MNunstructured")){
            spAr <- param[names(param)=="Abb"]
            spArs <- vector("list", 2)

               Ar.sds <- exp(spAr[1:(p+ncol(spdr)-1)])
               spAr <- spAr[-c(1:(p+ncol(spdr)-1))]
              
              spArs[[1]] <- diag(Ar.sds[1:ncol(spdr)], ncol(spdr))
              spArs[[2]] <- diag(c(.3,Ar.sds[-c(1:ncol(spdr))]))
              if(sp.Ar.struc == "MNunstructured" && ncol(spdr)>1){
              # row covariance
                    for(d in 1:(ncol(spdr)-1)){
                      for(r in (d+1):ncol(spdr)){
                        spArs[[1]][r,d] = spAr[1];
                        spAr <- spAr[-1]
                      }}
                }
                # column covariance
                
                  sp = 0;

                  for(cb in 1:length(blocks[-1])){
                    for (j in 1:Abranks[cb]){
                      for (r in (j+1):blocksp[cb]){
                        if(j<r && r<=blocksp[cb]){
                        spArs[[2]][r+sp,j+sp]=spAr[1];
                        spAr <- spAr[-1]
                        }
                      }
                    }
                  sp = sp +blocksp[cb]
                }

                spArs[[1]] <- spArs[[1]]%*%t(spArs[[1]])
                spArs[[2]] <- spArs[[2]]%*%t(spArs[[2]])
          }else if(sp.Ar.struc %in% c("diagonalCL1", "CL1")){
            Ar.sds <- exp(spAr[1:(sum(blocksp*ncol(spdr))+p-length(blocksp))])
            spAr <- spAr[-c(1:(sum(blocksp*ncol(spdr))+p-length(blocksp)))]
            SArmbs <- list()
            SArmC <- list()
            sp <- 1
            for(cb in 1:length(blocks[-1])){
              # build block diagonal matrices
              for(j in 1:blocksp[cb]){
                SArmbs[[sp]] <- diag(Ar.sds[1:ncol(spdr)], ncol(spdr))
                Ar.sds <- Ar.sds[-c(1:ncol(spdr))]
                if(sp.Ar.struc == "CL1" && ncol(spdr)>1){
                  for(d in 1:(ncol(spdr)-1)){
                    for(r in (d+1):ncol(spdr)){
                      SArmbs[[sp]][r,d] = spAr[1];
                      spAr <- spAr[-1]
                    }}
                }
                sp <- sp+1
              }
              # build second matrix, first diagonal entry fixed for identifiability
              SArmC[[cb]] = diag(c(1, if(blocksp[cb]>1)Ar.sds[1:(blocksp[cb]-1)]))
              Ar.sds <- Ar.sds[-c(1:(blocksp[cb]-1))]
              
              for (j in 1:Abranks[cb]){
                for (r in (j+1):blocksp[cb]){
                  if(j<r && r<=blocksp[cb]){
                    SArmC[[cb]][r,j]=spAr[1];
                    spAr <- spAr[-1]
                  }
                }
              }
              SArmC[[cb]] <- cov2cor(SArmC[[cb]]%*%t(SArmC[[cb]]))
            }
            spArs <- Matrix::bdiag(SArmbs)%*%kronecker(Matrix::bdiag(SArmC),diag(ncol(spdr)))%*%Matrix::t(Matrix::bdiag(SArmbs))
          }else if(sp.Ar.struc %in% c("CL2")){
            Ar.sds <- exp(spAr[1:(ncol(spdr)-1+p*ncol(spdr))])
            spAr <- spAr[-c(1:(ncol(spdr)-1+p*ncol(spdr)))]
            SArmR <- diag(c(1,Ar.sds[1:(ncol(spdr)-1)]), ncol(spdr))
            Ar.sds <- Ar.sds[-c(1:(ncol(spdr)-1))]
            if(sp.Ar.struc=="CL2"){
              for(d in 1:(ncol(spdr)-1)){
                for(r in (d+1):ncol(spdr)){
                  SArmR[r,d] = spAr[1];
                  spAr <- spAr[-1]
                }}
            }
            SArmPs <- vector("list", length(blocks)-1)
            for(cb in 1:length(blocks[-1])){
              SArmPs[[cb]] <- vector("list", ncol(spdr))
              for(d in 1:ncol(spdr)){
                if(blocksp[cb]>1){
                SArmPs[[cb]][[d]] <- diag(Ar.sds[1:blocksp[cb]], blocksp[cb])
                Ar.sds <- Ar.sds[-c(1:blocksp[cb])]
                
                  for (j in 1:Abranks[cb]){
                    for (r in (j+1):blocksp[cb]){
                      if(j<r && r<=blocksp[cb]){
                        SArmPs[[cb]][[d]][r,j]=spAr[1];
                        spAr <- spAr[-1]
                      }
                    }
                  }}else{
                    SArmPs[[cb]][[d]]<-diag(1,1)
                  }
              }
            }
            # need to bind covariate-wise
            spArs <- vector("list", length=ncol(spdr))
            for(d in 1:ncol(spdr)){
              spArs[[d]] <- Matrix::bdiag(sapply(SArmPs, "[[", d,simplify=FALSE)) # get every dth element for each block
            }
            spArs <- Matrix::bdiag(spArs)%*%kronecker(cov2cor(SArmR%*%t(SArmR)),diag(p))%*%Matrix::t(Matrix::bdiag(spArs))
          }else if(sp.Ar.struc == "diagonalCL2"){
            spArs <- vector("list", 1)
            for(d in 1:ncol(spdr)){
              Ar.sds <- exp((spAr)[1:p])
              spAr <- spAr[-c(1:p)]
              spArs[[d]] <- diag(Ar.sds)
            }
            for(d in 1:ncol(spdr)){
                sp = 0;
                for(cb in 1:length(blocks[-1])){
                    for (j in 1:Abranks[cb]){
                      for (r in (j+1):blocksp[cb]){
                        if(j<r && r<=blocksp[cb]){
                        spArs[[d]][r+sp,j+sp]=spAr[1];
                        spAr <- spAr[-1]
                        }
                      }
                    }
                  sp = sp +blocksp[cb]
                }
              spArs[[d]] <- spArs[[d]]%*%t(spArs[[d]])
            }
          }else if(sp.Ar.struc == "unstructured"){
            spArs <- vector("list", 1)
            
            Ar.sds <- exp((spAr)[1:(p*ncol(spdr))])
            spAr <- spAr[-c(1:(p*ncol(spdr)))]

                for(cb in 1:length(blocks[-1])){
                  spArs[[cb]] <- diag(Ar.sds[1:(blocksp[cb]*ncol(spdr))])
                  
                for(j in 1:Abranks[cb]){
                  for(r in (j+1):(blocksp[cb]*ncol(spdr))){
                    if(j<r && r<=(blocksp[cb]*ncol(spdr))){
                    spArs[[cb]][r,j] = spAr[1];
                    spAr <- spAr[-1]
                    }
                  }} 
                }
                spArs <- list(as.matrix(Matrix::bdiag(spArs)))
            spArs[[1]] <- spArs[[1]]%*%t(spArs[[1]])
          }
          out$Ab <- spArs
        }
      }
    }else{
      objrFinal <- list()
      optrFinal <- list()
    }
  
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
  
  
  
  return(out)
}

