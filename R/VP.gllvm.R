#' @title Calculate variance partitioning
#' @description Calculates variance partitioning for gllvm object with function \code{varPartitioning()}.
#' 
#' @param object an object of class 'gllvm'.
#' @param group a vector of integers identifying grouping of X covariates, the default is to use model terms formula and lv.formula.
#' @param groupnames a vector of strings given as names for the groups defined in group
#' @param adj.cov logical, whether or not to adjust co-variation within the group
#' @param grouplvs logical, whether or not to group latent variables to one group
#' 
#' @details
#' 
#' Variance for the linear predictor for response j can be calculated as 
#' 
#'  \deqn{Var(\eta_j) = \sum_k \beta_{jk}^2*var(z_{.k}) + 2 \sum_{(k1=1,...,K-1)} \sum_{(k2=k1+1,...,K)}  \beta_{j(k1)}\beta_{j(k2)} Cov(Z_{.k1},Z_{.k2}) , }
#'  
#' where \eqn{z_{.k}} is a vector consisting of predictor/latent variable/row effect etc values for all sampling units i.
#' If \eqn{z_{.k}}s are not correlated, covariance term is 0 and thus the variance explained of a response j for predictor \eqn{z_{.k}} is given as \eqn{\beta_{jk}^2*var(z_{.k})/Var(\eta_j)}.
#' 
#' In case of correlated predictors, it is advised to group them into a same group. The variance explained is calculated for the correlated group of predictors together and adjusted with the covariance term.
#' 
#' 
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#' @examples
#'# Extract subset of the microbial data to be used as an example
#'data(microbialdata)
#'X <- microbialdata$Xenv
#'y <- microbialdata$Y[, order(colMeans(microbialdata$Y > 0), 
#'                      decreasing = TRUE)[21:40]]
#'fit <- gllvm(y, X[,1:3], formula = ~ pH + Phosp, family = poisson(), 
#'              studyDesign = X[,4:5], row.eff = ~(1|Site))
#'VP <- varPartitioning(fit)
#'plotVarPartitioning(VP)
#'
#'\dontrun{
#'# Plot the result of  variance partitioning
#'plotVP(VP, col = palette(hcl.colors(5, "Roma")))
#'}
#'
#'@aliases varPartitioning VP varPartitioning.gllvm
#'@export
#'@export varPartitioning.gllvm
varPartitioning.gllvm <- function(object, group = NULL, groupnames=NULL, adj.cov = TRUE, grouplvs=FALSE, ...) {
  if (!any(class(object) == "gllvm"))
    stop("Class of the object isn't 'gllvm'.")
  if(!is.null(object$lv.X) && is.null(object$lv.X.design))object$lv.X.design <- object$lv.X #for backward compatibility
  Z <- NULL
  CoefMat <- NULL
  
  groupnamesF <- groupF <- NULL
  num.lv = object$num.lv
  
  r0 <- NULL
  p <- ncol(object$y)
  if(object$family == "betaH") p <- p*2
  n <- nrow(object$y)

  if (!is.null(object$X)) {
    formula <- formula(terms(object))
  } else {
    formula <- NULL
  }


  if (is.null(colnames(object$y))) {
    colnames(object$y) <- paste("y", 1:p, sep = "")
  }

  groupF <- groupnamesF <- NULL
  # Basec env model:
  if (!is.null(object$X) && is.null(object$TR)) {
    groupnamesF <- labels(terms(formula))
    groupF <- attr(model.matrix(formula, data = object$X), "assign")
    groupF <- groupF[groupF!=0]
    B <- object$params$Xcoef
    Z  <- cbind(Z, object$X.design)
    CoefMat <- rbind(CoefMat, t(B))
    # eta <- eta + X.d %*% t(B)
  }
  
  # Fourth corner model
  if (!is.null(object$X) && !is.null(object$TR)) {
    # stop(paste("VP for Fourth-corner model not implemented yet"))
    # groupnamesF <- labels(terms(formula))
    # groupnamesF <- groupnamesF[groupnamesF %in% colnames(object$X)]
    # groupF <- attr(model.matrix(formula, data = X.d), "assign")
    # groupF <- groupF[groupF!=0]
    # strsplit(colnames(X.d), split = ":")
    # warning("Note: Automatic grouping based on formula in fourth corner model not yet implemented, group manually using 'group' and 'groupnames' arguments.")
    
    X.d <- object$X #X.d[1:n,rownames(object$fourth.corner)]
    x_in_model = colnames(X.d)
    
    # colnames(object$fourth.corner)
    # Species specific effects total:
    BTR <- matrix(0, nrow = ncol(X.d), ncol = p)
    rownames(BTR) = x_in_model
    # Main effects for X
    BTR[,] = object$params$B[x_in_model]
    # Fourth corner elements X %*% B %*% TR
    BTR[rownames(object$fourth.corner),] <- BTR[rownames(object$fourth.corner),] + object$fourth.corner %*% t(object$TR)
    # Species specific random effects
    if(is.null(object$randomX)){
      BTR[rownames(object$params$Br),] <- BTR[rownames(object$params$Br),] + object$params$Br
    }
      
    Z  <- cbind(Z, X.d)
    CoefMat <- rbind(CoefMat, BTR)
  }
  
  # Phylogen/random 
  if(is.null(object$col.eff$col.eff))object$col.eff$col.eff <- FALSE # backward compatibility
  if (object$col.eff$col.eff == "random" ) {
    groupnamesF <- labels(terms(subbars1(reformulate(sprintf("(%s)", sapply(findbars1(model1$col.eff$col.eff.formula), deparse1))))))
    groupF <- attr(model.matrix(subbars1(reformulate(sprintf("(%s)", sapply(findbars1(model1$col.eff$col.eff.formula), deparse1)))), data = object$col.eff$Xt), "assign")
    groupF <- groupF[groupF!=0]

    X.d <- object$col.eff$spdr[, colnames(object$col.eff$spdr)!= "Intercept"]
    x_in_model = colnames(X.d)
    
    # Species specific effects total:
    BBr <- matrix(0, nrow = ncol(X.d), ncol = p)
    rownames(BBr) = x_in_model
    # Main effects for X
    if(length(object$params$B)>0) BBr[x_in_model,] = object$params$B[x_in_model]
    # Species specific random effects
    BBr[rownames(object$params$Br[x_in_model,]),] <- BBr[rownames(object$params$Br[x_in_model,]),] + object$params$Br[x_in_model,]
    
    Z  <- cbind(Z, as.matrix(X.d))
    CoefMat <- rbind(CoefMat, as.matrix(BBr))
  }
  
  LVgroups = NULL
  # Inclusion of lvs
    if (object$num.lv > 0 | (object$num.lv.c + object$num.RR) > 0) {
      
      if(!is.null(object$lvs) && inherits(object$lvCor,"formula")){
        if(nrow(object$lvs)!=n) object$lvs = as.matrix(object$TMBfn$env$data$dLV%*%object$lvs) # !!!
      }
      
      theta <- (object$params$theta[, 1:(object$num.lv + (object$num.lv.c + object$num.RR)), drop = F])

      # X:s for Constrained lvs/RR
      lv.X <- object$lv.X.design
      # lv.X variances Separated only if quaratic = FALSE
      if ((object$num.lv.c + object$num.RR) > 0 & object$quadratic == FALSE) {
        if(is.null(object$params$corsLvXcoef)){
        groupnamesF <- c(groupnamesF, paste("CLV:",labels(terms(object$lv.formula)), sep = ""))
        groupF <- c(groupF, attr(model.matrix(object$lv.formula, data = object$lv.X), "assign")[-1] + max(groupF,0))
        }else{
          # cannot use formula for this..
          groupnamesF <- c(groupnamesF, paste("CLV:",colnames(object$lv.X.design), sep = ""))
          groupF <- c(groupF, attr(object$lv.X.design, "assign")[-1] + max(groupF,0))
        }
        Z <- cbind(Z, lv.X)
        Bt <- object$params$LvXcoef %*% t((theta[, 1:(object$num.lv.c + object$num.RR), drop = F]))
        rownames(Bt) <- paste("CLV:",rownames(Bt), sep = "")
        CoefMat <- rbind(CoefMat, Bt)
      } else if((object$num.lv.c + object$num.RR) > 0) {
        LVCLV = object$lv.X.design %*% object$params$LvXcoef
        Z <- cbind(Z, LVCLV)
        Bt <- t((theta[, 1:(object$num.lv.c + object$num.RR), drop = F]))
        rownames(Bt) <- paste(rownames(Bt),":X", sep = "")
        CoefMat <- rbind(CoefMat, Bt)
        LVgroups = c(LVgroups, 1:(object$num.lv.c + object$num.RR))
      }
      
      # lvs/unconstr. part
      if (object$num.RR == 0) {
        lvs <- t(t(object$lvs) * object$params$sigma.lv)
        Z <- cbind(Z, lvs)
        CoefMat <- rbind(CoefMat, t(theta))
        LVgroups = c(1:ncol(lvs))
      } else {
        # To be checked:
        if (object$num.lv.c > 0) {
            lvs <- cbind(t(t(object$lvs[, 1:object$num.lv.c, drop =FALSE]) * 
                             object$params$sigma.lv[1:object$num.lv.c]), 
                         matrix(0, ncol = object$num.RR, nrow = n), 
                         t(t(object$lvs[, -c(1:object$num.lv.c), drop =FALSE]) * 
                             object$params$sigma.lv[-(1:object$num.lv.c)]))
        } else if (object$num.lv > 0 & object$num.lv.c == 0) {
            lvs <- cbind(matrix(0, ncol = object$num.RR, 
                        nrow = n), t(t(object$lvs) * object$params$sigma.lv))
        } else {
            lvs <- matrix(0, ncol = object$num.RR, nrow = n)
        }
        
        lv0 = !(colSums(lvs==0)==n)
        Z <- cbind(Z, lvs[, lv0, drop=FALSE])
        CoefMat <- rbind(CoefMat, t(theta[, lv0, drop=FALSE]))
        LVgroups = c(LVgroups, (1:ncol(lvs))[lv0])
      }
      # eta <- eta + lvs %*% t(theta)
      
      # Quadratic
      if (object$quadratic != FALSE) {
        # stop(paste("VP for quadratic model not implemented yet"))
        
        if(object$num.lv>0){
          theta2 <- (object$params$theta[, -c(1:(object$num.lv.c + object$num.RR+object$num.lv)), drop = F])
          theta2 <- (theta2[, (object$num.lv.c + object$num.RR+1):ncol(theta2), drop = F])
          Z <- cbind(Z, (lvs[,(ncol(lvs)-object$num.lv+1):ncol(lvs), drop = F])^2)
          CoefMat <- rbind(CoefMat, t(theta2))
          
          LVgroups = c(LVgroups, (1:ncol(lvs))[(ncol(lvs)-object$num.lv+1):ncol(lvs)])
        }
        if ((object$num.lv.c + object$num.RR) > 0) {
          theta2 <- object$params$theta[,-c(1:(object$num.lv+object$num.lv.c+object$num.RR))]
          theta2C <- abs(theta2[, 1:(object$num.lv.c +
                                       object$num.RR), drop = F])
          lvs <- lvs[,1:(object$num.lv.c+object$num.RR)] + lv.X%*%object$params$LvXcoef
          Z <- cbind(Z, (lvs)^2)
          CoefMat <- rbind(CoefMat, t(-theta2C))
          
          LVgroups = c(LVgroups, (1:ncol(lvs)))
          
          # for (j in 1:p) {
          #   eta[, j] <- eta[, j] - lvs^2%*%theta2C[j, ]
          # }
        }
      }
      if(!all(sort(unique(LVgroups)) == 1:length(unique(LVgroups)))){
        LVgroups = as.numeric(factor(LVgroups))
      }
      groupF <- c(groupF, LVgroups+max(groupF,0))
      
    }

  if (!isFALSE(object$row.eff) && is.null(r0)) {
    if(!is.null(object$params$row.params.random)){
      rnams <- unique(names(object$params$row.params.random))
      for (rn in rnams) {
        r0 <- cbind(r0, as.matrix(object$TMBfn$env$data$dr0[,names(object$params$row.params.random)==rn]%*%object$params$row.params.random[names(object$params$row.params.random)==rn]) )
        CoefMat <- rbind(CoefMat, rep(1,p))
        rownames(CoefMat)[nrow(CoefMat)] = paste("Random effect:",rn)
      }
    } 
    if (!is.null(object$params$row.params.fixed)){
      if(nrow(object$TMBfn$env$data$xr)!=nrow(object$y)){
        r0 <- cbind(r0, as.matrix(object$params$row.params.fixed))
      }else{
        r0 <- cbind(r0, object$TMBfn$env$data$xr%*%object$params$row.params.fixed)
      }
      CoefMat <- rbind(CoefMat, rep(1,p))
      rownames(CoefMat)[nrow(CoefMat)] = "Fixed row effect"        
    }
    Z <- cbind(Z, r0)
  }
  
  
  
#---------------
  
  
  
  # Cov calculations
  
  
  
  # Covariances of variables
  Zcov <- (cov(Z))
  
  # For higly correlated covariates variance partitioning for
  # each of those separately may not be wise, so they could be grouped together.
  if(is.null(group)) {
    if(!is.null(groupF)){
      if(length(groupF) < ncol(Zcov)){
        group <- c(groupF, max(groupF,0) + 1:(ncol(Zcov)-length(groupF)))
      } else {
        group <- groupF[1:ncol(Zcov)]
      }
    } else {
      group <- 1:ncol(Zcov)
    }
  } else if(length(group) < ncol(Zcov)){
    if(num.lv>0) {
      if(grouplvs){
        group = c(group,max(group) + rep(1,length(LVgroups)))
        # group = c(group,max(group) + rep(1,(num.lv+object$num.lv.c)))
      } else {
        group = c(group,max(group) + LVgroups)
        # group = c(group,max(group) + 1:(1+(num.lv+object$num.lv.c)-1))
      }
    }
    if(!isFALSE(object$row.eff)){
      group = c(group,max(group) + 1:ncol(r0))
    }
  }
  
  if(is.null(groupnames)) {
    if(!is.null(groupnamesF)){
      groupnames <- groupnamesF
    }
    
    if(max(group)>length(groupnames)){
      for (k in (length(groupnames)+1):max(group)) {
        groupnames[k] <- paste(rownames(CoefMat)[group==k],collapse="/")
      }
    }
      
  } else if(max(group)>length(groupnames)){
    for (k in (length(groupnames)+1):max(group)) {
      groupnames[k] <- paste(rownames(CoefMat)[group==k],collapse="/")
    }
  }
  
  # Not grouped but adjusted cov
  totCV1 = sum((CoefMat[,1])%*%(Zcov)%*%((CoefMat[,1]))) # Total variance for species 1 with covariances between x's
  totV2 = sum((CoefMat[,1])^2*diag(Zcov)) # Total variance for species 1 without covariance terms between x's

  # With covariances and grouped:
  # Var(Lj) = sum_k{\beta*_jk^2*var(Z_k)} + 2 sum_(k1=1,...,K-1) { sum_(k2=k1+1,...,K) { \beta*_j(k1)\beta*_j(k2) Cov(Z_k1,Z_k2) }} (Eq1)
  LVpartitGr = matrix(NA, nrow = max(group),ncol = ncol(CoefMat))
  rownames(LVpartitGr) = groupnames
  colnames(LVpartitGr) = colnames(CoefMat)
  for (i in 1:max(group)) {
    Zcovi = Zcov[group==i,group==i, drop=FALSE]
    CoefMati = CoefMat[group==i,, drop=FALSE]
    if(nrow(Zcovi)>1){ # If grouped
      for (j in 1:ncol(CoefMati)) {
        if(adj.cov){
          # Adjust grouped variables with covariances
          LVpartitGr[i,j] = sum(t(t(CoefMati[,j]))%*%t(CoefMati[,j])*Zcovi)
        } else {
          LVpartitGr[i,j] = sum(diag(t(t(CoefMati[,j]))%*%t(CoefMati[,j])*Zcovi))
        }
      }
    } else { # For unique
      LVpartitGr[i,] <- c(t(CoefMati^2)%*%Zcovi)
    }
    # rownames(LVpartitGr)[i] = (rownames(CoefMat)[group==i])[1]
  }
  
  # Variances
  # Species specific total variances
  LtotV <- colSums(LVpartitGr)
  # Species specific and covariate specific variances
  LVpartit <- LVpartitGr
  
  # Proportions for explained variance (in linear predictors)
  # Species specific variance partitioning:
  PropExplainedVarSp <- t(LVpartit)/LtotV
  
  if(object$family == "betaH"){
    return(list(PropExplainedVarSp=PropExplainedVarSp[1:(p/2),, drop=FALSE],PropExplainedVarHurdleSp=PropExplainedVarSp[-(1:(p/2)),, drop=FALSE], LtotV=LtotV, LVpartit=t(LVpartit), group = group, groupnames, family = object$family))
  } else {
    return(list(PropExplainedVarSp=PropExplainedVarSp, LtotV=LtotV, LVpartit=t(LVpartit), group = group, groupnames, family = object$family))
  }
  
}


#'@export varPartitioning
varPartitioning <- function(object, ...)
{
  UseMethod(generic="varPartitioning")
}

#'@export VP
VP <- function(object, ...)
{
  varPartitioning.gllvm(object, ...)
}
