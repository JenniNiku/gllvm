#' @title Plot latent variables from gllvm model
#' @description Plots latent variables and their corresponding coefficients (biplot).
#'
#' @param object   an object of class 'gllvm'.
#' @param biplot   \code{TRUE} if both latent variables and their coefficients are plotted, \code{FALSE} if only latent variables.
#' @param ind.spp  the number of response variables (usually, species) to include on the biplot. The default is none, or all if \code{biplot = TRUE}.
#' @param alpha    a numeric scalar between 0 and 1 that is used to control the relative scaling of the latent variables and their coefficients, when constructing a biplot.
#' @param main  main title.
#' @param which.lvs indices of two latent variables to be plotted if number of the latent variables is more than 2. A vector with length of two. Defaults to \code{c(1,2)}. 
#' @param jitter   if \code{TRUE}, jittering is applied on points.
#' @param jitter.amount   numeric, positive value indicating an amount of jittering for each point, defaults to 0.2 (jitter range).
#' @param s.colors colors for sites
#' @param symbols logical, if \code{TRUE} sites are plotted using symbols, if \code{FALSE} (default) site numbers are used
#' @param cex.spp size of species labels in biplot
#' @param cex.env size of labels for arrows in constrianed ordination
#' @param spp.colors colors for sites, defaults to \code{"blue"}
#' @param spp.arrows plot species scores as arrows if outside of the range of the plot? Defaults to \code{FALSE} for linear response models and \code{TRUE} for quadratic response models.
#' @param lab.dist distance between label and arrow heads. Value between 0 and 1
#' @param arrow.scale positive value, to scale arrows
#' @param arrow.ci represent statistical uncertainty for arrows in constrained ordinatioon using confidence interval? Defaults to \code{TRUE}
#' @param predict.region logical, if \code{TRUE} prediction regions for the predicted latent variables are plotted, defaults to \code{FALSE}.
#' @param level level for prediction regions.
#' @param lty.ellips line type for prediction ellipses. See graphical parameter lty.
#' @param lwd.ellips line width for prediction ellipses. See graphical parameter lwd.
#' @param col.ellips colors for prediction ellipses.
#' @param ...	additional graphical arguments.
#'
#' @details
#' Function constructs a scatter plot of two latent variables, i.e. an ordination plot. 
#' Latent variables are re-rotated to their principal direction using singular value decomposition,
#' so that the first plotted latent variable does not have to be the first latent variable in the model.
#' If only one latent variable is in the fitted model, latent variables are plotted against their corresponding row indices.
#' The latent variables are labeled using the row index of the response matrix y.
#'
#' Coefficients related to latent variables are plotted in the same figure with the latent
#' variables if \code{biplot = TRUE}. They are labelled using the column names of y. The number
#' of latent variable coefficients to be plotted can be controlled by ind.spp. An argument alpha
#' is used to control the relative scaling of the latent variables and their coefficients.
#' If \code{alpha = 0.5}, the latent variables and their coefficients are on the same scale.
#' For details for constructing a biplot, see Gabriel (1971).
#' 
#' Latent variable scores are always scaled by their estimated standard deviations, for plotting.
#' 
#' For a quadratic response model, species optima are plotted. Any species scores that are outside the range 
#' of the predicted site scores are not directly plotted, but their main direction is indicated with arrows instead.
#' This ensures that the plot remains on a reasonable scale.
#' 
#' Effects of environmental variables in constrained ordination are indicated with arrows.
#' If any of the arrows exceeds the range of the plot, arrows are scaled to 80% of the plot range,
#' but so that the relative contribution of predictors is maintained.
#' If standard errors are available in the provided model, the slopes of environmental variables
#' for which the 95% confidence intervals do not include zero are shown as red, while others 
#' are slightly less intensely coloured.
#' 
#' 
#' @note 
#' - If error is occurred when using \code{ordiplot()}, try full name of the function \code{ordiplot.gllvm()} as functions named 'ordiplot' might be found in other packages as well.
#' 
#' @references 
#' Gabriel, K. R. (1971). The biplot graphic display of matrices with application to principal component analysis. Biometrika, 58, 453-467.
#' 
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Francis K.C. Hui
#'
#' @examples
#' #'# Extract subset of the microbial data to be used as an example
#'data(microbialdata)
#'y <- microbialdata$Y[, order(colMeans(microbialdata$Y > 0), 
#'                      decreasing = TRUE)[21:40]]
#'fit <- gllvm(y, family = poisson())
#'fit$logL
#'ordiplot(fit, predict.region = TRUE)
#' \dontrun{
#' #'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'fit <- gllvm(y, family = poisson())
#'# Ordination plot:
#'ordiplot(fit)
#'# Biplot with 10 species
#'ordiplot(fit, biplot = TRUE, ind.spp = 10)
#'}
#'@aliases ordiplot ordiplot.gllvm
#'@export
#'@export ordiplot.gllvm
ordiplot.gllvm <- function(object, biplot = FALSE, ind.spp = NULL, alpha = 0.5, main = NULL, which.lvs = c(1, 2), predict.region = FALSE, level =0.95,
                           jitter = FALSE, jitter.amount = 0.2, s.colors = 1, symbols = FALSE, cex.spp = 0.7, spp.colors = "blue", arrow.scale = 0.8, arrow.ci = TRUE, spp.arrows = NULL, cex.env = 0.7, lab.dist = 0.1, lwd.ellips = 0.5, col.ellips = 4, lty.ellips = 1, ...) {
  if (!any(class(object) %in% "gllvm"))
    stop("Class of the object isn't 'gllvm'.")
  if(is.null(spp.arrows)){
    if(object$quadratic!=FALSE){
      spp.arrows <- TRUE
    }else{
      spp.arrows <- FALSE
    }
  }
  arrow.scale <- abs(arrow.scale)
  a <- jitter.amount
  n <- NROW(object$y)
  p <- NCOL(object$y)
  num.lv <- object$num.lv
  num.lv.c <- object$num.lv.c 
  num.RR <- object$num.RR
  quadratic <- object$quadratic
  if (!is.null(ind.spp)) {
    ind.spp <- min(c(p, ind.spp))
  } else {
    ind.spp <- p
  }
  if(length(spp.colors)==1){
    spp.colors <- rep(spp.colors,p)
  }else if(length(spp.colors)!=p){
    stop("spp.colors needs to be of length p or 1.")
  }
  if ((num.lv+(num.lv.c+num.RR)) == 0)
    stop("No latent variables to plot.")
  
  if (is.null(rownames(object$params$theta)))
    rownames(object$params$theta) = paste("V", 1:p)
  
  if((num.lv.c+num.RR)>0){
    type <- "constrained"
  }else{
    type <- "scaled"
  }
  
  lv <- getLV(object, type = type)
  
  if ((num.lv+(num.lv.c+num.RR)) == 1) {
    if(num.lv==1)plot(1:n, lv, ylab = "LV1", xlab = "Row index")
    if((num.lv.c+num.RR)==1)plot(1:n, lv, ylab = "CLV1", xlab = "Row index")
  }
  
  if ((num.lv+(num.lv.c+num.RR)) > 1) {
    if((num.lv.c+num.RR)>0){
      do_svd <- svd(lv)
    }else{
      # This gives the right rotation for unconstrained ordination
      do_svd <- svd(object$lvs)
    }
    # do_svd <- svd(lv)
    # do_svd <- svd(object$lvs)
    svd_rotmat_sites <- do_svd$v
    svd_rotmat_species <- do_svd$v

    choose.lvs <- lv
    if(quadratic == FALSE){choose.lv.coefs <- object$params$theta}else{choose.lv.coefs<-optima(object,sd.errors=F)}  
    
    #A check if species scores are within the range of the LV
    ##If spp.arrows=TRUE plots those that are not in range as arrows
    if(spp.arrows){
      lvth <- max(abs(choose.lvs))
      idx <- choose.lv.coefs>(-lvth)&choose.lv.coefs<lvth
    }else{
      idx <- matrix(TRUE,ncol=num.lv+num.lv.c+num.RR,nrow=p)
    }
    
    bothnorms <- vector("numeric",ncol(choose.lv.coefs))
    for(i in 1:ncol(choose.lv.coefs)){
      bothnorms[i] <- sqrt(sum(choose.lvs[,i]^2)) * sqrt(sum(choose.lv.coefs[idx[,i],i]^2))
    }
    
    # bothnorms <- sqrt(colSums(choose.lvs^2)) * sqrt(colSums(choose.lv.coefs^2)) 
    ## Standardize both to unit norm then scale using bothnorms. Note alpha = 0.5 so both have same norm. Otherwise "significance" becomes scale dependent
    scaled_cw_sites <- t(t(choose.lvs) / sqrt(colSums(choose.lvs^2)) * (bothnorms^alpha)) 
    # scaled_cw_species <- t(t(choose.lv.coefs) / sqrt(colSums(choose.lv.coefs^2)) * (bothnorms^(1-alpha))) 
    scaled_cw_species <- choose.lv.coefs
    for(i in 1:ncol(scaled_cw_species)){
      scaled_cw_species[,i] <- choose.lv.coefs[,i] / sqrt(sum(choose.lv.coefs[idx[,i],i]^2)) * (bothnorms[i]^(1-alpha)) 
    }
    
    choose.lvs <- scaled_cw_sites%*%svd_rotmat_sites
    choose.lv.coefs <- scaled_cw_species%*%svd_rotmat_species
    # 
    # if(spp.arrows){
    #   idx <- choose.lv.coefs>matrix(apply(choose.lvs,2,min),ncol=ncol(choose.lv.coefs),nrow=nrow(choose.lv.coefs),byrow=T)&choose.lv.coefs<matrix(apply(choose.lvs,2,max),ncol=ncol(choose.lv.coefs),nrow=nrow(choose.lv.coefs),byrow=T)
    # }else{
    #   idx <- matrix(TRUE,ncol=num.lv+num.lv.c+num.RR,nrow=p)
    # }
    # 
    B<-(diag((bothnorms^alpha)/sqrt(colSums(object$lvs^2)))%*%svd_rotmat_sites)
    
    
    # testcov <- object$lvs %*% t(object$params$theta)
    # do.svd <- svd(testcov, num.lv, num.lv)
    # choose.lvs <- do.svd$u * matrix( do.svd$d[1:num.lv] ^ alpha,
    #     nrow = n, ncol = num.lv, byrow = TRUE )
    # choose.lv.coefs <- do.svd$v * matrix(do.svd$d[1:num.lv] ^ (1 - alpha),
    #     nrow = p, ncol = num.lv, byrow = TRUE )
    
    
    if (!biplot) {
      plot(choose.lvs[, which.lvs],
           xlab = paste("Latent variable", which.lvs[1]), 
           ylab = paste("Latent variable", which.lvs[2]),
           main = main , type = "n", ... )
      
      if (predict.region) {
        if(length(col.ellips)!=n){ col.ellips =rep(col.ellips,n)}
        if (object$method == "LA") {
          for (i in 1:n) {
            covm <- (t(B)%*%object$prediction.errors$lvs[i,,]%*%B)[which.lvs,which.lvs];
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=num.lv+num.lv.c+num.RR)), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
          }        
          } else {
          
          sdb<-CMSEPf(object)$A
          if(object$row.eff=="random" && dim(object$A)[2]>dim(object$lvs)[2]){
            object$A<- object$A[,-1,-1]
          }
          if(num.RR>0){
            #variational covariances but add 0s for RRR
            A <- array(0,dim=c(n,num.lv.c+num.RR+num.lv,num.lv.c+num.RR+num.lv))
            A[,-c((num.lv.c+1):(num.lv.c+num.RR)),-c((num.lv.c+1):(num.lv.c+num.RR))] <- object$A
          }else{A<-object$A}
          object$A<-sdb+A
          
          r=0
          for (i in 1:n) {
            if(!object$TMB && object$Lambda.struc == "diagonal"){
              covm <- (t(B)%*%diag(object$A[i,1:num.lv+r])%*%B)[which.lvs,which.lvs];
              # covm <- diag(object$A[i,which.lvs+r]);
            } else {
              covm <- (t(B)%*%object$A[i,1:(num.lv+num.lv.c+num.RR)+r,1:(num.lv+num.RR+num.lv.c)+r]%*%B)[which.lvs,which.lvs];
              # covm <- object$A[i,which.lvs+r,which.lvs+r];
            }
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=num.lv+num.RR+num.lv.c)), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
          }
        }
      }
      
      
      if (!jitter)
        if (symbols) {
          points(choose.lvs[, which.lvs], col = s.colors, ...)
        } else {
          text(choose.lvs[, which.lvs], label = 1:n, cex = 1.2, col = s.colors)
        }
      if (jitter)
        if (symbols) {
          points(choose.lvs[, which.lvs][, 1] + runif(n,-a,a), choose.lvs[, which.lvs][, 2] + runif(n,-a,a), col =
                   s.colors, ...)
        } else {
          text(
            (choose.lvs[, which.lvs][, 1] + runif(n,-a,a)),
            (choose.lvs[, which.lvs][, 2] + runif(n,-a,a)),
            label = 1:n, cex = 1.2, col = s.colors )
        }
    }
    
    if (biplot) {
      if(quadratic==F)largest.lnorms <- order(apply(object$params$theta ^ 2, 1, sum), decreasing = TRUE)[1:ind.spp]
      if(quadratic!=F)largest.lnorms <- order(apply(object$params$theta[,1:((num.lv.c+num.RR)+num.lv)] ^ 2, 1, sum)+2*apply(object$params$theta[,-c(1:((num.lv.c+num.RR)+num.lv))] ^ 2, 1, sum) , decreasing = TRUE)[1:ind.spp]
      
      plot(
        rbind(choose.lvs[, which.lvs], choose.lv.coefs[apply(idx,1,all), which.lvs]),
             xlab = paste("Latent variable", which.lvs[1]), 
             ylab = paste("Latent variable", which.lvs[2]),
            main = main, type = "n", ... )
      
      if (predict.region) {
        if(length(col.ellips)!=n){ col.ellips =rep(col.ellips,n)}
        if (object$method == "LA") {
          for (i in 1:n) {
            covm <- (t(B)%*%object$prediction.errors$lvs[i,,]%*%B)[which.lvs,which.lvs];
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=num.lv+num.lv.c+num.RR)), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
          }
        } else {
          
          sdb<-CMSEPf(object)$A
          if(object$row.eff=="random" && dim(object$A)[2]>dim(object$lvs)[2]){
            object$A<- object$A[,-1,-1]
          }
          if(num.RR>0){
            A <- array(0,dim=c(n,num.lv.c+num.RR+num.lv,num.lv.c+num.RR+num.lv))
            A[,-c((num.lv.c+1):(num.lv.c+num.RR)),-c((num.lv.c+1):(num.lv.c+num.RR))] <- object$A
          }else{A<-object$A}
          object$A<-sdb+A
          r=0
          for (i in 1:n) {
            if(!object$TMB && object$Lambda.struc == "diagonal"){
              covm <- (t(B)%*%diag(object$A[i,1:num.lv+r])%*%B)[which.lvs,which.lvs];
              # covm <- diag(object$A[i,which.lvs+r]);
            } else {
              covm <- (t(B)%*%object$A[i,1:(num.lv+num.lv.c+num.RR)+r,1:(num.lv+num.RR+num.lv.c)+r]%*%B)[which.lvs,which.lvs];
              # covm <- object$A[i,which.lvs+r,which.lvs+r];
            }
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=num.lv+num.RR+num.lv.c)), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
          }
        }
      }

      if (!jitter){
        if (symbols) {
          points(choose.lvs[, which.lvs], col = s.colors, ...)
        } else {
          text(choose.lvs[, which.lvs], label = 1:n, cex = 1.2, col = s.colors)
        }
    
        spp.colors <- spp.colors[largest.lnorms][1:ind.spp]
        text(
          matrix(choose.lv.coefs[largest.lnorms[1:ind.spp],which.lvs][apply(idx[largest.lnorms[1:ind.spp],which.lvs],1,function(x)all(x)),], nrow = sum(apply(!idx[largest.lnorms[1:ind.spp],which.lvs],1,function(x)!any(x)))),
          label = rownames(object$params$theta)[largest.lnorms[1:ind.spp]][apply(idx[largest.lnorms[1:ind.spp],which.lvs],1,function(x)all(x))],
          col = spp.colors, cex = cex.spp )
      }
      if (jitter){
        if (symbols) {
          points(choose.lvs[, which.lvs[1]] + runif(n,-a,a), (choose.lvs[, which.lvs[2]] + runif(n,-a,a)), col =
                   s.colors, ...)
        } else {
          text(
            (choose.lvs[, which.lvs[1]] + runif(n,-a,a)),
            (choose.lvs[, which.lvs[2]] + runif(n,-a,a)),
            label = 1:n, cex = 1.2, col = s.colors )
        }
        spp.colors <- spp.colors[apply(!idx[largest.lnorms[1:ind.spp],which.lvs],1,function(x)!any(x))]
        text(
          (matrix(choose.lv.coefs[largest.lnorms,which.lvs][apply(!idx[largest.lnorms[1:ind.spp],which.lvs],1,function(x)!any(x)),], nrow =sum(apply(!idx[largest.lnorms[1:ind.spp],which.lvs],1,function(x)!any(x)))) + runif(2*sum(apply(!idx[largest.lnorms[1:ind.spp],which.lvs],1,function(x)!any(x))),-a,a)),
          label = rownames(object$params$theta)[largest.lnorms][apply(!idx[largest.lnorms[1:ind.spp],which.lvs],1,function(x)!any(x))],
          col = spp.colors, cex = cex.spp )
      }
      
      ##Here add arrows for species with optima outside of range LV
      # if(quadratic!=FALSE){
      if(spp.arrows){
        marg<-par("usr")
        Xlength<-sum(abs(marg[1:2]))/2
        Ylength<-sum(abs(marg[3:4]))/2
        origin<- c(mean(marg[1:2]),mean(marg[3:4]))
        #scores_to_plot <- choose.lv.coefs[largest.lnorms[!apply(idx[largest.lnorms,which.lvs,drop=F],1,all)],which.lvs,drop=F]
        scores_to_plot <- choose.lv.coefs[largest.lnorms[1:ind.spp],which.lvs][apply(idx[largest.lnorms[1:ind.spp],which.lvs],1,function(x)!all(x)),,drop=F]
        if(nrow(scores_to_plot)>0){
        ends <- t(t(t(t(scores_to_plot)-origin)/sqrt((scores_to_plot[,1]-origin[1])^2+(scores_to_plot[,2]-origin[2])^2)*min(Xlength,Ylength)))*arrow.scale
        for(i in 1:nrow(scores_to_plot)){
          arrows(origin[1],origin[2],ends[i,1]+origin[1],ends[i,2]+origin[2],col=spp.colors[i], cex = cex.spp, length = 0.2)
          text(x=ends[i,1]*(1+lab.dist)+origin[1],y=ends[i,2]*(1+lab.dist)+origin[2],labels = row.names(scores_to_plot)[i],col=spp.colors[i], cex = cex.spp)
        }
        }
      # }
      }
    }
    
    #Only draw arrows when no unconstrained LVs are present currently: diffcult otherwise due to rotation
    #Could alternatively post-hoc regress unconstrained LVs..but then harder to distinguish which is post-hoc in the plot..
    if(num.lv==0&(num.lv.c+num.RR)>0){
      LVcoef <- (object$params$LvXcoef%*%svd_rotmat_sites)[,which.lvs]

      if(!is.logical(object$sd)&arrow.ci){
        covB <- object$Hess$cov.mat.mod
        colnames(covB) <- row.names(covB) <- names(object$TMBfn$par)[object$Hess$incl]
        covB <- covB[row.names(covB)=="b_lv",colnames(covB)=="b_lv"]
        rotSD <- matrix(0,ncol=num.RR+num.lv.c,nrow=ncol(object$lv.X)) 
        #using svd_rotmat_sites instead of B so that uncertainty of the predictors is not affected by the scaling using alpha and sigma.lv
        for(i in 1:ncol(object$lv.X)){
          rotSD[i,] <- sqrt(abs(diag(t(svd_rotmat_sites[1:(num.lv.c+num.RR),1:(num.lv.c+num.RR)])%*%covB[seq(i,(num.RR+num.lv.c)*ncol(object$lv.X),by=ncol(object$lv.X)),seq(i,(num.RR+num.lv.c)*ncol(object$lv.X),by=ncol(object$lv.X))]%*%svd_rotmat_sites[1:(num.lv.c+num.RR),1:(num.lv.c+num.RR)])))
        }
        rotSD <- rotSD[,which.lvs]
        cilow <- LVcoef+qnorm( (1 - 0.95) / 2)*rotSD[,which.lvs]
        ciup <-LVcoef+qnorm(1- (1 - 0.95) / 2)*rotSD[,which.lvs]
        lty <- rep("solid",ncol(object$lv.X))
        col <- rep("red", ncol(object$lv.X))
        lty[sign(cilow[,1])!=sign(ciup[,1])|sign(cilow[,2])!=sign(ciup[,2])] <- "solid"
        col[sign(cilow[,1])!=sign(ciup[,1])|sign(cilow[,2])!=sign(ciup[,2])] <- hcl(0, 100, 80)#rgb(1,0,0,alpha=0.3)
        
      }else{
        lty <- rep("solid",ncol(object$lv.X))
        col<-rep("red",ncol(object$lv.X))
      }
      
      #account for variance of the predictors
      LVcoef <- LVcoef/apply(object$lv.X,2,sd)
      marg<-par("usr")
    
      origin<- c(mean(marg[1:2]),mean(marg[3:4]))
      Xlength<-sum(abs(marg[1:2]))/2
      Ylength<-sum(abs(marg[3:4]))/2

      ends <- LVcoef/max(abs(LVcoef))*min(Xlength,Ylength)*arrow.scale
      for(i in 1:nrow(LVcoef)){
        # arrows(x0=origin[1],y0=origin[2],x1=((LVcoef[i,1])/max(abs(LVcoef[,1])*Xlength*0.8)-origin[1]),y1=((LVcoef[i,2])/max(abs(LVcoef[,2])*Ylength*0.8)-origin[2]),col=col[i],lty=lty[i])
        tryCatch({arrows(x0=origin[1],y0=origin[2],x1=ends[i,1]+origin[1],y1=ends[i,2]+origin[2],col=col[i],length=0.2,lty=lty[i]);
        text(x=origin[1]+ends[i,1]*(1+lab.dist),y=origin[2]+ends[i,2]*(1+lab.dist),labels = row.names(LVcoef)[i],col=col[i], cex = cex.env)},
        warning=function(w){print(paste("The effect for", row.names(LVcoef)[i],"was too small too plot an arrow."))}
        )
      }
    }
   
  }
}


#'@export ordiplot
ordiplot <- function(object, ...)
{
  UseMethod(generic = "ordiplot")
}
