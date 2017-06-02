## Starting values by fitting marginal glm's

start.values.gllvm <- function(y, X = NULL, T=NULL, family, offset= NULL, trial.size = 1, num.lv = 0, start.lvs = NULL, seed = NULL,power=NULL) {
  if(!is.null(seed)) set.seed(seed)
       if(!is.null(T)) X=NULL;
       N <- nrow(y); p <- ncol(y); y <- as.matrix(y)
       if(!is.numeric(y))
              stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")


       if(!(family %in% c("poisson","negative.binomial","binomial","ordinal","tweedie")))
              stop("inputed family not allowed...sorry =(")

       if(num.lv > 0) {
              unique.ind <- which(!duplicated(y))
              if(is.null(start.lvs)) {
                index <- mvtnorm::rmvnorm(N, rep(0, num.lv));
                     unique.index <- as.matrix(index[unique.ind,])
                     }

              if(!is.null(start.lvs)) {
                     index <- as.matrix(start.lvs)
                     unique.index <- as.matrix(index[unique.ind,])
                     }
              }

       if(num.lv == 0) { index <- NULL }

       y <- as.matrix(y)

       if(family == "ordinal") {
              max.levels <- apply(y,2,max);
              if(any(max.levels == 1) || all(max.levels == 2)) stop("Ordinal data requires all columns to have at least has two levels. If al columns only have two levels, please use family == binomial instead. Thanks")
              }

       if(is.null(rownames(y))) rownames(y) <- paste("row",1:N,sep="")
       if(is.null(colnames(y))) colnames(y) <- paste("col",1:p,sep="")

       options(warn = -1)

       if(family != "ordinal" && family != "tweedie") { ## Using logistic instead of prbit regession here for binomial, but whatever...
              if(!is.null(X) & num.lv > 0) fit.mva <- mvabund::manyglm(y ~ X + index, family = family, K = trial.size)
              if(is.null(X) & num.lv > 0) fit.mva <- mvabund::manyglm(y ~ index, family = family, K = trial.size)
              if(!is.null(X) & num.lv == 0) fit.mva <- mvabund::manyglm(y ~ X, family = family, K = trial.size)
              if(is.null(X) & num.lv == 0) fit.mva <- mvabund::manyglm(y ~ 1, family = family, K = trial.size)
              params <- t(fit.mva$coef)
              }


       if(family == "negative.binomial") { phi <- fit.mva$phi  + 1e-5} else { phi <- NULL }

       if(family == "tweedie") {
         indeX <- cbind(index, X)
         if(num.lv == 0) indeX <- X#cbind(index = rep(1,N), X)

         phi0<-NULL
         coefs0<-NULL
         for(j in 1:p){
           fitj<-try(glm( y[,j] ~ indeX, family=statmod::tweedie(var.power=power, link.power=0) ),silent = TRUE)
           if(!inherits(fitj, "try-error")){
             coefs0<-rbind(coefs0,fitj$coefficients)
             phi0<-c(phi0,median(fitj$weights))
           } else {coefs0<-rbind(coefs0,rnorm(num.lv+dim(X)[2]+1)); phi0<-c(phi0,runif(1,0.2,3));}
         }

         params=as.matrix(coefs0)
         phi=phi0
       }

       if(family == "ordinal") {
              max.levels <- max(y)
              params <- matrix(NA,p,ncol(cbind(1,X,index)))
              zeta <- matrix(NA,p,max.levels - 1)
              zeta[,1] <- 0 ## polr parameterizes as no intercepts and all cutoffs vary freely. Change this to free intercept and first cutoff to zero
              for(j in 1:p) {
                     y.fac <- factor(y[,j])
                     if(length(levels(y.fac)) > 2) {
                            if(!is.null(X) & num.lv > 0) cw.fit <- MASS::polr(y.fac ~ X+index, method = "probit")
                            if(is.null(X) & num.lv > 0) cw.fit <- MASS::polr(y.fac ~ index, method = "probit")
                            if(!is.null(X) & num.lv == 0) cw.fit <- MASS::polr(y.fac ~ X, method = "probit")
                            if(is.null(X) & num.lv == 0) cw.fit <- MASS::polr(y.fac ~ 1, method = "probit")
                            params[j,] <- c(cw.fit$zeta[1],-cw.fit$coefficients)
                            zeta[j,2:length(cw.fit$zeta)] <- cw.fit$zeta[-1]-cw.fit$zeta[1]
                            }
                     if(length(levels(y.fac)) == 2) {
                            if(!is.null(X) & num.lv > 0) cw.fit <- glm(y.fac ~ X+index, family = binomial(link = "probit"))
                            if(is.null(X) & num.lv > 0) cw.fit <- glm(y.fac ~ index, family = binomial(link = "probit"))
                            if(!is.null(X) & num.lv == 0) cw.fit <- glm(y.fac ~ X, family = binomial(link = "probit"))
                            if(is.null(X) & num.lv == 0) cw.fit <- glm(y.fac ~ 1, family = binomial(link = "probit"))
                            params[j,] <- cw.fit$coef
                            }
                     }
              }

       out <- list(params=params,phi=phi)
       if(num.lv > 0) out$index <- index
       if(family == "ordinal") out$zeta <- zeta
       options(warn = 0)

       return(out)
       }

## Calculates the commonly used (1/2) theta'_j Lambda_i theta_j
## lambda is a n x num.lv x nmu.lv array, theta is p x num.lv matrix
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



## Calculates information matrix based on scores of VA log-likelihood

calc.infomat <- function(theta = NULL, beta0=NULL, env = NULL, row.params = NULL, vameans = NULL, lambda=NULL, phi = NULL, zeta = NULL, num.lv, family, Lambda.struc, row.eff, y, X = NULL,T = NULL, inter = NULL,trait=NULL, offset=0,get.fourth=TRUE,get.trait=TRUE) {
       n <- nrow(y); p <- ncol(y)
       if(!is.null(X)) num.X <- ncol(X)
       if(!is.null(T)) num.T <- ncol(T)
       trial.size <- 1
       if(family == "ordinal") {
              max.levels <- max(y)
              if(any(max.levels == 1) || all(max.levels == 2)) stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks")
              }
       if (is.null(offset))  offset <- matrix(0, nrow = n, ncol = p)

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
                     if(is.null(T)) { new.env <- matrix(x2[1:(p * num.X)],p,num.X); x2 <- x2[-(1:(p * num.X))] }
                     else { new.env <- x2[1:num.X]; x2 <- x2[-(1:num.X)]; new.env=matrix(rep(new.env,each=p),nrow=p) }
                     }
              new.inter <- NULL; if(!is.null(T) && get.fourth) { new.inter <- x2[1:(num.X*num.T)]; x2 <- x2[-(1:(num.X * num.T))] }
              new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
              new.trait <- NULL; if(!is.null(T) && get.trait) { new.trait <- x2[1:num.T]; x2 <- x2[-(1:num.T)]; new.trait=matrix(rep(new.trait,each=n),nrow=n)}
              new.phi <- NULL; if(family == "negative.binomial") { new.phi <- x2[1:p]; x2 <- x2[-(1:p)] }
              new.zeta <- NULL; if(family == "ordinal") {
                     new.zeta <- matrix(NA,p,max.levels - 2)
                     for(j in 1:p) { if(max(y[,j]) > 2) { new.zeta[j,1:(max(y[,j])-2)] <- x2[1:(max(y[,j])-2)]; x2 <- x2[-(1:(max(y[,j])-2))] } }
                     new.zeta <- cbind(0,new.zeta);
                     }

              eta.mat <- matrix(new.beta0,n,p,byrow=TRUE)  + offset
              if(!is.null(X)) eta.mat <- eta.mat + X%*%t(new.env)
              if(!is.null(T) && get.fourth) eta.mat <- eta.mat + matrix((kronecker(T,X) %*% new.inter),n,p)
              if(!is.null(T) && get.trait) eta.mat <- eta.mat +  new.trait %*% t(T)
              if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta)
              if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)

              grad.theta <- grad.beta0 <- grad.env <- grad.inter <- grad.phi <- grad.trait <- grad.row.params <- grad.zeta <- NULL;

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
                     if(!is.null(X)) {
                            for (l in 1:num.X) { sum1 <- sweep(y-exp(eta.mat),1,X[,l],"*")
                                   if(is.null(T)) { grad.env <- c(grad.env,colSums(sum1))} else { grad.env <- c(grad.env,sum(sum1)) }
                                   }
                     }
                     if(!is.null(T) && get.fourth) {
                       sum1 <- c(y-exp(eta.mat)) * kronecker(T,X)
                       grad.inter <- c(grad.inter,colSums(sum1))
                     }
                     if(!is.null(T) && get.trait) {
                       grad.trait<-colSums(y-exp(eta.mat))%*%T
                     }
                     if(row.eff) { grad.row.params <- rowSums(y - exp(eta.mat)) }
                     }

              if(family=="negative.binomial") {
                     phi.mat <- matrix(new.phi,n,p,byrow=TRUE)
                     eta.mat <- matrix(new.beta0,n,p,byrow=TRUE)  + offset
                     if(!is.null(X)) eta.mat <- eta.mat + X%*%t(new.env)
                     if(!is.null(T) && get.fourth) eta.mat <- eta.mat + matrix((kronecker(T,X) %*% new.inter),n,p)
                     if(!is.null(T) && get.trait) eta.mat <- eta.mat +  new.trait %*% t(T)
                     if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
                     if(num.lv > 0) mu.mat <- eta.mat + new.vameans %*% t(new.theta) - calc.quad(new.lambda,new.theta,Lambda.struc)$mat
                     if(num.lv > 0) mu.mat.noquad <- eta.mat + new.vameans %*% t(new.theta)
                     if(num.lv == 0) { mu.mat <- eta.mat;  mu.mat.noquad=eta.mat }
                     grad.beta0 <- colSums(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(mu.mat) + phi.mat)))

                     if(num.lv > 0) {
                            for(l in 1:num.lv) {
                                   if(Lambda.struc == "unstructured") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep(t(Lambda.theta[,,l]),1,-new.vameans[,l],"+") * (-(y + phi.mat) * (phi.mat / (exp(mu.mat) + phi.mat)))
                                   if(Lambda.struc == "diagonal") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep((new.lambda[,l] %*% t(new.theta[,l])),1,-new.vameans[,l],"+") * (-(y + phi.mat) * (phi.mat / (exp(mu.mat) + phi.mat)))
                                   grad.theta <- c(grad.theta,colSums(sum1))
                                   }
                            }
                            if(!is.null(X)) {
                                   for(l in 1:num.X) { sum1 <- sweep(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(mu.mat) + phi.mat)),1,X[,l],"*")
                                          if(is.null(T)) { grad.env <- c(grad.env,colSums(sum1)) } else { grad.env <- c(grad.env,sum(sum1)) }
                                          }
                            }
                            if(!is.null(T) && get.fourth) {
                              sum1 <- c(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(eta.mat) + phi.mat))) * kronecker(T,X)
                              grad.inter <- c(grad.inter,colSums(sum1))
                            }
                            if(row.eff) { grad.row.params <- rowSums(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(mu.mat) + phi.mat))) }
                            if(!is.null(T) && get.trait) {
                              grad.trait<-colSums(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(eta.mat) + phi.mat)))%*%T
                            }
                            grad.phi <- colSums(-mu.mat.noquad + 1 + log(phi.mat) - digamma(phi.mat) - log(1 + phi.mat / exp(mu.mat)) - (y + phi.mat)/(phi.mat + exp(mu.mat)) + digamma(y + phi.mat))
                            }

                     if(family=="binomial") {
                            grad.beta0 <- colSums(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat))),na.rm=TRUE)

                            if(num.lv > 0) {
                                   for(l in 1:num.lv) {
                                          if(Lambda.struc == "unstructured") sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat))/(pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
                                          if(Lambda.struc == "diagonal") sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,new.vameans[,l],"*") - (new.lambda[,l] %*% t(new.theta[,l]))
                                          grad.theta <- c(grad.theta,colSums(sum1,na.rm=TRUE)) }
                                          }

                            if(!is.null(X)) {
                                   for(l in 1:num.X) {
                                          sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,X[,l],"*")
                                          if(is.null(T)) { grad.env <- c(grad.env,colSums(sum1)) } else { grad.env <- c(grad.env,sum(sum1))}
                                          }
                            }
                            if(!is.null(T) && get.fourth) {
                              sum1 <- c(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10)) * kronecker(T,X)
                              grad.inter <- c(grad.inter,colSums(sum1))
                            }
                            if(!is.null(T) && get.trait) {
                              grad.trait<-colSums(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10))%*%T
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
                            if(!is.null(X)) {
                                   for(l in 1:num.X) { sum1 <- sweep(deriv.trunnorm,1,X[,l],"*");
                                          if(is.null(T)) { grad.env <- c(grad.env,colSums(sum1,na.rm=TRUE)) } else { grad.env <- c(grad.env,sum(sum1,na.rm=TRUE)) }
                                          }
                            }
                            if(!is.null(T) && get.fourth) {
                              sum1 <- c(deriv.trunnorm) * kronecker(T,X)
                              grad.inter <- c(grad.inter,colSums(sum1))
                            }
                            if(!is.null(T) && get.trait) {
                              grad.trait<-colSums(deriv.trunnorm)%*%T
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
       return(c(grad.theta, grad.beta0, grad.env, grad.inter, grad.row.params,grad.trait, grad.phi, grad.zeta))
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
              if(is.null(T)) {new.env <- matrix(x2[1:(p * num.X)],p,num.X); x2 <- x2[-(1:(p * num.X))]
              } else { new.env <- x2[1:num.X]; x2 <- x2[-(1:num.X)]; new.env <- matrix(rep(new.env,each=p),nrow=p) }
              }
       new.inter <- NULL; if(!is.null(T) && get.fourth) { new.inter <- x2[1:(num.X*num.T)]; x2 <- x2[-(1:(num.X * num.T))] }
       new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
       new.trait <- NULL; if(!is.null(T) && get.trait) { new.trait <- x2[1:num.T]; x2 <- x2[-(1:num.T)]; new.trait=matrix(rep(new.trait,each=n),nrow=n)}
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
       if(!is.null(X)) eta.mat <- eta.mat + X%*%t(new.env)
       if(!is.null(T) && get.fourth) eta.mat <- eta.mat + matrix((kronecker(T,X) %*% new.inter),n,p)
       if(!is.null(T) && get.trait) eta.mat <- eta.mat +  new.trait %*% t(T)
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
                     if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.env)
                     if(!is.null(T) && get.fourth) eta.mat <- eta.mat + matrix((kronecker(T,X) %*% new.inter),n,p)
                     if(!is.null(T) && get.trait) eta.mat <- eta.mat +  new.trait %*% t(T)
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
              if(is.null(T)) { new.env <- matrix(x2[1:(p * num.X)],p,num.X); x2 <- x2[-(1:(p * num.X))] }
              else { new.env <- x2[1:num.X]; x2 <- x2[-(1:num.X)]; new.env <- matrix(rep(new.env,each=p),nrow=p) }
              }
       new.inter <- NULL; if(!is.null(T) && get.fourth) { new.inter <- x2[1:(num.X*num.T)]; x2 <- x2[-(1:(num.X * num.T))] }
       new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
       new.trait <- NULL; if(!is.null(T) && get.trait) { new.trait <- x2[1:num.T]; x2 <- x2[-(1:num.T)]; new.trait=matrix(rep(new.trait,each=n),nrow=n)}
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
       if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.env)
       if(!is.null(T) && get.fourth) eta.mat <- eta.mat + matrix((kronecker(T,X) %*% new.inter),n,p)
       if(!is.null(T) && get.trait) eta.mat <- eta.mat +  new.trait %*% t(T)
       if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
       if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta)
       grad.theta <- grad.beta0 <- grad.env <- grad.inter <- grad.phi <- grad.row.params <-grad.trait <- grad.zeta <- NULL

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
                            if(Lambda.struc == "diagonal") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep((new.lambda[,l] %*% t(new.theta[,l])),1,new.vameans[,l],"+") * exp(eta.mat)
                            grad.theta <- c(grad.theta,colSums(sum1))
                            }
              }
              if(!is.null(X)) {
                     for(l in 1:num.X) { sum1 <- sweep(y-exp(eta.mat),1,X[,l],"*")
                            if(is.null(T)) { grad.env <- c(grad.env,colSums(sum1))} else { grad.env <- c(grad.env,sum(sum1)) }
                            }
              }
              if(!is.null(T) && get.fourth) {
                sum1 <- c(y-exp(eta.mat)) * kronecker(T,X)
                grad.inter <- c(grad.inter,colSums(sum1))
              }
              if(!is.null(T) && get.trait) {
                grad.trait<-colSums(y-exp(eta.mat))%*%T
              }
              if(row.eff) { grad.row.params <- rowSums(y-exp(eta.mat)) }
              }

       if(family=="negative.binomial") {
              phi.mat <- matrix(new.phi,n,p,byrow=TRUE)
              eta.mat <- matrix(new.beta0,n,p,byrow=TRUE)  + offset
              if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.env)
              if(!is.null(T) && get.fourth) eta.mat <- eta.mat + matrix((kronecker(T,X) %*% new.inter),n,p)
              if(!is.null(T) && get.trait) eta.mat <- eta.mat +  new.trait %*% t(T)
              if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
              if(num.lv > 0) mu.mat <- eta.mat + new.vameans %*% t(new.theta) - calc.quad(new.lambda,new.theta,Lambda.struc)$mat
              if(num.lv > 0) mu.mat.noquad <- eta.mat + new.vameans %*% t(new.theta)
              if(num.lv == 0) mu.mat <- eta.mat
              grad.beta0 <- colSums(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(mu.mat) + phi.mat)))

              if(num.lv > 0) {
                     for(l in 1:num.lv) {
                            if(Lambda.struc == "unstructured") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep(t(Lambda.theta[,,l]),1,-new.vameans[,l],"+") * (-(y + phi.mat) * (phi.mat / (exp(mu.mat) + phi.mat)))
                            if(Lambda.struc == "diagonal") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep((new.lambda[,l] %*% t(new.theta[,l])),1,-new.vameans[,l],"+") * (-(y + phi.mat) * (phi.mat / (exp(mu.mat) + phi.mat)))
                            grad.theta <- c(grad.theta,colSums(sum1))
                            }
                     }
              if(!is.null(X)) {
                     for(l in 1:num.X) { sum1 <- sweep(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(mu.mat) + phi.mat)),1,X[,l],"*")
                            if(is.null(T)) { grad.env <- c(grad.env,colSums(sum1))} else { grad.env <- c(grad.env,sum(sum1)) }
                            }
              }
              if(!is.null(T) && get.fourth) {
                sum1 <- c(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(eta.mat) + phi.mat))) * kronecker(T,X)
                grad.inter <- c(grad.inter,colSums(sum1))
              }
              if(!is.null(T) && get.trait) {
                grad.trait<-colSums(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(eta.mat) + phi.mat)))%*%T
              }
              if(row.eff) { grad.row.params <- rowSums(-phi.mat - (y + phi.mat) * (-phi.mat / (exp(mu.mat) + phi.mat))) }

              grad.phi <- colSums(-mu.mat.noquad + 1 + log(phi.mat) - digamma(phi.mat) - log(1 + phi.mat / exp(mu.mat)) - (y + phi.mat) / (phi.mat + exp(mu.mat)) + digamma(y + phi.mat))
              }

       if(family=="binomial") {
              grad.beta0 <- colSums(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),na.rm=TRUE)
              if(num.lv > 0) {
                     for(l in 1:num.lv) {
                     if(Lambda.struc == "unstructured") sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
                     if(Lambda.struc == "diagonal") sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,new.vameans[,l],"*") - (new.lambda[,l] %*% t(new.theta[,l]))
                     grad.theta <- c(grad.theta,colSums(sum1,na.rm=TRUE)) }
                     }
              if(!is.null(X)) {
                     for(l in 1:num.X) {
                            sum1 <- sweep(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10),1,X[,l],"*")
                            if(is.null(T)) { grad.env <- c(grad.env,colSums(sum1))} else { grad.env <- c(grad.env,sum(sum1)) }
                            }
              }
              if(!is.null(T) && get.fourth) {
                sum1 <- c(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10)) * kronecker(T,X)
                grad.inter <- c(grad.inter,colSums(sum1))
              }
              if(!is.null(T) && get.trait) {
                grad.trait<-colSums(dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10))%*%T
              }
              if(row.eff) { grad.row.params <- rowSums((dnorm(eta.mat) * (y - trial.size * pnorm(eta.mat)) / (pnorm(eta.mat) * (1 - pnorm(eta.mat)) + 1e-10)),na.rm=TRUE) }
              }

       if(family=="ordinal") {
              grad.zeta <- matrix(0,p,ncol(new.zeta))
              deriv.trunnorm <- matrix(0,n,p)
              for(j in 1:p) {
                     deriv.trunnorm[y[,j] == 1,j] <- -dnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j]) / pnorm(new.zeta[j,1] - eta.mat[y[,j] == 1,j])
                     deriv.trunnorm[y[,j] == max(y[,j]),j] <- dnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]) / (1 - pnorm(new.zeta[j,max(y[,j]) - 1] - eta.mat[y[,j] == max(y[,j]),j]))
                     j.levels <- (1:max(y[,j]))[-c(1,max(y[,j]))]
                     for(k in j.levels) {
                            deriv.trunnorm[y[,j] == k,j] <- (-dnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) + dnorm(new.zeta[j,k - 1] - eta.mat[y[,j] == k,j])) / (pnorm(new.zeta[j,k] - eta.mat[y[,j] == k,j]) - pnorm(new.zeta[j,k - 1] - eta.mat[y[,j] == k,j]))
                            }
                     }
              deriv.trunnorm[!is.finite(deriv.trunnorm)] <- 0
              grad.beta0 <- colSums(deriv.trunnorm,na.rm=TRUE)
              if(num.lv > 0) {
                     for(l in 1:num.lv) {
                            if(Lambda.struc == "unstructured") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
                            if(Lambda.struc == "diagonal") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - (new.lambda[,l]%*%t(new.theta[,l]))
                            grad.theta <- c(grad.theta,colSums(sum1,na.rm=TRUE))
                            }
                     }

              if(!is.null(X)) {
                     for(l in 1:num.X) { sum1 <- sweep(deriv.trunnorm,1,X[,l],"*")
                            if(is.null(T)) { grad.env <- c(grad.env,colSums(sum1,na.rm=TRUE))} else { grad.env <- c(grad.env,sum(sum1,na.rm=TRUE)) }
                            }
              }
              if(!is.null(T) && get.fourth) {
                sum1 <- c(deriv.trunnorm) * kronecker(T,X)
                grad.inter <- c(grad.inter,colSums(sum1))
              }
              if(!is.null(T) && get.trait) {
                grad.trait<-colSums(deriv.trunnorm)%*%T
              }
              if(row.eff) { grad.row.params <- rowSums(deriv.trunnorm,na.rm=TRUE) }

              for(j in 1:p) {
                     zeta0 <- new.zeta[j,]
                     grad.zeta[j,max(y[,j]) - 1] <- grad.zeta[j,max(y[,j]) - 1] - sum(dnorm(zeta0[max(y[,j]) - 1] - eta.mat[which(y[,j] == max(y[,j])),j]) / (1 - pnorm(zeta0[max(y[,j]) - 1] - eta.mat[which(y[,j] == max(y[,j])),j])))
                     j.levels <- (1:max(y[,j]))[-c(1,max(y[,j]))]
                     for(k in j.levels) {
                            grad.zeta[j,k] <- grad.zeta[j,k] + sum(dnorm(zeta0[k] - eta.mat[y[,j] == k,j]) / (pnorm(zeta0[k] - eta.mat[y[,j] == k,j]) - pnorm(zeta0[k - 1] - eta.mat[y[,j] == k,j])))
                            grad.zeta[j,k - 1] <- grad.zeta[j,k - 1] - sum(dnorm(zeta0[k - 1]-eta.mat[y[,j] == k,j]) / (pnorm(zeta0[k] - eta.mat[y[,j] == k,j]) - pnorm(zeta0[k - 1] - eta.mat[y[,j] == k,j])))
                            }
                     }

              grad.zeta <- grad.zeta[,-1]
              if(length(unique(apply(y,2,max))) > 1) { grad.zeta <- c(grad.zeta[-which(grad.zeta == 0)]) }
              }

  return(c(grad.theta, grad.beta0, grad.env, grad.inter, grad.row.params,grad.trait, grad.phi, grad.zeta))
}

       zero.cons <- which(theta==0)
       if(family == "ordinal") {
              zeta2 <- as.vector(t(zeta[,-1]))
              find.na.zeta <- which(is.na(zeta2)); if(length(find.na.zeta) > 0) zeta2 <- zeta2[-find.na.zeta]
              }
       if(family != "ordinal") { zeta2 <- NULL }

       A.mat <- -nd2(x0=c(theta,beta0,env,inter,row.params,trait,phi,zeta2), f = grad.all.mod, vameans = vameans, lambda = lambda)
       if(length(zero.cons) > 0) { A.mat <- A.mat[-zero.cons,]; A.mat <- A.mat[,-zero.cons] }

       w1=FALSE
       if(num.lv > 0) {
              D.mat.fun <- function(i3) {
                     if(Lambda.struc == "unstructured") { lambda.i <- lambda[i3,,][lower.tri(lambda[i3,,],diag=TRUE)] }
                     if(Lambda.struc == "diagonal") { lambda.i <- lambda[i3,]; }

              return(-nd2(f = grad.all.var, x0=c(vameans[i3,],lambda.i), vameans = vameans, lambda = lambda, mod.x = c(theta,beta0,env,inter,row.params,trait,phi,zeta2), sel.i = i3))
                     }

              unique.ind <- which(!duplicated(y))
              D.mat <- vector("list",n);
              for(k in 1:length(unique.ind)) {
                     subD <- D.mat.fun(i3=unique.ind[k]); match.seq <- which(apply(y,1,identical,y[unique.ind[k],]) == 1) ## Unfortunately replace() only works if you want to replace elements with numerics
                     if( any(is.nan(subD)) ){ logi<-is.nan(subD); subD[logi]=diag(dim(subD)[1])[logi]; w1=TRUE}
                     for(k2 in match.seq) { D.mat[[k2]] <- subD }
                     }
              D.mat <- as.matrix(Matrix::bdiag(D.mat)); rm(subD)
              B.mat <- vector("list",n)
              for(i3 in 1:length(unique.ind)) {
                     #cat("Onto row",i3,"\n")
                     if(Lambda.struc == "unstructured") { lambda.i <- lambda[unique.ind[i3],,][lower.tri(lambda[unique.ind[i3],,],diag=TRUE)] }
                     if(Lambda.struc == "diagonal") { lambda.i <- lambda[unique.ind[i3],]; }

                     subB.mat <- -nd2(f = grad.all.cross, x0=c(vameans[unique.ind[i3],],lambda.i), vameans = vameans, lambda = lambda, mod.x = c(theta,beta0,env,inter,row.params,trait,phi,zeta2), sel.i = unique.ind[i3])
                     if(length(zero.cons) > 0) { subB.mat <- subB.mat[-zero.cons,] }
                     if( any(is.nan(subB.mat)) ){ logi<-is.nan(subB.mat); subB.mat[logi]=0; w1=TRUE}
                     match.seq <- which(apply(y,1,identical,y[unique.ind[i3],]) == 1)
                     for(k2 in match.seq) { B.mat[[k2]] <- subB.mat }
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
              if(is.null(T)){beta.errors <- matrix(sd.errors[1:(p*ncol(X))],p,ncol(X));
              rownames(beta.errors) <- colnames(y); colnames(beta.errors) <- colnames(X);
              sd.errors <- sd.errors[-(1:(p*ncol(X)))]
              sdout$Xcoef <- beta.errors
              } else {
              env.errors <- sd.errors[1:ncol(X)]
              names(env.errors) <- colnames(X)
              sd.errors <- sd.errors[-(1:(ncol(X)))]
              sdout$Xcoef <- env.errors

              if(get.fourth) {
                inter.errors <- matrix(sd.errors[1:(ncol(T)*ncol(X))],ncol(T),ncol(X));
                sd.errors <- sd.errors[-(1:(ncol(T)*ncol(X)))]
                rownames(inter.errors) <- colnames(T);  colnames(inter.errors) <- colnames(X);
                sdout$fourth <- inter.errors
                }
              }
       }
       if(!is.null(T) && get.trait){
         trait.errors <- sd.errors[1:ncol(T)];
         names(trait.errors) <- colnames(T)
         sd.errors <- sd.errors[-(1:ncol(T))];
         sdout$Tcoef <- trait.errors
       }

       row.params.errors <- NULL;
       if(row.eff) {
         row.params.errors <- sd.errors[1:n]; sd.errors <- sd.errors[-(1:n)];
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

## A function to compute highly accurate first-order derivatives
## From Fornberg and Sloan (Acta Numerica, 1994, p. 203-267; Table 1, page 213)
## Adapted by Scott Foster from code nicked off the net 2007
nd2.la <- function(x0, func, m = NULL, D.accur = 2, ...) {
  D.n <- length(x0)
  if(is.null(m)) { D.f0 <- func(x0, ...); m <- length(D.f0) }
  if(D.accur == 2) { D.w <- tcrossprod(rep(1,m), c(-1/2,1/2)); D.co <- c(-1, 1) }
  else {
    D.w <- tcrossprod(rep(1,m), c(1/12, -2/3, 2/3, -1/12))
    D.co <- c(-2, -1, 1, 2)
  }
  D.n.c <- length(D.co)
  macheps <- .Machine$double.eps
  D.h <- c(macheps^(1/3)*abs(x0))
  D.deriv <- matrix(NA, nrow = m, ncol = D.n)

  for (ii in 1:D.n) {
    D.temp.f <- matrix(0, m, D.n.c)
    for (jj in 1:D.n.c) {
      D.xd <- x0 + D.h[ii]*D.co[jj]*(1:D.n == ii)
      D.temp.f[,jj] <- func( D.xd, ...)
    }
    D.deriv[,ii] <- rowSums(D.w*D.temp.f)/D.h[ii]
  }
  return(D.deriv)
}
