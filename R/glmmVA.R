#' Fitting Generalized Linear Mixed-Effects Models with Variational Approximation
#'
#' Wraps the \code{"\link{gllvm}"} function to fit a (univaroate) generalized linear mixed-effects model (GLMM) in a more familiar syntax. 
#' Both fixed and random effects are specified in \code{lme4}-style via the formula argument.
#'
#' @param formula a formula object describing both the fixed- and random-effects part of the model. A response should be present on the left-hand side of the operator, and otherwise passed as a named 'y' argument to the function. Random effects are written in \code{lme4}-style. Some structured random effects are supported, see \code{"\link{gllvm}"} for more information.
#' @param data an optional data frame containing the variables named in formula.
#' @param family a family as supported by the \code{"\link{gllvm}"} function.
#' @param ... other arguments passed onto the \code{"\link{gllvm}"} function.
#' 
#' @return An object of class "glmmVA" that inherits from the "gllvm" class.
#' 
#' @seealso \code{\link{gllvm}}
#' 
#' @author Bert van der Veen
#'
#' @examples
#' data(eSpider)
#' data <- cbind(data.frame(y = c(eSpider$abund[eSpider$nonNA,]),
#'                          species = factor(rep(1:ncol(eSpider$abund), each = length(eSpider$nonNA))),
#'                          site = factor(rep(1:length(eSpider$nonNA), ncol(eSpider$abund)))),
#'                          do.call(rbind, replicate(ncol(eSpider$abund), scale(eSpider$X[eSpider$nonNA,]), simplify = FALSE)))
#' # Example 1: crossed random slope effects                          
#' model <- glmmVA(y~species + ConWate + ConHumu + (0+ConHumu|species) + (0+ConWate|species), family = "poisson", data = data)
#'
#' # Example 2: correlated random slopes
#' model <- glmmVA(y~species + ConWate + ConHumu + (0+ConWate+ConHumu|species), family = "poisson", data = data)
#' 
#' # Example 3: 
#' @export
glmmVA <- function(formula, data, family, ...) {
  
  args  <- list(...)
    
  mf <- model.frame(subbars1(formula), data = data)
  y <- as.matrix(model.response(mf))
  offset <- model.offset(mf)
  formula <- no.offset(formula)
    
  if(is.null(y) && !"y" %in% names(args)){
    stop("Response variable not found")
  }else if("y" %in% names(args)){
    y <- args[[names(args)]]
    args[[-which("y"%in%names(args))]]
  }else if(!is.null(y)){
    row.eff.formula <- formula[-2]
  }
  
  X = model.frame(subbars1(row.eff.formula), data)
  
  # Facilitate y to be passed as a 2-column matrix for binomial responses
  if(is.matrix(y) && family %in% c("binomial", "ZIB", "ZNIB") && !"Ntrials" %in% names(args)){
    Ntrials <- rowSums(y)
    y <- y[, 1,drop=FALSE]
  }else if(is.matrix(y) && family %in% c("binomial", "ZIB", "ZNIB") && !"Ntrials"){
    stop("Cannot both pass 'y' as a matrix and provide the 'Ntrials' argument.")
  }
  
  object <- do.call(gllvm, c(list(y = y, studyDesign = X, family = family, num.lv = 0, row.eff = row.eff.formula, offset = offset), args))
  
  object$call <-  match.call()
  class(object) <-  c("glmmVA", "gllvm")
  
  return(object)
}
