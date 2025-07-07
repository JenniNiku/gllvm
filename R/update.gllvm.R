#' @title Update and Re-fit a gllvm Model Call
#' @description Models of class gllvm have three formula interfaces: formula, lv.formula, and row.eff. Each can include fixed and random effects.
#'
#'
#' @param object an object of class 'gllvm'.
#' @param formula interface for column-specific effects.
#' @param lv.formula interface for latent variable (reduced-rank) effects.
#' @param row.eff interface for row-specific (i.e., same for all columns) effects.
#' @param eval default to TRUE, in which case the reformulated model is fitted directly. if FALSE returns the new model call.
#' @param ...	not used.
#' @details
#' Updates the formula of a gllvm object.
#'
#' @author Bert van der Veen
#' 
#' @aliases update update.gllvm
#' @method update gllvm
#' @importFrom stats update
#' 
#' @export
#' @export update.gllvm
update.gllvm <- function(object, formula = NULL, lv.formula = NULL, row.eff = NULL, eval = TRUE, ...){
  if(inherits(object, "gllvm")){
    call <- getCall(object)  
  }else{
    call <- object
  }
  
  match.call(expand.dots = FALSE, call = call)
  
  args <- list(...)
  
  if(!missing(formula)){
    formula.new <- formula
    formula_ <- call$formula
    if(!is.null(formula_) && !is.null(formula.new)){
      call$formula <- update.formula(formula_, formula.new)
    }else{
      call$formula <- formula.new
    }
  }
  if(!missing(lv.formula)){
    lv.formula.new <- lv.formula
    lv.formula_ <- call$lv.formula
    
    if(!is.null(lv.formula_)  && !is.null(lv.formula.new)){
      call$lv.formula <- update.formula(lv.formula_, lv.formula.new)
    }else{
      call$lv.formula <- lv.formula.new
    }
  }
  if(!missing(row.eff) && is.language(row.eff)){
    row.eff.new <- row.eff
    row.eff_ <- call$row.eff
    
    if(!is.null(row.eff_) && !is.null(row.eff.new) && is.formula(row.eff)){
      call$row.eff <- update.formula(row.eff_, row.eff.new)
    }else{
      call$row.eff <- row.eff.new
    }
  }else if(!is.null(row.eff)){
    call[["row.eff"]] <- row.eff
  }
  
  # Add or override other arguments from ...
  for (arg_name in names(args)) {
    if(exists(arg_name, parent.frame())){
      call[[arg_name]] <- as.name(arg_name)
    }else {
      call[[arg_name]] <- args[[arg_name]]
    }
  }

  if(eval){
    eval(call, parent.frame())
  }else{
    call
  }
}