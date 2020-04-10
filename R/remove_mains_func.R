#' @title Main Effects Reduction
#'
#' @description Additional variable reduction after the genetic algorithm has finished. Internal use only.
#'
#' @keywords internal
#'
#' @param best_ind Vector contain the model.
#' @param all_dataframe Data frame with all variables.
#' @param Kvar Maximum number of variables allowed in the final model
#' @param Kint Maximum number of interactions allowed in the final model
#' @param k Basis dimension for nonparametric terms estimated using cubic splines.
#' @param family Specifies the family for the gam (see ?family and ?family.mgcv). Default is gaussian().
#' @param method Specifies the metric for smoothing parameter selection (see ?gam). Default is "REML".
#' @param optimizer Specifies the numerical optimization algorithm for the gam (see ?gam). Default is c("outer","newton").
#'
#' @return Reduced model.
#'
#' @import utils
#' @import stats
#' @import mgcv
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import doMC
remove_mains_func <- function(best_ind,all_dataframe,Kvar,Kint,k,family,method,optimizer){
  `%nin%` <- Negate(`%in%`)
  best_rmse <- biccalc(best_ind,all_dataframe,Kvar,Kint,k,family,method,optimizer)
  if (Kint == 0){
    to_remove <- NULL
    for (i in 1:Kvar){
      if (best_ind[i] == "0"){
        next
      }
      best_ind_n <- best_ind
      best_ind_n[i] <- "0"
      best_ind_n[as.numeric(Kvar+Kint+i)] <- "0"
      bic_n <- biccalc(best_ind_n,all_dataframe,Kvar,Kint,k,family,method,optimizer)
      if (as.numeric(bic_n-best_rmse) < 7){
        to_remove <- c(to_remove,i)
      }

    }
    if (length(to_remove) > 0){
      best_ind[to_remove] <- "0"
      best_ind[to_remove+Kint+Kvar] <- "0"
    }
    return(best_ind)
  } else {
    #remove mains
    to_remove <- NULL
    for (i in 1:Kvar){
      if (best_ind[i] == "0"){
        next
      }
      best_ind_n <- best_ind
      best_ind_n[i] <- "0"
      best_ind_n[as.numeric(Kvar+Kint+i)] <- "0"

      #remove illegal ints
      mains <- best_ind_n[1:Kvar]
      ints <- best_ind_n[as.numeric(Kvar+1):as.numeric(Kvar+Kint)]

      if (length(as.numeric(substr(mains[which(mains != "0")],2,100)))>1){
        combin <- combn(sort(as.numeric(substr(mains[which(mains != "0")],2,100))),2)
        legalints <- paste0("x",combin[1,],":x",combin[2,])
      } else {
        legalints <- NULL
      }

      legalints_plus_zero <- c(legalints,"0")

      pos <- which(ints %nin% legalints_plus_zero)

      best_ind_n[as.numeric(Kvar+pos)] <- "0"

      best_ind_n[as.numeric(Kvar+Kint+Kvar+pos)] <- "0"

      bic_n <- biccalc(best_ind_n,all_dataframe,Kvar,Kint,k,family,method,optimizer)

      if (as.numeric(bic_n-best_rmse) < 7){
        to_remove <- c(to_remove,i)
      }

    }

    if (length(to_remove) > 0){
      #remove mains and nonpars on mains
      best_ind[to_remove] <- "0"
      best_ind[as.numeric(Kvar+Kint+to_remove)] <- "0"

      #remove illegal ints and nonpars on illegal ints

      mains <- best_ind[1:Kvar]
      ints <- best_ind[as.numeric(Kvar+1):as.numeric(Kvar+Kint)]

      if (length(as.numeric(substr(mains[which(mains != "0")],2,100)))>1){
        combin <- combn(sort(as.numeric(substr(mains[which(mains != "0")],2,100))),2)
        legalints <- paste0("x",combin[1,],":x",combin[2,])
      } else {
        legalints <- NULL
      }

      legalints_plus_zero <- c(legalints,"0")

      pos <- which(ints %nin% legalints_plus_zero)

      best_ind[as.numeric(Kvar+pos)] <- "0"

      best_ind[as.numeric(Kvar+Kint+Kvar+pos)] <- "0"
    }

    best_rmse <- biccalc(best_ind,all_dataframe,Kvar,Kint,k,family,method,optimizer)

    #remove ints
    to_remove <- NULL
    for (i in 1:Kint){
      if (best_ind[as.numeric(Kvar+i)] == "0"){
        next
      }
      best_ind_n <- best_ind
      best_ind_n[as.numeric(Kvar+i)] <- "0"
      best_ind_n[as.numeric(Kvar+Kint+Kvar+i)] <- "0"

      bic_n <- biccalc(best_ind_n,all_dataframe,Kvar,Kint,k,family,method,optimizer)
      if (as.numeric(bic_n - best_rmse) < 7){
        to_remove <- c(to_remove,i)
      }
    }
    if (length(to_remove) > 0){
      best_ind[as.numeric(Kvar+to_remove)] <- "0"
      best_ind[as.numeric(Kvar+Kint+Kvar+to_remove)] <- "0"
    }
    return(best_ind)
  }
}
