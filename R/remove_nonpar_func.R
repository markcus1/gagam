#' @title Reduction of Nonparametric Terms
#'
#' @description Additional reduction of nonparametric terms after the genetic algorithm has finished. Internal use only.
#'
#' @keywords internal
#'
#' @param best_ind Vector contain the model.
#' @param all_dataframe Data frame with all variables.
#' @param Kvar Maximum number of variables allowed in the final model
#' @param Kint Maximum number of interactions allowed in the final model
#' @param k Basis dimension for nonparametric terms estimated using cubic splines.
#' @param bs Spline basis.
#' @param family Specifies the family for the gam (see ?family and ?family.mgcv). Default is gaussian().
#' @param method Specifies the metric for smoothing parameter selection (see ?gam). Default is "REML".
#' @param optimizer Specifies the numerical optimization algorithm for the gam (see ?gam). Default is c("outer","newton").
#' @param always_par Specifies which variables are always estimated parametrically
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
remove_nonpar_func <- function(best_ind,all_dataframe,Kvar,Kint,k,bs,family,method,optimizer,always_par){
  `%nin%` <- Negate(`%in%`)
  best_rmse <- biccalc(best_ind,all_dataframe,Kvar,Kint,k,bs,family,method,optimizer,always_par)
  if (Kint == 0){
    to_remove <- NULL
    for (i in 1:Kvar){
      if (best_ind[i] == "0"){
        next
      }

      best_ind_n <- best_ind
      best_ind_n[as.numeric(Kvar+Kint+i)] <- "0"

      bic_n <- biccalc(best_ind_n,all_dataframe,Kvar,Kint,k,bs,family,method,optimizer,always_par)

      if (as.numeric(bic_n-best_rmse) < 7){
        to_remove <- c(to_remove,i)
      }
    }
    if (length(to_remove) > 0){
      best_ind[as.numeric(Kvar+Kint+to_remove)] <- "0"
    }
    return(best_ind)
  } else {
    to_remove1 <- NULL
    to_remove2 <- NULL

    for (i in 1:Kvar){
      if (best_ind[i] == "0"){
        next
      }

      best_ind_n <- best_ind
      best_ind_n[as.numeric(Kvar+Kint+i)] <- "0"

      bic_n <- biccalc(best_ind_n,all_dataframe,Kvar,Kint,k,bs,family,method,optimizer,always_par)

      if (as.numeric(bic_n-best_rmse) < 7){
        to_remove1 <- c(to_remove1,i)
      }
    }
    for (i in 1:Kint){
      if (best_ind[as.numeric(Kvar+Kint+Kvar+i)] == "0"){
        next
      }

      best_ind_n <- best_ind
      best_ind_n[as.numeric(Kvar+Kint+Kvar+i)] <- "0"

      bic_n <- biccalc(best_ind_n,all_dataframe,Kvar,Kint,k,bs,family,method,optimizer,always_par)

      if (as.numeric(bic_n-best_rmse) < 7){
        to_remove2 <- c(to_remove2,i)
      }
    }
    if (length(to_remove1) > 0){
      best_ind[as.numeric(Kvar+Kint+to_remove1)] <- "0"
    }
    if (length(to_remove2) > 0){
      best_ind_n[as.numeric(Kvar+Kint+Kvar+to_remove2)] <- "0"
    }
    return(best_ind)
  }
}
