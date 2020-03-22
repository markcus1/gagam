#' @title Extract Terms from Model
#'
#' @description Extracts linear and nonparametric variables from a model. Internal use only.
#'
#' @param individual Vector containing a model
#' @param Kvar Maximum number of variables allowed in the final model
#' @param Kint Maximum number of interactions allowed in the final model
#'
#' @return A list specifying which variables are included linearly and which are included nonparametrically.
#'
#' @import utils
#' @import stats
#' @import mgcv
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import doMC
extract_terms <- function(individual,Kvar,Kint){
  `%nin%` = Negate(`%in%`)

  best_ind_mains <- individual[1:Kvar]
  best_ind_nonpars <- individual[as.numeric(Kvar+Kint+1):as.numeric(Kvar+Kint+Kvar)]

  nonlinear_mains_final <- best_ind_mains[which(best_ind_nonpars == "1")]
  linear_mains_final <- best_ind_mains[which(best_ind_mains %nin% c(nonlinear_mains_final,0))]

  if (Kint > 0){
    best_ind_ints <- individual[as.numeric(Kvar+1):as.numeric(Kvar+Kint)]
    best_ind_ints_nonpars <- individual[as.numeric(Kvar+Kint+Kvar+1):as.numeric(Kvar+Kint+Kvar+Kint)]

    nonlinear_ints_final <- best_ind_ints[which(best_ind_ints_nonpars == "1")]
    linear_ints_final <- best_ind_ints[which(best_ind_ints %nin% c(nonlinear_ints_final,0))]
  }

  res_list <- as.list(NULL)

  res_list$linear_mains <- linear_mains_final
  res_list$nonlinear_mains <- nonlinear_mains_final

  if (Kint > 0){
    res_list$linear_ints <- linear_ints_final
    res_list$nonlinear_ints <- nonlinear_ints_final
  }

  return(res_list)

}
