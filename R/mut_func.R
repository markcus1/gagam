#' @title Model Mutation
#'
#' @description Mutates both the variables included and the model structure. Internal use only.
#'
#' @param mutation_row Vector containing a model.
#' @param Kvar Maximum number of variables allowed in the final model
#' @param Kint Maximum number of interactions allowed in the final model
#' @param pm Mutation rate for variables.
#' @param pint Mutation rate for interactions of variables.
#' @param pnonpar Mutation rate for model structure of main terms.
#' @param pintnonpar Mutation rate for model structure of interaction terms.
#' @param regressors Total number of variables considered.
#'
#' @return Vector containing the mutated model.
#'
#' @examples
#' mutation_row <- c("x1","0","x4","x3","0","x8","1","0","0","0","0","1")
#' Kvar <- 10
#' Kint <- 0
#' pm <- 0.05
#' pint <- 0.1
#' pnonpar <- 0.1
#' pintnonpar <- 0.1
#' regressors <- 20
#'
#' \dontrun{
#' mut_func(mutation_row,Kvar,Kint,pm,pint,pnonpar,pintnonpar,regressors)
#' }
#'
#' @import utils
#' @import stats
#' @import mgcv
mut_func <- function(mutation_row,Kvar,Kint,pm,pint,pnonpar,pintnonpar,regressors){
  `%nin%` <- Negate(`%in%`)
  mains <- mutation_row[1:Kvar]
  nonpar <- mutation_row[as.numeric(Kvar+Kint+1):as.numeric(Kvar+Kint+Kvar)]

  #mains
  N <- length(which(mains!="0"))

  if (N == 0){
    mutant <- mutation_row
    return(mutant)
  }

  all_regressors <- paste0("x",1:regressors)

  basket_mains <- all_regressors[all_regressors %nin% mains]

  mask <- sample(c(1,0),Kvar,replace=TRUE,prob = c(pm,1-pm))

  mask_mains_not_zero <- as.numeric(mains!="0")

  mask_muts_nonzeros <- as.numeric((mask+mask_mains_not_zero)==2)

  mask_muts_zeros <- mask-mask_muts_nonzeros

  mask_for_nonpar <- abs(mask-1)
  mask_for_nonpar <- mask_for_nonpar - as.numeric(mains=="0")
  mask_for_nonpar[mask_for_nonpar %nin% c(0,1)] <- 0

  mains[which(mask_muts_nonzeros==1)] <- "0"
  nonpar[which(mask_muts_nonzeros==1)] <- "0"

  if (length(basket_mains) >= sum(mask_muts_zeros)){
    mains[which(mask_muts_zeros==1)] <- sample(basket_mains,sum(mask_muts_zeros),replace=FALSE)
    nonpar[which(mask_muts_zeros==1)] <- sample(c("1","0"),sum(mask_muts_zeros),replace=TRUE,prob=c(0.5,0.5))
  } else {
    len_basket <- length(basket_mains)
    len_diff <- sum(mask_muts_zeros) - len_basket
    sampled_mains <- c(basket_mains,replicate(len_diff,c(0)))
    sampled_nonpars <- sample(c("1","0"),len_basket,replace=TRUE,prob=c(0.5,0.5))
    sampled_nonpars <- c(sampled_nonpars,replicate(len_diff,c(0)))
    mains[which(mask_muts_zeros==1)] <- sampled_mains
    nonpar[which(mask_muts_zeros==1)] <- sampled_nonpars
  }

  curr_nonpar <- nonpar[which(mask_for_nonpar==1)]

  mask_nonpar <- sample(c(1,0),length(curr_nonpar),replace=TRUE,prob = c(as.numeric(pnonpar),as.numeric(1-pnonpar)))

  nonpar[which(mask_for_nonpar==1)] <- abs(as.numeric(curr_nonpar)-mask_nonpar)

  #ints

  if (Kint > 0){

    ints <- mutation_row[as.numeric(Kvar+1):as.numeric(Kvar+Kint)]
    ints_nonpar <- mutation_row[as.numeric(Kvar+Kint+Kvar+1):as.numeric(Kvar+Kint+Kvar+Kint)]

    if (length(as.numeric(substr(mains[which(mains != "0")],2,100)))>1){
      combin <- combn(sort(as.numeric(substr(mains[which(mains != "0")],2,100))),2)
      legalints <- paste0("x",combin[1,],":x",combin[2,])
    } else {
      legalints <- NULL
    }

    legalints_plus_zero <- c(legalints,"0")

    pos <- which(ints %nin% legalints_plus_zero)

    ints[pos] <- "0"
    ints_nonpar[pos] <- "0"

    basket_ints <- legalints[legalints %nin% ints]

    mask <- sample(c(1,0),Kint,replace=TRUE,prob=c(pint,1-pint))

    mask_ints_not_zero <- as.numeric(ints!="0")

    mask_muts_nonzeros <- as.numeric((mask+mask_ints_not_zero)==2)

    mask_muts_zeros <- mask-mask_muts_nonzeros

    mask_for_nonpar <- abs(mask-1)
    mask_for_nonpar <- mask_for_nonpar - as.numeric(ints=="0")
    mask_for_nonpar[mask_for_nonpar %nin% c(0,1)] <- 0

    ints[which(mask_muts_nonzeros==1)] <- "0"
    ints_nonpar[which(mask_muts_nonzeros==1)] <- "0"

    if (length(basket_ints) >= sum(mask_muts_zeros)){
      ints[which(mask_muts_zeros==1)] <- sample(basket_ints,sum(mask_muts_zeros),replace=FALSE)
      ints_nonpar[which(mask_muts_zeros==1)] <- sample(c("1","0"),sum(mask_muts_zeros),replace=TRUE,prob=c(0.5,0.5))
    } else {
      len_basket <- length(basket_ints)
      len_diff <- sum(mask_muts_zeros) - len_basket
      sampled_ints <- c(basket_ints,replicate(len_diff,c(0)))
      sampled_nonpars <- sample(c("1","0"),len_basket,replace=TRUE,prob=c(0.5,0.5))
      sampled_nonpars <- c(sampled_nonpars,replicate(len_diff,c(0)))
      ints[which(mask_muts_zeros==1)] <- sampled_ints
      ints_nonpar[which(mask_muts_zeros==1)] <- sampled_nonpars
    }

    curr_nonpar <- ints_nonpar[which(mask_for_nonpar==1)]

    mask_nonpar <- sample(c(1,0),length(curr_nonpar),replace=TRUE,prob = c(as.numeric(pintnonpar),as.numeric(1-pintnonpar)))

    ints_nonpar[which(mask_for_nonpar==1)] <- abs(as.numeric(curr_nonpar)-mask_nonpar)

    mutant <- c(mains,ints,nonpar,ints_nonpar)
  } else {
    mutant <- c(mains,nonpar)
  }

  return(mutant)
}
