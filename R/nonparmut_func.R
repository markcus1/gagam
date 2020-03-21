#' @title Mutate Model Structure
#'
#' @description Mutates the model structure, leaving the included variables unchanged. Internal use only.
#'
#' @param mutation_row Vector containing a model.
#' @param Kvar Maximum number of variables allowed in the final model
#' @param Kint Maximum number of interactions allowed in the final model
#' @param pnonpar Mutation rate for model structure of main terms.
#' @param pintnonpar Mutation rate for model structure of interaction terms.
#'
#' @return Vector containing the mutated model.
#'
#' @examples
#' mutation_row <- c("x1","0","x4","x3","0","x8","1","0","0","0","0","1")
#' Kvar <- 10
#' Kint <- 0
#' pnonpar <- 0.1
#' pintnonpar <- 0.1
#'
#' \dontrun{
#' nonparmut_func(mutation_row,Kvar,Kint,pnonpar,pintnonpar)
#' }
#'
#' @import utils
#' @import stats
#' @import mgcv
nonparmut_func <- function(mutation_row,Kvar,Kint,pnonpar,pintnonpar){
  `%nin%` <- Negate(`%in%`)
  mains <- mutation_row[1:Kvar] #Extract mains
  nonpar <- mutation_row[as.numeric(Kvar+Kint+1):as.numeric(Kvar+Kint+Kvar)] #Extract nonpars

  mask <- sample(c(1,0),length(which(mains!="0")),replace=TRUE,prob = c(pnonpar,1-pnonpar)) #Create mutation mask
  nonpar[which(mains!="0")] <- abs(as.numeric(nonpar[which(mains!="0")]) - mask) #Mutate according to the mask

  if (Kint > 0){
    ints <- mutation_row[as.numeric(Kvar+1):as.numeric(Kvar+Kint)]
    ints_nonpar <- mutation_row[as.numeric(Kvar+Kint+Kvar+1):as.numeric(Kvar+Kint+Kvar+Kint)]

    mask <- sample(c(1,0),length(which(ints!="0")),replace=TRUE,prob=c(pintnonpar,1-pintnonpar))
    ints_nonpar[which(ints!="0")] <- abs(as.numeric(ints_nonpar[which(ints!="0")]) - mask)

    mutant <- c(mains,ints,nonpar,ints_nonpar)
  } else {
    mutant <- c(mains,nonpar)
  }

  return(mutant)
}
