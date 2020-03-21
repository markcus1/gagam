#' @title Crossover Function
#'
#' @description Performs crossover between two models. Internal use only.
#'
#' @param crossover_matrix_row Vector containing both parent models.
#' @param Kvar Maximum number of variables allowed in the final model
#' @param Kint Maximum number of interactions allowed in the final model
#' @param pc Legacy variable (equal to 0.5)
#'
#' @return Vector containing two children models.
#'
#' @examples
#' p1 <- c("x1","0","x4","x3","0","x8","1","0","0","0","0","1")
#' p2 <- c("x1","x2","x4","x7","0","x8","1","0","1","1","0","0")
#' crossover_matrix_row <- c(p1,p2)
#' Kvar <- 6
#' Kint <- 0
#' pc <- 0.5
#'
#' \dontrun{
#' crossover_func(crossover_matrix_row,Kvar,Kint,pc)
#' }
#'
#' @import utils
#' @import stats
#' @import mgcv
crossover_func <- function(crossover_matrix_row,Kvar,Kint,pc){
  `%nin%` <- Negate(`%in%`)
  p1 <- crossover_matrix_row[1:as.numeric(2*Kvar+2*Kint)] #extract parent 1
  p2 <- crossover_matrix_row[as.numeric(2*Kvar+2*Kint+1):as.numeric(4*Kvar+4*Kint)] #extract parent 2

  p1_mains <- p1[1:Kvar] #extract mains from parent 1
  p2_mains <- p2[1:Kvar] #extract mains from parent 2

  if (Kint > 0){
    p1_ints <- p1[as.numeric(Kvar+1):as.numeric(Kvar+Kint)] #extract ints from parent 1
    p1_nonpar_i <- p1[as.numeric(Kvar+Kint+Kvar+1):as.numeric(Kvar+Kint+Kvar+Kint)] #extract the nonpars on ints from parent 1
    p2_ints <- p2[as.numeric(Kvar+1):as.numeric(Kvar+Kint)] #extract ints from parent 2
    p2_nonpar_i <- p2[as.numeric(Kvar+Kint+Kvar+1):as.numeric(Kvar+Kint+Kvar+Kint)] #extract the nonpars on ints from parent 2
  }

  p1_nonpar <- p1[as.numeric(Kvar+Kint+1):as.numeric(Kvar+Kint+Kvar)] #extract the nonpars on the mains from parent 1
  p2_nonpar <- p2[as.numeric(Kvar+Kint+1):as.numeric(Kvar+Kint+Kvar)] #extract the nonpars on the mains from parent 2

  all_common_mains <- intersect(p1_mains[which(p1_mains!="0")],p2_mains[which(p2_mains!="0")]) #find which mains are found in both models (including those not equal in nonpar terms)
  have_common_nonpar <- p1_nonpar[which(p1_mains %in% all_common_mains)] == p2_nonpar[which(p2_mains %in% all_common_mains)] #find which common mains also have common nonpar terms

  child1_mains <- NULL #initialize children
  child2_mains <- NULL
  child1_nonpar <- NULL
  child2_nonpar <- NULL
  if (length(all_common_mains) > 0){

    if (sum(have_common_nonpar) > 0){
      child1_mains <- all_common_mains[have_common_nonpar] #put the mains that are in both models and have the same nonpar status into both children
      child2_mains <- all_common_mains[have_common_nonpar] #put the mains that are in both models and have the same nonpar status into both children

      pos <- c(replicate(length(child1_mains),99))
      for (i in 1:length(child1_mains)){
        pos[i] <- which(child1_mains[i]==p1_mains) #error
      }

      child1_nonpar <- p1_nonpar[pos] #put the nonpars on variables which are in both models and have the same nonpar into the children nonpar
      child2_nonpar <- p1_nonpar[pos] #put the nonpars on variables which are in both models and have the same nonpar into the children nonpar
    }

    common_mains_wo_nopar <- all_common_mains[which(all_common_mains %nin% child1_mains)] #find which mains are in both models but with different nonpars
    crossed_nonpar1 <- NULL #initialize
    crossed_nonpar2 <- NULL

    if (length(common_mains_wo_nopar) > 0){
      pos1 <- c(replicate(length(common_mains_wo_nopar),99)) #initialize
      pos2 <- c(replicate(length(common_mains_wo_nopar),99))

      for (i in 1:length(common_mains_wo_nopar)){
        pos1[i] <- which(common_mains_wo_nopar[i]==p1_mains) #find positions (in parent 1) of vars that both models have in common but differ in nonpar
        pos2[i] <- which(common_mains_wo_nopar[i]==p2_mains) #find positions (in parent 2) of vars that both models have in common but differ in nonpar
      }

      #the followign lines perform crossover of nonpars on the common mains in vectorized form
      nonpar_mask1 <- sample(c(1,0),length(pos1),replace=TRUE,prob = c(pc,1-pc))
      nonpar_mask2 <- abs(nonpar_mask1-1)

      crossed_nonpar1 <- as.numeric(p1_nonpar)[pos1]*nonpar_mask1 + as.numeric(p2_nonpar)[pos2]*nonpar_mask2
      crossed_nonpar2 <- as.numeric(p1_nonpar)[pos1]*nonpar_mask2 + as.numeric(p2_nonpar)[pos2]*nonpar_mask1

    }

    child1_mains <- c(child1_mains,common_mains_wo_nopar) #put all common mains into the children
    child2_mains <- c(child2_mains,common_mains_wo_nopar)

    child1_nonpar <- c(child1_nonpar,crossed_nonpar1) #put the nonpars into the children
    child2_nonpar <- c(child2_nonpar,crossed_nonpar2)
  }

  p1_mains_wo_common <- p1_mains[p1_mains %nin% all_common_mains] #mains in parent 1 that were not in parent 2
  p2_mains_wo_common <- p2_mains[p2_mains %nin% all_common_mains] #mains in parent 2 that were not in parent 1

  p1_nonpar_wo_common <- p1_nonpar[p1_mains %nin% all_common_mains] #nonpars corresponding to the mains in parent 1 that were not in parent 2
  p2_nonpar_wo_common <- p2_nonpar[p2_mains %nin% all_common_mains] #nonpars corresponding to the mains in parent 2 that were not in parent 1

  p1_mains_wo_common_nonpar <- paste0(p1_mains_wo_common,"_",p1_nonpar_wo_common) #put the mains together with the nonpars
  p2_mains_wo_common_nonpar <- paste0(p2_mains_wo_common,"_",p2_nonpar_wo_common) #put the mains together with the nonpars

  length_temp <- length(p1_mains_wo_common_nonpar) #number of vars that are not common to both models

  mask <- sample(c(1,0),length_temp,replace=TRUE,prob = c(pc,1-pc)) #the following lines perform crossover in vectorized form

  child1_mains_second_part <- c(p1_mains_wo_common_nonpar[which(mask==1)],p2_mains_wo_common_nonpar[which(mask==0)])
  child2_mains_second_part <- c(p2_mains_wo_common_nonpar[which(mask==1)],p1_mains_wo_common_nonpar[which(mask==0)])

  child1_mains_second_part_mains <- substr(child1_mains_second_part,1,sapply(gregexpr("_",child1_mains_second_part),"[",1)-1)
  child1_mains_second_part_nonpar <- substr(child1_mains_second_part,sapply(gregexpr("_",child1_mains_second_part),"[",1)+1,sapply(gregexpr("_",child1_mains_second_part),"[",1)+1)

  child2_mains_second_part_mains <- substr(child2_mains_second_part,1,sapply(gregexpr("_",child2_mains_second_part),"[",1)-1)
  child2_mains_second_part_nonpar <- substr(child2_mains_second_part,sapply(gregexpr("_",child2_mains_second_part),"[",1)+1,sapply(gregexpr("_",child2_mains_second_part),"[",1)+1)

  child1_mains <- c(child1_mains,child1_mains_second_part_mains)
  child2_mains <- c(child2_mains,child2_mains_second_part_mains)

  child1_nonpar <- c(child1_nonpar,child1_mains_second_part_nonpar)
  child2_nonpar <- c(child2_nonpar,child2_mains_second_part_nonpar)

  sort1 <- sample(c(1:Kvar))
  sort2 <- sample(c(1:Kvar))

  child1_mains <- child1_mains[sort1] #randomize order of vars in both children
  child1_nonpar <- child1_nonpar[sort1]

  child2_mains <- child2_mains[sort2]
  child2_nonpar <- child2_nonpar[sort2]

  #Interactions
  if (Kint > 0){

    p1_c1_int <- p1_ints
    p2_c1_int <- p2_ints
    p1_c1_int_nonpar <- p1_nonpar_i
    p2_c1_int_nonpar <- p2_nonpar_i

    p1_c2_int <- p1_ints
    p2_c2_int <- p2_ints
    p1_c2_int_nonpar <- p1_nonpar_i
    p2_c2_int_nonpar <- p2_nonpar_i

    mains_in_c1 <- sort(as.numeric(substr(child1_mains[-which(child1_mains=="0")],2,100)))
    mains_in_c2 <- sort(as.numeric(substr(child2_mains[-which(child2_mains=="0")],2,100)))

    if (length(mains_in_c1) > 1){
      comb1 <- combn(mains_in_c1,2)
      legalint1 <- paste0("x",comb1[1,],":x",comb1[2,])
    } else {
      legalint1 <- NULL
    }

    if (length(mains_in_c2) > 1){
      comb2 <- combn(mains_in_c2,2)
      legalint2 <- paste0("x",comb2[1,],":x",comb2[2,])
    } else {
      legalint2 <- NULL
    }

    p1_c1_int_nonpar[p1_c1_int %nin% legalint1] <- "0" #from both parents, we set the nonpar status of ints illegal in child 1 to 0
    p2_c1_int_nonpar[p2_c1_int %nin% legalint1] <- "0"

    p1_c2_int_nonpar[p1_c2_int %nin% legalint2] <- "0" #from both parents, we set the nonpar status of ints illegal in child 2 to 0
    p2_c2_int_nonpar[p2_c2_int %nin% legalint2] <- "0"

    p1_c1_int[p1_c1_int %nin% legalint1] <- "0" #from both parents, we remove the interactions that are illegal in child 1
    p2_c1_int[p2_c1_int %nin% legalint1] <- "0"

    p1_c2_int[p1_c2_int %nin% legalint2] <- "0" #from both parents we remove interactions that are illegal in child 2
    p2_c2_int[p2_c2_int %nin% legalint2] <- "0"

    #child 1
    all_common_ints <- intersect(p1_c1_int[which(p1_c1_int != "0")],p2_c1_int[which(p2_c1_int != "0")]) #find which interactions are found in both parents (including those not equal in nonpar terms)
    have_common_nonpar <- p1_c1_int_nonpar[which(p1_c1_int %in% all_common_ints)] == p2_c1_int_nonpar[which(p2_c1_int %in% all_common_ints)] #find which common ints also have common nonpar terms

    c1_int <- NULL #initialize child
    c1_int_nonpar <- NULL
    if (length(all_common_ints) > 0){
      if (sum(have_common_nonpar) > 0){
        c1_int <- all_common_ints[have_common_nonpar] #put the ints that are in both models and have same nonpar into the child

        pos <- c(replicate(length(c1_int),99))
        for (i in 1:length(c1_int)){
          pos[i] <- which(c1_int[i]==p1_c1_int)
        }

        c1_int_nonpar <- p1_c1_int_nonpar[pos]

      }
    }

    common_ints_wo_nonpar <- all_common_ints[which(all_common_ints %nin% c1_int)]

    crossed_nonpar <- NULL

    if (length(common_ints_wo_nonpar) > 0){
      pos1 <- c(replicate(length(common_ints_wo_nonpar),99))
      pos2 <- c(replicate(length(common_ints_wo_nonpar),99))

      for (i in 1:length(common_ints_wo_nonpar)){
        pos1[i] <- which(common_ints_wo_nonpar[i]==p1_c1_int)
        pos2[i] <- which(common_ints_wo_nonpar[i]==p2_c1_int)
      }

      nonpar_mask1 <- sample(c(1,0),length(pos1),replace = TRUE,prob=c(pc,1-pc))
      nonpar_mask2 <- abs(nonpar_mask1-1)

      crossed_nonpar <- as.numeric(p1_c1_int_nonpar)[pos1]*nonpar_mask1 + as.numeric(p2_c1_int_nonpar)[pos2]*nonpar_mask2

    }

    c1_int <- c(c1_int,common_ints_wo_nonpar)
    c1_int_nonpar <- c(c1_int_nonpar,crossed_nonpar)

    p1_c1_int_wo_common <- p1_c1_int[p1_c1_int %nin% all_common_ints]
    p2_c1_int_wo_common <- p2_c1_int[p2_c1_int %nin% all_common_ints]

    p1_c1_int_nonpar_wo_common <- p1_c1_int_nonpar[p1_c1_int %nin% all_common_ints]
    p2_c1_int_nonpar_wo_common <- p2_c1_int_nonpar[p2_c1_int %nin% all_common_ints]

    p1_c1_wo_common_ints_and_nonpar <- paste0(p1_c1_int_wo_common,"_",p1_c1_int_nonpar_wo_common)
    p2_c1_wo_common_ints_and_nonpar <- paste0(p2_c1_int_wo_common,"_",p2_c1_int_nonpar_wo_common)

    length_temp <- length(p1_c1_wo_common_ints_and_nonpar)

    mask <- sample(c(1,0),length_temp,replace=TRUE,prob=c(pc,1-pc))

    c1_ints_second_part <- c(p1_c1_wo_common_ints_and_nonpar[which(mask==1)],p2_c1_wo_common_ints_and_nonpar[which(mask==0)])

    c1_ints_second_part_ints <- substr(c1_ints_second_part,1,sapply(gregexpr("_",c1_ints_second_part),"[",1)-1)
    c1_ints_second_part_nonpar <- substr(c1_ints_second_part,sapply(gregexpr("_",c1_ints_second_part),"[",1)+1,sapply(gregexpr("_",c1_ints_second_part),"[",1)+1)

    c1_int <- c(c1_int,c1_ints_second_part_ints)
    c1_int_nonpar <- c(c1_int_nonpar,c1_ints_second_part_nonpar)

    sort1 <- sample(c(1:Kint))

    c1_int <- c1_int[sort1]
    c1_int_nonpar <- c1_int_nonpar[sort1]

    #child 2
    all_common_ints <- intersect(p1_c2_int[which(p1_c2_int != "0")],p2_c2_int[which(p2_c2_int != "0")]) #find which interactions are found in both parents (including those not equal in nonpar terms)
    have_common_nonpar <- p1_c2_int_nonpar[which(p1_c2_int %in% all_common_ints)] == p2_c2_int_nonpar[which(p2_c2_int %in% all_common_ints)] #find which common ints also have common nonpar terms

    c2_int <- NULL #initialize child
    c2_int_nonpar <- NULL
    if (length(all_common_ints) > 0){
      if (sum(have_common_nonpar) > 0){
        c2_int <- all_common_ints[have_common_nonpar] #put the ints that are in both models and have same nonpar into the child

        pos <- c(replicate(length(c2_int),99))
        for (i in 1:length(c2_int)){
          pos[i] <- which(c2_int[i]==p1_c2_int)
        }

        c2_int_nonpar <- p1_c2_int_nonpar[pos]

      }
    }

    common_ints_wo_nonpar <- all_common_ints[which(all_common_ints %nin% c2_int)]

    crossed_nonpar <- NULL

    if (length(common_ints_wo_nonpar) > 0){
      pos1 <- c(replicate(length(common_ints_wo_nonpar),99))
      pos2 <- c(replicate(length(common_ints_wo_nonpar),99))

      for (i in 1:length(common_ints_wo_nonpar)){
        pos1[i] <- which(common_ints_wo_nonpar[i]==p1_c2_int)
        pos2[i] <- which(common_ints_wo_nonpar[i]==p2_c2_int)
      }

      nonpar_mask1 <- sample(c(1,0),length(pos1),replace = TRUE,prob=c(pc,1-pc))
      nonpar_mask2 <- abs(nonpar_mask1-1)

      crossed_nonpar <- as.numeric(p1_c2_int_nonpar)[pos1]*nonpar_mask1 + as.numeric(p2_c2_int_nonpar)[pos2]*nonpar_mask2

    }

    c2_int <- c(c2_int,common_ints_wo_nonpar)
    c2_int_nonpar <- c(c2_int_nonpar,crossed_nonpar)

    p1_c2_int_wo_common <- p1_c2_int[p1_c2_int %nin% all_common_ints]
    p2_c2_int_wo_common <- p2_c2_int[p2_c2_int %nin% all_common_ints]

    p1_c2_int_nonpar_wo_common <- p1_c2_int_nonpar[p1_c2_int %nin% all_common_ints]
    p2_c2_int_nonpar_wo_common <- p2_c2_int_nonpar[p2_c2_int %nin% all_common_ints]

    p1_c2_wo_common_ints_and_nonpar <- paste0(p1_c2_int_wo_common,"_",p1_c2_int_nonpar_wo_common)
    p2_c2_wo_common_ints_and_nonpar <- paste0(p2_c2_int_wo_common,"_",p2_c2_int_nonpar_wo_common)

    length_temp <- length(p1_c2_wo_common_ints_and_nonpar)

    mask <- sample(c(1,0),length_temp,replace=TRUE,prob=c(pc,1-pc))

    c2_ints_second_part <- c(p1_c2_wo_common_ints_and_nonpar[which(mask==1)],p2_c2_wo_common_ints_and_nonpar[which(mask==0)])

    c2_ints_second_part_ints <- substr(c2_ints_second_part,1,sapply(gregexpr("_",c2_ints_second_part),"[",1)-1)
    c2_ints_second_part_nonpar <- substr(c2_ints_second_part,sapply(gregexpr("_",c2_ints_second_part),"[",1)+1,sapply(gregexpr("_",c2_ints_second_part),"[",1)+1)

    c2_int <- c(c2_int,c2_ints_second_part_ints)
    c2_int_nonpar <- c(c2_int_nonpar,c2_ints_second_part_nonpar)

    sort1 <- sample(c(1:Kint))

    c2_int <- c2_int[sort1]
    c2_int_nonpar <- c2_int_nonpar[sort1]

    child1 <- c(child1_mains,c1_int,child1_nonpar,c1_int_nonpar)
    child2 <- c(child2_mains,c2_int,child2_nonpar,c2_int_nonpar)

  } else {

    child1 <- c(child1_mains,child1_nonpar)
    child2 <- c(child2_mains,child2_nonpar)

  }

  child <- c(child1,child2)
  return(child)
}
