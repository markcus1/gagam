#' @title BIC Calculator 2
#'
#' @description Alternative function to calculate the BIC of the model. Internal use only.
#'
#' @param pop_matrix_row Vector containing a model
#' @param all_dataframe Dataframe with all variables
#' @param Kvar Maximum number of variables allowed in the final model
#' @param Kint Maximum number of interactions allowed in the final model
#' @param k Basis dimension for nonparametric terms estimated using cubic splines.
#'
#' @return BIC value of the model
#'
#' @import utils
#' @import stats
#' @import mgcv
biccalc2 <- function(pop_matrix_row,all_dataframe,Kvar,Kint,k){ #gonna take each element in the list (using mclapply)
  `%nin%` <- Negate(`%in%`)
  mains <- pop_matrix_row[1:Kvar] #extract mains into a vector
  if (length(mains) == 0){ #if there are no mains then exit (return a huge number)
    return(99999)
  }
  if (sum(mains == "0") == Kvar){
    return(99999)
  }
  if (Kint > 0){
    ints <- pop_matrix_row[as.numeric(Kvar+1):as.numeric(Kvar+Kint)] #extract interactions into a vector
    nonparint <- pop_matrix_row[as.numeric(Kvar+Kint+Kvar+1):as.numeric(Kvar+Kint+Kvar+Kint)] #extract nonpars on interactions into a vector

    nonparint <- ints[which(nonparint!="0")] #stores nonparametric interactions
    ints <- ints[which(ints!="0")] #remove zeros in ints
    ints <- ints[which(ints %nin% nonparint)] #stores parametric interactions

    if (length(nonparint) > 0){
      colon_pos <- sapply(gregexpr(":",nonparint),"[",1) #extract the position of the colon in every int

      part1 <- as.matrix(substr(nonparint,1,as.numeric(colon_pos-1))) #matrix of first vars in each int
      part2 <- as.matrix(substr(nonparint,as.numeric(colon_pos+1),100)) #matrix of second vars in each int

      nonparint <- paste0("ti(",part1,",",part2,",bs='cr',k=",k,")") #put the ints together along with the wrapper

    } else {
      nonparint <- NULL
    }

  } else {
    ints <- NULL
    nonparint <- NULL
  }
  nonparmain <- pop_matrix_row[as.numeric(Kvar+Kint+1):as.numeric(Kvar+Kint+Kvar)] #extract nonpars on mains into a vector
  nonparmain <- mains[which(nonparmain!="0")]  #stores nonparametric mains
  mains <- mains[which(mains!="0")] #removes zeros in mains
  if (length(mains) == 0){
    mains <- NULL
  }
  mains <- mains[which(mains %nin% nonparmain)] #stores parametric mains
  if (length(nonparmain)>0){
    nonparmain <- paste0("s(",nonparmain,",bs='cr',k=",k,")") #adds the wrapper to the nonpar mains
  } else {
    nonparmain <- NULL
  }
  all_ind <- c(mains,ints,nonparmain,nonparint) #combines all terms together
  formula_f <- as.formula(paste0("y~",paste(all_ind,collapse = "+"))) #creates the formula

  gamfit <- gam(formula_f,data = all_dataframe,family=gaussian(),method = "REML",optimizer=c("outer","newton")) #estimate gam

  if (gamfit$converged == FALSE){ #sometimes the model won't converge on the first try but will on the second
    gamfit <- gam(formula_f,data = all_dataframe,family=gaussian(),method = "REML",optimizer=c("outer","newton")) #so we give it another try if it doesn't work
    if (gamfit$converged == FALSE){
      return(99999) #and if it still doesn't work we give up
    }
  }

  return(gamfit)
}
