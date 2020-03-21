#' @title Model Sort
#'
#' @description Sorts the variables in the model and summarizes the variables included and whether they enter linearly or not in a string. For internal use in memoization.
#'
#' @param individual_row Vector containing the model
#' @param Kvar Maximum number of variables allowed in the final model
#' @param Kint Maximum number of interactions allowed in the final model
#'
#' @return A string summarizing the model.
#'
#' @examples
#' individual_row <- c("x1","0","x4","x3","0","x8","1","0","0","0","0","1")
#' Kvar <- 6
#' Kint <- 0
#'
#' \dontrun{
#' sorting_func(individual_row,Kvar,Kint)
#' }
#'
#' @import utils
#' @import stats
#' @import mgcv
sorting_func <- function(individual_row,Kvar,Kint){
  `%nin%` <- Negate(`%in%`)
  ind_to_sort <- individual_row
  mains <- ind_to_sort[1:Kvar] #extract main effects into a vector
  nonpar_m <- ind_to_sort[as.numeric(Kvar+Kint+1):as.numeric(Kvar+Kint+Kvar)] #extract the nonpar on mains into a vector

  mains <- substr(mains,2,100) #remove x's
  mains[which(mains=="")] <- 0 #put zeros in empty places
  mains <- as.numeric(mains) #make numeric

  sort_temp <- sort(mains,index.return=TRUE) #sort mains

  mains <- sort_temp$x #apply mains sort
  nonpar_m <- nonpar_m[sort_temp$ix] #sort nonpars on mains according to the same order as we used for the mains

  mains_sorted <- paste0(mains,collapse = "_") #put mains together into a string
  nonpar_m_sorted <- paste0(nonpar_m,collapse = "_") #put nonpars on mains together into a string

  if (Kint > 0){
    ints <- ind_to_sort[as.numeric(Kvar+1):as.numeric(Kvar+Kint)] #extract interactions into a vector
    nonpar_i <- ind_to_sort[as.numeric(Kvar+Kint+Kvar+1):as.numeric(Kvar+Kint+Kvar+Kint)] #extract nonpars interactions into a vector

    int_no_zeros <- ints[ints %nin% c(0)] #remove zeros
    if (length(int_no_zeros) == 0){ #if there were no interactions we just put zeros everywhere in both the ints vector and the nonpar on ints vector and we're done
      int_sorted <- paste(replicate(Kint,c(0)),collapse = "_")
      nonpar_i_sorted <- paste(replicate(Kint,c(0)),collapse = "_")
      sorted_ind <- paste(mains_sorted,int_sorted,nonpar_m_sorted,nonpar_i_sorted,sep = "-")
    } else {
      leng_int <- length(int_no_zeros) #how many interactions are there

      int_colon_pos <- sapply(gregexpr(":",int_no_zeros),"[",1) #extract the position of the colon in every int

      int_part1 <- as.matrix(as.numeric(substr(int_no_zeros,2,as.numeric(int_colon_pos-1)))) #extract the first part of every int
      int_part2 <- as.matrix(as.numeric(substr(int_no_zeros,as.numeric(int_colon_pos+2),100))) #extract the second part of every int

      int_matr <- cbind(int_part1,int_part2) #put those together into a matrix

      int_matr <- t(apply(int_matr,1,sort)) #sort so that the first variable in an interaction always has a smaller index, e.g. x1:x2 not x2:x1
      int_matr <- int_matr[order(int_matr[,1],int_matr[,2],decreasing=FALSE),] #sort so that the first interaction has a smaller second variable index, e.g. x1:x2 then x1:x3 not x1:x3 then x1:x2

      if (NCOL(int_matr) == 1){
        int_matr <- t(int_matr) #the matrix does form properly if there's only one int - this corrects that
      }

      #Sort the nonpar of ints
      ints_only_num_sorted <- ints #now we'll keep the order of ints the same as originally and just sort the indexes (just sort to make sure x1:x2 not x2:x1)

      int_colon_pos <- sapply(gregexpr(":",ints_only_num_sorted),"[",1)
      int_part1 <- substr(ints_only_num_sorted,2,as.numeric(int_colon_pos-1))
      int_part1[which(int_part1=="")] <- "0"
      int_part2 <- substr(ints_only_num_sorted,as.numeric(int_colon_pos+2),100)
      int_part1 <- t(as.numeric(int_part1))
      int_part2 <- t(as.numeric(int_part2))
      int_together <- rbind(int_part1,int_part2)
      int_part1 <- apply(int_together,2,min)
      int_part2 <- apply(int_together,2,max)
      ints_only_num_sorted <- paste0("x",int_part1,":x",int_part2)
      ints_only_num_sorted[which(ints_only_num_sorted=="x0:x0")] <- "0"

      int_sorted_w_xs <- paste0("x",int_matr[,1],":","x",int_matr[,2]) #write the sorted interactions in the same form as they were originally, i.e. like x1:x2 (part 1)
      int_sorted_w_xs <- c(as.numeric(replicate(as.numeric(Kint-leng_int),c(0))),int_sorted_w_xs) #part 2

      matches <- match(int_sorted_w_xs,ints_only_num_sorted) #find the positions of the ints in the sorted vector in the original vector

      nonpar_i_sorted <- nonpar_i[matches] #sort the nonpars on ints appropriately (i.e. according to the new order of ints)
      nonpar_i_sorted <- paste0(nonpar_i_sorted,collapse = "_") #put the nonpars on ints together into a string

      #

      int_sorted <- paste0(int_matr[,1],":",int_matr[,2]) #merge the interactions, e.g. 1 and 2 into 1:2

      int_sorted <- c(as.numeric(replicate(as.numeric(Kint-leng_int),c(0))),int_sorted) #add the appropriate number of zeros

      int_sorted <-  paste0(int_sorted,collapse = "_") #put them together into a string

      sorted_ind <- paste(mains_sorted,int_sorted,nonpar_m_sorted,nonpar_i_sorted,sep = "-") #put everything together into one string
    }
  } else {
    sorted_ind <- paste(mains_sorted,nonpar_m_sorted,sep = "-") #put everything together into one string
  }
  return(sorted_ind)
}
