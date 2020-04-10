#' @title Extract Interaction Names
#'
#' @description Names the interactions in the output correctly. Internal use only.
#'
#' @keywords internal
#'
#' @param ints Vector like temp_extract_vars$linear_ints
#' @param xcolnames Vector of column names of matrix x.
#'
#' @return Vector of correctly named interactions
#'
#' @import utils
#' @import stats
ints_names_func <- function(ints,xcolnames){
  temp_ints <- ints
  if (length(temp_ints) > 0){
    for (i in 1:length(temp_ints)){
      temp_int <- temp_ints[i]
      if (temp_int == "0"){
        next
      }
      col_pos <- sapply(gregexpr(":",temp_int),"[",1)
      temp_int_1 <- as.numeric(substr(temp_int,2,as.numeric(col_pos-1)))
      temp_int_2 <- as.numeric(substr(temp_int,as.numeric(col_pos+2),10))
      temp_ints[i] <- paste0(xcolnames[temp_int_1],":",xcolnames[temp_int_2])
    }
    return(temp_ints)
  } else {
    return(temp_ints)
  }
}
