#' @title Summarize GAGAM
#'
#' @description Summary method for \code{gagam} objects.
#'
#' @param object Fitted \code{gagam} object
#' @param reduc If you wish to summarize one of the reduced models, specify the index of reduction. Specify only one index at a time. Default is NULL.
#' @param ... Any other arguments to pass to \code{\link[mgcv]{summary.gam}}.
#'
#' @return Summary of the \code{gam} object.
#' @rdname summary.gagam
#' @export
#'
#' @import utils
#' @import stats
#' @import mgcv
summary.gagam <- function(object,..., reduc=NULL){
  gagam <- object
  if(is.null(reduc)){
    gam_obj <- gagam$best_gam
    return(summary.gam(gam_obj))
  } else {
    if (is.null(gagam[[paste0("best_gam_red",reduc)]])){
      stop("gagam object doesn't contain the specified reduced model")
    }
    gam_obj <- gagam[[paste0("best_gam_red",reduc)]]
    return(summary.gam(gam_obj))
  }
}
