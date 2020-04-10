#' @title Plot GAGAM
#'
#' @description Plot method for \code{gagam} objects.
#'
#' @param x Fitted \code{gagam} object
#' @param reduc If you wish to summarize one of the reduced models, specify the index of reduction. Specify only one index at a time. Default is NULL.
#' @param ... Any other arguments to pass to \code{\link[mgcv]{plot.gam}}.
#'
#' @rdname plot.gagam
#' @export
#'
#' @import utils
#' @import stats
#' @import mgcv
plot.gagam <- function(x,..., reduc=NULL){
  gagam <- x
  if(is.null(reduc)){
    gam_obj <- gagam$best_gam
    plot.gam(gam_obj,...)
  } else {
    if (is.null(gagam[[paste0("best_gam_red",reduc)]])){
      stop("gagam object doesn't contain the specified reduced model")
    }
    gam_obj <- gagam[[paste0("best_gam_red",reduc)]]
    plot.gam(gam_obj,...)
  }
}
