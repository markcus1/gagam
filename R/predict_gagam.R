#' @title Prediction from GAGAM
#'
#' @description Takes the result of variable and structure discovery search using a genetic algorithm (\code{gagam()} function) and produces predictions given a set of values for the explanatory variables.
#'
#' @param gagam \code{gagam} object (result of the \code{gagam()} function)
#' @param newdata Matrix or data frame with explanatory variables. The order of columns must be the same as that used to construct the gagam object. If missing, the function extracts in-sample predictions.
#' @param reduc If prediction should be from one of the reduced models, specify the index of reduction. Specify only one index at a time. Default is NULL.
#'
#' @return A vector of predictions.
#' @export
#'
#' @examples
#' N <- 1000
#' set.seed(123)
#' xdat2 <- matrix(rnorm(N*10,0,1),nrow=N,ncol=10)
#' ydat2 <- 4*xdat2[,1]+5*xdat2[,2]+6*xdat2[,3]+(xdat2[,4])^2 + 4 + rnorm(N,0,0.25)
#'
#' xdat <- xdat2[1:500,]
#' ydat <- ydat2[1:500]
#'
#' \dontrun{
#' example_gagam <- gagam(ydat,xdat,Kvar = 6,no_gen = 50)
#' }
#'
#' xdat <- xdat2[501:1000,]
#' ydat <- ydat2[501:1000]
#'
#' newdata <- as.data.frame(cbind(as.matrix(ydat),xdat))
#'
#' \dontrun{
#' predict_gagam(example_gagam,newdata)
#' }
#'
#' @import utils
#' @import stats
#' @import mgcv
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import doMC
predict_gagam <- function(gagam,newdata,reduc=NULL){
  if (is.null(reduc)){
    gam_obj <- gagam$best_gam

    if (missing(newdata)){
      predict.gam(gam_obj)
    } else {
      ndata <- as.data.frame(newdata)
      colnames(ndata) <- c("y",paste0("x",1:as.numeric(NCOL(ndata)-1)))

      predict.gam(gam_obj,newdata=ndata)
    }
  } else {
    if (is.null(gagam[[paste0("best_gam_red",reduc)]])){
      stop("gagam doesn't contain the specified reduced model")
    }
    if (length())
    gam_obj <- gagam[[paste0("best_gam_red",reduc)]]

    if (missing(newdata)){
      predict.gam(gam_obj)
    } else {
      ndata <- as.data.frame(newdata)
      colnames(ndata) <- c("y",paste0("x",1:as.numeric(NCOL(ndata)-1)))

      predict.gam(gam_obj,newdata=ndata)
    }

  }

}
