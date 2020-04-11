#' @title Genetic Algorithm for Generalized Additive Models
#'
#' @description Implements the genetic algorithm for simultaneous variable selection and structure discovery in generalized additive models. For a given dependent variable and a set of explanatory variables, the genetic algorithm determines which regressors should be included linearly, which nonparametrically, and which should be excluded from the regression equation. The aim is to minimize the Bayesian Information Criterion value of the model.
#'
#' @param y Vector, matrix, data frame, or factor containing observations of the dependent variable.
#' @param x Matrix or data frame containing all considered explanatory variables. If the columns have names those will be used for variable names in the final output.
#' @param pop_size Size of the population (needs to be a multiple of 500). Default is 500.
#' @param Kvar Maximum number of variables allowed in the final model. Default is 15.
#' @param Kint Maximum number of interactions allowed in the final model. Default is 0.
#' @param no_gen Number of generations until convergence. Default is 100.
#' @param p_m Mutation rate for variables. Default is 0.05.
#' @param p_int Mutation rate for interactions of variables. Default is 0.1.
#' @param p_nonpar Mutation rate for the linear/nonparametric indicators for variables. Default is 0.1.
#' @param p_int_nonpar Mutation rate for the linear/nonparametric indicators for interactions. Default is 0.1.
#' @param multicore Whether to use multiple cores in computation. Strongly recommended but may not work on Windows. Default is TRUE.
#' @param cores Number of cores to use with multicore. Default (\code{NULL}) uses all cores.
#' @param k Basis dimension for nonparametric terms estimated using splines. Default is 10.
#' @param bs Spline basis for nonparametric terms. Specified as a two letter character string. Default is the natural cubic spline, bs="cr". See \code{\link[mgcv]{smooth.terms}} for an overview of what is available.
#' @param family Specifies the family for the gam (see \code{\link[stats]{family}} and \code{\link[mgcv]{family.mgcv}}). Default is gaussian().
#' @param method Specifies the metric for smoothing parameter selection (see \code{\link[mgcv]{gam}}). Default is "REML".
#' @param optimizer Specifies the numerical optimization algorithm for the gam (see \code{\link[mgcv]{gam}}). Default is c("outer","newton").
#' @param reduc Implements additional variable elimination methods at the end of the run of the genetic algorithm. User can choose between methods 1, 2, and 3. Multiple methods can be chosen. E.g. reduc=c(1) or reduc=c(1,3). See the \href{https://github.com/markcus1/gagam/blob/master/GAGAMpaper.pdf}{GAGAM paper} for an explanation of the methods. Default is NULL.
#' @param always_par Vector of the column numbers (in x) of the variables always estimated parametrically (for noncontinuous predictors).
#'
#' @return A list containing: \code{\link[mgcv]{gam}} object (fitted best model), vector of indexes or names of variables included linearly, vector of indexes or names of variables included nonparametrically (and the same lists for interactions if Kint is greater than 0).
#' @export
#'
#' @examples
#' N <- 500
#' set.seed(123)
#' xdat <- matrix(rnorm(N*10,0,1),nrow=N,ncol=10)
#' ydat <- 4*xdat[,1]+5*xdat[,2]+6*xdat[,3]+(xdat[,4])^2 + 4 + rnorm(N,0,0.25)
#'
#' \dontrun{
#' example_gagam <- gagam(ydat,xdat,Kvar = 6,no_gen = 50)
#' }
#'
#' @import utils
#' @import stats
#' @import mgcv
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import doMC
#'
#' @references Cus, Mark. 2020. "Simultaneous Variable Selection And Structure Discovery In Generalized Additive Models". https://github.com/markcus1/gagam/blob/master/GAGAMpaper.pdf.
gagam <- function(y,x,pop_size = 500,Kvar = 15,Kint = 0,no_gen = 100,p_m = 0.05,p_int = 0.1,p_nonpar=0.1,p_int_nonpar=0.1,multicore=TRUE,cores=NULL,k=10,bs="cr",family=gaussian(),method="REML",optimizer=c("outer","newton"),reduc=NULL,always_par=NULL){
  if(missing(y)){
    stop("Dependent variable is missing.")
  }
  if(missing(x)){
    stop("Matrix of explanatory variables is missing.")
  }
  if (is.matrix(y) == FALSE & is.data.frame(y) == FALSE & is.vector(y) == FALSE & is.factor(y) == FALSE){
    warning("y is not a vector, matrix, data frame, or factor - is it specified correctly?")
  }
  if (is.matrix(x) == FALSE & is.data.frame(x) == FALSE){
    stop("x must be a matrix or a data frame.")
  }
  y <- as.data.frame(y)
  if (NCOL(y) > 1){
    stop("The dependent variable must be univariate.")
  }
  xmatr <- as.data.frame(x)
  if (NROW(y) != NROW(xmatr)){
    stop("The lenghts of x and y do not match.")
  }
  if (Kvar > NCOL(xmatr)){
    Kvar <- NCOL(xmatr)
    warning("The maximum number of variables in the model (Kvar) was greater than the number of variables in x - Kvar is now set to NCOL(x).")
  }
  if (Kvar == 1){
    stop("Kvar should be more than 1.")
  }
  if (Kint > 0 & !is.null(always_par)){
    stop("Always parametric terms currently not supported when interactions are enabled.")
  }
  if (!is.null(colnames(x))){
    xcolnames <- colnames(x)
  } else {
    xcolnames <- NULL
  }
  x_names_m <- paste0("x",c(1:NCOL(xmatr)))
  all_data <- cbind(y,xmatr)
  colnames(all_data) <- c("y",x_names_m)
  if (pop_size %% 500 != 0){
    stop("Population size is not a multiple of 500.")
  }

  #########Additional specification
  regressors <- NCOL(xmatr)

  if (multicore == TRUE){
    if (!is.null(cores)){
      no_cores <- cores
    } else {
      no_cores <- detectCores()
    }
    registerDoParallel(cores=no_cores)
    registerDoMC(cores=no_cores)
  }

  `%nin%` <- Negate(`%in%`)

  percent_parents <- 0.5
  var10 <- 0.02 #elitism (e.g. top 10 out of 500)
  var240 <- 0.48 #top 240 out of 500 - don't receive crossover (well, they do but they're also kept separate)
  var40 <- 0.08 #top 40 out of 500
  var190 <- 19/24 #190 randomly out of 240
  var10_by_pop_size <- as.numeric(var10*pop_size)
  var240_by_pop_size <- as.numeric(var240*pop_size)
  var40_by_pop_size <- as.numeric(var40*pop_size)
  var190_by240_by_pop_size <- 190*(pop_size/500)
  pm <- p_m
  pint <- p_int
  pnonpar <- p_nonpar
  pintnonpar <- p_int_nonpar
  pc <- 0.5

  if (!is.null(always_par)){
    always_par <- paste0("x",always_par)
  }

  #########Algorithm
  bic_dict <- as.list(NULL)

  ########Generate initial population
  allvars <- paste0("x",1:regressors)
  allvars <- c(allvars,replicate(length(allvars),c(0)))

  temp_main <- sample(allvars,Kvar,replace=FALSE)

  temp_nonpar <- c(replicate(Kvar,c(0)))
  temp_nonpar[which(temp_main!="0")] <- sample(c(1,0),length(which(temp_main!="0")),replace = TRUE,prob = c(0.5,0.5))

  if (Kint > 0){
    if (length(temp_main[-which(temp_main=="0")]) >= 2){
      temp_main_no_zeros <- temp_main[-which(temp_main=="0")]
      temp_main_no_zeros <- as.numeric(substr(temp_main_no_zeros,2,100))
      temp_main_no_zeros <- sort(temp_main_no_zeros)
      comb <- combn(temp_main_no_zeros,2)
      allowed_int <- paste0("x",comb[1,],":x",comb[2,])
      allowed_int <- c(allowed_int,replicate(length(allowed_int),c(0)))
      if (length(allowed_int) < Kint){
        diff <- Kint - length(allowed_int)
        allowed_int <- c(allowed_int,replicate(diff,c(0)))
      }
      temp_int <- sample(allowed_int,Kint,replace=FALSE)

      no_nonzero_ints_here <- sum(temp_int!="0")
      pos_nonzero_ints_here <- which(temp_int!="0")

      if (no_nonzero_ints_here > 0){
        temp_int_nonpar <- replicate(Kint,c("0"))
        temp_int_nonpar[pos_nonzero_ints_here] <- sample(c(1,0),no_nonzero_ints_here,replace=TRUE,prob=c(0.5,0.5))
      } else {
        temp_int_nonpar <- replicate(Kint,c("0"))
      }

    } else {
      temp_int <- replicate(Kint,c(0))
      temp_int_nonpar <- replicate(Kint,c("0"))
    }

    initpop <- c(temp_main,temp_int,temp_nonpar,temp_int_nonpar)
  } else {
    initpop <- c(temp_main,temp_nonpar)
  }

  i <- 2

  while (i <= pop_size){
    allvars <- paste0("x",1:regressors)
    allvars <- c(allvars,replicate(length(allvars),c(0)))
    temp_main <- sample(allvars,Kvar,replace=FALSE)

    temp_nonpar <- c(replicate(Kvar,c(0)))
    temp_nonpar[which(temp_main!="0")] <- sample(c(1,0),length(which(temp_main!="0")),replace = TRUE,prob = c(0.5,0.5))

    if (Kint > 0){

      if (length(temp_main[-which(temp_main=="0")]) >= 2){
        temp_main_no_zeros <- temp_main[-which(temp_main=="0")]
        temp_main_no_zeros <- as.numeric(substr(temp_main_no_zeros,2,100))
        temp_main_no_zeros <- sort(temp_main_no_zeros)
        comb <- combn(temp_main_no_zeros,2)
        allowed_int <- paste0("x",comb[1,],":x",comb[2,])
        allowed_int <- c(allowed_int,replicate(length(allowed_int),c(0)))
        if (length(allowed_int) < Kint){
          diff <- Kint - length(allowed_int)
          allowed_int <- c(allowed_int,replicate(diff,c(0)))
        }
        temp_int <- sample(allowed_int,Kint,replace=FALSE)

        no_nonzero_ints_here <- sum(temp_int!="0")
        pos_nonzero_ints_here <- which(temp_int!="0")

        if (no_nonzero_ints_here > 0){
          temp_int_nonpar <- replicate(Kint,c("0"))
          temp_int_nonpar[pos_nonzero_ints_here] <- sample(c(1,0),no_nonzero_ints_here,replace=TRUE,prob=c(0.5,0.5))
        } else {
          temp_int_nonpar <- replicate(Kint,c("0"))
        }

      } else {
        temp_int <- replicate(Kint,c(0))
        temp_int_nonpar <- replicate(Kint,c("0"))
      }

      initpop <- rbind(initpop,c(temp_main,temp_int,temp_nonpar,temp_int_nonpar))

    } else {

      initpop <- rbind(initpop,c(temp_main,temp_nonpar))

    }

    i <- i + 1
  }

  if (Kint > 0){
    rm(allvars,temp_main,comb,allowed_int,temp_int,temp_nonpar,temp_int_nonpar,no_nonzero_ints_here,pos_nonzero_ints_here,temp_main_no_zeros)
  } else {
    rm(allvars,temp_main,temp_nonpar)
  }

  best_ind <- NULL
  best_rmse <- 10^6

  ########Generations
  gen <- 1
  while (gen <= no_gen){
    initpop_split <- lapply(seq_len(NROW(initpop)),function(x) initpop[x,])

    #######Calculate fitnesses
    ######Memoization
    if (multicore==TRUE){
      initpop_sort_for_dict <- mclapply(initpop_split,sorting_func,Kvar,Kint,mc.preschedule = TRUE,mc.cores=getOption("mc.cores",no_cores))
    } else {
      initpop_sort_for_dict <- lapply(initpop_split,sorting_func,Kvar,Kint)
    }
    initpop_sort_for_dict <- sapply(initpop_sort_for_dict,"[",1)

    temp <- bic_dict[initpop_sort_for_dict]

    types <- sapply(temp,typeof)

    in_dict <- as.numeric(which(types=="double"))
    not_in_dict <- as.numeric(which(types!="double"))

    rm(temp,types)

    rmse <- as.list(NULL) #not really rmse but bic (uses legacy notation)

    rmse[in_dict] <- bic_dict[initpop_sort_for_dict[in_dict]]

    if (length(not_in_dict) > 0){
      initpop2 <- initpop[not_in_dict,]

      if (NCOL(initpop2) != 1){
        initpop_split2 <- lapply(seq_len(NROW(initpop2)),function(x) initpop2[x,])
      } else {
        initpop_split2 <- as.list(NULL)
        initpop_split2[[1]] <- initpop2
      }

      if (multicore == TRUE){
        rmse2 <- mclapply(initpop_split2,biccalc,all_data,Kvar,Kint,k,bs,family,method,optimizer,always_par,mc.preschedule = TRUE,mc.cores=getOption("mc.cores",no_cores))
      } else {
        rmse2 <- lapply(initpop_split2,biccalc,all_data,Kvar,Kint,k,bs,family,method,optimizer,always_par)
      }

      rmse[not_in_dict] <- sapply(rmse2,"[",1)

      err <- as.numeric(which(sapply(sapply(rmse2,"[",1),typeof) != "double"))

      if (length(err) > 0){
        not_in_dict <- not_in_dict[-err]

        bic_dict[initpop_sort_for_dict[not_in_dict]] <- sapply(rmse2,"[",1)[-err]
      } else {
        bic_dict[initpop_sort_for_dict[not_in_dict]] <- sapply(rmse2,"[",1)
      }

    }

    errors <- NULL

    errors <- as.numeric(which(sapply(sapply(rmse,"[",1),typeof) != "double"))

    if (length(errors) == 0){
      errors <- NULL
    }

    for (i in errors){

      rmse[[i]] <- tryCatch({
        return(biccalc(initpop_split[[i]],all_data,Kvar,Kint,k,bs,family,method,optimizer,always_par))
      }, warning = function(w) {
      }, error = function(e) {
        return(99999)
      }, finally ={}
      )

    }

    rmse_num <- as.numeric(rmse)

    rmses_and_inds <- cbind(as.matrix(rmse_num),as.matrix(c(1:pop_size)))

    rmses_and_inds <- rmses_and_inds[order(rmses_and_inds[,1]),]

    #Remember best individual
    if (min(rmse_num) < best_rmse){
      best_rmse <- min(rmse_num)
      best_ind <- initpop[which(rmse_num==min(rmse_num))[1],]
    }

    #Misc
    best10indexes <- rmses_and_inds[1:var10_by_pop_size,2]
    best240indexes <- rmses_and_inds[1:var240_by_pop_size,2]
    best40indexes <- rmses_and_inds[1:var40_by_pop_size,2]

    #Start forming the sample of individuals to be mutatated
    pre_mut_pop <- initpop[best240indexes,]

    #Start forming new initpop
    newpop <- initpop[best10indexes,]

    #Crossover top 10 within themselves; add those crosses to the new pop
    best10split <- initpop[best10indexes,]
    best10split <- best10split[sample(NROW(best10split)),]
    best10split <- cbind(best10split[1:as.numeric(NROW(best10split)/2),],best10split[as.numeric(NROW(best10split)/2+1):as.numeric(NROW(best10split)),])
    best10split <- lapply(seq_len(NROW(best10split)),function(x) best10split[x,])

    best10crossed_over <- lapply(best10split,crossover_func,Kvar,Kint,pc)

    best10crossed_over <- as.matrix(cbind(sapply(best10crossed_over,"[",1:as.numeric(Kvar+Kint+Kvar+Kint)),sapply(best10crossed_over,"[",as.numeric(Kvar+Kint+Kvar+Kint+1):as.numeric(Kvar+Kint+Kvar+Kint+Kvar+Kint+Kvar+Kint))))
    best10crossed_over <- t(best10crossed_over)

    newpop <- rbind(newpop,best10crossed_over)

    #Mutate only the nonparametric bits of the top 10; add those mutants to the new pop
    best10split <- initpop[best10indexes,]
    best10split <- lapply(seq_len(NROW(best10split)),function(x) best10split[x,])

    best10mutated_nonpar <- lapply(best10split,nonparmut_func,Kvar,Kint,pnonpar,pintnonpar)

    best10mutated_nonpar <- as.matrix(sapply(best10mutated_nonpar,"["))
    best10mutated_nonpar <- t(best10mutated_nonpar)

    newpop <- rbind(newpop,best10mutated_nonpar)

    #From the sample for crossover
    random_crosses <- sample(best240indexes,var190_by240_by_pop_size,replace=TRUE) #var190*var240*pop_size
    crossover_pop <- initpop[random_crosses,]
    crossover_pop <- cbind(crossover_pop[1:as.numeric(NROW(crossover_pop)/2),],crossover_pop[as.numeric(NROW(crossover_pop)/2+1):as.numeric(NROW(crossover_pop)),])

    best40_crosses <- initpop[best40indexes,]
    best40_crosses <- best40_crosses[sample(NROW(best40_crosses)),]
    best40_crosses <- cbind(best40_crosses[1:as.numeric(NROW(best40_crosses)/2),],best40_crosses[as.numeric(NROW(best40_crosses)/2+1):as.numeric(NROW(best40_crosses)),])

    crossover_pop <- rbind(crossover_pop,best40_crosses)

    #Perform crossover
    crossover_pop_split <- lapply(seq_len(NROW(crossover_pop)),function(x) crossover_pop[x,])

    if (multicore == TRUE){
      crossed_over <- mclapply(crossover_pop_split,crossover_func,Kvar,Kint,pc,mc.preschedule = TRUE,mc.cores=getOption("mc.cores",no_cores))
    } else {
      crossed_over <- lapply(crossover_pop_split,crossover_func,Kvar,Kint,pc)
    }

    pre_mut_pop <- rbind(pre_mut_pop,t(sapply(crossed_over,"[",1:as.numeric(Kvar+Kint+Kvar+Kint))),t(sapply(crossed_over,"[",as.numeric(Kvar+Kint+Kvar+Kint+1):as.numeric(Kvar+Kint+Kvar+Kint+Kvar+Kint+Kvar+Kint))))

    #Mutation
    pre_mut_pop_split <- lapply(seq_len(NROW(pre_mut_pop)),function(x) pre_mut_pop[x,])

    if (multicore == TRUE){
      mutants <- mclapply(pre_mut_pop_split,mut_func,Kvar,Kint,pm,pint,pnonpar,pintnonpar,regressors,mc.preschedule = TRUE,mc.cores=getOption("mc.cores",no_cores))
    } else {
      mutants <- lapply(pre_mut_pop_split,mut_func,Kvar,Kint,pm,pint,pnonpar,pintnonpar,regressors)
    }

    mutants <- t(sapply(mutants,"["))

    #Form new population
    newpop <- rbind(newpop,mutants)

    initpop <- newpop

    #Go to next gen
    print(paste0("Generation: ",gen," of ",no_gen))
    gen <- gen + 1
  }
  initpop_split <- lapply(seq_len(NROW(initpop)),function(x) initpop[x,])

  if (multicore == TRUE){
    rmse_fin_pop <- mclapply(initpop_split,biccalc,all_data,Kvar,Kint,k,bs,family,method,optimizer,always_par,mc.preschedule = TRUE,mc.cores=getOption("mc.cores",no_cores))
  } else {
    rmse_fin_pop <- lapply(initpop_split,biccalc,all_data,Kvar,Kint,k,bs,family,method,optimizer,always_par)
  }

  rmse_num <- as.numeric(rmse_fin_pop)

  #Remember best individual
  if (min(rmse_num) < best_rmse){
    best_rmse <- min(rmse_num)
    best_ind <- initpop[which(rmse_num==min(rmse_num))[1],]
  }

  if (Kint == 0){
    best_ind_mains <- best_ind[1:Kvar]
    best_ind_nonpar <- best_ind[as.numeric(Kvar+1):as.numeric(Kvar+Kvar)]
    if (!is.null(always_par)){
      best_ind_nonpar[best_ind_mains %in% always_par] <- "0"
      best_ind <- c(best_ind_mains,best_ind_nonpar)
    }
  }


  #######Prepare results
  ##Extract included vars
  temp_extract_vars <- extract_terms(best_ind,Kvar,Kint)

  results_list <- as.list(NULL)

  if (is.null(xcolnames)){
    results_list$best_gam <- biccalc2(best_ind,all_data,Kvar,Kint,k,bs,family,method,optimizer,xcolnames,always_par)

    results_list$linear_mains <- temp_extract_vars$linear_mains
    results_list$nonparametric_mains <- temp_extract_vars$nonlinear_mains

    if (Kint > 0){
      results_list$linear_interactions <- temp_extract_vars$linear_ints
      results_list$nonparametric_interactions <- temp_extract_vars$nonlinear_ints
    }

  } else {
    results_list$best_gam <- biccalc2(best_ind,all_data,Kvar,Kint,k,bs,family,method,optimizer,xcolnames,always_par)

    results_list$linear_mains <- xcolnames[as.numeric(substr(temp_extract_vars$linear_mains,2,10))]
    results_list$nonparametric_mains <- xcolnames[as.numeric(substr(temp_extract_vars$nonlinear_mains,2,10))]

    if (Kint > 0){
      results_list$linear_interactions <- ints_names_func(temp_extract_vars$linear_ints,xcolnames)
      results_list$nonparametric_interactions <- ints_names_func(temp_extract_vars$nonlinear_ints,xcolnames)
    }

  }


  ##Reduction
  if (!is.null(reduc)){
    for (elem in reduc){
      if (elem == 1){
        best_ind_red1 <- remove_mains_func(best_ind,all_data,Kvar,Kint,k,bs,family,method,optimizer,always_par)
        temp_extract_vars <- extract_terms(best_ind_red1,Kvar,Kint)
        if (is.null(xcolnames)){
          results_list$best_gam_red1 <- biccalc2(best_ind_red1,all_data,Kvar,Kint,k,bs,family,method,optimizer,xcolnames,always_par)
          results_list$linear_mains_red1 <- temp_extract_vars$linear_mains
          results_list$nonparametric_mains_red1 <- temp_extract_vars$nonlinear_mains
          if (Kint > 0){
            results_list$linear_interactions_red1 <- temp_extract_vars$linear_ints
            results_list$nonparametric_interactions_red1 <- temp_extract_vars$nonlinear_ints
          }
        } else {
          results_list$best_gam_red1 <- biccalc2(best_ind_red1,all_data,Kvar,Kint,k,bs,family,method,optimizer,xcolnames,always_par)
          results_list$linear_mains_red1 <- xcolnames[as.numeric(substr(temp_extract_vars$linear_mains,2,10))]
          results_list$nonparametric_mains_red1 <- xcolnames[as.numeric(substr(temp_extract_vars$nonlinear_mains,2,10))]
          if (Kint > 0){
            results_list$linear_interactions_red1 <- ints_names_func(temp_extract_vars$linear_ints,xcolnames)
            results_list$nonparametric_interactions_red1 <- ints_names_func(temp_extract_vars$nonlinear_ints,xcolnames)
          }

        }
      } else if (elem == 2){
        best_ind_red2 <- remove_nonpar_func(remove_mains_func(best_ind,all_data,Kvar,Kint,k,bs,family,method,optimizer,always_par),all_data,Kvar,Kint,k,bs,family,method,optimizer,always_par)
        temp_extract_vars <- extract_terms(best_ind_red2,Kvar,Kint)
        if (is.null(xcolnames)){
          results_list$best_gam_red2 <- biccalc2(best_ind_red2,all_data,Kvar,Kint,k,bs,family,method,optimizer,xcolnames,always_par)
          results_list$linear_mains_red2 <- temp_extract_vars$linear_mains
          results_list$nonparametric_mains_red2 <- temp_extract_vars$nonlinear_mains
          if (Kint > 0){
            results_list$linear_interactions_red2 <- temp_extract_vars$linear_ints
            results_list$nonparametric_interactions_red2 <- temp_extract_vars$nonlinear_ints
          }
        } else {
          results_list$best_gam_red2 <- biccalc2(best_ind_red2,all_data,Kvar,Kint,k,bs,family,method,optimizer,xcolnames,always_par)
          results_list$linear_mains_red2 <- xcolnames[as.numeric(substr(temp_extract_vars$linear_mains,2,10))]
          results_list$nonparametric_mains_red2 <- xcolnames[as.numeric(substr(temp_extract_vars$nonlinear_mains,2,10))]
          if (Kint > 0){
            results_list$linear_interactions_red2 <- ints_names_func(temp_extract_vars$linear_ints,xcolnames)
            results_list$nonparametric_interactions_red2 <- ints_names_func(temp_extract_vars$nonlinear_ints,xcolnames)
          }
        }
      } else if (elem == 3){
        best_ind_red3 <- remove_nonpar_func(best_ind,all_data,Kvar,Kint,k,bs,family,method,optimizer,always_par)
        temp_extract_vars <- extract_terms(best_ind_red3,Kvar,Kint)
        if (is.null(xcolnames)){
          results_list$best_gam_red3 <- biccalc2(best_ind_red3,all_data,Kvar,Kint,k,bs,family,method,optimizer,xcolnames,always_par)
          results_list$linear_mains_red3 <- temp_extract_vars$linear_mains
          results_list$nonparametric_mains_red3 <- temp_extract_vars$nonlinear_mains
          if (Kint > 0){
            results_list$linear_interactions_red3 <- temp_extract_vars$linear_ints
            results_list$nonparametric_interactions_red3 <- temp_extract_vars$nonlinear_ints
          }
        } else {
          results_list$best_gam_red3 <- biccalc2(best_ind_red3,all_data,Kvar,Kint,k,bs,family,method,optimizer,xcolnames,always_par)
          results_list$linear_mains_red3 <- xcolnames[as.numeric(substr(temp_extract_vars$linear_mains,2,10))]
          results_list$nonparametric_mains_red3 <- xcolnames[as.numeric(substr(temp_extract_vars$nonlinear_mains,2,10))]
          if (Kint > 0){
            results_list$linear_interactions_red3 <- ints_names_func(temp_extract_vars$linear_ints,xcolnames)
            results_list$nonparametric_interactions_red3 <- ints_names_func(temp_extract_vars$nonlinear_ints,xcolnames)
          }
        }

      }
    }
  }

  to_return <- structure(results_list,class="gagam")

  return(to_return)

}
