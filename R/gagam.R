#' @title Genetic Algorithm for Generalized Additive Models
#'
#' @description Implements the genetic algorithm for simultaneous variable selection and structure discovery in generalized additive models. For a given dependent variable and a set of explanatory variables, the genetic algorithm determines which regressors should be included linearly, which nonparametrically, and which should be excluded from the regression equation. The aim is to minimize the Bayesian Information Criterion value of the model.
#'
#' @param y Vector, matrix, or data frame containing observations of the dependent variable.
#' @param x Matrix or data frame containing all considered explanatory variables.
#' @param pop_size Size of the population (in multiples of 500). Default is 500.
#' @param Kvar Maximum number of variables allowed in the final model. Default is 15.
#' @param Kint Maximum number of interactions allowed in the final model. Default is 0.
#' @param no_gen Number of generations until convergence. Default is 100.
#' @param p_m Mutation rate for variables. Default is 0.05.
#' @param p_int Mutation rate of interactions of variables. Default is 0.1.
#' @param p_nonpar Mutation rate for the linear/nonparametric indicators for variables. Default is 0.1.
#' @param p_int_nonpar Mutation rate for the linear/nonparametric indicators for interactions. Default is 0.1.
#' @param multicore Whether to use multiple cores in computation. Strongly recommended but may not work on Windows. Default is TRUE.
#' @param k Basis dimension for nonparametric terms estimated using cubic splines. Default is 10.
#'
#' @return A list containing: gam object from the mgcv package (estimated best model), vector of indexes of variables included linearly, vector of indexes of variables included nonparametrically (and the same lists for interactions if Kint is greater than 0).
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
gagam <- function(y,x,pop_size = 500,Kvar = 15,Kint = 0,no_gen = 100,p_m = 0.05,p_int = 0.1,p_nonpar=0.1,p_int_nonpar=0.1,multicore=TRUE,k=10){
  if(missing(y)){
    stop("Dependent variable (y) is missing.")
  }
  if(missing(x)){
    stop("Matrix of explanatory variables is missing.")
  }
  if (is.matrix(y) == FALSE & is.data.frame(y) == FALSE & is.vector(y) == FALSE){
    stop("y must be a vector, a matrix, or a data frame.")
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
    stop("The maximum number of variables included in the model (Kvar) is greater than the number of variables in x.")
  }
  if (Kvar == 1){
    stop("Kvar should be more than 1.")
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
    no_cores <- detectCores()
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
        rmse2 <- mclapply(initpop_split2,biccalc,all_data,Kvar,Kint,k,mc.preschedule = TRUE,mc.cores=getOption("mc.cores",no_cores))
      } else {
        rmse2 <- lapply(initpop_split2,biccalc,all_data,Kvar,Kint,k)
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
        return(biccalc(initpop_split[[i]],all_data))
      }, warning = function(w) {
      }, error = function(e) {
        return(10000)
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
    rmse_fin_pop <- mclapply(initpop_split,biccalc,all_data,Kvar,Kint,k,mc.preschedule = TRUE,mc.cores=getOption("mc.cores",no_cores))
  } else {
    rmse_fin_pop <- lapply(initpop_split,biccalc,all_data,Kvar,Kint,k)
  }

  rmse_num <- as.numeric(rmse_fin_pop)

  #Remember best individual
  if (min(rmse_num) < best_rmse){
    best_rmse <- min(rmse_num)
    best_ind <- initpop[which(rmse_num==min(rmse_num))[1],]
  }

  best_ind_copy <- best_ind

  #######Prepare results
  ##Extract included vars
  best_ind_mains <- best_ind[1:Kvar]
  best_ind_nonpars <- best_ind[as.numeric(Kvar+Kint+1):as.numeric(Kvar+Kint+Kvar)]

  nonlinear_mains_final <- best_ind_mains[which(best_ind_nonpars == "1")]
  linear_mains_final <- best_ind_mains[which(best_ind_mains %nin% c(nonlinear_mains_final,0))]

  if (Kint > 0){
    best_ind_ints <- best_ind[as.numeric(Kvar+1):as.numeric(Kvar+Kint)]
    best_ind_ints_nonpars <- best_ind[as.numeric(Kvar+Kint+Kvar+1):as.numeric(Kvar+Kint+Kvar+Kint)]

    nonlinear_ints_final <- best_ind_ints[which(best_ind_ints_nonpars == "1")]
    linear_ints_final <- best_ind_ints[which(best_ind_ints %nin% c(nonlinear_ints_final,0))]
  }

  results_list <- as.list(NULL)

  results_list$best_gam <- biccalc2(best_ind,all_data,Kvar,Kint,k)

  results_list$linear_mains <- linear_mains_final
  results_list$nonparametric_mains <- nonlinear_mains_final

  if (Kint > 0){
    results_list$linear_interactions <- linear_ints_final
    results_list$nonparametric_interactions <- nonlinear_ints_final
  }

  return(results_list)

}
