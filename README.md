
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Genetic Algorithm for Generalized Additive Models

<!-- badges: start -->

<!-- badges: end -->

This package implements a genetic algorithm which performs simulatenous
variable selection and structure discovery in generalized additive
models. For a given dependent variable and a set of explanatory
variables, the genetic algorithm determines which regressors should be
included linearly, which nonparametrically, and which should be excluded
from the regression equation. The aim is to minimize the Bayesian
Information Criterion value of the model.

## Installation

You can install the released version of gagam from GitHub with:

``` r
library(devtools)
install_github("markcus1/gagam")
```

## Example

This is a basic example which shows you how to use the package for
variable selection and structure discovery:

``` r
library(gagam)
#> Loading required package: mgcv
#> Loading required package: nlme
#> This is mgcv 1.8-31. For overview type 'help("mgcv-package")'.
#> Loading required package: parallel
#> Loading required package: foreach
#> Loading required package: doParallel
#> Loading required package: iterators
#> Loading required package: doMC
## basic example code
N <- 500
set.seed(123)
xdat <- matrix(rnorm(N*10,0,1),nrow=N,ncol=10)
ydat <- 4*xdat[,1]+5*xdat[,2]+6*xdat[,3]+(xdat[,4])^2 + 4 + rnorm(N,0,0.25)
example_gagam <- gagam(ydat,xdat,Kvar = 6,no_gen = 50)
#> [1] "Generation: 1 of 50"
#> [1] "Generation: 2 of 50"
#> [1] "Generation: 3 of 50"
#> [1] "Generation: 4 of 50"
#> [1] "Generation: 5 of 50"
#> [1] "Generation: 6 of 50"
#> [1] "Generation: 7 of 50"
#> [1] "Generation: 8 of 50"
#> [1] "Generation: 9 of 50"
#> [1] "Generation: 10 of 50"
#> [1] "Generation: 11 of 50"
#> [1] "Generation: 12 of 50"
#> [1] "Generation: 13 of 50"
#> [1] "Generation: 14 of 50"
#> [1] "Generation: 15 of 50"
#> [1] "Generation: 16 of 50"
#> [1] "Generation: 17 of 50"
#> [1] "Generation: 18 of 50"
#> [1] "Generation: 19 of 50"
#> [1] "Generation: 20 of 50"
#> [1] "Generation: 21 of 50"
#> [1] "Generation: 22 of 50"
#> [1] "Generation: 23 of 50"
#> [1] "Generation: 24 of 50"
#> [1] "Generation: 25 of 50"
#> [1] "Generation: 26 of 50"
#> [1] "Generation: 27 of 50"
#> [1] "Generation: 28 of 50"
#> [1] "Generation: 29 of 50"
#> [1] "Generation: 30 of 50"
#> [1] "Generation: 31 of 50"
#> [1] "Generation: 32 of 50"
#> [1] "Generation: 33 of 50"
#> [1] "Generation: 34 of 50"
#> [1] "Generation: 35 of 50"
#> [1] "Generation: 36 of 50"
#> [1] "Generation: 37 of 50"
#> [1] "Generation: 38 of 50"
#> [1] "Generation: 39 of 50"
#> [1] "Generation: 40 of 50"
#> [1] "Generation: 41 of 50"
#> [1] "Generation: 42 of 50"
#> [1] "Generation: 43 of 50"
#> [1] "Generation: 44 of 50"
#> [1] "Generation: 45 of 50"
#> [1] "Generation: 46 of 50"
#> [1] "Generation: 47 of 50"
#> [1] "Generation: 48 of 50"
#> [1] "Generation: 49 of 50"
#> [1] "Generation: 50 of 50"
example_gagam
#> $best_gam
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> y ~ x1 + x2 + x3 + s(x4, bs = "cr", k = 10)
#> 
#> Estimated degrees of freedom:
#> 8.37  total = 12.37 
#> 
#> REML score: 36.30185     
#> 
#> $linear_mains
#> [1] "x1" "x2" "x3"
#> 
#> $nonparametrics_mains
#> [1] "x4"
```

## Further Information

The help pages and the package manual contain more information on the
parameters one can specify.

For details on the algorithm see PAPER.
