
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Genetic Algorithm for Generalized Additive Models

<!-- badges: start -->

<!-- badges: end -->

This package implements a genetic algorithm which performs simulatenous
variable selection and structure discovery in generalized additive
models of the form:
$$g(\mathbb{E}(y_i)) = \alpha + \sum_{j\in \mathcal{L}} \beta_j x_{ij} + \sum_{k \in \mathcal{N}} f_k(x_{ik}) + \varepsilon_i$$
For a given dependent variable and a set of explanatory variables, the
genetic algorithm determines which regressors should be included
linearly (set $\mathcal{L}$), which nonparametrically (set
$\mathcal{N}$), and which should be excluded from the regression
equation. The aim is to minimize the Bayesian Information Criterion
value of the model.

For more details on the motivation behind GAGAM and how it works see
the paper: Cus, Mark. 2020. “Simultaneous Variable Selection And
Structure Discovery In Generalized Additive Models”.
<https://github.com/markcus1/gagam/blob/master/GAGAMpaper.pdf>.

A pdf `gagam` package manual is available here:
<https://github.com/markcus1/gagam/blob/master/gagam_package_manual.pdf>.

## Getting Started

### Installation

`gagam` can be installed and loaded with:

``` r
#library(devtools)
#install_github("markcus1/gagam")
library(gagam)
#> Loading required package: mgcv
#> Loading required package: nlme
#> This is mgcv 1.8-31. For overview type 'help("mgcv-package")'.
#> Loading required package: parallel
#> Loading required package: foreach
#> Loading required package: doParallel
#> Loading required package: iterators
#> Loading required package: doMC
```

### Simple Example

Suppose we have a dependent variable, $y$, ten explanatory variables,
$x_1,x_2,\ldots,x_{10}$, and the true model only contains the first
four. Moreover, $x_1$, $x_2$, $x_3$ have a linear effect on $y$ while
$x_4$ has a nonlinear relationship with the dependent variable. We now
demonstrate how to use `gagam()` to extract the true model.

Begin by generating the data:

``` r
N <- 500 #number of observations
set.seed <- 125
xdat <- matrix(rnorm(N*10,0,1),nrow=N,ncol=10) #matrix of explanatory variables (each drawn from N(0,1))
ydat <- 4*xdat[,1]+5*xdat[,2]+6*xdat[,3]+(xdat[,4])^2 + 4 + rnorm(N,0,0.25) #generate y according to the true model
```

We can now run GAGAM using:

``` r
example_gagam <- gagam(y=ydat,x=xdat,Kvar=6,no_gen=50,multicore=FALSE)
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
#> y ~ x3 + x1 + x2 + s(x4, bs = "cr", k = 10)
#> 
#> Estimated degrees of freedom:
#> 8.1  total = 12.1 
#> 
#> REML score: 52.69499     
#> 
#> $linear_mains
#> [1] "x3" "x1" "x2"
#> 
#> $nonparametric_mains
#> [1] "x4"
#> 
#> attr(,"class")
#> [1] "gagam"
```

The output is a list. The first element contains the optimal model as
determined by GAGAM as a fitted `mgcv::gam` object. The second and third
elements show which variables the algorithm included in the model
linearly and nonparametrically, respectively.

We can see that GAGAM correctly extracted the true model -
$x_1$, $x_2$, $x_3$ are included linearly, $x_4$ is included
nonparametrically, and the other variables are excluded.

## Arguments

As the example above demonstrates, the main function implementing GAGAM
is `gagam()`. This section explains the options available to the user
when implementing `gagam()`.

### Data

`gagam()` has only two required arguments: `y` and `x`.

`y` should contain the values of the dependent variable. It can take the
form of a vector, matrix, data frame, or factor. `y` should always be
univariate - even if a matrix or a data frame is used, it should always
have only one column.

`x` should contain the observations of all explanatory variables. It can
be a matrix or a data frame.

Note that if the columns of `x` are named then the explanatory variables
will inherit those names. Otherwise, the explanatory variables will
automatically be named `x1,x2,...,xp`.

### Algorithm Specification

The user is allowed to specify a number of algorithm parameters\*:

`no_gen` sets the number of iterations (generations) of the algorithm.
By default, this is 100 and that is usually sufficient. For smaller
problems `no_gen` might be set to a lower number, e.g. 50.

`pop_size` is the population size. By default, this is set to 500 but
the user can increase it in steps of 500 (e.g. to 1000, 1500, 2000, …).
`pop_size` should always be a multiple of 500.

`Kvar` sets the maximum number of variables allowed in the model. In the
example above we had 10 explanatory variables but we decided that we
will only choose among models with six or fewer predictors. This
argument is especially important if we have more explanatory variables
than observations ($p>n$). Then, `Kvar` should always be set to a
number lower than $n$, otherwise the algorithm might arrive at an
overdetermined model and return an error.

`Kint` extends the model from the introduction to also allow first-order
interactions. The argument specifies the maximum number of interactions
allowed in the model. By default this is set to 0, i.e. no interactions
allowed. `gagam` with `Kint` set to more than 0 is usually very (\!)
slow and is thus not recommended.

`p_m` is the mutation rate for the main effects (predictors).

`p_nonpar` is the mutation rate for the functional forms.

`p_int` is the mutation rate for the interactions.

`p_int_nonpar` is the mutation rate for functional forms of
interactions.

`multicore` is a logical which specifies whether multiple CPU cores
should be used in fitting `gagam`. Currently implemented using
`mclapply` and so will not work on Windows. By default this is set to
`TRUE`. Sometimes you might get an error “Error in mcfork: unable to
fork, possible reason: Resource temporarily unavailable”. If this
happens, try restarting the R session.

`cores` sets the number of CPU cores to use with `multicore`. By default
`gagam` will use all available cores.

\*Refer to the paper mentioned at the top for more details on these
parameters.

### Model and Estimation Specification

`gagam` is essentially a wrapper for the `mgcv` package and particularly
the `mgcv::gam` function. Thus, the user can pass arguments to `gam()`
which determine how the model should be specified and estimated.

`bs` is the spline basis used to estimate nonparametric terms. By
default, we use natural cubic splines but other bases might be more
appropriate in certain applications. See `mgcv::smooth.terms` for an
overview of what is available.

`k` is the basis dimension. This is set to 10 by default. Currently, all
nonparametric terms are estimated using the same basis and basis
dimension. Also, note that due to the way generalized additive models
are estimated, `Kvar` times `k` should be less than the number of
observations, otherwise `mgcv::gam` might return the error “Model has
more coefficients than observations”.

`family` specifies the distribution of the dependent variable and the
link function. Gaussian() is the default. Binomial(link=“logit”) might
be used for classification. Any family supported by `mgcv::gam` is
allowed.

`method` specifies the criterion for choosing the smoothing parameter.
By default, we use the restricted maximum likelihood (REML). Any
`method` supported by `mgcv::gam` is allowed.

`optimizer` specifies the optimization method to optimize the criterion
given by `method`. See `mgcv::gam` for more details.

`always_par` allows the user to force some terms to always enter the
model parametrically. This is useful if we have noncontinuous predictors
such as dummy variables. In the extreme, one might make all variables
parametric and use GAGAM just for variable selection. `always_par`
should be specified as a numerical vector where the elements give the
column numbers in `x` of variables to be estimated parametrically, e.g.
`always_par=c(1,2,3)` will ensure that the variables in the first three
columns of `x` are always parametric.

### Model Reduction

As mentioned in the introduction and in the paper, GAGAM uses the BIC to
judge the quality of the models. In small samples some noise variables
can reduce the BIC and get included in the recommended model. To counter
this, `gagam` implements three variable elimination methods (chosen
using the argument `reduc`):

`reduc=c(1)` loops through the variables in the model returned by
`gagam`, removing one at a time and noting the changes in BIC. In the
end, all variables whose removal increases the BIC by less than 7 are
eliminated from the model.

`reduc=c(3)` is similar to `reduc=c(1)` but tries making nonparametric
variables linear.

`reduc=c(2)` applies `reduc=c(1)` and then `reduc=c(3)`.

More than one reduction can be specified, e.g. `reduc=c(1,2,3)`. `gagam`
will return the original recommended model and the reduced models.

`reduc` should generally only be used when the number of explanatory
variables is very large and in small samples.

## Summary, Plot, Predict

Summary, plot, and predict methods are available for objects returned by
`gagam`. These essentially just take the fitted `mgcv::gam` object
contained in the `gagam` output and pass it to
summary.gam/plot.gam/predict.gam. If one of the `reduc` options was
used, one can specify which model to use (`NULL` for the original model,
`reduc = 1` for the model reduced using `reduc=c(1)` etc.). We also
allow the user to pass through additional arguments for
`plot.gam`/`predict.gam`.

``` r
summary(example_gagam)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> y ~ x3 + x1 + x2 + s(x4, bs = "cr", k = 10)
#> 
#> Parametric coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  4.98238    0.01136   438.5   <2e-16 ***
#> x3           5.99686    0.01095   547.6   <2e-16 ***
#> x1           4.01328    0.01122   357.7   <2e-16 ***
#> x2           5.00179    0.01137   439.9   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>         edf Ref.df    F p-value    
#> s(x4) 8.104  8.706 1410  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> R-sq.(adj) =  0.999   Deviance explained = 99.9%
#> -REML = 52.695  Scale est. = 0.064403  n = 500
plot(example_gagam)
```

## Boston Housing Example

The following code shows how to apply GAGAM to the Boston Housing data
set.

``` r
library(mlbench)
data(BostonHousing)

y <- BostonHousing[,14]
x <- BostonHousing[,1:13]

#boston_gagam <- gagam(y,x,Kvar = 13,always_par = c(2,4,9))

#summary(boston_gam)
```
