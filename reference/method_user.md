# User-Defined Functions for Estimating Weights

This page explains the details of estimating weights using a
user-defined function. The function must take in arguments that are
passed to it by
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
and return a vector of weights or a list containing the weights.

To supply a user-defined function, the function object should be entered
directly to `method`; for example, for a function `fun`, `method = fun`.

### Point Treatments

The following arguments are automatically passed to the user-defined
function, which should have named parameters corresponding to them:

- `treat`: a vector of treatment status for each unit. This comes
  directly from the left hand side of the formula passed to
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  and so will have it's type (e.g., numeric, factor, etc.), which may
  need to be converted.

- `covs`: a data frame of covariate values for each unit. This comes
  directly from the right hand side of the formula passed to
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md).
  The covariates are processed so that all columns are numeric; all
  factor variables are split into dummies and all interactions are
  evaluated. All levels of factor variables are given dummies, so the
  matrix of the covariates is not full rank. Users can use
  [`make_full_rank()`](https://ngreifer.github.io/WeightIt/reference/make_full_rank.md),
  which accepts a numeric matrix or data frame and removes columns to
  make it full rank, if a full rank covariate matrix is desired.

- `s.weights`: a numeric vector of sampling weights, one for each unit.

- `ps`: a numeric vector of propensity scores.

- `subset`: a logical vector the same length as `treat` that is `TRUE`
  for units to be included in the estimation and `FALSE` otherwise. This
  is used to subset the input objects when `exact` is used. `treat`,
  `covs`, `s.weights`, and `ps`, if supplied, will already have been
  subsetted by `subset`.

- `estimand`: a character vector of length 1 containing the desired
  estimand. The characters will have been converted to uppercase. If
  "ATC" was supplied to estimand,
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  sets `focal` to the control level (usually 0 or the lowest level of
  `treat`) and sets `estimand` to "ATT".

- `focal`: a character vector of length 1 containing the focal level of
  the treatment when the estimand is the ATT (or the ATC as detailed
  above).
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  ensures the value of focal is a level of `treat`.

- `stabilize`: a logical vector of length 1. It is not processed by
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  before it reaches the fitting function.

None of these parameters are required to be in the fitting function.
These are simply those that are automatically available.

In addition, any additional arguments supplied to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
will be passed on to the fitting function.
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
ensures the arguments correspond to the parameters of the fitting
function and throws an error if an incorrectly named argument is
supplied and the fitting function doesn't include `...` as a parameter.

The fitting function must output either a numeric vector of weights or a
list (or list-like object) with an entry named wither "w" or "weights".
If a list, the list can contain other named entries, but only entries
named "w", "weights", "ps", and "fit.obj" will be processed. "ps" is a
vector of propensity scores and "fit.obj" should be an object used in
the fitting process that a user may want to examine and that is included
in the `weightit` output object as "obj" when `include.obj = TRUE`. The
"ps" and "fit.obj" components are optional, but "weights" or "w" is
required.

### Longitudinal Treatments

Longitudinal treatments can be handled either by running the fitting
function for point treatments for each time point and multiplying the
resulting weights together or by running a method that accommodates
multiple time points and outputs a single set of weights. For the
former,
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
can be used with the user-defined function just as it is with
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md).
The latter method is not yet accommodated by
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md),
but will be someday, maybe.

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)

## Examples

``` r
data("lalonde", package = "cobalt")

#A user-defined version of method = "ps"
my.ps <- function(treat, covs, estimand, focal = NULL, ...) {
  covs <- make_full_rank(covs)
  d <- data.frame(treat, covs)
  f <- formula(d)
  ps <- glm(f, data = d, family = "binomial")$fitted
  w <- get_w_from_ps(ps, treat = treat, estimand = estimand,
                     focal = focal)

  list(w = w, ps = ps)
}

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = my.ps, estimand = "ATT"))
#> A weightit object
#>  - method: "my.ps" (a user-defined method)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74
summary(W1)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                 Max
#> treated 1.                  ||              1.   
#> control 0.022 |---------------------------| 2.044
#> 
#> - Units with the 5 most extreme weights by group:
#>                                    
#>             1     2   3     4     5
#>  treated    1     1   1     1     1
#>           410   226 224   111    84
#>  control 1.33 1.437 1.5 1.637 2.044
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.000 0.000    0.00       0
#> control       0.823 0.701    0.33       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    255.99     185
cobalt::bal.tab(W1)
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.0199
#> age         Contin.   0.0459
#> educ        Contin.  -0.0360
#> married      Binary   0.0044
#> nodegree     Binary   0.0080
#> re74        Contin.  -0.0275
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    255.99     185

data("msmdata")
(W2 <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                        A_2 ~ X1_1 + X2_1 +
                          A_1 + X1_0 + X2_0,
                        A_3 ~ X1_2 + X2_2 +
                          A_2 + X1_1 + X2_1 +
                          A_1 + X1_0 + X2_0),
                   data = msmdata,
                   method = my.ps))
#> A weightitMSM object
#>  - method: "my.ps" (a user-defined method)
#>  - number of obs.: 7500
#>  - sampling weights: none
#>  - number of time points: 3 (A_1, A_2, A_3)
#>  - treatment:
#>     + time 1: 2-category
#>     + time 2: 2-category
#>     + time 3: 2-category
#>  - covariates:
#>     + baseline: X1_0, X2_0
#>     + after time 1: X1_1, X2_1, A_1, X1_0, X2_0
#>     + after time 2: X1_2, X2_2, A_2, X1_1, X2_1, A_1, X1_0, X2_0

summary(W2)
#>                         Time 1                        
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                   Max
#> treated 1.079 |---------------------------| 403.483
#> control 1.276 |-------------------|         284.764
#> 
#> - Units with the 5 most extreme weights by group:
#>                                                 
#>             3172    3065    2025    1938     731
#>  treated 166.992 170.555 196.414 213.193 403.483
#>             2301    1275    1145    1121     832
#>  control 155.625 168.964  172.42 245.882 284.764
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       1.914 0.816   0.649       0
#> control       1.706 0.862   0.670       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted 3306.    4194. 
#> Weighted    845.79   899.4
#> 
#>                         Time 2                        
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                   Max
#> treated 1.079 |---------------------------| 403.483
#> control 1.276 |----------------|            245.882
#> 
#> - Units with the 5 most extreme weights by group:
#>                                                 
#>             2902    1869    1779    1509    1313
#>  treated 168.964 170.555 196.414 284.764 403.483
#>             2684    2549    1250     911     620
#>  control 155.625 166.992  172.42 213.193 245.882
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       1.892 0.819   0.652       0
#> control       1.748 0.869   0.686       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted 3701.   3799.  
#> Weighted    912.87  829.87
#> 
#>                         Time 3                        
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                   Max
#> treated 1.079 |---------------------------| 403.483
#> control 1.276 |---------|                   148.155
#> 
#> - Units with the 5 most extreme weights by group:
#>                                                 
#>             1991    1254     893     668     468
#>  treated 196.414 213.193 245.882 284.764 403.483
#>             4479    4021    2455    2427     112
#>  control  88.072  97.827 104.623 121.845 148.155
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       1.832 0.975   0.785       0
#> control       1.254 0.683   0.412       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted 4886.   2614.  
#> Weighted   1900.26  600.12
#> 
cobalt::bal.tab(W2)
#> Balance summary across all time points
#>        Times    Type Max.Diff.Adj
#> X1_0 1, 2, 3 Contin.       0.0342
#> X2_0 1, 2, 3  Binary       0.0299
#> X1_1    2, 3 Contin.       0.0657
#> X2_1    2, 3  Binary       0.0299
#> A_1     2, 3  Binary       0.0262
#> X1_2       3 Contin.       0.0643
#> X2_2       3  Binary       0.0096
#> A_2        3  Binary       0.0054
#> 
#> Effective sample sizes
#>  - Time 1
#>            Control Treated
#> Unadjusted 3306.    4194. 
#> Adjusted    845.79   899.4
#>  - Time 2
#>            Control Treated
#> Unadjusted 3701.   3799.  
#> Adjusted    912.87  829.87
#>  - Time 3
#>            Control Treated
#> Unadjusted 4886.   2614.  
#> Adjusted   1900.26  600.12

# Kernel balancing using the `kbal` package, available
# using `pak::pak_install("chadhazlett/KBAL")`.
# Only the ATT and ATC are available.

if (FALSE) { # \dontrun{
  kbal.fun <- function(treat, covs, estimand, focal, verbose, ...) {
    args <- list(...)

    if (!estimand %in% c("ATT", "ATC")) {
      stop('`estimand` must be "ATT" or "ATC".', call. = FALSE)
    }

    treat <- as.numeric(treat == focal)

    args <- args[names(args) %in% names(formals(kbal::kbal))]
    args$allx <- covs
    args$treatment <- treat
    args$printprogress <- verbose

    cat_cols <- apply(covs, 2L, function(x) length(unique(x)) <= 2)

    if (all(cat_cols)) {
      args$cat_data <- TRUE
      args$mixed_data <- FALSE
      args$scale_data <- FALSE
      args$linkernel <- FALSE
      args$drop_MC <- FALSE
    }
    else if (any(cat_cols)) {
      args$cat_data <- FALSE
      args$mixed_data <- TRUE
      args$cat_columns <- colnames(covs)[cat_cols]
      args$allx[,!cat_cols] <- scale(args$allx[,!cat_cols])
      args$cont_scale <- 1
    }
    else {
      args$cat_data <- FALSE
      args$mixed_data <- FALSE
    }

    k.out <- do.call(kbal::kbal, args)
    w <- k.out$w

    list(w = w, fit.obj = k.out)
  }

  (Wk <- weightit(treat ~ age + educ + married +
                    nodegree + re74, data = lalonde,
                  method = kbal.fun, estimand = "ATT",
                  include.obj = TRUE))
  summary(Wk)
  cobalt::bal.tab(Wk, stats = c("m", "ks"))
} # }
```
