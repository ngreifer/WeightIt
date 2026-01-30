# Weighting methods

`.weightit_methods` is a list containing the allowable weighting methods
that can be supplied by name to the `method` argument of
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md),
and
[`weightit.fit()`](https://ngreifer.github.io/WeightIt/reference/weightit.fit.md).
Each entry corresponds to an allowed method and contains information
about what options are and are not allowed for each method. While this
list is primarily for internal use by checking functions in WeightIt, it
might be of use for package authors that want to support different
weighting methods.

## Usage

``` r
.weightit_methods
```

## Format

An object of class `list` of length 11.

## Details

Each component is itself a list containing the following components:

- `treat_type`:

  at least one of `"binary"`, `"multinomial"`, or `"continuous"`
  indicating which treatment types are available for this method.

- `estimand`:

  which estimands are available for this method. All methods that
  support binary and multi-category treatments accept `"ATE"`, `"ATT"`,
  and `"ATC"`, as well as some other estimands depending on the method.
  See
  [`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
  for more details about what each estimand means.

- `alias`:

  a character vector of aliases for the method. When an alias is
  supplied, the corresponding method will still be dispatched. For
  example, the canonical method to request entropy balancing is
  `"ebal"`, but `"ebalance"` and `"entropy"` also work. The first value
  is the canonical name.

- `description`:

  a string containing the description of the name in English.

- `ps`:

  a logical for whether propensity scores are returned by the method for
  binary treatments. Propensity scores are never returned for
  multi-category or continuous treatments.

- `msm_valid`:

  a logical for whether the method can be validly used with longitudinal
  treatments.

- `msm_method_available`:

  a logical for whether a version of the method can be used that
  estimates weights using a single model rather than multiplying the
  weights across time points. This is related to the `is.MSM.method`
  argument of
  [`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).

- `subclass_ok`:

  a logical for whether `subclass` can be supplied to compute
  subclassification weights from the propensity scores.

- `packages_needed`:

  a character vector of the minimal packages required to use the method.
  Some methods may require additional packages for certain options.

- `package_versions_needed`:

  a named character vector of the minimal package versions required to
  use the method. Only present when `packages_needed` is not empty.

- `s.weights_ok`:

  a logical for whether sampling weights can be used with the method.

- `missing`:

  a character vector of the allowed options that can be supplied to
  `missing` when missing data is present. All methods accept `"ind"` for
  the missingness indicator approach; some other methods accept
  additional values.

- `moments_int_ok`:

  a logical for whether `moments`, `int`, and `quantile` can be used
  with the method.

- `moments_default`:

  when `moments_int_ok` is `TRUE`, the default value of `moments` used
  with the method. For most methods, this is 1.

- `density_ok`:

  a logical for whether arguments that control the density can be used
  with the method when used with a continuous treatment.

- `stabilize_ok`:

  a logical for whether the `stabilize` argument (and `num.formula` for
  longitudinal treatments) can be used with the method.

- `plot.weightit_ok`:

  a logical for whether
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) can be used
  on the `weightit` output with the method.

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
and
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
for how the methods are used. Also see the individual methods pages for
information on whether and how each option can be used.

## Examples

``` r
# Get all acceptable names
names(.weightit_methods)
#>  [1] "glm"       "bart"      "cbps"      "ebal"      "energy"    "gbm"      
#>  [7] "ipt"       "npcbps"    "optweight" "super"     "cfd"      

# Get all acceptable names and aliases
lapply(.weightit_methods, `[[`, "alias")
#> $glm
#> [1] "glm" "ps" 
#> 
#> $bart
#> [1] "bart"
#> 
#> $cbps
#> [1] "cbps"  "cbgps"
#> 
#> $ebal
#> [1] "ebal"     "ebalance" "entropy" 
#> 
#> $energy
#> [1] "energy" "dcows" 
#> 
#> $gbm
#> [1] "gbm" "gbr"
#> 
#> $ipt
#> [1] "ipt"
#> 
#> $npcbps
#> [1] "npcbps"  "npcbgps"
#> 
#> $optweight
#> [1] "optweight" "sbw"      
#> 
#> $super
#> [1] "super"        "superlearner"
#> 
#> $cfd
#> [1] "cfd"    "kernel"
#> 

# Which estimands are allowed with `method = "bart"`
.weightit_methods[["bart"]]$estimand
#> [1] "ATE"  "ATT"  "ATC"  "ATO"  "ATM"  "ATOS"

# All methods that support continuous treatments
supp <- sapply(.weightit_methods, function(x) {
  "continuous" %in% x$treat_type
})
names(.weightit_methods)[supp]
#> [1] "glm"       "bart"      "cbps"      "ebal"      "energy"    "gbm"      
#> [7] "npcbps"    "optweight" "super"    

# All methods that return propensity scores (for
# binary treatments only)
supp <- sapply(.weightit_methods, `[[`, "ps")
names(.weightit_methods)[supp]
#> [1] "glm"   "bart"  "cbps"  "gbm"   "ipt"   "super"
```
