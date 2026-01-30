# Plot information about the weight estimation process

`plot.weightit()` plots information about the weights depending on how
they were estimated. Currently, only weighting using `method = "gbm"` or
`"optweight"` are supported. To plot the distribution of weights, see
[`plot.summary.weightit()`](https://ngreifer.github.io/WeightIt/reference/summary.weightit.md).

## Usage

``` r
# S3 method for class 'weightit'
plot(x, ...)
```

## Arguments

- x:

  a `weightit` object; the output of a call to
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md).

- ...:

  unused.

## Value

A `ggplot` object.

## Details

### `method = "gbm"`

After weighting with generalized boosted modeling,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) displays the
results of the tuning process used to find the optimal number of trees
(and tuning parameter values, if modified) that are used in the final
weights. The plot produced has the number of trees on the x-axis and the
value of the criterion on the y-axis with a diamond at the optimal
point. When multiple parameters are selected by tuning, a separate line
is displayed on the plot for each combination of tuning parameters. When
`by` is used in the call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
the plot is faceted by the `by` variable. See
[`method_gbm`](https://ngreifer.github.io/WeightIt/reference/method_gbm.md)
for more information on selecting tuning parameters.

### `method = "optweight"`

After estimating stable balancing weights,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) displays the
values of the dual variables for each balance constraint in a bar graph.
Large values of the dual variables indicate the covariates for which the
balance constraint is causing increases in the variability of the
weights, i.e., the covariates for which relaxing the imbalance tolerance
would yield the greatest gains in effective sample size. For continuous
treatments, the dual variables are split into those for the target
(i.e., ensuring the mean of each covariate after weighting is equal to
its unweighted mean) and those for balance (i.e., ensuring the
treatment-covariate correlations are no larger than the imbalance
tolerance). This is essentially a wrapper for
[`optweight::plot.optweight()`](https://ngreifer.github.io/optweight/reference/plot.optweight.html)
. See
[`method_optweight`](https://ngreifer.github.io/WeightIt/reference/method_optweight.md)
for details.

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`plot.summary.weightit()`](https://ngreifer.github.io/WeightIt/reference/summary.weightit.md)

## Examples

``` r
# See example at the corresponding methods page
```
