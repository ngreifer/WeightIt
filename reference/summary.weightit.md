# Print and Summarize Output

[`summary()`](https://rdrr.io/r/base/summary.html) generates a summary
of the `weightit` or `weightitMSM` object to evaluate the properties of
the estimated weights.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) plots the
distribution of the weights.
[`nobs()`](https://rdrr.io/r/stats/nobs.html) extracts the number of
observations.

## Usage

``` r
# S3 method for class 'weightit'
summary(object, top = 5L, ignore.s.weights = FALSE, weight.range = TRUE, ...)

# S3 method for class 'summary.weightit'
plot(x, binwidth = NULL, bins = NULL, ...)

# S3 method for class 'weightitMSM'
summary(object, top = 5L, ignore.s.weights = FALSE, weight.range = TRUE, ...)

# S3 method for class 'summary.weightitMSM'
plot(x, binwidth = NULL, bins = NULL, time = 1L, ...)
```

## Arguments

- object:

  a `weightit` or `weightitMSM` object; the output of a call to
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  or
  [`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).

- top:

  how many of the largest and smallest weights to display. Default is 5.
  Ignored when `weight.range = FALSE`.

- ignore.s.weights:

  `logical`; whether or not to ignore sampling weights when computing
  the weight summary. If `FALSE`, the default, the estimated weights
  will be multiplied by the sampling weights (if any) before values are
  computed.

- weight.range:

  `logical`; whether to display statistics about the range of weights
  and the highest and lowest weights for each group. Default is `TRUE`.

- ...:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  additional arguments passed to
  [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html) to
  determine the number of bins, though
  [`ggplot2::geom_histogram()`](https://ggplot2.tidyverse.org/reference/geom_histogram.html)
  is actually used to create the plot.

- x:

  a `summary.weightit` or `summary.weightitMSM` object; the output of a
  call to `summary.weightit()` or `summary.weightitMSM()`.

- binwidth, bins:

  arguments passed to
  [`ggplot2::geom_histogram()`](https://ggplot2.tidyverse.org/reference/geom_histogram.html)
  to control the size and/or number of bins.

- time:

  `numeric`; the time point for which to display the distribution of
  weights. Default is to plot the distribution for the first time
  points.

## Value

For point treatments (i.e., `weightit` objects),
[`summary()`](https://rdrr.io/r/base/summary.html) returns a
`summary.weightit` object with the following elements:

- weight.range:

  The range (minimum and maximum) weight for each treatment group.

- weight.top:

  The units with the greatest weights in each treatment group; how many
  are included is determined by `top`.

- coef.of.var (Coef of Var):

  The coefficient of variation (standard deviation divided by mean) of
  the weights in each treatment group and overall.

- scaled.mad (MAD):

  The mean absolute deviation of the weights in each treatment group and
  overall divided by the mean of the weights in the corresponding group.

- negative entropy (Entropy):

  The negative entropy (\\\sum w log(w)\\) of the weights in each
  treatment group and overall divided by the mean of the weights in the
  corresponding group.

- num.zeros:

  The number of weights equal to zero.

- effective.sample.size:

  The effective sample size for each treatment group before and after
  weighting. See
  [`ESS()`](https://ngreifer.github.io/WeightIt/reference/ESS.md).

For longitudinal treatments (i.e., `weightitMSM` objects),
[`summary()`](https://rdrr.io/r/base/summary.html) returns a list of the
above elements for each treatment period.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
`ggplot` object with a histogram displaying the distribution of the
estimated weights. If the estimand is the ATT or ATC, only the weights
for the non-focal group(s) will be displayed (since the weights for the
focal group are all 1). A dotted line is displayed at the mean of the
weights.

[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns a single number.
Note that even units with `weights` or `s.weights` of 0 are included.

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md),
[`summary()`](https://rdrr.io/r/base/summary.html)

## Examples

``` r
# See example at ?weightit or ?weightitMSM
```
