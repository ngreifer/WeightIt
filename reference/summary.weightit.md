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
summary(
  object,
  top = 5L,
  ignore.s.weights = FALSE,
  weight.range = TRUE,
  which.time,
  ...
)

# S3 method for class 'summary.weightitMSM'
plot(x, binwidth = NULL, bins = NULL, which.time = 1L, ..., time)
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

- which.time:

  for [`summary()`](https://rdrr.io/r/base/summary.html), which models
  to summarize, given as a vector of positions in `formula.list` or of
  treatment or censoring variable names; omit it to summarize all of
  them. For [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  which single model to display the distribution of weights for, given
  as one position or one name; default is the first. This mirrors the
  argument of the same name in
  [`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
  , except that a value matching no model is an error rather than a
  warning, since here it decides what is computed. Note that when
  censoring is modeled, the censoring models occupy positions of their
  own, so positions do not count treatment time points alone; naming the
  variable avoids having to count.

- time:

  for [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the
  former name of `which.time`, which still works.

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

  The negative entropy (\\\frac{1}{n}\sum w \log(w)\\) of the weights in
  each treatment group and overall, after dividing the weights by their
  mean in the corresponding group.

- weight.mean (Mean of Weights):

  The mean of the weights in each treatment group and overall. Only
  included when the weights are stabilized.

- num.zeros:

  The number of weights equal to zero.

- effective.sample.size:

  The effective sample size for each treatment group before and after
  weighting. See
  [`ESS()`](https://ngreifer.github.io/WeightIt/reference/ESS.md).

For longitudinal treatments (i.e., `weightitMSM` objects),
[`summary()`](https://rdrr.io/r/base/summary.html) returns a list of the
above elements for each treatment period. When censoring is modeled (see
[`.cens()`](https://ngreifer.github.io/WeightIt/reference/dot-cens.md)),
each censoring model gets an entry of its own, placed among the
treatment entries in the order the models were fit and named for its
censoring indicator. Every entry summarizes the same weights – the
product across all the models – and differs only in the sample it
summarizes them over: a treatment entry splits by treatment group, and a
censoring entry covers the units still under observation when that model
was fit, since the censored units have a weight of exactly 0 and are not
part of the weighted sample.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
`ggplot` object with a histogram displaying the distribution of the
estimated weights. For a censoring model, only the weights of the units
still under observation are displayed, matching what
[`summary()`](https://rdrr.io/r/base/summary.html) reports. If the
estimand is the ATT or ATC, only the weights for the non-focal group(s)
will be displayed (since the weights for the focal group are all 1). A
dotted line is displayed at the mean of the weights.

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
