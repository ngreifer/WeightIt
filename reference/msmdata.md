# Simulated data for a 3 time point sequential study

This is a simulated dataset of 7500 units with covariates and treatment
measured three times and the outcome measured at the end from a
hypothetical observational study examining the effect of treatment
delivered at each time point on an adverse event.

The data were generated using a simple simulation mechanism. For further
details on how the dataset was built, see the code at
[data-raw/msmdata.R](https://github.com/ngreifer/Weightit/blob/master/data-raw/msmdata.R).

The dataset is provided to illustrate the features of
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
and is not based on a realistic data-generating process, so it should
not be used as a benchmark.

For simulating realistic data with a known data-generating mechanism,
consider using the simcausal package.

## Usage

``` r
msmdata
```

## Format

A data frame with 7500 observations on the following 10 variables.

- `X1_0`:

  a count covariate measured at baseline

- `X2_0`:

  a binary covariate measured at baseline

- `A_1`:

  a binary indicator of treatment status at the first time point

- `X1_1`:

  a count covariate measured at the first time point (after the first
  treatment)

- `X2_1`:

  a binary covariate measured at the first time point (after the first
  treatment)

- `A_2`:

  a binary indicator of treatment status at the second time point

- `X1_2`:

  a count covariate measured at the second time point (after the second
  treatment)

- `X2_2`:

  a binary covariate measured at the first time point (after the first
  treatment)

- `A_3`:

  a binary indicator of treatment status at the third time point

- `Y_B`:

  a binary indicator of the outcome event (e.g., death)

## Examples

``` r
data("msmdata")
```
