# Make a design matrix full rank

When writing [user-defined
methods](https://ngreifer.github.io/WeightIt/reference/method_user.md)
for use with
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
it may be necessary to take the potentially non-full rank `covs` data
frame and make it full rank for use in a downstream function. This
function performs that operation.

## Usage

``` r
make_full_rank(mat, with.intercept = TRUE)
```

## Arguments

- mat:

  a numeric matrix or data frame to be transformed. Typically this
  contains covariates. `NA`s are not allowed.

- with.intercept:

  `logical`; whether an intercept (i.e., a vector of 1s) should be added
  to `mat` before making it full rank. If `TRUE`, the intercept will be
  used in determining whether a column is linearly dependent on others.
  Regardless, no intercept will be included in the output.

## Value

An object of the same type as `mat` containing only linearly independent
columns.

## Details

`make_full_rank()` calls [`qr()`](https://rdrr.io/r/base/qr.html) to
find the rank and linearly independent columns of `mat`, which are
retained while others are dropped. If `with.intercept` is set to `TRUE`,
an intercept column is added to the matrix before calling
[`qr()`](https://rdrr.io/r/base/qr.html). Note that dependent columns
that appear later in `mat` will be dropped first.

See example at
[`method_user`](https://ngreifer.github.io/WeightIt/reference/method_user.md).

## Note

Older versions would drop all columns that only had one value. With
`with.intercept = FALSE`, if only one column has only one value, it will
not be removed, and it will function as though there was an intercept
present; if more than only column has only one value, only the first one
will remain.

## See also

[`method_user`](https://ngreifer.github.io/WeightIt/reference/method_user.md),
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)

## Examples

``` r
set.seed(1234)
n <- 20
c1 <- rbinom(n, 1, .4)
c2 <- 1 - c1
c3 <- rnorm(n)
c4 <- 10 * c3
mat <- data.frame(c1, c2, c3, c4)

make_full_rank(mat) #leaves c2 and c4
#>    c1          c3
#> 1   0 -0.47719270
#> 2   1 -0.99838644
#> 3   1 -0.77625389
#> 4   1  0.06445882
#> 5   1  0.95949406
#> 6   1 -0.11028549
#> 7   0 -0.51100951
#> 8   0 -0.91119542
#> 9   1 -0.83717168
#> 10  0  2.41583518
#> 11  1  0.13408822
#> 12  0 -0.49068590
#> 13  0 -0.44054787
#> 14  1  0.45958944
#> 15  0 -0.69372025
#> 16  1 -1.44820491
#> 17  0  0.57475572
#> 18  0 -1.02365572
#> 19  0 -0.01513830
#> 20  0 -0.93594860

make_full_rank(mat, with.intercept = FALSE) #leaves c1, c2, and c4
#>    c1 c2          c3
#> 1   0  1 -0.47719270
#> 2   1  0 -0.99838644
#> 3   1  0 -0.77625389
#> 4   1  0  0.06445882
#> 5   1  0  0.95949406
#> 6   1  0 -0.11028549
#> 7   0  1 -0.51100951
#> 8   0  1 -0.91119542
#> 9   1  0 -0.83717168
#> 10  0  1  2.41583518
#> 11  1  0  0.13408822
#> 12  0  1 -0.49068590
#> 13  0  1 -0.44054787
#> 14  1  0  0.45958944
#> 15  0  1 -0.69372025
#> 16  1  0 -1.44820491
#> 17  0  1  0.57475572
#> 18  0  1 -1.02365572
#> 19  0  1 -0.01513830
#> 20  0  1 -0.93594860
```
