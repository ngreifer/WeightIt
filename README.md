
<!-- README.md is generated from README.Rmd. Please edit that file -->
WeightIt
========

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/WeightIt)](https://cran.r-project.org/package=WeightIt)

`WeightIt` is a one-stop package to generate balancing weights for point treatments in observational studies. Contained within `WeightIt` are methods that call on other R packages to estimate weights. The value of `WeightIt` is in its unified and familiar syntax used to generate the weights, as each of these other packages have their own, often challenging to navigate, syntax. `WeightIt` extends the capabilities of these packages to generate weights used to estimate the ATE, ATT, and ATC for binary or multinomial treatments, and treatment effects for continuous treatments when available. In these ways, `WeightIt` does for weighting what `MatchIt` has done for matching, and `MatchIt` users will find the syntax familiar.

To install and load `WeightIt`, use the code below:

``` r
install.packages("WeightIt")
library("WeightIt")
```

The workhorse function of `WeightIt` is `weightit()`, which generates weights from a given formula and data input according to methods and other parameters sepcified by the user. Below is an example of the use of `weightit()` to generate weights for estimating the ATE:

``` r
data("lalonde", package = "cobalt")
W <- weightit(treat ~ age + educ + nodegree + married + race + re74 + re75, 
    data = lalonde, method = "ps", estimand = "ATE")
print(W)
```

    A weightit object
     - method: "ps" (propensity score weighting)
     - number of obs.: 614
     - sampling weights: none
     - treatment: 2-category
     - estimand: ATE
     - covariates: age, educ, race, married, nodegree, re74, re75

Evaluating weights has two components: evaluating the covariate balance produces by the weights, and evaluating whether the weights will allow for sufficient precision in the eventual effect estimate. For the first goal, functions in the `cobalt` package, which are fully compatible with `WeightIt`, can be used, as demonstrated below:

``` r
library("cobalt")
bal.tab(W, un = TRUE)
```

    Balance Measures:
                    Type Diff.Un Diff.Adj
    prop.score  Distance  1.7569   0.1360
    age          Contin. -0.2419  -0.1676
    educ         Contin.  0.0448   0.1296
    race_black    Binary  0.6404   0.0499
    race_hispan   Binary -0.0827   0.0047
    race_white    Binary -0.5577  -0.0546
    married       Binary -0.3236  -0.0944
    nodegree      Binary  0.1114  -0.0547
    re74         Contin. -0.5958  -0.2740
    re75         Contin. -0.2870  -0.1579

    Effective sample sizes:
               Control Treated
    Unadjusted     429  185.00
    Adjusted       329   58.33

For the second goal, qualities of the distributions of weights can be assessed using `summary()`, as demonstrated below.

``` r
summary(W)
```

    Summary of weights:

    - Weight ranges:
               Min X.............................     Max
    treated 1.1721  |---------------------------| 40.0773
    control 1.0092  |-|                            4.7432

    - Units with 5 greatest weights by group:
                                                    
                 137     124     116      68      10
     treated 13.5451 15.9884 23.2967 23.3891 40.0773
                 412     388     226     196     118
     control  4.0301  4.0592  4.2397  4.5231  4.7432

              Ratio Coef.of.Var
    treated 34.1921      1.4777
    control  4.7002      0.5519
    overall 39.7134      1.3709

    - Effective Sample Sizes:
               Control Treated
    Unweighted 429.000 185.000
    Weighted   329.008  58.327

Desirable qualities include ratios close to 1, coefficients of variation close to 0, and large effective sample sizes.

The table below contains the available methods in `WeightIt` for estimating weights for binary, multinomial, and continuous treatments using various methods and functions from various packages.

<table style="width:74%;">
<colgroup>
<col width="20%" />
<col width="29%" />
<col width="12%" />
<col width="11%" />
</colgroup>
<thead>
<tr class="header">
<th>Treatment type</th>
<th>Method (<code>method =</code>)</th>
<th>Function</th>
<th>Package</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><strong>Binary</strong></td>
<td>Binary regression PS (<code>&quot;ps&quot;</code>)</td>
<td><code>glm()</code></td>
<td><code>base</code></td>
</tr>
<tr class="even">
<td>-</td>
<td>Generalized boosted modeling PS (<code>&quot;pgbm&quot;</code>)</td>
<td><code>ps()</code></td>
<td><code>twang</code></td>
</tr>
<tr class="odd">
<td>-</td>
<td>Covariate Balancing PS (<code>&quot;cbps&quot;</code>)</td>
<td><code>CBPS()</code></td>
<td><code>CBPS</code></td>
</tr>
<tr class="even">
<td>-</td>
<td>Non-Parametric Covariate Balancing PS (<code>&quot;npcbps&quot;</code>)</td>
<td><code>npCBPS()</code></td>
<td><code>CBPS</code></td>
</tr>
<tr class="odd">
<td>-</td>
<td>Entropy Balancing (<code>&quot;ebal&quot;</code>)</td>
<td><code>ebalance()</code></td>
<td><code>ebal</code></td>
</tr>
<tr class="even">
<td>-</td>
<td>Stabilized Balancing Weights (<code>&quot;sbw&quot;</code>)</td>
<td><code>sbw()</code></td>
<td><code>sbw</code></td>
</tr>
<tr class="odd">
<td>-</td>
<td>Empirical Balancing Calibration Weights (<code>&quot;ebcw&quot;</code>)</td>
<td><code>ATE()</code></td>
<td><code>ATE</code></td>
</tr>
<tr class="even">
<td><strong>Multinomial</strong></td>
<td>Multiple binary regression PS (<code>&quot;ps&quot;</code>)</td>
<td><code>glm()</code></td>
<td><code>base</code></td>
</tr>
<tr class="odd">
<td>-</td>
<td>Multinomial regression PS (<code>&quot;ps&quot;</code>)</td>
<td><code>mlogit()</code></td>
<td><code>mlogit</code></td>
</tr>
<tr class="even">
<td>-</td>
<td>Bayesian multinomial regression PS (<code>&quot;ps&quot;, link = &quot;bayes.probit&quot;</code>)</td>
<td><code>MNP()</code></td>
<td><code>MNP</code></td>
</tr>
<tr class="odd">
<td>-</td>
<td>Generalized boosted modeling PS (<code>&quot;gbm&quot;</code>)</td>
<td><code>ps()</code></td>
<td><code>twang</code></td>
</tr>
<tr class="even">
<td>-</td>
<td>Covariate Balancing PS (<code>&quot;cbps&quot;</code>)</td>
<td><code>CBPS()</code></td>
<td><code>CBPS</code></td>
</tr>
<tr class="odd">
<td>-</td>
<td>Non-Parametric Covariate Balancing PS (<code>&quot;npcbps&quot;</code>)</td>
<td><code>npCBPS()</code></td>
<td><code>CBPS</code></td>
</tr>
<tr class="even">
<td>-</td>
<td>Entropy Balancing (<code>&quot;ebal&quot;</code>)</td>
<td><code>ebalance()</code></td>
<td><code>ebal</code></td>
</tr>
<tr class="odd">
<td>-</td>
<td>Stabilized Balancing Weights (<code>&quot;sbw&quot;</code>)</td>
<td><code>sbw()</code></td>
<td><code>sbw</code></td>
</tr>
<tr class="even">
<td>-</td>
<td>Empirical Balancing Calibration Weights (<code>&quot;ebcw&quot;</code>)</td>
<td><code>ATE()</code></td>
<td><code>ATE</code></td>
</tr>
<tr class="odd">
<td><strong>Continuous</strong></td>
<td>Generalized linear model PS (<code>&quot;ps&quot;</code>)</td>
<td><code>glm()</code></td>
<td><code>base</code></td>
</tr>
<tr class="even">
<td>-</td>
<td>Covariate Balancing PS (<code>&quot;cbps&quot;</code>)</td>
<td><code>CBPS()</code></td>
<td><code>CBPS</code></td>
</tr>
<tr class="odd">
<td>-</td>
<td>Non-Parametric Covariate Balancing PS (<code>&quot;npcbps&quot;</code>)</td>
<td><code>npCBPS()</code></td>
<td><code>CBPS</code></td>
</tr>
</tbody>
</table>

If you would like to see your package integrated into `WeightIt`, or for any other questions or comments about `WeightIt`, please contact Noah Greifer at <noah@unc.edu>.
