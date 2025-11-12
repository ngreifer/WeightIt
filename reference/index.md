# Package index

## Estimate weights

- [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  : Estimate Balancing Weights
- [`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
  : Generate Balancing Weights for Longitudinal Treatments
- [`summary(`*`<weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/summary.weightit.md)
  [`plot(`*`<summary.weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/summary.weightit.md)
  [`summary(`*`<weightitMSM>`*`)`](https://ngreifer.github.io/WeightIt/reference/summary.weightit.md)
  [`plot(`*`<summary.weightitMSM>`*`)`](https://ngreifer.github.io/WeightIt/reference/summary.weightit.md)
  : Print and Summarize Output

## Specific estimation methods

- [`method_bart`](https://ngreifer.github.io/WeightIt/reference/method_bart.md)
  : Propensity Score Weighting Using BART
- [`method_cbps`](https://ngreifer.github.io/WeightIt/reference/method_cbps.md)
  : Covariate Balancing Propensity Score Weighting
- [`method_ebal`](https://ngreifer.github.io/WeightIt/reference/method_ebal.md)
  [`method_entropy`](https://ngreifer.github.io/WeightIt/reference/method_ebal.md)
  : Entropy Balancing
- [`method_energy`](https://ngreifer.github.io/WeightIt/reference/method_energy.md)
  : Energy Balancing
- [`method_gbm`](https://ngreifer.github.io/WeightIt/reference/method_gbm.md)
  : Propensity Score Weighting Using Generalized Boosted Models
- [`method_glm`](https://ngreifer.github.io/WeightIt/reference/method_glm.md)
  : Propensity Score Weighting Using Generalized Linear Models
- [`method_ipt`](https://ngreifer.github.io/WeightIt/reference/method_ipt.md)
  : Inverse Probability Tilting
- [`method_npcbps`](https://ngreifer.github.io/WeightIt/reference/method_npcbps.md)
  : Nonparametric Covariate Balancing Propensity Score Weighting
- [`method_optweight`](https://ngreifer.github.io/WeightIt/reference/method_optweight.md)
  [`method_sbw`](https://ngreifer.github.io/WeightIt/reference/method_optweight.md)
  : Stable Balancing Weights
- [`method_super`](https://ngreifer.github.io/WeightIt/reference/method_super.md)
  : Propensity Score Weighting Using SuperLearner
- [`method_user`](https://ngreifer.github.io/WeightIt/reference/method_user.md)
  : User-Defined Functions for Estimating Weights

## Fit weighted regression models

- [`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
  [`lm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
  [`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
  [`multinom_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
  [`coxph_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
  : Fitting Weighted Generalized Linear Models

- [`predict(`*`<glm_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/predict.glm_weightit.md)
  [`predict(`*`<ordinal_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/predict.glm_weightit.md)
  [`predict(`*`<multinom_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/predict.glm_weightit.md)
  :

  Predictions for `glm_weightit` objects

- [`anova(`*`<glm_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/anova.glm_weightit.md)
  :

  Methods for
  [`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
  objects

- [`summary(`*`<glm_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/glm_weightit-methods.md)
  [`summary(`*`<multinom_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/glm_weightit-methods.md)
  [`summary(`*`<ordinal_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/glm_weightit-methods.md)
  [`summary(`*`<coxph_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/glm_weightit-methods.md)
  [`print(`*`<glm_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/glm_weightit-methods.md)
  [`vcov(`*`<glm_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/glm_weightit-methods.md)
  [`estfun(`*`<glm_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/glm_weightit-methods.md)
  [`update(`*`<glm_weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/glm_weightit-methods.md)
  :

  Methods for
  [`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
  objects

## Modify weights

- [`trim()`](https://ngreifer.github.io/WeightIt/reference/trim.md) :
  Trim (Winsorize) Large Weights
- [`calibrate()`](https://ngreifer.github.io/WeightIt/reference/calibrate.md)
  : Calibrate Propensity Score Weights
- [`sbps()`](https://ngreifer.github.io/WeightIt/reference/sbps.md) :
  Subgroup Balancing Propensity Score

## Other functions

- [`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
  : Compute weights from propensity scores

- [`ESS()`](https://ngreifer.github.io/WeightIt/reference/ESS.md) :
  Compute effective sample size of weighted sample

- [`weightit.fit()`](https://ngreifer.github.io/WeightIt/reference/weightit.fit.md)
  : Generate Balancing Weights with Minimal Input Processing

- [`as.weightit()`](https://ngreifer.github.io/WeightIt/reference/as.weightit.md)
  [`as.weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/as.weightit.md)
  :

  Create a `weightit` object manually

- [`make_full_rank()`](https://ngreifer.github.io/WeightIt/reference/make_full_rank.md)
  : Make a design matrix full rank

- [`plot(`*`<weightit>`*`)`](https://ngreifer.github.io/WeightIt/reference/plot.weightit.md)
  : Plot information about the weight estimation process

- [`.weightit_methods`](https://ngreifer.github.io/WeightIt/reference/dot-weightit_methods.md)
  : Weighting methods

## Datasets

- [`msmdata`](https://ngreifer.github.io/WeightIt/reference/msmdata.md)
  : Simulated data for a 3 time point sequential study
