# Notes on `vignette("longitudinal-treatments")` for review

Written 2026-10-08 alongside the first draft of `vignettes/longitudinal-treatments.Rmd`. Each item is a decision I made without asking, or something I could not settle from the package and want you to look at.

## Decisions made

1. **Structure.** Introduction; a conceptual section (estimand, time-varying confounding, assumptions, the weights and stabilization, balance); a full analysis of `msmdata` (initial imbalance, unstabilized then stabilized weights, balance at each time point with a love plot, a CBPS alternative, effect estimation with a saturated MSM, and per-time-point effects recovered from it with `avg_comparisons()`); a shorter censoring analysis (estimand, simulated dropout, weights, balance, effect); and an "Other Considerations" section in place of a reporting section. The dataset is `msmdata`, as in the other vignettes, so the three vignettes agree on the example.

2. **Stabilized weights are the main analysis.** The unstabilized weights are shown first to motivate stabilization (maximum weight about 403, effective sample sizes a quarter of the group sizes). The effect is estimated with the stabilized logistic regression weights because they support M-estimation.

3. **Prior-treatment imbalance with stabilized weights is explained rather than hidden.** With `stabilize = TRUE`, `bal.tab()` shows standardized mean differences of about .13 and .15 for `A_1` at time 2 and `A_2` at time 3, roughly their unadjusted values. The vignette explains that the saturated numerator preserves the association among the treatments, that this requires the MSM to include the treatment history, and that balance within levels of the treatment history is what the weights guarantee. This is the point most likely to confuse readers who compare the two balance tables, so it may deserve more space, or a demonstration of balance within treatment-history strata (cobalt's `cluster` argument on the treatment history could do it).

4. **The CBPS alternative is shown but not used for the effect.** `method = "cbps"` with `is.MSM.method = TRUE` gives exact mean balance at every time point, prior treatments included, with lower effective sample sizes, and no M-estimation. I present it as the natural next step when balance is inadequate and then proceed with the GLM weights. Adding covariate interactions to the GLM treatment models, the other obvious refinement, made the time-3 balance on `X1_2` worse (.14 against .10) and lowered the ESS, so it is mentioned as an option but not shown.

5. **No parsimonious MSM is fit, and no effects of individual treatments are reported** (both removed at your request). The estimands are the comparisons between regimes. An earlier version recovered the average effect of treatment at each time point with `avg_comparisons()`, but that averages over the observed distribution of the other treatments, which is a reference distribution with no causal meaning for an MSM; the equal-weight factorial main effect would be the defensible alternative if this is ever added back. The conceptual section still says that a parsimonious MSM is a necessity only with many time points.

6. **Censoring is simulated after the time-3 covariates and before `A_3`,** depending on `X1_2` and `A_2` (`plogis(-4 + .35 * X1_2 + .8 * A_2)`, about 18% censored), with `A_3` and `Y_B` set to missing. This matches the timing in cobalt's longitudinal vignette. An earlier draft placed the censoring before the time-3 covariates (as in `vignette("WeightIt")`), which set `X1_2` and `X2_2` to missing and produced `X1_2:<NA>` rows in the time-4 balance table; see item 11.

7. **No complete-case comparison** (removed at your request). The prose reports the censoring-weighted estimate next to the full-data estimate from the first analysis and explains in words why dropping the censored units would bias the estimate (the lower-risk remaining population and the induced association between `A_2` and `X1_2`).

8. **No Rd markup, no dashes, headings as noun phrases**, and the style greps from the r-doc-style skill are clean. Variable names (`A_1`, `X1_2`, `msmdata`) appear as sentence subjects; the argument-subject rule is about arguments, so I left those.

9. **Added to `_pkgdown.yml`** after estimating-effects and before installing-packages, and a NEWS item was added below the weighting-methods one. Reorder if you prefer it elsewhere.

10. **Cross-references (2026-10-08).** Both new vignettes are now pointed to from the README (the vignette list, the methods table, and the censoring paragraph), from `weightit()` (Details and a MatchIt-style `@seealso` list of all five vignettes), `weightitMSM()`, `.cens()`, `msmdata`, `.weightit_methods`, `get_w_from_ps()`, `ESS()`, `trim()`, `calibrate()`, and every `method_*` page (one `@seealso` line each), and from `vignette("WeightIt")`, `vignette("estimating-effects")`, `vignette("installing-packages")`, and `vignette("weighting-methods")`. README.md was regenerated from README.Rmd; the diff is the added prose only.

## Things to look at

10. **The censoring example is seed-sensitive.** I ran five seeds under two censoring designs (about 18% and about 34% censored) and two model specifications. With the stronger design, the final weights are very variable, the censoring-weighted estimate of the always-versus-never contrast ranges from about -.20 to -.34 across seeds (full-data estimate -.265), and the complete-case estimate is further from the truth in most but not all seeds. With the milder design used in the vignette, time-3 balance is good (maximum standardized difference .08 to .16) and the censoring-weighted estimate is close to the full-data estimate; the complete-case estimate, no longer shown, was consistently about .02 to .04 further from it. Seed 1234 under either design produces one enormous weight that destroys the time-3 balance (maximum standardized difference .32 to .59) and makes both estimates equally biased; I chose seed 7, which is representative of the other seeds. If the data-generating process for `msmdata` ever changes, the prose around the censoring numbers should be re-read.

11. **Product-of-weights balance after censoring.** In the pathological seed, the time-3 treatment weights alone balanced the covariates among the units at risk (standardized differences below .04), but the product with the censoring weights did not (`X1_2` at .32). The censoring weights upweight high-`X1_2` units, where the main-effects logistic model for `A_3` is misspecified (the true model has `X1_2:X2_2` and `X1_1:X2_1` interactions), and the misspecification is amplified. The vignette now ends the censoring section with a sentence saying that balance after a censoring time point deserves particular care. Two things might be worth considering in the package: whether the treatment models after a censoring point should be fit with the accumulated censoring weights as sampling weights (so they target balance in the IPCW pseudo-population under misspecification), and whether `summary.weightitMSM()` or `bal.tab()` should flag a single unit carrying a large share of the weight.

12. **`bal.tab()` with missing post-censoring covariates.** When covariates measured after the censoring point are set to `NA` for censored units, `bal.tab()` on the `weightitMSM` object prints `X1_2:<NA>` and `X2_2:<NA>` rows (all zeros) in the later treatment's table and warns about missing values, even though those units are not in the at-risk set for that model. The vignette avoids this by setting only `A_3` and `Y_B` to missing. Filed as a GitHub issue on cobalt (2026-10-08) with a reproducible example against cobalt 5.0.0 and WeightIt 2.1.0.

13. **The marginal balance check with stabilized weights** (item 3) might merit a sentence in cobalt's longitudinal vignette as well, since `bal.tab()` there is run on stabilized weights without comment on the prior-treatment rows.

14. **`love.plot()` call.** I used `love.plot(W, stats = "m", binary = "std", abs = TRUE, thresholds = .1, which.time = .all)`. I extracted the rendered figure and checked it: three panels, one per time point, with the adjusted points at or below the .1 line except for the prior treatments, as the prose says.

## Bibliography

Ten entries were added: eight from Zotero through Better BibTeX (`robinsNewApproachCausal1986`, `robinsMarginalStructuralModels2000a`, `hernanEstimatingCausalEffects2006`, `imaiRobustEstimationInverse2015`, `jacksonDiagnosticsConfoundingTimevarying2016`, `lefebvreImpactMisspecificationTreatment2008`, `vanderweeleCausalInferenceLongitudinal2016`, `westreichInvitedCommentaryPositivity2010`) and two copied from cobalt's `references.bib` (`hernanMarginalStructuralModels2000`, `robinsAnalysisSemiparametricRegression1995`). Entries already in the bib that are cited: `robinsMarginalStructuralModels2000`, `coleConstructingInverseProbability2008`, `thoemmesPrimerInverseProbability2016`, `hernanCausalInferenceWhat2020`, `imaiCovariateBalancingPropensity2014`.

Added after your review (2026-10-08): `stallworthyInvestigatingCausalQuestions2026` from Zotero, cited in the introduction, in the estimand section on timing and dosage questions, and in Other Considerations for the *devMSMs* package; `luClearerDefinitionSelection2022` from the updated Zotero record; and eight entries written from Crossref metadata (`danielMethodsDealingTimedependent2013`, `mansourniaHandlingTimeVarying2017`, `robinsCorrectingNoncomplianceDependent2000`, `hernanStructuralApproachSelection2004`, `howeLimitationInverseProbabilityofcensoring2011`, `plattPositivityAssumptionMarginal2012`, `seamanReviewInverseProbability2013`, `williamsonMarginalStructuralModels2017`). For Daniel et al., Platt et al., and Seaman and White, Crossref's BibTeX gives the online-first year (2012, 2011, 2011); the entries use the print year and volume, which is how the papers are usually cited. Fewell et al. (2004) and Blackwell (2013) were skipped at your request. None of the Crossref entries is in Zotero, so they will not be found by a citation-key search there.

## Claims to double-check

- The saturated MSM equation and the definitions of sequential exchangeability and positivity follow Robins, Hernán, and Brumback (2000) and Hernán and Robins (2020); I wrote them from those sources' standard statements, not from a specific page.
- "Stabilized weights have a mean of 1" (with a saturated numerator) and "balance within levels of the treatment history" are standard; the second is my phrasing of Cole and Hernán's (2008) point that the MSM must include the numerator's variables.
- The statement that `glm_weightit()`'s M-estimation standard errors include the estimation of the stabilization factors comes from `vignette("estimating-effects")`.
- The description of the CBPS MSM option ("estimates all the treatment models at once so that the product of the weights exactly balances the covariate means at every time point") comes from `?method_cbps` and the balance output, which shows zeros at every time point.
- The collider explanation for the complete-case bias (dropout depends on `A_2` and `X1_2`; restricting to the uncensored induces an association between them that the treatment weights estimated among them do not remove) is my reading of the data-generating process, not a cited result.
- "Misspecification of the treatment models biases the weights and the estimate" is attributed to Lefebvre, Delaney, and Platt (2008), whose abstract I read.

## Things left out on purpose

Continuous and multi-category treatments at some time points (mentioned only as possible), sampling weights, multiple imputation, `by`, random effects, more than one censoring time point, survival outcomes, and the `num.formula` mechanics beyond a pointer. Bootstrap standard errors are mentioned but not run, to keep the build time short. A complete-case comparison and a parsimonious (additive) MSM were in the first draft and removed at your request.

## Checks run

- The vignette renders with `rmarkdown::render()` under `devtools::load_all()`, with every citation resolving and no warnings from the code.
- No Rd markup in the markdown; headings, counts, clefts, dashes, and excluded-word greps are clean.
- `pkgdown::check_pkgdown()` passes with the new article entry.
