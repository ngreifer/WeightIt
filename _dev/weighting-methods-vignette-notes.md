# Notes on `vignette("weighting-methods")` for review

Written 2026-10-08 alongside the first draft of `vignettes/weighting-methods.Rmd`. Each item is a decision I made without asking, or something I could not settle from the package and want you to look at.

## Decisions made

1. **`"ps"` was read as `method = "glm"`.** The request listed `ps`, `cbps`, and `ipt` as the parametric methods. The canonical name for propensity score weighting with a GLM is `"glm"`, with `"ps"` kept as an alias, and `?method_ps` now documents supplying your own propensity scores through the `ps` argument. The vignette's parametric section covers `method = "glm"` and notes the old name; supplying propensity scores is covered as one of the Additional Options.

2. **Section structure follows the MatchIt vignette, plus an Additional Options section.** The request named the opening material, the three method families, and the closing table and flowchart. I added a short Additional Options section (moments/int/quantile, stabilize, trim(), calibrate(), subclass, ps) because the MatchIt vignette has a parallel "Customizing the Matching Specification" section and the method sections kept referring to these. Cut it if the vignette should stay strictly to what was asked.

3. **The flowchart is a hand-written SVG** at `vignettes/weighting-methods-flowchart.svg`, included with markdown image syntax and embedded as a data URI by html_vignette. The alternatives (DiagrammeR, a mermaid block, drawing it with grid in a chunk) would add a dependency or put drawing code in a prose-only vignette. The boxes use light fills with dark text so they read on both the light and dark pkgdown themes. Check how it looks in dark mode and at phone width.

4. **The summary table (revised 2026-10-08 after your comments) is one wide table inside a scrolling container with the method column fixed.** The pipe table sits in a pandoc fenced div (`::: {.scroll-table}`), and a raw `{=html}` style block at the top of the file gives the div `overflow-x: auto`, makes the first column `position: sticky` with a background of `var(--bs-body-bg, #fff)` (so it follows pkgdown's dark theme and falls back to white in a plain html_vignette), and sets per-column minimum widths so the long text columns wrap at about 18em. I confirmed the rendered HTML has the div wrapping the table and that pandoc emits no column widths, but I could not take a screenshot (the sandbox blocks the browser pane on local files and headless Chrome hung), so please check the scroll and the frozen column in the browser, in both themes. If the sticky column's background ever shows the wrong color, the `--bs-body-bg` variable is the thing to look at.

   The earlier version had two static tables. They could instead be generated in an `echo = FALSE` chunk from `.weightit_methods`, which would keep the capabilities table in sync with the package automatically. I kept them static because the request said no code should run, but a metadata-driven table is not data analysis, so this is worth reconsidering. The second table (strengths, limitations, when to consider) has to be hand-written either way.

5. **Not added to `_pkgdown.yml`.** The file has no `articles:` section, so pkgdown lists vignettes alphabetically and the new one lands between "installing-packages" and "WeightIt". MatchIt's `_pkgdown.yml` orders its articles explicitly. If you want "Weighting Methods" to appear second, after the getting-started vignette, an `articles:` block like MatchIt's is the way to do it.

6. **`TASKS.md` was created at the project root** and added to `.Rbuildignore`, per the global instructions. The repository never had one; delete it if you do not want it here.

7. **A NEWS item was added** at the top of the development version section.

8. **Length.** The vignette runs to about 10,600 words from the introduction to the references, against about 8,600 for MatchIt's matching methods vignette. The tables and the Additional Options section account for most of the difference. Trim there first if it should be shorter.

9. **The ATO, ATM, and ATOS paragraph was shortened** at your request (from about 260 words to about 170). Dropped: the statement that overlap weights minimize the asymptotic variance among tilting functions under homoscedasticity, the explicit tilting function in prose (it is in the table), and the sentence on model misspecification. Those are in the Li, Morgan, and Zaslavsky (2018) and Zhou et al. (2020) citations that remain.

## Bibliography

Sixty-five entries were added to `vignettes/references.bib`. Most were exported from Zotero through Better BibTeX with the bulky fields (abstract, file, keywords) stripped, seven were copied from MatchIt's `references.bib` (`mccaffrey2004`, `desai2017`, `rubin2001`, `stuart2010`, `austin2009`, `rosenbaum1983`, `delosangelesresaDirectStableWeight2020`), and `firthBiasReductionMaximum1993` was written by hand from `?method_glm` because it is not in Zotero.

Updated after your review (2026-10-08):

- `gutmanImprovingInverseProbability2024` now comes from Zotero, with the pages (473–480) from Crossref.
- `vanderlaanStabilizedInverseProbability2025` is the Zotero record (2025, arXiv v3); the vignette cites it as 2025.
- `shook-saPowerSampleSize2022` replaces the 2020 early-online entry with the published details from Crossref: *Biometrics*, 78(1), 388–398, 2022. The Zotero record still has "biom.13405" as its pages, so this entry will drift from Zotero until that record is updated.
- `hulingIndependenceWeightsCausal2024` is now the single entry for the independence weights paper (JASA, 119(546), 1657–1670, 2024, from Crossref); `huling2023` and the arXiv entry were removed. The same update was made in `?method_energy` (the reference and the three in-text years) and in the 2.x NEWS item that introduced energy balancing for continuous treatments. The Zotero record is dated 2023 with "0(0), 1–14" as its pages, so its citation key will change if you update it.
- `ben-michaelBalancingActCausal2021`: see the note in the summary message. I could not find a published version through Crossref or arXiv's journal-ref field, so the entry is unchanged pending the venue.

Entries still worth a look:

- `hiranoPropensityScoreContinuous2005`: the Zotero record gives the year as 2005 and the book title as "Wiley Series in Probability and Statistics"; the chapter is usually cited as 2004, in *Applied Bayesian Modeling and Causal Inference from Incomplete-Data Perspectives* (Gelman and Meng, eds.).
- `greiferChoosingEstimandWhen2021` is still the arXiv preprint, as elsewhere in the package.
- `maoPropensityScoreWeighting2018`: I removed the early-online page number from the Zotero record; the published version is in *Statistical Methods in Medical Research* 28(8).

## Claims to double-check

These are statements in the vignette that go beyond what the help pages say, with where each came from.

- Overlap weights "minimize the asymptotic variance of the weighted estimate among all tilting functions when the outcome variance is constant" (Estimands section). Li, Morgan, and Zaslavsky (2018) prove this under homoscedasticity.
- CBPS "has been found to balance the covariates better and yield less biased estimates than maximum likelihood when the propensity score model is misspecified" (CBPS section). This is from Wyss et al. (2014), whose abstract I read in full; the claim matches it.
- The IPT double robustness statement and the contrast with entropy balancing and just-identified CBPS are paraphrased from `?method_ipt`.
- "@chanGloballyEfficientNonparametric2016 show that a broad class of them ... can attain the semiparametric efficiency bound for the ATE when the set of balanced functions grows with the sample size" (Optimization-Based Methods). Paraphrased from their abstract (sieve estimator, global efficiency).
- The energy distance formula is the standard one; I did not re-derive the weighted version used in the package.
- The description of the characteristic function distance as "closely related to the maximum mean discrepancy" follows the wording on `?method_cfd`. I did not read Santra, Chen, and Park (2026) beyond its abstract.
- The ESS formula is given per treatment group, as `summary.weightit()` reports it.
- The statement that `method = "glm"` does not support M-estimation with a kernel density estimate is from `?method_glm`.
- The remark that evidence for BART-estimated propensity scores "is thinner than for GBM" is my reading of the note on `?method_bart` that most BART causal inference references concern the outcome model.
- "Stable balancing weights ... scale well for moderate numbers of constraints" (large datasets paragraph) is my characterization, not a cited result.

## Things left out on purpose

Longitudinal treatments, censoring weights, sampling weights, clustering, missing data handling, the `by` argument, and `weightit.fit()`, per the request. Random effects are mentioned in one sentence each under `"glm"` and `"bart"` because those are options of the methods themselves.

## Checks run

- No Rd markup in the markdown files (`grep -n '\\[A-Za-z]\+{' NEWS.md README.md vignettes/*.Rmd` returns only math).
- Headings carry no question words, negatives, or counts.
- No announced counts, clefts, or arguments as bare sentence subjects (greps from the r-doc-style skill).
- The vignette renders with `rmarkdown::render()` and every citation key resolves.
