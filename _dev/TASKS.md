# Tasks

## In Progress

## To Do

- [ ] Review the notes in `_dev/longitudinal-treatments-vignette-notes.md` and decide on the open questions about `vignette("longitudinal-treatments")`.

- [ ] Review the notes in `_dev/weighting-methods-vignette-notes.md` and decide on the open questions about `vignette("weighting-methods")`.

## Done

- [x] Remove the effects of individual treatments from `vignette("longitudinal-treatments")`, so its estimands are the comparisons between regimes only (2026-10-09).

- [x] Break the saturated MSM equation in `vignette("longitudinal-treatments")` over two lines so it fits the pkgdown margins, and replace the worked code in the longitudinal section of `vignette("estimating-effects")` with a pointer to the new vignette (2026-10-09).

- [x] Add cross-references to `vignette("weighting-methods")` and `vignette("longitudinal-treatments")` in the README, the help pages (`weightit()`, `weightitMSM()`, `.cens()`, `msmdata`, `.weightit_methods`, `get_w_from_ps()`, `ESS()`, `trim()`, `calibrate()`, and the eleven method pages), and the other vignettes (2026-10-08).

- [x] Write `vignette("longitudinal-treatments")`, a worked analysis of a longitudinal treatment with `weightitMSM()`, with a second analysis adding censoring weights (2026-10-08).

- [x] Write `vignette("weighting-methods")`, a prose guide to the weighting methods modeled on MatchIt's `vignette("matching-methods")`, with a summary table and a selection flowchart (2026-10-08).
