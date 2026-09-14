# WeightIt 2.1.0

## Test environments

* local: macOS Tahoe 26.5.2, aarch64-apple-darwin23, R 4.6.1 (2026-06-24)
* GitHub Actions: see `.github/workflows/`

## R CMD check results

`R CMD check --as-cran` on the built tarball:

```
Status: OK
```

0 errors | 0 warnings | 0 notes

The full test suite was run separately with `NOT_CRAN=true`: 36 test files, 0
failures and 0 errors. Most tests are behind `skip_on_cran()` because they fit many
models across the optional weighting-method packages, so the suite CRAN runs is
deliberately a small subset, reporting
`[ FAIL 0 | WARN 0 | SKIP 145 | PASS 46 ]` in about 6 seconds.

## Reverse dependencies

All 23 reverse dependencies were checked with `revdepcheck` against both the CRAN
and the development version. One was newly broken:

* *cobalt* 5.0.0 -- two test failures, caused by a *WeightIt* bug that this release
  fixes: a censoring marker written as `WeightIt::.cens(C)` was named after the whole
  call rather than after the indicator, and *cobalt*'s fixtures use that spelling.
  The revdepcheck run above predates the fix. Re-checking *cobalt* 5.0.0 against the
  fixed version gives `Status: OK`, with `checking tests ... OK`.

Nothing else regressed. *jointVIP* fails identically against both versions, so its
error is unrelated to this release, and *mvGPS* loses a warning. The remaining 20 are
unchanged.

## Notes on this release

The `.cens()` function, added in 2.0.0, is now re-exported from *cobalt* rather
than defined here. *cobalt* 5.0.0 added its own `.cens()` so that
`cobalt::bal.tab()` could recognize a censoring model, and two identically named
exports meant one masked the other whenever both packages were attached. There is
now a single function, and `Imports:` requires `cobalt (>= 5.0.0)` accordingly.
The behavior of `.cens()` is unchanged.
