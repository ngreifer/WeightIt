# WeightIt 2.2.0

## Test environments

* local: macOS 27.0.1, aarch64-apple-darwin23, R 4.6.2 (2026-10-29)
* GitHub Actions: see `.github/workflows/`

## R CMD check results

`R CMD check --as-cran` on the built tarball:

```
Status: OK
```

0 errors | 0 warnings | 0 notes

The full test suite was run separately with `NOT_CRAN=true`: 39 test files, 0
failures and 0 errors (`[ FAIL 0 | WARN 0 | SKIP 0 | PASS 8894 ]`). Most tests are
behind `skip_on_cran()` because they fit many models across the optional
weighting-method packages, so the suite CRAN runs is deliberately a small subset,
reporting `[ FAIL 0 | WARN 0 | SKIP 183 | PASS 78 ]` in about 4 seconds.

## Reverse dependencies

All 24 reverse dependencies were checked with `revdepcheck` against both the CRAN
and the development version. There were no new problems, and no package failed to
be checked.

## Notes on this release

`method = "bart"` now fits its models with the *bartisan* package (on CRAN since
version 0.1.0) instead of *dbarts* and *stan4bart*, which are no longer in
`Suggests:`.
