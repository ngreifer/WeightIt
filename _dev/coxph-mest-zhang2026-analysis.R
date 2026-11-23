# Summarizes the runs of _dev/coxph-mest-zhang2026.R (see that file)
root <- "/Users/NoahGreifer/Dropbox/Research/R/WeightIt"
rd <- function(v) readRDS(file.path(root, "_dev", sprintf("coxph-mest-zhang2026-%s.rds", v)))
r <- Reduce(function(a, b) merge(a, b, by = c("rep", "spec", "method")),
            list(setNames(rd("v1709")[c("rep", "spec", "method", "est", "se")],
                          c("rep", "spec", "method", "est_1709", "se_1709")),
                 setNames(rd("cran")[c("rep", "spec", "method", "est", "se")],
                          c("rep", "spec", "method", "est_cran", "se_cran")),
                 setNames(rd("dev")[c("rep", "spec", "method", "est", "se", "vcov_type")],
                          c("rep", "spec", "method", "est", "se_new", "vcov_type"))))
stopifnot(nrow(r) == 80000, all(r$vcov_type == "asympt"))
cat("max |est diff| 1709 vs dev:", max(abs(r$est_1709 - r$est)),
    " cran vs dev:", max(abs(r$est_cran - r$est)), "\n")
cat("max |SE diff| 1709 vs 2.1.0:", max(abs(r$se_1709 - r$se_cran)), "\n")

z <- qnorm(.975)
r$rej_old <- abs(r$est / r$se_1709) > z
r$rej_new <- abs(r$est / r$se_new) > z
r$rel <- r$se_new / r$se_1709 - 1

s <- do.call(rbind, lapply(split(r, list(r$spec, r$method)), function(d) {
  esd <- sd(d$est)
  data.frame(spec = d$spec[1], method = d$method[1],
             ESD = esd,
             meanSE_old = mean(d$se_1709), meanSE_new = mean(d$se_new),
             SEratio_old = mean(d$se_1709) / esd, SEratio_new = mean(d$se_new) / esd,
             relchange_mean_pct = 100 * mean(d$rel),
             relchange_maxabs_pct = 100 * max(abs(d$rel)),
             relchange_p99_pct = 100 * quantile(abs(d$rel), .99),
             RR_old = 100 * mean(d$rej_old), RR_new = 100 * mean(d$rej_new),
             n_flips = sum(d$rej_old != d$rej_new))
}))
s <- s[order(s$method, match(s$spec, c("oracle", "over", "omit", "wrong"))), ]
rownames(s) <- NULL
print(format(s, digits = 4), row.names = FALSE)
saveRDS(s, file.path(root, "_dev", "coxph-mest-zhang2026-summary.rds"))
