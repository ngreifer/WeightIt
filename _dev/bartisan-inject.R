# Load the development version of WeightIt and replace the dbarts-based
# `weightit2bart*()` functions in its namespace with the bartisan-based ones from
# `_dev/weightit2bart_bartisan.R`, patching `.weightit_methods` to match. Sourced by
# `_dev/test-weightit2bart_bartisan.R` and other scripts in `_dev/` that exercise
# the new method. Run from the package root.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

stopifnot(requireNamespace("bartisan", quietly = TRUE))

# --- Inject the new implementation ---------------------------------------------

# The namespace is locked, so the new helper functions cannot be added to it.
# Instead, the new functions live in an environment whose parent is the namespace,
# so they can see both the helpers and WeightIt's internals, and only the four
# `weightit2bart*()` functions, which already exist, are replaced in the namespace.
ns <- asNamespace("WeightIt")

new_env <- new.env(parent = ns)
sys.source("_dev/weightit2bart_bartisan.R", envir = new_env)

for (fn in ls(new_env, all.names = TRUE)) {
  f <- get(fn, envir = new_env)
  if (!is.function(f)) next
  environment(f) <- new_env
  assign(fn, f, envir = new_env)
  if (exists(fn, envir = ns, inherits = FALSE)) {
    assignInNamespace(fn, f, ns = "WeightIt")
  }
}

methods <- ns$.weightit_methods
methods$bart$packages_needed <- "bartisan"
methods$bart$package_versions_needed <- NULL
methods$bart$s.weights_ok <- TRUE
methods$bart$density_ok <- FALSE
assignInNamespace(".weightit_methods", methods, ns = "WeightIt")

# Keep the old implementation around for the comparison at the end
old_env <- new.env()
sys.source("R/weightit2bart.R", envir = old_env)
for (fn in ls(old_env, all.names = TRUE)) {
  f <- get(fn, envir = old_env)
  if (is.function(f)) environment(f) <- ns
}

