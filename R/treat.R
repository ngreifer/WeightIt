#Treatment class

#The class every processed treatment carries. It is cobalt's (see
#`cobalt::treat-class`), and so is the `[` method that preserves the attributes across
#subsetting: cobalt registers it and WeightIt does not, so there is one method for one
#class and no chance of the two packages overwriting each other's.
#
#`treat` is the shared contract and `cobalt.treat` a transitional alias, which cobalt
#introduced because WeightIt 2.0.0 registered a competing `[.treat`. Both are set, so the
#object is indistinguishable from one cobalt processed itself and dispatch finds the
#method under whichever name a given cobalt registers it: 5.0.0 has only the alias, later
#versions have both. The alias can go once a cobalt registering `[.treat` is the minimum
#in `Imports:`. It goes first because a multi-category treatment is a factor underneath,
#and `[.factor` would otherwise win and drop every attribute.
.treat_classes <- c("cobalt.treat", "treat")

.set_treat_class <- function(x) {
  class(x) <- unique(c(.treat_classes, class(x)))

  x
}

as.treat <- function(x, process = NULL, censoring = NULL) {
  if (is_null(process)) {
    process <- !inherits(x, "treat")
  }

  arg::arg_flag(process)

  if (process || !has_treat_type(x)) {
    #Multi-category treatments are passed through `factor()`, here and inside
    #`assign_treat_type()`, which drops every attribute; the treatment's name is
    #carried across by hand so that it survives processing as it does for the other
    #treatment types. `weightitMSM()` names its `treat.list` from it.
    treat.name <- .attr(x, "treat.name")

    x <- assign_treat_type(x, censoring = censoring)
    treat.type <- get_treat_type(x)

    if (treat.type %in% c("multinomial", "multi-category")) {
      x <- assign_treat_type(factor(x))
    }

    attr(x, "treat.name") <- treat.name
  }

  .set_treat_class(x)
}
