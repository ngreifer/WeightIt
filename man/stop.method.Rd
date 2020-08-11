\name{stop.method}
\alias{stop.method}

\title{
Balance criteria for tuning
}
\description{
Some methods involve tuning of parameters and use some measure of covariate balance as the crierion to select the optimal parameter values. For example with \code{\link[=method_gbm]{method = "gbm"}}, a balance criterion can be used to select the optimal number of trees to use. In \pkg{WeightIt}, the argument \code{stop.method} controls which criterion is optimized, using the naming conventions originally used in \pkg{twang}. This page documents the available arguments to \code{stop.method}. There has been some research on which balance criteria perform better than others; see the References section for some articles. Of course, how each method fares in any given dataset depends on a variety of factors, and multiple methods should be tried and compared before moving forward with a set of weights.

Each treatment type has its own balance criteria available.

\subsection{Binary Treatments}{
\describe{
\item{\code{"es.mean"}, \code{"es.max"}, \code{"es.rms"}}{The average, maximum, or root mean squared absolute standardized mean difference (ASMD) among the covariates, respectively. The ASMD is computed using \code{\link[cobalt:balance.summary]{col_w_smd()}} in \pkg{cobalt}. All covariates are standardized, including binary covariates (note that in \pkg{cobalt}, raw mean differences are used for binary variables by default). The standardization factor (i.e., \code{s.d.denom}) depends on the estimand requested and follows the conventions used in \pkg{cobalt}.}
\item{\code{"ks.mean"}, \code{"ks.max"}, \code{"ks.rms"}}{The average, maximum, or root mean squared Kolmogorov-Smirnov (KS) statistic among the covariates, respectively. The KS statistic is computed using \code{\link[cobalt:balance.summary]{col_w_ks()}} in \pkg{cobalt}.}
\item{\code{"mahalanobis"}}{The sample Mahalanobis distace in the weighted sample. This is similar to \code{"es.rms"} except that variables that are redundant with each other will be downweighted. The sample Mahalanobis distance is computed using the (generalized) inverse of the unweighted covariance matrix computed in the focal group when the estimand is the ATT or ATC and using the (generalized) inverse of the average of the unweighted covariance matrices computed within each treatment group (analogous the the ASMD, which uses the average of the group variances in its denominator).}
\item{\code{"energy.dist"}}{The energy distance between the weighted samples, as described by Huling & Mak (2020). The "improved" energy distance is used when the estimand is not the ATT. (Note that weights directly minimizing the energy distance can be found using \code{\link[=method_energy]{method = "energy"}}.)}
\item{\code{"r2"}}{The pseudo-R2 of a logistic regression of the treatment on the covariates with the weights applied. Franklin et al. (2014) consider a similar metric, the post-matching C-statistic, but the pseduo-R2 accomplishes the same goal without requiring a decision boundary. The pseudo-R2 used is the McKelvy & Zavoina pseudo-R2.}
}
}
\subsection{Multinomial Treatments}{
\describe{
\item{\code{"es.mean"}, \code{"es.max"}, \code{"es.rms"}}{The average, maximum, or root mean squared absolute standardized mean difference (ASMD) among the covariates, respectively, across all pairs of treatments. The ASMD is computed using \code{\link[cobalt:balance.summary]{col_w_smd()}} in \pkg{cobalt}. All covariates are standardized, including binary covariates (note that in \pkg{cobalt}, raw mean differences are used for binary variables by default). The standardization factor (i.e., \code{s.d.denom}) depends on the estimand requested and follows the conventions used in \pkg{cobalt}. The same standardization factor is used across treatment pairs for each covariate. When the estimand is the ATT, only differences between the focal group and each other group are computed.}
\item{\code{"ks.mean"}, \code{"ks.max"}, \code{"ks.rms"}}{The average, maximum, or root mean squared Kolmogorov-Smirnov (KS) statistic among the covariates, respectively, across all pairs of treatments. The KS statistic is computed using \code{\link[cobalt:balance.summary]{col_w_ks()}} in \pkg{cobalt}. When the estimand is the ATT, only differences between the focal group and each other group are computed.}
\item{\code{"energy.dist"}}{The total energy distance among the weighted samples, as described by Huling & Mak (2020). The "improved" energy distance is used when the estimand is not the ATT. (Note that weights directly minimizing the energy distance can be found using \code{\link[=method_energy]{method = "energy"}}.)}
}
}
\subsection{Continuous Treatments}{
\describe{
\item{\code{"p.mean"}, \code{"p.max"}, \code{"p.rms"}}{The average, maximum, or root mean squared absolute Pearson correlation between the treatment and covariates, respectively. The Pearson correlation is computed using \code{\link[cobalt:balance.summary]{col_w_cov()}} in \pkg{cobalt}. The correlation uses the unweighted standard deviations of the treatment and covariates in the denominator.}
\item{\code{"s.mean"}, \code{"s.max"}, \code{"s.rms"}}{The average, maximum, or root mean squared absolute Spearman correlation between the treatment and covariates, respectively. The Spearman correlation is computed using \code{\link[cobalt:balance.summary]{col_w_cov()}} in \pkg{cobalt}. The correlation uses the unweighted standard deviations of the rank-transformed treatment and covariates in the denominator.}
\item{\code{"r2"}}{The model R2 of a linear regression of the treatment on the covariates with the weights applied. The standard R2 is used.}
}
}
}

\references{
Ali, M. S., Groenwold, R. H. H., Pestman, W. R., Belitser, S. V., Roes, K. C. B., Hoes, A. W., de Boer, A., & Klungel, O. H. (2014). Propensity score balance measures in pharmacoepidemiology: A simulation study. Pharmacoepidemiology and Drug Safety, 23(8), 802–811. \doi{10.1002/pds.3574}

Belitser, S. V., Martens, E. P., Pestman, W. R., Groenwold, R. H. H., de Boer, A., & Klungel, O. H. (2011). Measuring balance and model selection in propensity score methods. Pharmacoepidemiology and Drug Safety, 20(11), 1115–1129. \doi{10.1002/pds.2188}

Franklin, J. M., Rassen, J. A., Ackermann, D., Bartels, D. B., & Schneeweiss, S. (2014). Metrics for covariate balance in cohort studies of causal effects. Statistics in Medicine, 33(10), 1685–1699. \doi{10.1002/sim.6058}

Huling, J. D., & Mak, S. (2020). Energy Balancing of Covariate Distributions. ArXiv:2004.13962 [Stat]. \url{https://arxiv.org/abs/2004.13962}

Griffin, B. A., McCaffrey, D. F., Almirall, D., Burgette, L. F., & Setodji, C. M. (2017). Chasing Balance and Other Recommendations for Improving Nonparametric Propensity Score Models. Journal of Causal Inference, 5(2). \doi{10.1515/jci-2015-0026}

}
\seealso{
\code{\link{method_gbm}}, \code{\link{method_super}}, which use these balance criteria.

\code{\link[cobalt:balance.summary]{balance summary}} in \pkg{cobalt} for details of some calculations.
}