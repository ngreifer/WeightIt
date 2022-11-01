\name{ESS}
\alias{ESS}
\title{
Compute effective sample size of weighted sample
}
\description{
Computes the effective sample size (ESS) of a weighted sample, which represents the size of an unweighted sample with approximately the same amount of precision as the weighted sample under consideration.}
\usage{
ESS(w)
}
\arguments{
  \item{w}{
a vector of weights
}
}
\details{
The ESS is calculated as \eqn{(\sum w)^2/\sum w^2}.
}

\references{
McCaffrey, D. F., Ridgeway, G., & Morral, A. R. (2004). Propensity Score Estimation With Boosted Regression for Evaluating Causal Effects in Observational Studies. Psychological Methods, 9(4), 403–425. \doi{10.1037/1082-989X.9.4.403}

Shook‐Sa, B. E., & Hudgens, M. G. (2020). Power and sample size for observational studies of point exposure effects. Biometrics, biom.13405. \doi{doi.org/10.1111/biom.13405}
}

\seealso{
\fun{summary.weightit}
}
\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "glm", estimand = "ATE"))
summary(W1)
ESS(W1$weights[W1$treat == 0])
ESS(W1$weights[W1$treat == 1])
}
