\name{ESS}
\alias{ESS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute effective sampole size of weighted sample
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
%%  ~~ If necessary, more details than the description above ~~
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (w)
{
    sum(w)^2/sum(w^2)
  }
}
