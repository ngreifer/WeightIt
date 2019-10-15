\name{make_full_rank}
\alias{make_full_rank}

\title{
Make a matrix full rank
}
\description{
When writing \link[=method_user]{user-defined methods} for use with \code{weightit}, it may be necessary to take the potentially non-full rank \code{covs} data frame and make it full rank for use in a downstream function. This function performs that operation.
}
\usage{
make_full_rank(mat,
               with.intercept = TRUE)

}
\arguments{
  \item{mat}{
a numeric matrix or data frame to be transformed. Typically this contains covariates. \code{NA}s are not allowed.
}
  \item{with.intercept}{
whether an intercept (i.e., a vector of 1s) should be added to \code{mat} before making it full rank. If \code{TRUE}, the intercept will be used in determining whether a column is linearly dependent on others. Regardless, no intercept will be included in the output.
}
}
\details{
\code{make_full_rank} makes a matrix full rank by removing columns one at a time and determining whether the rank of the matrix changes. If it does not, that column is deleted. First, all columns that only contain one value are deleted. Then, if \code{with.intercept} is set to \code{TRUE}, an intercept column is added to the matrix. After determining which columns can be removed without changing the rank of the matrix, a matrix is returned with only those columns (and not the added intercept).

See example at \code{\link{method_user}}.
}
\value{
An object of the same type as \code{mat} containing only linearly independent columns.
}

\author{
Noah Greifer
}

\seealso{
\code{\link{method_user}}, \code{\link{model.matrix}}
}
\examples{
set.seed(1000)
c1 <- rbinom(10, 1, .4)
c2 <- 1-c1
c3 <- rnorm(10)
c4 <- 10*c3
mat <- data.frame(c1, c2, c3, c4)

make_full_rank(mat) #leaves c2 and c4

make_full_rank(mat, with.intercept = FALSE) #leaves c1, c2, and c4
}
