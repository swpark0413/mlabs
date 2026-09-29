\encoding{UTF-8}
\name{print.basis_analysis}
\alias{print.basis_analysis}
\title{Print method for MLABS basis-structure analyses}
\usage{
\method{print}{basis_analysis}(x, top_n = 10, digits = 4, ...)
}
\arguments{
\item{x}{An object from \strong{basis_analysis}.}

\item{top_n}{The maximum number of predictors and of predictor pairs displayed.}

\item{digits}{The number of decimal places used when printing.}

\item{...}{Additional arguments passed to or from other methods.}
}
\value{
The input object \code{x}, invisibly. The summary is printed as a side effect.
}
\description{
Prints the posterior summary of the selected-order basis fraction, the
main-effect basis counts by predictor, the leading pairwise co-occurrence
counts, and the distribution of the B-spline factor degrees.
}
\details{
The output has four blocks. The first reports the posterior mean, median,
standard deviation, and 95\% credible interval of the fraction of active basis
functions with \eqn{K_j \ge} \code{min_Kj}. The second lists predictors by
their mean number of univariate basis functions per draw. The third lists the
predictor pairs that most often occur inside the same basis function, with the
count over all basis functions (\code{all}) and over selected-order ones
(\code{selected}). The fourth gives the percentage of univariate factors at
each degree.

Only the leading \code{top_n} rows of the second and third blocks are shown;
the full tables are in the \code{main_effect} and \code{interaction} components
of the object.
}
\examples{

\dontrun{

set.seed(413)
n <- 250
X <- matrix(runif(5 * n), ncol = 5)
y <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
     10 * X[, 4] + 5 * X[, 5] + rnorm(n)

# fit MLABS with quadratic B-spline factors and two-way interactions
fit <- mlabs::mlabs_collapsed(y, X, nmcmc = 100000, nburn = 50000, nthin = 10,
                              maxInt = 2, allowed_degrees = 2, aJ = 10,
                              exp_ratio = 1, shrink_c = 5)

ba <- mlabs::basis_analysis(fit)
ba
print(ba, top_n = 20)
}

}
\seealso{
\code{basis_analysis}, \code{plot_basis_analysis}
}
\author{
Sewon Park
}
\keyword{print}
