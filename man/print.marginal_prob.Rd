\encoding{UTF-8}
\name{print.marginal_prob}
\alias{print.marginal_prob}
\title{Print method for MLABS predictor-use summaries}
\usage{
\method{print}{marginal_prob}(x, ...)
}
\arguments{
\item{x}{An object from \strong{marginal_var_prob}.}

\item{...}{Additional arguments passed to or from other methods.}
}
\value{
The input object \code{x}, invisibly. The table is printed as a side effect.
}
\description{
Prints the posterior mean MIP, ABR, and ICR of every predictor, with a flag for
predictors whose MIP exceeds the stored threshold.
}
\details{
Predictors are listed in decreasing order of MIP. The header reports the number
of predictors, the number of retained draws, and the threshold, and a legend
restates the definition of each measure. An ICR entry is shown as \code{NA}
when no retained draw contains a basis function of interaction order two or
higher. Credible intervals for the three measures are in the \code{post_MIP},
\code{post_ABR}, and \code{post_ICR} components of the object.
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

mp <- mlabs::marginal_var_prob(fit)
mp
}

}
\seealso{
\code{marginal_var_prob}, \code{plot_marginal_var_prob}
}
\author{
Sewon Park
}
\keyword{print}
