\encoding{UTF-8}
\name{plot_fitted}
\alias{plot_fitted}
\title{Plot fitted values against observed responses}
\usage{
plot_fitted(fit, Y_train = .training_y(fit), title = "Fitted vs Observed")
}
\arguments{
\item{fit}{An object from \strong{mlabs_collapsed}.}

\item{Y_train}{A numeric vector of observed training responses on the original
response scale, with one value per training observation. The default
reconstructs them from the fitted object.}

\item{title}{The plot title.}
}
\value{
A \pkg{ggplot2} object.
}
\description{
Plots the posterior mean fitted values against the observed training responses,
with the identity line added, as a global check of how well the fitted
expansion reproduces the data.
}
\details{
The fitted values are the posterior means \eqn{E[f(\mathbf{x}_i) \mid
\mathbf{Y}]} obtained by evaluating the regression function in every retained
draw and averaging, so the display summarizes the posterior mean surface rather
than any single draw. Both axes use the original response scale, and the panel
is drawn with a fixed aspect ratio and common limits, so the identity line has
unit slope and departures from it can be read directly.

The comparison uses the training observations, so it describes the fit to the
data used for estimation rather than predictive performance; the latter is
assessed on held-out data with \code{\link{predict.mlabs}}.

The function requires the \pkg{ggplot2} package.
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

mlabs::plot_fitted(fit)
}

}
\seealso{
\code{plot_residuals}, \code{predict.mlabs}, \code{diagnosis}
}
\author{
Sewon Park
}
\keyword{hplot}
