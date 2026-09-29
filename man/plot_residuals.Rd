\encoding{UTF-8}
\name{plot_residuals}
\alias{plot_residuals}
\title{Residual diagnostics for an MLABS fit}
\usage{
plot_residuals(fit, Y_train = .training_y(fit), title = "Residual Diagnostics")
}
\arguments{
\item{fit}{An object from \strong{mlabs_collapsed}.}

\item{Y_train}{A numeric vector of observed training responses on the original
response scale, with one value per training observation. The default
reconstructs them from the fitted object.}

\item{title}{The overall title of the combined figure.}
}
\value{
A \pkg{patchwork} object containing three \pkg{ggplot2} panels.
}
\description{
Displays the residuals of an MLABS fit against the posterior mean fitted
values, a normal quantile-quantile plot of the standardized residuals, and a
histogram of the residuals, as checks of the Gaussian observation model.
}
\details{
Residuals are formed as \eqn{y_i - E[f(\mathbf{x}_i) \mid \mathbf{Y}]} on the
original response scale, using the posterior mean surface rather than any single
draw. The three panels display them in the three standard ways used to examine
the Gaussian observation model \eqn{\varepsilon_i \sim N(0, \sigma^2)} of
\code{\link{mlabs_collapsed}}.

The first panel plots the residuals against the posterior mean fitted values,
with a horizontal reference line at zero and a LOESS trend. The second is a
normal quantile-quantile plot of the residuals after centring and scaling, with
the reference line of a standard normal sample. The third is a histogram of the
residuals on the original scale.

Because the residuals are computed at the training observations, they describe
the fit to the data used for estimation. The posterior of \eqn{\sigma^2} from
\code{\link{summary.mlabs}} gives the noise level that the model itself infers.

The function requires the \pkg{ggplot2} and \pkg{patchwork} packages.
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

mlabs::plot_residuals(fit)
}

}
\seealso{
\code{plot_fitted}, \code{summary.mlabs}, \code{diagnosis}
}
\author{
Sewon Park
}
\keyword{hplot}
