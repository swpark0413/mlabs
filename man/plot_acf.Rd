\encoding{UTF-8}
\name{plot_acf}
\alias{plot_acf}
\title{Autocorrelation plots for an MLABS fit}
\usage{
plot_acf(
  fit,
  pars = c("sigma2", "M", "pi"),
  which_pi = NULL,
  fitted_index = NULL,
  lag_max = 50L,
  title = "MLABS Autocorrelation Diagnostics"
)
}
\arguments{
\item{fit}{An object from \strong{mlabs_collapsed}.}

\item{pars}{A character vector of distinct quantities to display, taken from
\code{"sigma2"}, \code{"M"}, \code{"pi"}, \code{"f"}, \code{"J"},
\code{"logLik"}, and \code{"mse"}. Panels are drawn in the given order. See
\sQuote{Details}.}

\item{which_pi}{Used when \code{pars} contains \code{"pi"}: the predictors
whose selection probabilities are shown, given as indices or names.
\code{NULL} uses the three predictors with the largest posterior mean
\eqn{\pi_\ell}.}

\item{fitted_index}{Used when \code{pars} contains \code{"f"}: the distinct
observation indices at which the regression function is shown. \code{NULL}
uses the observations whose responses are closest to the lower quartile, the
median, and the upper quartile of the training response.}

\item{lag_max}{The maximum lag at which the autocorrelation is computed. It is
truncated internally at one less than the number of retained draws.}

\item{title}{The overall title of the combined figure.}
}
\value{
A \pkg{patchwork} object containing one \pkg{ggplot2} panel per chain, arranged
in at most three columns.
}
\description{
Plots the autocorrelation functions of selected scalar chains of an MLABS fit
over the retained draws.
}
\details{
The available chains, all defined in every retained draw \eqn{s}, are
\describe{
  \item{\code{"sigma2"}}{the noise variance \eqn{(\sigma^{(s)})^2};}
  \item{\code{"M"}}{the total-mass parameter \eqn{M^{(s)}};}
  \item{\code{"pi"}}{the predictor-selection probabilities
    \eqn{\pi_\ell^{(s)}}, one panel per predictor in \code{which_pi};}
  \item{\code{"f"}}{the regression function \eqn{f^{(s)}(\mathbf{x}_i)} at the
    observations in \code{fitted_index}, one panel each;}
  \item{\code{"J"}}{the number of active basis functions \eqn{J^{(s)}};}
  \item{\code{"logLik"}}{the log-likelihood;}
  \item{\code{"mse"}}{the training mean squared residual.}
}
The noise variance, the regression-function values, the log-likelihood, and the
mean squared residual are reported on the original response scale after a fit
with \code{normalize = TRUE}; the log-likelihood is reconstructed from the
saved mean squared residuals and noise standard deviations when it is not
stored.

The dashed lines mark the approximate white-noise bounds
\eqn{\pm z_{0.975}/\sqrt{S}} for \eqn{S} retained draws. Autocorrelation is
computed on the thinned chain, so lag one corresponds to \code{nthin} sampler
iterations. A panel reports that the ACF is unavailable when its chain has
fewer than two values, contains non-finite values, or is constant; the last
case occurs for \eqn{\pi} when \code{control$pi_warmup} exceeds \code{nburn}.

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

mlabs::plot_acf(fit)
mlabs::plot_acf(fit, pars = c("sigma2", "f"), fitted_index = c(1, 25, 50))
mlabs::plot_acf(fit, pars = "pi", which_pi = c("x1", "x2"), lag_max = 100L)
}

}
\seealso{
\code{diagnosis}, \code{summary.mlabs}
}
\author{
Sewon Park
}
\keyword{hplot}
\keyword{ts}
