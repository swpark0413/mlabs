\encoding{UTF-8}
\name{summary.mlabs}
\alias{summary.mlabs}
\alias{print.summary.mlabs}
\title{Posterior summaries of the fixed-dimensional quantities of an MLABS fit}
\usage{
\method{summary}{mlabs}(
  object,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  include_pi = TRUE,
  include_fitted = TRUE,
  fitted_index = NULL,
  ...
)

\method{print}{summary.mlabs}(
  x,
  digits = max(3L, getOption("digits") - 3L),
  max_fitted = 20L,
  ...
)
}
\arguments{
\item{object}{An object from \strong{mlabs_collapsed}.}

\item{probs}{A numeric vector of distinct probabilities in \eqn{(0,1)} giving
the posterior quantiles to report.}

\item{include_pi}{A logical value indicating whether the saved
predictor-selection probabilities \code{pi_var} are summarized, one row per
predictor.}

\item{include_fitted}{A logical value indicating whether posterior summaries of
the regression function at the observed predictor vectors are computed.}

\item{fitted_index}{An optional vector of distinct observation indices at which
the regression function is summarized. \code{NULL} uses all observations.
Restricting the index set is useful for large \eqn{n}, since the regression
function must be re-evaluated for every retained draw.}

\item{...}{Additional arguments passed to or from other methods.}

\item{x}{An object from \strong{summary.mlabs}.}

\item{digits}{The number of significant digits used when printing.}

\item{max_fitted}{The maximum number of regression-function rows printed. All
computed rows remain available in the returned matrices.}
}
\value{
An object of class \code{"summary.mlabs"} containing:
\describe{
  \item{\code{statistics}}{A matrix of empirical posterior means and variances,
    one row per summarized quantity.}
  \item{\code{quantiles}}{A matrix of the requested posterior quantiles, with
    the same rows.}
  \item{\code{n_draws}}{The number of retained posterior draws.}
  \item{\code{saved_iterations}}{The iteration indices of the retained draws,
    or their sequence when the indices are not stored.}
  \item{\code{include_pi}, \code{include_fitted}, \code{probs}}{The options
    used to build the summary.}
  \item{\code{fitted_index}, \code{fitted_names}}{The observation indices and
    row names of the summarized regression-function values.}
}
}
\description{
Computes empirical posterior means, variances, and quantiles for the
fixed-dimensional quantities of an MLABS fit, namely the number of basis
functions, the Lévy total-mass parameter, the noise variance, the training mean
squared residual, the log-likelihood, the predictor-selection probabilities,
and the regression function evaluated at observed predictor vectors.
}
\details{
MLABS is a trans-dimensional model: the number of basis functions changes along
the chain, and a basis function in one draw has no counterpart in another.
Componentwise summaries of the coefficients \eqn{\beta_j}, the degrees
\eqn{c_l^{(j)}}, or the knots \eqn{\xi_l^{(j)}} are therefore not meaningful,
and they are deliberately excluded here. Two kinds of quantities remain
well defined across draws and are reported instead: global parameters whose
dimension is fixed, and functionals of the fitted regression function. For
structure-level inference use \code{\link{marginal_var_prob}} and
\code{\link{basis_analysis}}, which summarize how predictors and
interaction orders are used rather than the identity of individual basis
functions.

The global block contains the basis count \eqn{J}, the total-mass parameter
\eqn{M}, the noise variance \eqn{\sigma^2}, the training mean squared residual,
the log-likelihood, and, when \code{include_pi = TRUE}, each component of
\eqn{\pi}. The posterior of \eqn{J} is the most direct summary of the fitted
model size; mass piling up at \code{NB_max} indicates that the upper bound is
binding and should be raised. When log-likelihood draws are not stored, they are
reconstructed from the saved mean squared residuals and noise standard
deviations under the Gaussian observation model,
\deqn{\log L = -\frac{n}{2}\left\{\log(2\pi\sigma^2)
+ \frac{\mathrm{MSE}}{\sigma^2}\right\}.}

When \code{include_fitted = TRUE}, rows named \code{f[1]}, \code{f[2]}, and so
on report posterior summaries of \eqn{f(\mathbf{x}_i)} at the selected observed
predictor vectors, obtained by evaluating the regression function separately in
each retained draw. Their quantiles give pointwise credible intervals for the
fitted surface.

If the model was fitted with \code{normalize = TRUE}, the noise variance, the
mean squared residual, and the regression-function values are reported on the
original response scale, and the log-likelihood is adjusted by
\eqn{-n\log(\mathrm{sd}(\mathbf{Y}))} accordingly. The print method displays at
most \code{max_fitted} regression-function rows to keep the output readable,
while the returned matrices contain every requested row.
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

summary(fit)

# global parameters only, with 90% intervals
summary(fit, probs = c(0.05, 0.5, 0.95), include_fitted = FALSE)

# pointwise credible intervals for the first ten fitted values
s <- summary(fit, fitted_index = 1:10, include_pi = FALSE)
s$quantiles
}

}
\references{
Park, S., Oh, J., and Lee, J. (2026+), "Multivariate Lévy adaptive B-spline
	regression", \emph{Manuscript}.
}
\seealso{
\code{mlabs_collapsed}, \code{marginal_var_prob},
\code{basis_analysis}, \code{diagnosis}
}
\author{
Sewon Park
}
\keyword{regression}
\keyword{nonparametric}
