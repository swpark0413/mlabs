\encoding{UTF-8}
\name{mlabs-package}
\alias{mlabs-package}
\alias{mlabs}
\docType{package}
\title{
\packageTitle{mlabs}
}
\description{
\packageDescription{mlabs}
}
\details{
The \pkg{mlabs} package implements Multivariate Lévy Adaptive B-Spline
Regression (MLABS), a Bayesian nonparametric regression model for functions
with nonlinear effects, low-order interactions, and locally varying features.
The regression function is represented as a random sum of low-order
tensor-product B-spline basis functions under a finite Lévy random measure
prior,
\deqn{f(\mathbf{x}) = \beta_0 + \sum_{j=1}^{J}\beta_j
\prod_{l=1}^{K_j} B_{c_l^{(j)}}\!\left(x_{\nu_l^{(j)}}; \xi_l^{(j)}\right),}
so that the number of basis terms, the interaction orders, the active
predictors, the spline degrees, and the knot locations are all inferred from
the data rather than fixed in advance.

Three consequences of this construction motivate the package. The Lévy prior
makes only a finite random subset of the candidate basis family active, so the
fitted expansion is sparse. B-spline factors have local support and free
degrees and knots, so smooth global trends and sharp localized features, including
discontinuities represented by degree-zero factors, can coexist in a single
fit. Finally, each basis is a product of at most \eqn{K_{\max}} univariate
factors, so interactions are represented explicitly while the number of basis
functions required for a target accuracy grows only polynomially in the ambient
dimension \eqn{p}, in contrast to a full \eqn{p}-dimensional tensor-product
spline space.

Model fitting is performed by \code{\link{mlabs_collapsed}}, whose rarely
changed prior hyperparameters and tuning constants are collected in
\code{\link{mlabs_control}}. Posterior computation uses a reversible-jump MCMC
sampler combining Bayesian backfitting, collapsed birth and death moves, joint
degree-knot refresh moves, and residual-guided structural proposals. The fitted
object of class \code{"mlabs"} supports the following workflow:
\describe{
  \item{Prediction}{\code{\link{predict.mlabs}}.}
  \item{Posterior summaries}{\code{\link{summary.mlabs}} for the basis count,
    the total mass, the noise variance, the log-likelihood, the
    predictor-selection probabilities, and pointwise summaries of the fitted
    regression function.}
  \item{Predictor importance}{\code{\link{marginal_var_prob}} and
    \code{\link{plot_marginal_var_prob}}.}
  \item{Basis structure}{\code{\link{basis_analysis}} and
    \code{\link{plot_basis_analysis}}.}
  \item{MCMC diagnostics}{\code{\link{diagnosis}} and
    \code{\link{plot_acf}}.}
  \item{Fit diagnostics}{\code{\link{plot_fitted}} and
    \code{\link{plot_residuals}}.}
}

Predictors must be numeric and complete; categorical predictors and missing
values are not handled automatically. Because inference is based on MCMC,
results depend on the random seed, and independent chains should be compared
before posterior quantities are reported.
}
\author{
\packageAuthor{mlabs}

Maintainer: \packageMaintainer{mlabs}
}
\references{
Park, S., Oh, J., and Lee, J. (2026+), "Multivariate Lévy adaptive B-spline
	regression", \emph{Manuscript}.

Park, S., Oh, H.-S., and Lee, J. (2023), "Lévy adaptive B-spline regression via
	overcomplete systems", \emph{Statistica Sinica}, 33(4), 2715–2737.

Wolpert, R. L., Clyde, M. A., and Tu, C. (2011), "Stochastic expansions using
	continuous dictionaries: Lévy adaptive regression kernels",
	\emph{The Annals of Statistics}, 39(4), 1916–1962.

Green, P. J. (1995), "Reversible jump Markov chain Monte Carlo computation and
	Bayesian model determination", \emph{Biometrika}, 82(4), 711–732.
}
\seealso{
\code{mlabs_collapsed}, \code{mlabs_control}
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
mlabs::marginal_var_prob(fit)
mlabs::diagnosis(fit)
}

}
\keyword{package}
\keyword{regression}
\keyword{nonparametric}
