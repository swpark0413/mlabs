\encoding{UTF-8}
\name{mlabs_collapsed}
\alias{mlabs_collapsed}
\title{Bayesian nonparametric regression with multivariate Lévy adaptive B-splines (MLABS)}
\usage{
mlabs_collapsed(
  y,
  X,
  nmcmc = 200000,
  nburn = 100000,
  nthin = 10,
  NB_max = 100,
  maxInt = 2,
  allowed_degrees = c(0, 1, 2),
  aJ = 5,
  shrink_c = 3.0,
  exp_ratio = 0.1,
  normalize = TRUE,
  lambda_J = 0.0,
  k_weights = NULL,
  verbose = TRUE,
  control = list()
)
}
\arguments{
\item{y}{A numeric response vector of length \eqn{n}.}

\item{X}{A numeric \eqn{n \times p} predictor matrix, or a data frame whose
columns are all numeric.}

\item{nmcmc}{The total number of MCMC iterations.}

\item{nburn}{The number of burn-in iterations to discard. It must be smaller
than \code{nmcmc}.}

\item{nthin}{The thinning interval for storing posterior samples. After
burn-in, iterations divisible by \code{nthin} are retained, which yields
\code{(nmcmc - nburn) \%/\% nthin} saved draws.}

\item{NB_max}{An integer giving the maximum number of active basis functions,
\eqn{N_{B,\max}}. The sampler restricts the state space to
\eqn{1 \le J \le N_{B,\max}}, so a birth move is always rejected at
\eqn{J = N_{B,\max}} and a death move at \eqn{J = 1}. It acts as a
computational safeguard rather than as a modelling assumption, and should be
set large enough that the posterior of \eqn{J} is not truncated.}

\item{maxInt}{An integer giving the maximum interaction order
\eqn{K_{\max}} of a single basis function. A basis with \eqn{K_j = 1}
represents a nonlinear main effect, and a basis with \eqn{K_j \ge 2}
represents an interaction among \eqn{K_j} distinct predictors. Small values
such as \eqn{K_{\max} \le 3} are recommended, since higher orders enlarge the
knot space and the cost of trans-dimensional updates.}

\item{allowed_degrees}{An integer vector giving the set \eqn{S} of permitted
univariate B-spline degrees. Degrees 0, 1, 2, and 3 are supported. Allowing
several degrees lets components with different smoothness coexist in one
expansion; including degree 0 permits locally constant factors and hence
discontinuous fits.}

\item{aJ}{The shape parameter \eqn{a_\gamma} of the Gamma prior on the
total-mass parameter \eqn{M}. Because \eqn{J \mid M \sim
\mathrm{Poi}(Me^{-\lambda_J})}, larger values favour larger expansions a
priori.}

\item{shrink_c}{The coefficient-scale multiplier \eqn{c_{\mathrm{shrink}}} in
the conditional coefficient prior; see \sQuote{Details}. Setting
\code{shrink_c = NULL} disables this scaling and uses the fixed prior standard
deviation \code{control$sigmab} instead, which no longer adapts to the expected
model size. The scale of the coefficient prior is controlled by whichever of
the two is in force: \code{shrink_c} in the default shrinkage parameterization,
and \code{control$sigmab} only when \code{shrink_c = NULL}.}

\item{exp_ratio}{The boundary-expansion ratio \eqn{E}. The knot domain of
predictor \eqn{\ell} is the observed range expanded by \eqn{E} times its width
at each end, so knots can be placed slightly outside the observed data and
basis functions can carry signal near the boundary.}

\item{normalize}{A logical value indicating whether the predictors are min-max
scaled to \eqn{[0,1]} and the response is standardized before fitting.
Predictions and posterior summaries are returned on the original response
scale.}

\item{lambda_J}{A nonnegative model-size penalty \eqn{\lambda_J}. The
conditional Poisson mean of the basis count is \eqn{\bar J_M = Me^{-\lambda_J}},
so larger values shrink the expansion towards fewer basis functions.}

\item{k_weights}{An optional nonnegative vector
\eqn{(w_1, \ldots, w_{K_{\max}})} of length \code{maxInt} giving the
categorical prior on the interaction order, \eqn{K_j \sim
\mathrm{Categorical}(w_1, \ldots, w_{K_{\max}})}. \code{NULL} uses equal
weights. For \code{maxInt = 2}, \code{c(0.8, 0.2)} favours main effects while
still allowing interactions, whereas \code{c(0, 1)} restricts every basis
function to a two-way interaction. The vector is normalized internally.}

\item{verbose}{A logical value indicating whether the resolved settings and
periodic progress messages are printed during sampling.}

\item{control}{A named list of sampler tuning parameters, or an empty list
(the default) for the standard settings. Entries that are not supplied keep
their defaults, and an unrecognised name raises an error listing the allowed
names. See \code{\link{mlabs_control}}.}
}
\value{
A list of class \code{"mlabs"} containing:
\describe{
  \item{\code{x}, \code{y}}{The training predictors and response on the
    fitting scale.}
  \item{\code{prior}}{The prior settings actually used, including the knot
    boundary matrix, the allowed degree set, and the resolved knot
    temperature.}
  \item{\code{control}}{The resolved list of tuning parameters, as returned by
    \code{\link{mlabs_control}}.}
  \item{\code{beta0}}{The fixed intercept on the fitting scale.}
  \item{\code{beta}}{A list of posterior samples of the coefficient vectors
    (length: number of retained iterations). The length of each element equals
    the basis count in that draw.}
  \item{\code{xi}, \code{C}, \code{nu}, \code{K}}{Lists of posterior samples of
    the knot sequences, spline degrees, zero-based active predictor indices,
    and interaction orders of every basis function.}
  \item{\code{sigma}}{A numeric vector of posterior samples of the noise
    standard deviation on the fitting scale.}
  \item{\code{J}, \code{M}}{Numeric vectors of posterior samples of the basis
    count and of the total-mass parameter.}
  \item{\code{pi_var}}{A matrix of posterior samples of the predictor-selection
    probabilities, with one row per retained draw.}
  \item{\code{mse}}{A numeric vector of training mean squared residuals on the
    fitting scale, one per retained draw.}
  \item{\code{scale}}{The normalization information used to map between the
    fitting and original scales, or \code{NULL} when
    \code{normalize = FALSE}.}
  \item{\code{type}}{The character string \code{"regression"}.}
  \item{\code{mcmc}}{The chain settings \code{nmcmc}, \code{nburn},
    \code{nthin} and the number of saved draws.}
  \item{\code{acceptance}}{The empirical birth and death acceptance rates.}
}
}
\description{
Provides posterior samples of the basis structure and of the global parameters
of the Multivariate Lévy Adaptive B-Spline Regression (MLABS) model, in which
the unknown regression function is a random sum of low-order tensor-product
B-spline basis functions under a finite Lévy random measure prior. The number
of basis functions, their interaction orders, active predictors, spline
degrees, and knot locations are inferred jointly from the data by a collapsed
reversible-jump MCMC algorithm.
}
\details{
Park, Oh, and Lee (2026+) proposed MLABS as a multivariate extension of Lévy
adaptive B-spline regression (LABS; Park, Oh, and Lee, 2023), which in turn
specializes the Lévy adaptive regression kernels framework (LARK; Wolpert,
Clyde, and Tu, 2011) to B-spline generating functions. For observations
\eqn{(\mathbf{x}_i, Y_i)_{i=1}^n} with \eqn{Y_i = f(\mathbf{x}_i)+\varepsilon_i}
and \eqn{\varepsilon_i \sim N(0,\sigma^2)}, the regression function is
represented as
\deqn{f(\mathbf{x}) = \beta_0 + \sum_{j=1}^{J} \beta_j
\mathbf{B}_j(\mathbf{x}; \psi_j), \qquad
\mathbf{B}_j(\mathbf{x}; \psi_j) = \prod_{l=1}^{K_j}
B_{c_l^{(j)}}\!\left(x_{\nu_l^{(j)}}; \xi_l^{(j)}\right),}
where \eqn{B_c(\cdot\,; \xi)} is a degree-\eqn{c} B-spline determined by an
ordered knot sequence \eqn{\xi = (\xi_1, \ldots, \xi_{c+2})}, and
\eqn{\psi_j = (K_j, \nu^{(j)}, c^{(j)}, \xi^{(j)})} collects the interaction
order, the active predictors, the degrees, and the knot sequences of the
\eqn{j}th basis function. Because each basis is a product of at most
\eqn{K_{\max}} univariate factors, interactions are represented explicitly
without constructing a full \eqn{p}-dimensional tensor-product expansion.

Writing \eqn{\bar J_M = Me^{-\lambda_J}} for the expected number of basis
terms, the hierarchical model is
\deqn{\sigma^2 \sim \mathrm{IG}(r/2, rR/2), \qquad
M \sim \mathrm{Ga}(a_\gamma, b_\gamma), \qquad
J \mid M, \lambda_J \sim \mathrm{Poi}(\bar J_M),}
\deqn{\beta_j \mid M \sim N(0, \phi_M^2), \qquad
\phi_M = c_{\mathrm{shrink}}\,\mathrm{sd}(\mathbf{Y})/\sqrt{\bar J_M},}
\deqn{\pi \sim \mathrm{Dir}(\alpha/p, \ldots, \alpha/p), \qquad
K_j \sim \mathrm{Categorical}(w_1, \ldots, w_{K_{\max}}), \qquad
\nu^{(j)} \mid K_j, \pi \sim \mathrm{WOR}_{K_j}(\pi),}
where \eqn{\mathrm{WOR}_{K}(\pi)} denotes weighted sampling of \eqn{K} distinct
predictors without replacement. Sampling the active predictors without
replacement prevents a predictor from appearing twice in one tensor-product
basis, so \eqn{K_j} retains its interpretation as an interaction order. The
coefficient scale \eqn{\phi_M} depends on \eqn{\bar J_M} rather than on the
realized \eqn{J}, so the prior signal variance does not grow with the expected
number of basis terms and the prior scale of an existing coefficient is
unchanged by birth and death moves. Each degree is uniform on the allowed set
\eqn{S}. Conditional on \eqn{M}, the construction is a finite Lévy random
measure on the joint structural space, so only a random finite subset of the
candidate basis family is active and the representation is sparse.

Since both the number and the structure of the basis functions are unknown,
posterior inference requires trans-dimensional updates. The sampler (Green,
1995) performs four kinds of updates in each sweep:
\itemize{
\item Bayesian backfitting (Hastie and Tibshirani, 2000) of the coefficients,
  together with local Metropolis-Hastings random-walk updates of the knot
  coordinates, with step size \code{control$xi_step_size};
\item joint degree-knot refresh moves for at most \code{control$n_refresh}
  existing basis functions, in which the degree and the knot sequence of a
  basis are proposed together because changing a degree also changes the
  dimension of its knot vector;
\item collapsed birth and death moves for the number of basis functions;
\item conjugate updates of \eqn{\sigma^2} and \eqn{M}, and a Metropolis-Hastings
  update of \eqn{\pi} with a DART-style Dirichlet independence proposal
  (Linero, 2018), with \eqn{\pi} held at the uniform vector during the first
  \code{control$pi_warmup} iterations.
}

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

# override two tuning parameters only; the rest keep their defaults
fit2 <- mlabs::mlabs_collapsed(y, X, nmcmc = 100000, nburn = 50000, nthin = 10,
                               control = list(sigmab = 2, alpha_dir = 10.0))

pred <- predict(fit, X)
}

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

Hastie, T. and Tibshirani, R. (2000), "Bayesian backfitting",
	\emph{Statistical Science}, 15(3), 196–223.

DiMatteo, I., Genovese, C. R., and Kass, R. E. (2001), "Bayesian curve-fitting
	with free-knot splines", \emph{Biometrika}, 88(4), 1055–1071.

Linero, A. R. (2018), "Bayesian regression trees for high-dimensional prediction
	and variable selection", \emph{Journal of the American Statistical
	Association}, 113(522), 626–636.
}
\seealso{
\code{mlabs_control}, \code{predict.mlabs}, \code{summary.mlabs},
\code{marginal_var_prob}, \code{basis_analysis}, \code{diagnosis}
}
\author{
Sewon Park
}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}
