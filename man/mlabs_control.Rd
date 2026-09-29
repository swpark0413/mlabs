\encoding{UTF-8}
\name{mlabs_control}
\alias{mlabs_control}
\title{Prior hyperparameters and tuning constants for the MLABS sampler}
\usage{
mlabs_control(
  bJ = 1.0,
  sigr = 0.01,
  sigR = 0.01,
  alpha_dir = 5.0,
  beta0 = NULL,
  sigmab = NULL,
  init_J = NULL,
  jump_prob = NULL,
  pi_warmup = NULL,
  ref_sigma = 1.0,
  knot_temperature = NULL,
  boost = 5.0,
  xi_step_size = 0.03,
  n_refresh = 5L,
  odd_anchor_uniform = 0.05
)
}
\arguments{
\item{bJ}{The rate parameter \eqn{b_\gamma} of the Gamma prior on the
total-mass parameter \eqn{M}, whose shape is the \code{aJ} argument of
\code{\link{mlabs_collapsed}}. It was fixed at 1 in the numerical studies of
the paper.}

\item{sigr}{The degrees-of-freedom parameter \eqn{r} in the inverse-Gamma prior
\eqn{\sigma^2 \sim \mathrm{IG}(r/2, rR/2)} on the noise variance.}

\item{sigR}{The scale parameter \eqn{R} in the same inverse-Gamma prior. Small
values of \eqn{r} and \eqn{R} give a weakly informative prior on the noise
level.}

\item{alpha_dir}{The total concentration \eqn{\alpha} of the symmetric
Dirichlet prior \eqn{\pi \sim \mathrm{Dir}(\alpha/p, \ldots, \alpha/p)} on the
predictor-selection probabilities. Smaller values concentrate the posterior on
a few predictors and are useful when many predictors are expected to be
irrelevant.}

\item{beta0}{The fixed intercept \eqn{\beta_0} on the fitting scale. It is held
fixed throughout sampling. \code{NULL} uses the sample mean of the response on
the fitting scale, which is zero when \code{normalize = TRUE}.}

\item{sigmab}{The fixed coefficient prior standard deviation on the fitting
scale. It takes effect only when \code{shrink_c = NULL} in
\code{\link{mlabs_collapsed}}; under the default shrinkage parameterization the
coefficient scale is \eqn{\phi_M}, determined by \code{shrink_c}, and this entry
is ignored. \code{NULL} uses the sample standard deviation of the response on
the fitting scale.}

\item{init_J}{The number of basis functions used to initialize the chain.
\code{NULL} derives a value from \eqn{p} and \code{NB_max}. A value at or above
\code{NB_max} blocks birth moves from the first iteration and triggers a
warning.}

\item{jump_prob}{A numeric vector \eqn{(p_b, p_d, 1-p_b-p_d)} of length three
giving the birth, death, and no-change probabilities of the trans-dimensional
move, with nonnegative entries summing to one. \code{NULL} uses
\code{c(0.4, 0.4, 0.2)}. At the boundary states \eqn{J = 1} and
\eqn{J = N_{B,\max}} a birth or a death is forced, and the reversible-jump
acceptance ratio accounts for the corresponding asymmetry.}

\item{pi_warmup}{The number of initial iterations during which the
predictor-selection probabilities \eqn{\pi} are held at the uniform vector.
This prevents premature concentration on a small number of predictors before
the basis structure has stabilized. \code{NULL} uses \code{nburn \%/\% 2}.}

\item{ref_sigma}{A positive reference noise scale used only to set the default
knot-proposal temperature, \code{knot_temperature = max(1, sd(y)/ref_sigma)},
evaluated on the fitting scale. It is ignored when \code{knot_temperature} is
supplied directly. Because the response is standardized when
\code{normalize = TRUE}, the default reduces to \eqn{\max(1, 1/}
\code{ref_sigma}\eqn{)}, so values below one flatten the residual-guided knot
proposal and values at or above one leave it at \eqn{T = 1}.}

\item{knot_temperature}{The temperature \eqn{T \ge 1} of the residual-guided
knot proposal. Bin \eqn{b} of predictor \eqn{\ell} receives weight proportional
to \eqn{(\overline{|r|}_{\ell b}+\varepsilon_r)^{1/T}}, where
\eqn{\overline{|r|}_{\ell b}} is the mean absolute residual in that bin, so
larger values flatten the proposal and \eqn{T = 1} follows the residual profile
most closely. \code{NULL} uses \code{max(1, sd(y)/ref_sigma)} on the fitting
scale.}

\item{boost}{The strength \eqn{\kappa \ge 0} of the residual guidance in the
predictor proposal, \eqn{\tilde\pi_\ell \propto \pi_\ell(1+\kappa\rho_\ell)},
where \eqn{\rho_\ell \in [0,1]} is the residual relevance score of predictor
\eqn{\ell}. Setting \code{boost = 0} recovers proposals from the current
\eqn{\pi} alone. The guidance affects the proposal only; the prior is
unchanged and the difference enters the Metropolis-Hastings ratio.}

\item{xi_step_size}{The standard deviation \eqn{s_\xi} of the local knot
random walk, expressed as a fraction of the knot-domain width; it must lie in
\eqn{(0,1)}. Proposals falling outside the knot domain or violating the
ordering constraint are rejected immediately, so larger values move knots
further per attempt at the cost of a lower acceptance rate.}

\item{n_refresh}{The maximum number \eqn{n_{\mathrm{ref}}} of existing basis
functions subjected to a joint degree-knot refresh move in each iteration. The
actual number is \eqn{\min(n_{\mathrm{ref}}, J)}, and the basis functions are
selected uniformly without replacement.}

\item{odd_anchor_uniform}{The uniform-mixture weight \eqn{\eta \in (0,1]} in
the residual-guided proposal
\eqn{g_\ell(u \mid r) = (1-\eta)\tilde g_\ell(u \mid r) + \eta U_\ell(u)}
for the continuous central knot of an odd-degree factor. The uniform component
keeps the proposal positive over the whole knot domain, which guarantees that
every region remains reachable. It controls the proposal only, not the
design-based target prior.}
}
\value{
A named list with one element per argument, in the order shown above, suitable
for the \code{control} argument of \code{\link{mlabs_collapsed}}.
}
\description{
Collects the prior hyperparameters and sampler tuning constants that are rarely
changed into a single validated list, so that the signature of
\code{\link{mlabs_collapsed}} keeps only the arguments that are usually set for
a given analysis.
}
\details{
The \code{control} argument of \code{\link{mlabs_collapsed}} accepts either the
result of this function or a plain named list containing any subset of these
names; entries that are not supplied keep their defaults. Passing an empty list
or \code{NULL} therefore means \dQuote{use all defaults}. Entries documented as
\code{NULL} are not constants but placeholders resolved from the data inside
\code{mlabs_collapsed}, which is why they cannot be given a fixed default here.

Each entry is validated on construction, and an unrecognised name produces an
error that lists the allowed names. A misspelled tuning parameter therefore
fails loudly instead of being silently absorbed and ignored.

The entries fall into three groups. The first group, \code{bJ}, \code{sigr},
\code{sigR}, \code{alpha_dir}, \code{beta0}, and \code{sigmab}, completes the
prior specification of \code{\link{mlabs_collapsed}}. The second group,
\code{init_J}, \code{jump_prob}, and \code{pi_warmup}, governs the chain
initialization and the move schedule. The third group, \code{ref_sigma},
\code{knot_temperature}, \code{boost}, \code{xi_step_size}, \code{n_refresh},
and \code{odd_anchor_uniform}, governs the residual-guided proposals and the
local moves; these entries change how efficiently the sampler explores the space of
basis configurations, but not the posterior it targets, because every
residual-guided proposal is corrected by the corresponding Metropolis-Hastings
ratio.
}
\examples{

str(mlabs_control())

# change two entries and keep everything else at its default
ctrl <- mlabs_control(boost = 10, alpha_dir = 1)

\dontrun{

set.seed(413)
n <- 250
X <- matrix(runif(5 * n), ncol = 5)
y <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
     10 * X[, 4] + 5 * X[, 5] + rnorm(n)

# fit MLABS with quadratic B-spline factors and two-way interactions
fit <- mlabs::mlabs_collapsed(y, X, nmcmc = 100000, nburn = 50000, nthin = 10,
                              maxInt = 2, allowed_degrees = 2, aJ = 10,
                              exp_ratio = 1, shrink_c = 5, control = ctrl)

# equivalently, pass the entries directly
fit <- mlabs::mlabs_collapsed(y, X, nmcmc = 100000, nburn = 50000, nthin = 10,
                              maxInt = 2, allowed_degrees = 2, aJ = 10,
                              exp_ratio = 1, shrink_c = 5,
                              control = list(boost = 10, alpha_dir = 1))
}

}
\references{
Park, S., Oh, J., and Lee, J. (2026+), "Multivariate Lévy adaptive B-spline
	regression", \emph{Manuscript}.

DiMatteo, I., Genovese, C. R., and Kass, R. E. (2001), "Bayesian curve-fitting
	with free-knot splines", \emph{Biometrika}, 88(4), 1055–1071.

Linero, A. R. (2018), "Bayesian regression trees for high-dimensional prediction
	and variable selection", \emph{Journal of the American Statistical
	Association}, 113(522), 626–636.
}
\seealso{
\code{mlabs_collapsed}
}
\author{
Sewon Park
}
\keyword{regression}
\keyword{nonparametric}
