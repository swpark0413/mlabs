\encoding{UTF-8}
\name{marginal_var_prob}
\alias{marginal_var_prob}
\title{Posterior measures of predictor use in an MLABS fit}
\usage{
marginal_var_prob(fit, X_train = fit$x, var_names = NULL, threshold = 0.5)
}
\arguments{
\item{fit}{An object from \strong{mlabs_collapsed}.}

\item{X_train}{The training design used to obtain the predictor names. It must
have the same number of columns as the fitted design.}

\item{var_names}{An optional character vector of distinct predictor names in
the column order of the fitted design. \code{NULL} uses the column names of
\code{X_train}, or \code{x1}, \code{x2}, and so on when these are absent.}

\item{threshold}{A value in \eqn{(0,1)}. Predictors whose MIP exceeds it are
flagged in the returned summary.}
}
\value{
An object of class \code{"marginal_prob"} containing:
\describe{
  \item{\code{summary}}{A data frame with one row per predictor, sorted by
    decreasing MIP, with columns \code{MIP}, \code{ABR}, \code{ICR}, and the
    flag \code{important}.}
  \item{\code{post_MIP}, \code{post_ABR}, \code{post_ICR}}{Data frames of
    posterior mean, standard deviation, median, and 2.5\% and 97.5\% quantiles
    of the three measures.}
  \item{\code{samples_MIP}, \code{samples_ABR}, \code{samples_ICR}}{Matrices of
    the draw-wise values, one row per retained draw and one column per
    predictor.}
  \item{\code{var_names}, \code{threshold}, \code{n_save}, \code{p}}{The
    predictor names, the threshold, the number of retained draws, and the
    number of predictors.}
}
}
\description{
Computes three draw-wise measures of how often each predictor is used by the
fitted basis expansion.
}
\details{
Let \eqn{J^{(s)}} be the number of active basis functions in draw \eqn{s}, let
\eqn{A_j^{(s)}} be the number of them whose active predictor set contains
predictor \eqn{j}, and let \eqn{S} be the number of retained draws. The three
measures are
\deqn{\mathrm{MIP}_j = \frac{1}{S}\sum_{s=1}^{S}
\mathbf{1}\{A_j^{(s)} > 0\}, \qquad
\mathrm{ABR}_j = \frac{1}{S}\sum_{s=1}^{S} \frac{A_j^{(s)}}{J^{(s)}},}
the marginal inclusion probability and the average basis ratio. The
interaction-conditional ratio \eqn{\mathrm{ICR}_j} is the same ratio computed
over the basis functions with interaction order \eqn{K_j \ge 2} only, averaged
over the draws that contain at least one such basis function; it is \code{NA}
when no retained draw does.

All three are computed draw by draw before averaging, so the \code{samples_*}
components describe their posterior distributions. They count how often a
predictor is used and not the magnitude of its effect.
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
mp$post_MIP
}

}
\seealso{
\code{plot_marginal_var_prob}, \code{basis_analysis},
\code{mlabs_collapsed}
}
\author{
Sewon Park
}
\keyword{regression}
