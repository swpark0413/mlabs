\encoding{UTF-8}
\name{basis_analysis}
\alias{basis_analysis}
\title{Posterior analysis of the basis structure of an MLABS fit}
\usage{
basis_analysis(fit, X_train = fit$x, var_names = NULL, min_Kj = 2L)
}
\arguments{
\item{fit}{An object from \strong{mlabs_collapsed}.}

\item{X_train}{The training design used to obtain the predictor names. It must
have the same number of columns as the fitted design.}

\item{var_names}{An optional character vector of distinct predictor names in
the column order of the fitted design. \code{NULL} uses the column names of
\code{X_train}, or \code{x1}, \code{x2}, and so on when these are absent.}

\item{min_Kj}{A positive integer giving the minimum interaction order used in
the selected-order summaries. The default \code{2L} counts every basis
function that involves an interaction.}
}
\value{
An object of class \code{"basis_analysis"} containing:
\describe{
  \item{\code{pi_j}}{A posterior summary (draws, mean, standard deviation,
    median, and 2.5\% and 97.5\% quantiles) of the fraction of active basis
    functions with \eqn{K_j \ge} \code{min_Kj}.}
  \item{\code{main_effect}}{A data frame of univariate (\eqn{K_j = 1}) basis
    counts by predictor, as a mean per draw and as a total.}
  \item{\code{interaction}}{A data frame of pairwise co-occurrence counts for
    every predictor pair, over all basis functions and over selected-order
    basis functions.}
  \item{\code{degree_dist}}{A data frame of counts and proportions of the
    univariate factor degrees used.}
  \item{\code{order_dist}}{A data frame of counts and proportions of the
    interaction orders \eqn{K_j} used.}
  \item{\code{Kj_count}}{A data frame with the selected-order and total basis
    counts in each retained draw.}
  \item{\code{inter_mat_all}, \code{inter_mat_Kj}}{Symmetric \eqn{p \times p}
    matrices of mean pairwise co-occurrence counts per draw.}
  \item{\code{p}, \code{var_names}, \code{min_Kj}, \code{n_save}}{The number of
    predictors, the predictor names, the order threshold, and the number of
    retained draws.}
}
}
\description{
Summarizes how the retained posterior draws use the tensor-product dictionary:
the fraction of basis functions that are interactions, the univariate basis
counts by predictor, the co-occurrence of predictor pairs within a basis
function, and the distribution of the B-spline degrees.
}
\details{
The structural parameter of the \eqn{j}th basis function is
\eqn{\psi_j = (K_j, \nu^{(j)}, c^{(j)}, \xi^{(j)})}, and every component is
sampled. This function aggregates those components over the retained draws.

A basis function with \eqn{K_j = 3} contributes to all three of its predictor
pairs, so the co-occurrence counts are pairwise marginals of the interaction
structure rather than a decomposition of it. Degrees are counted once per
univariate factor. All quantities are counts of basis functions, not effect
sizes.
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

ba <- mlabs::basis_analysis(fit, min_Kj = 2L)
ba
head(ba$interaction)
ba$order_dist
}

}
\seealso{
\code{plot_basis_analysis}, \code{marginal_var_prob},
\code{mlabs_collapsed}
}
\author{
Sewon Park
}
\keyword{regression}
