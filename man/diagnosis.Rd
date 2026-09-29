\encoding{UTF-8}
\name{diagnosis}
\alias{diagnosis}
\title{Trace plots for assessing MCMC convergence of an MLABS fit}
\usage{
diagnosis(
  fit,
  title = "MLABS MCMC Trace Diagnostics",
  smooth = TRUE,
  span = 0.25
)
}
\arguments{
\item{fit}{An object from \strong{mlabs_collapsed}.}

\item{title}{The overall title of the combined figure.}

\item{smooth}{A logical value indicating whether a locally weighted trend is
overlaid on each trace.}

\item{span}{A numeric value in \eqn{(0,1)} giving the fraction of draws used in
each local fit. Smaller values give a more responsive trend.}
}
\value{
A \pkg{patchwork} object containing three vertically arranged \pkg{ggplot2}
panels. It is displayed automatically when the function is called at the
console.
}
\description{
Plots the traces of the number of basis functions, the log-likelihood, and the
noise variance over the retained draws, so that mixing and convergence of the
reversible-jump sampler can be assessed.
}
\details{
The dimension of the MLABS state changes along the chain, so a fixed parameter
vector cannot be traced. The three scalar summaries plotted here are defined in
every draw \eqn{s}:
\itemize{
\item the number of active basis functions \eqn{J^{(s)}};
\item the log-likelihood
  \deqn{\log L^{(s)} = -\frac{n}{2}\left\{\log\!\left(2\pi(\sigma^{(s)})^2\right)
  + \frac{1}{(\sigma^{(s)})^2}\frac{1}{n}\sum_{i=1}^{n}
  \left(y_i - f^{(s)}(\mathbf{x}_i)\right)^2\right\};}
\item the noise variance \eqn{(\sigma^{(s)})^2}.
}


The function requires the \pkg{ggplot2} and \pkg{patchwork} packages.
\code{\link{plot_acf}} shows the autocorrelation of the same three chains, and
\code{fit$acceptance} records the empirical birth and death acceptance rates.
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

mlabs::diagnosis(fit)
mlabs::diagnosis(fit, smooth = FALSE)
fit$acceptance
}

}
\references{
Green, P. J. (1995), "Reversible jump Markov chain Monte Carlo computation and
	Bayesian model determination", \emph{Biometrika}, 82(4), 711–732.

Park, S., Oh, J., and Lee, J. (2026+), "Multivariate Lévy adaptive B-spline
	regression", \emph{Manuscript}.
}
\seealso{
\code{plot_acf}, \code{summary.mlabs}, \code{plot_residuals}
}
\author{
Sewon Park
}
\keyword{hplot}
