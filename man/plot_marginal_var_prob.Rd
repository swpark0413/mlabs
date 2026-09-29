\encoding{UTF-8}
\name{plot_marginal_var_prob}
\alias{plot_marginal_var_prob}
\title{Plot posterior measures of predictor use}
\usage{
plot_marginal_var_prob(
  mp,
  threshold = NULL,
  title = "Marginal Variable Inclusion Probability",
  save_path = NULL,
  width = 1600,
  height = 1200
)
}
\arguments{
\item{mp}{An object from \strong{marginal_var_prob}.}

\item{threshold}{The value at which the dashed reference line is drawn in every
panel. \code{NULL} uses the threshold stored in \code{mp}.}

\item{title}{The overall title of the combined figure.}

\item{save_path}{An optional file path for a PNG copy of the figure.
\code{NULL} prints the figure on the current graphics device.}

\item{width}{The width in pixels of the PNG file, used only when
\code{save_path} is supplied.}

\item{height}{The height in pixels of the PNG file, used only when
\code{save_path} is supplied.}
}
\value{
The input object \code{mp}, invisibly. The combined figure is printed or written
to \code{save_path} as a side effect.
}
\description{
Displays the posterior means and the draw-wise distributions of the marginal
inclusion probability (MIP), the average basis ratio (ABR), and the
interaction-conditional ratio (ICR) of every predictor in a single six-panel
figure.
}
\details{
The figure shows MIP, ABR, and ICR in separate rows. For each measure, the left
panel displays posterior means and the right panel shows the corresponding
draw-wise distributions. Predictors are ordered by MIP in all panels.

The dashed line indicates \code{threshold}. It is used as an inclusion cutoff
for MIP and as a visual reference for ABR and ICR. Undefined ICR values are
displayed as \code{NA}.

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

mp <- mlabs::marginal_var_prob(fit)

mlabs::plot_marginal_var_prob(mp)
mlabs::plot_marginal_var_prob(mp, threshold = 0.8,
                              save_path = "predictor_use.png")
}

}
\seealso{
\code{marginal_var_prob}, \code{plot_basis_analysis}
}
\author{
Sewon Park
}
\keyword{hplot}
