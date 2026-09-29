\encoding{UTF-8}
\name{plot_basis_analysis}
\alias{plot_basis_analysis}
\title{Plot the posterior basis structure of an MLABS fit}
\usage{
plot_basis_analysis(
  ba,
  top_n = 15,
  title = "MLABS Basis Analysis",
  save_path = NULL,
  width = 1600,
  height = 1200
)
}
\arguments{
\item{ba}{An object from \strong{basis_analysis}.}

\item{top_n}{The maximum number of predictor pairs displayed in the
co-occurrence panel.}

\item{title}{The overall title of the combined figure.}

\item{save_path}{An optional file path for a PNG copy of the figure.
\code{NULL} prints the figure on the current graphics device.}

\item{width}{The width in pixels of the PNG file, used only when
\code{save_path} is supplied.}

\item{height}{The height in pixels of the PNG file, used only when
\code{save_path} is supplied.}
}
\value{
The input object \code{ba}, invisibly. The figure is printed or written to
\code{save_path} as a side effect.
}
\description{
Displays the components of an \code{\link{basis_analysis}} object in a
four-panel figure.
}
\details{
The first panel is a bar chart of the fraction of all active basis functions at
each interaction order \eqn{K_j}. The second is a horizontal bar chart of the
mean number of univariate basis functions per draw for each predictor. The
third shows the \code{top_n} predictor pairs with the highest co-occurrence,
each bar split into the selected-order and the remaining contributions; when no
interaction basis is present it is replaced by an explanatory message. The
fourth is a bar chart of the factor-degree distribution, as percentages of all
univariate factors.

The function requires the \pkg{ggplot2}, \pkg{patchwork}, and \pkg{scales}
packages.
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

ba <- mlabs::basis_analysis(fit)

mlabs::plot_basis_analysis(ba)
mlabs::plot_basis_analysis(ba, top_n = 10, save_path = "basis_structure.png")
}

}
\seealso{
\code{basis_analysis}, \code{plot_marginal_var_prob}
}
\author{
Sewon Park
}
\keyword{hplot}
