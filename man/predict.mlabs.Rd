\encoding{UTF-8}
\name{predict.mlabs}
\alias{predict.mlabs}
\title{Posterior and posterior predictive quantities from an MLABS fit}
\usage{
\method{predict}{mlabs}(
  object,
  newdata,
  type = c("mean", "function", "predictive"),
  draws = FALSE,
  level = 0.95,
  ...
)
}
\arguments{
\item{object}{An object from \strong{mlabs_collapsed}.}

\item{newdata}{A numeric matrix of new predictor values, with the same columns
in the same order as the design used for fitting. Values are supplied on the
original predictor scale.}

\item{type}{The quantity to return. \code{"mean"} gives the posterior mean of
the regression function, \code{"function"} its full posterior, and
\code{"predictive"} the posterior predictive distribution of a new
observation. See \sQuote{Details}.}

\item{draws}{A logical value used when \code{type} is \code{"function"} or
\code{"predictive"}. If \code{TRUE}, the individual draws are returned instead
of their summary.}

\item{level}{The probability content of the reported interval, used when
\code{draws = FALSE} and \code{type} is not \code{"mean"}.}

\item{...}{Additional arguments passed to or from other methods.}
}
\value{
For \code{type = "mean"}, a numeric vector with one posterior mean prediction
per row of \code{newdata}.

Otherwise, with \code{draws = FALSE}, a matrix with one row per row of
\code{newdata} and columns \code{mean}, \code{sd}, \code{lower}, and
\code{upper}, the last two being the equal-tailed interval of probability
\code{level}. With \code{draws = TRUE}, a matrix with one row per retained
posterior draw and one column per row of \code{newdata}.

All quantities are on the original response scale.
}
\description{
Evaluates the fitted MLABS model at new predictor values, returning either the
posterior mean of the regression function, its full posterior, or the posterior
predictive distribution of a new observation.
}
\details{
Each retained draw stores a complete basis configuration: the interaction
orders \eqn{K_j}, the active predictors \eqn{\nu^{(j)}}, the spline degrees
\eqn{c^{(j)}}, the knot sequences \eqn{\xi^{(j)}}, and the coefficients
\eqn{\beta_j}. For draw \eqn{s}, the corresponding tensor-product design matrix
is reconstructed at the new predictor values and the regression function is
evaluated as
\deqn{f^{(s)}(\mathbf{x}) = \beta_0 + \sum_{j=1}^{J^{(s)}} \beta_j^{(s)}
\prod_{l=1}^{K_j^{(s)}} B_{c_l^{(j,s)}}\!\left(x_{\nu_l^{(j,s)}};
\xi_l^{(j,s)}\right).}
Because quantities are formed from the function values and not from the basis
parameters, the trans-dimensional nature of the draws causes no difficulty:
basis functions need not be matched across draws.

The three values of \code{type} answer different questions.
\itemize{
\item \code{"mean"} returns the Monte Carlo average
  \eqn{S^{-1}\sum_{s} f^{(s)}(\mathbf{x})}, which estimates the posterior mean
  \eqn{E[f(\mathbf{x}) \mid \mathbf{Y}]}. It is a point estimate and carries no
  uncertainty.
\item \code{"function"} returns the posterior of \eqn{f(\mathbf{x})} itself,
  represented by the draws \eqn{f^{(1)}(\mathbf{x}), \ldots, f^{(S)}(\mathbf{x})}.
  Its interval is a pointwise credible interval for the mean surface, and it
  narrows as the sample size grows.
\item \code{"predictive"} returns the posterior predictive distribution of a
  new observation at \eqn{\mathbf{x}}, that is, of
  \eqn{Y^\star = f(\mathbf{x}) + \varepsilon} with
  \eqn{\varepsilon \sim N(0,\sigma^2)}. Each draw is obtained by adding
  observation noise from the same iteration,
  \deqn{\tilde y^{(s)} = f^{(s)}(\mathbf{x}) + \sigma^{(s)} z^{(s)},
  \qquad z^{(s)} \sim N(0,1),}
  so that the result is an exact Monte Carlo sample from
  \eqn{p(Y^\star \mid \mathbf{Y})}. Its interval is a prediction interval and
  does not shrink to zero width, since it retains the irreducible noise
  \eqn{\sigma^2}.
}


If the model was fitted with \code{normalize = TRUE}, the new predictors are
mapped with the stored training minima and ranges before evaluation, and both
the function values and the noise standard deviations are transformed back to
the original response scale. Predictor values outside the training range are
still evaluated, but they may fall outside the knot domain, where every
B-spline factor vanishes and the fitted function reverts to the intercept;
extrapolation beyond the observed range should therefore be interpreted with
care.

With \code{draws = TRUE} the returned matrix has \eqn{S \times} \code{nrow(newdata)}
entries, which can be large for long chains and many prediction points; the
default summary avoids storing it.
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

# point estimate of the regression function
fitted_values <- predict(fit, X)
plot(y, fitted_values); abline(0, 1)

# pointwise 95% credible intervals for f(x)
head(predict(fit, X, type = "function"))

# 95% prediction intervals for a new observation
pi95 <- predict(fit, X, type = "predictive")
mean(y >= pi95[, "lower"] & y <= pi95[, "upper"])   # empirical coverage

# keep the predictive sample itself
ydraws <- predict(fit, X[1:5, ], type = "predictive", draws = TRUE)
dim(ydraws)
hist(ydraws[, 1])
}

}
\references{
Park, S., Oh, J., and Lee, J. (2026+), "Multivariate Lévy adaptive B-spline
	regression", \emph{Manuscript}.
}
\seealso{
\code{mlabs_collapsed}, \code{summary.mlabs}, \code{plot_fitted}
}
\author{
Sewon Park
}
\keyword{regression}
