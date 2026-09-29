# Deterministic in-memory data generation for the surface and Friedman studies.
# Generator version generated_v1: new random realizations, not archived CSVs.
# No mlabs installation is needed to source this file or generate data.
SIM_DATA_VERSION <- "generated_v1"
SIM_DATA_BASE_SEED <- 0413

simulation_seed <- function(group, function_id, rsnr, rep_idx) {
  if (length(group) != 1L || !group %in% c("surface", "friedman"))
    stop("group must be surface or friedman.")
  max_id <- if (group == "surface") 4L else 3L
  if (length(function_id) != 1L || !is.numeric(function_id) ||
      !is.finite(function_id) || !function_id %in% seq_len(max_id))
    stop("Invalid function_id.")
  if (length(rsnr) != 1L || !is.numeric(rsnr) || !rsnr %in% c(1, 5))
    stop("rsnr must be 1 or 5.")
  if (length(rep_idx) != 1L || !is.numeric(rep_idx) || !is.finite(rep_idx) ||
      rep_idx != floor(rep_idx) || rep_idx < 1L || rep_idx > 9999L)
    stop("rep_idx must be an integer in 1:9999.")
  as.integer(SIM_DATA_BASE_SEED + match(group, c("surface", "friedman")) *
               1000000L + function_id * 100000L +
               match(rsnr, c(1, 5)) * 10000L + rep_idx)
}

# Preserve the caller's RNG state so generating data cannot change MCMC seeds.
with_simulation_seed <- function(seed, code) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  force(code)
}

f_radial <- function(X) {
  r2 <- (X[, 1] - 0.5)^2 + (X[, 2] - 0.5)^2
  24.234 * r2 * (0.75 - r2)
}
f_mexican_hat <- function(X) {
  tau <- 0.12
  r2 <- (X[, 1] - 0.5)^2 + (X[, 2] - 0.5)^2
  (1 - r2 / tau^2) * exp(-r2 / (2 * tau^2))
}
f_genz <- function(X) {
  ifelse(X[, 1] > 0.4 | X[, 2] > 0.6, 0, exp(4 * X[, 1] + 4 * X[, 2]))
}
f_oblique <- function(X) {
  x1 <- X[, 1]; x2 <- X[, 2]
  upper <- x2 >= -0.6 * x1 + 0.75
  out <- numeric(nrow(X))
  out[upper] <- 0.2 + x1[upper]^2 + 0.1 * x2[upper]
  out[!upper] <- 0.7 + 0.01 * abs(4 * x1[!upper] + 10 * x2[!upper] - 9)^1.5
  out
}

# source/source_id are retained only for the existing MCMC seed schedule.
SURFACE_SPECS <- list(
  radial = list(paper_id=1L, label="Radial", source="A", source_id=1L,
                truth=f_radial),
  mexican_hat = list(paper_id=2L, label="Mexican hat", source="B", source_id=2L,
                     truth=f_mexican_hat),
  genz = list(paper_id=3L, label="Genz", source="B", source_id=3L, truth=f_genz),
  oblique = list(paper_id=4L, label="Oblique", source="A", source_id=3L,
                 truth=f_oblique)
)
SURFACE_KEYS <- names(SURFACE_SPECS)

load_surface_rep <- function(surface_key, rsnr, rep_idx) {
  if (length(surface_key) != 1L || !surface_key %in% SURFACE_KEYS)
    stop("Unknown surface_key.")
  spec <- SURFACE_SPECS[[surface_key]]
  seed <- simulation_seed("surface", spec$paper_id, rsnr, rep_idx)
  with_simulation_seed(seed, {
    g <- seq(0, 1, length.out = 30L)
    X_train <- as.matrix(expand.grid(x1=g, x2=g))
    X_test <- matrix(runif(2500L * 2L), nrow=2500L, ncol=2L)
    colnames(X_test) <- colnames(X_train)
    mu_train <- spec$truth(X_train)
    sigma <- sd(mu_train) / rsnr
    Y_train <- mu_train + rnorm(nrow(X_train), sd=sigma)
    list(X_train=X_train, Y_train=Y_train, X_test=X_test,
         mu_test=spec$truth(X_test), mu_train=mu_train, sigma=sigma,
         data_seed=seed, data_version=SIM_DATA_VERSION)
  })
}

friedman_truth <- function(X, fn_id) {
  if (fn_id == 1L) {
    10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
      10 * X[, 4] + 5 * X[, 5]
  } else if (fn_id == 2L) {
    sqrt(X[, 1]^2 + (X[, 2] * X[, 3] - 1 / (X[, 2] * X[, 4]))^2)
  } else if (fn_id == 3L) {
    atan((X[, 2] * X[, 3] - 1 / (X[, 2] * X[, 4])) / X[, 1])
  } else stop("fn_id must be 1, 2, or 3.")
}

load_friedman_rep <- function(fn_id, rsnr, rep_idx) {
  seed <- simulation_seed("friedman", fn_id, rsnr, rep_idx)
  with_simulation_seed(seed, {
    p <- if (fn_id == 1L) 10L else 4L
    make_X <- function(n) {
      X <- matrix(runif(n * p), nrow=n, ncol=p)
      if (fn_id != 1L) {
        X[, 1] <- 100 * X[, 1]
        X[, 2] <- 40 * pi + 520 * pi * X[, 2]
        X[, 4] <- 1 + 10 * X[, 4]
      }
      colnames(X) <- paste0("x", seq_len(p))
      X
    }
    X_train <- make_X(250L)
    X_test <- make_X(1000L)
    mu_train <- friedman_truth(X_train, fn_id)
    sigma <- sd(mu_train) / rsnr
    Y_train <- mu_train + rnorm(nrow(X_train), sd=sigma)
    list(X_train=X_train, Y_train=Y_train, X_test=X_test,
         mu_test=friedman_truth(X_test, fn_id), mu_train=mu_train, sigma=sigma,
         data_seed=seed, data_version=SIM_DATA_VERSION)
  })
}
