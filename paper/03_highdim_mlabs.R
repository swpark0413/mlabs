# =============================================================================
# Reproduce the MLABS high-dimensional Friedman-1 experiment
# =============================================================================
# Paper design:
#   * Friedman-1 signal depends on the first five predictors
#   * p varies on an approximately logarithmic grid up to 1,000
#   * n_train = 250, n_test = 1,000
#   * noise variance = 1 or 10
#   * 10 replications per setting
#
# This public script fits MLABS only. Competing methods and tuning/search code
# are intentionally excluded.
# The dimension-dependent MLABS settings used in this experiment are encoded
# directly below because they are part of the high-dimensional design rule.
# =============================================================================

THIS_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
})
source(file.path(THIS_DIR, "common.R"))

set.seed(0413)
P_VEC     <- floor(exp(seq(log(5), log(1000), length.out = 20))) + 1L
NREP      <- 10L
SEEDS     <- sample(1:1e6, NREP)
SIGMA_VEC <- c(1, sqrt(10))
N_CORES   <- get_n_cores()
RESULT_DIR <- file.path(dirname(THIS_DIR), "results")
RAW_CSV   <- file.path(RESULT_DIR, "highdim_mlabs_raw.csv")
SUM_CSV   <- file.path(RESULT_DIR, "highdim_mlabs_summary.csv")

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

friedman1 <- function(X) {
  10 * sin(pi * X[, 1] * X[, 2]) +
    20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] + 5 * X[, 5]
}

simulate_highdim <- function(p, sigma, seed) {
  set.seed(seed)
  n_train <- 250L
  n_test <- 1000L
  X_train <- matrix(runif(n_train * p), nrow = n_train)
  X_test  <- matrix(runif(n_test * p), nrow = n_test)
  mu_test <- friedman1(X_test)
  Y_train <- friedman1(X_train) + rnorm(n_train, mean = 0, sd = sigma)
  list(X_train = X_train, Y_train = Y_train, X_test = X_test, mu_test = mu_test)
}

# Dimension-dependent setting used in the supplied final high-dimensional code.
# This is written directly in the mlabs 1.0.0 interface.
highdim_spec <- function(p, y_train) {
  dim_lvl <- if (p <= 30) 1L else if (p <= 200) 2L else 3L

  nmcmc <- if (p <= 300) 100000L else 140000L
  nburn <- nmcmc %/% 2L
  NB_max <- if (p <= 100) 100L else if (p <= 300) 200L else 300L

  list(
    nmcmc = nmcmc,
    nburn = nburn,
    nthin = 10L,
    NB_max = NB_max,
    maxInt = 2L,
    allowed_degrees = 2L,
    aJ = 10,
    shrink_c = c(3.0, 3.0, 2.0)[dim_lvl],
    exp_ratio = 1.0,
    normalize = FALSE,
    lambda_J = 0,
    k_weights = NULL,
    control = list(
      bJ = 1,
      sigr = 0.01,
      sigR = 0.01,
      alpha_dir = 5.0,
      init_J = 15L,
      jump_prob = rep(1 / 3, 3),
      # normalize=FALSE, so the historical ref_sigma_multiplier=1 maps to sd(y).
      ref_sigma = stats::sd(y_train),
      boost = 5.0,
      xi_step_size = 0.05,
      n_refresh = 5L
    )
  )
}

fit_one <- function(p, sigma, rep_idx) {
  dat <- simulate_highdim(p, sigma, SEEDS[[rep_idx]])
  params <- highdim_spec(p, dat$Y_train)
  validate_mlabs_spec(params,
                      sprintf("high-dimensional p=%d, sigma=%.6g", p, sigma))

  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch({
    fit <- fit_mlabs_spec(
      y = dat$Y_train,
      X = dat$X_train,
      spec = params,
      context = sprintf("high-dimensional p=%d, sigma=%.6g", p, sigma)
    )
    pred <- predict_mlabs_mean(fit, dat$X_test)
    list(RMSE = rmse(pred, dat$mu_test), Mean_J = safe_mean_j(fit), error = NA_character_)
  }, error = function(e) {
    list(RMSE = NA_real_, Mean_J = NA_real_, error = conditionMessage(e))
  })

  data.frame(
    Method = "MLABS",
    p = p,
    rep = rep_idx,
    sigma = sigma,
    sigma2 = sigma^2,
    RMSE = out$RMSE,
    Mean_J = out$Mean_J,
    sec = as.numeric(proc.time()[["elapsed"]] - t0),
    error = out$error,
    stringsAsFactors = FALSE
  )
}

existing <- read_results(RAW_CSV)
done <- if (is.null(existing)) {
  character(0)
} else {
  existing_key <- paste(existing$p, existing$rep,
                        sprintf("%.12g", existing$sigma), sep = "|")
  completed_latest_keys(existing, existing_key)
}

jobs <- expand.grid(
  p = P_VEC,
  rep = seq_len(NREP),
  sigma = SIGMA_VEC,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
jobs$key <- paste(jobs$p, jobs$rep, sprintf("%.12g", jobs$sigma), sep = "|")
jobs <- jobs[!(jobs$key %in% done), , drop = FALSE]

cat(sprintf("MLABS high-dimensional reproduction: %d jobs remaining; workers=%d\n",
            nrow(jobs), N_CORES))

highdim_error_row <- function(z, message) {
  sigma <- as.numeric(z$sigma[[1L]])
  data.frame(
    Method = "MLABS", p = z$p[[1L]], rep = z$rep[[1L]], sigma = sigma,
    sigma2 = sigma^2, RMSE = NA_real_, Mean_J = NA_real_, sec = NA_real_,
    error = as.character(message), stringsAsFactors = FALSE
  )
}

run_job_batches(
  jobs,
  FUN = function(z) fit_one(z$p[[1L]], z$sigma[[1L]], z$rep[[1L]]),
  out_csv = RAW_CSV,
  n_cores = N_CORES,
  label = "High-dimensional Friedman-1",
  error_row = highdim_error_row,
  progress = function(z, completed, total, elapsed) {
    cat(sprintf(
      "  [%3d/%d] p=%4d sigma2=%g rep=%2d RMSE=%s sec=%.1f\n",
      completed, total, z$p[[1L]], z$sigma2[[1L]], z$rep[[1L]],
      ifelse(is.na(z$RMSE[[1L]]), "NA", sprintf("%.4f", z$RMSE[[1L]])),
      z$sec[[1L]]
    ))
  }
)

raw <- read_results(RAW_CSV)
if (is.null(raw) || nrow(raw) == 0L) stop("No high-dimensional MLABS results were produced.")
raw_key <- paste(raw$p, raw$rep, sprintf("%.12g", raw$sigma), sep = "|")
raw <- latest_by_key(raw, raw_key)
raw <- raw[!is.na(raw$RMSE), , drop = FALSE]
if (nrow(raw) == 0L) stop("All high-dimensional MLABS fits failed; no summary can be produced.")

summary_df <- do.call(rbind, lapply(
  split(raw, list(raw$p, raw$sigma), drop = TRUE),
  function(z) data.frame(
    Method = "MLABS",
    p = z$p[[1L]],
    sigma = z$sigma[[1L]],
    sigma2 = z$sigma2[[1L]],
    Mean_RMSE = mean(z$RMSE),
    SD_RMSE = sd_or_na(z$RMSE),
    Mean_J = mean_or_na(z$Mean_J),
    Mean_sec = mean_or_na(z$sec),
    N = nrow(z)
  )
))
summary_df <- summary_df[order(summary_df$sigma2, summary_df$p), ]
rownames(summary_df) <- NULL
write.csv(summary_df, SUM_CSV, row.names = FALSE)
print(summary_df, row.names = FALSE)
cat("Saved:", SUM_CSV, "\n")
