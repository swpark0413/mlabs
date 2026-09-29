# =============================================================================
# Reproduce the MLABS results for Friedman's benchmark functions
# =============================================================================
# Paper design:
#   * Friedman 1, 2, and 3
#   * RSNR = 1 and 5
#   * 250 training observations, 1,000 test observations (generated in memory)
#   * 100 independent replications per setting
#
# This public script fits MLABS only. Hyperparameter search and competing-method
# code are intentionally excluded. Selected paper settings are read from
# parameters.R.
# =============================================================================

THIS_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
})
source(file.path(THIS_DIR, "common.R"))
source_parameters(THIS_DIR)
source(file.path(THIS_DIR, "simulation_data.R"))

FN_IDS    <- 1:3
RSNR_VEC  <- c(1, 5)
NREP      <- 100L
RESULT_DIR <- file.path(dirname(THIS_DIR), "results", SIM_DATA_VERSION)
RAW_CSV   <- file.path(RESULT_DIR, "friedman_mlabs_raw.csv")
SUM_CSV   <- file.path(RESULT_DIR, "friedman_mlabs_summary.csv")
N_CORES   <- get_n_cores()

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

# Selected settings are validated against the mlabs 1.0.0 interface in common.R.

fit_one <- function(fn_id, rsnr, rep_idx) {
  dat <- load_friedman_rep(fn_id, rsnr, rep_idx)
  params <- get_friedman_params(fn_id, rsnr)
  validate_mlabs_spec(params, sprintf("Friedman %d / RSNR %g", fn_id, rsnr))

  set.seed(800000L + fn_id * 10000L + as.integer(rsnr) * 1000L + rep_idx)
  t0 <- proc.time()[["elapsed"]]

  out <- tryCatch({
    fit <- fit_mlabs_spec(
      y = dat$Y_train,
      X = dat$X_train,
      spec = params,
      context = sprintf("Friedman %d / RSNR %g", fn_id, rsnr)
    )
    pred <- predict_mlabs_mean(fit, dat$X_test)
    list(RMSE = rmse(pred, dat$mu_test), Mean_J = safe_mean_j(fit), error = NA_character_)
  }, error = function(e) {
    list(RMSE = NA_real_, Mean_J = NA_real_, error = conditionMessage(e))
  })

  data.frame(
    Method = "MLABS",
    function_id = fn_id,
    RSNR = rsnr,
    rep = rep_idx,
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
  existing_key <- paste(existing$function_id, existing$RSNR, existing$rep, sep = "|")
  completed_latest_keys(existing, existing_key)
}

jobs <- do.call(rbind, lapply(FN_IDS, function(fn_id) {
  do.call(rbind, lapply(RSNR_VEC, function(rsnr) {
    data.frame(fn_id = fn_id, rsnr = rsnr, rep_idx = seq_len(NREP))
  }))
}))
jobs$key <- paste(jobs$fn_id, jobs$rsnr, jobs$rep_idx, sep = "|")
jobs <- jobs[!(jobs$key %in% done), , drop = FALSE]

cat(sprintf("MLABS Friedman reproduction: %d jobs remaining; workers=%d\n",
            nrow(jobs), N_CORES))

friedman_error_row <- function(z, message) {
  data.frame(
    Method = "MLABS", function_id = z$fn_id[[1L]], RSNR = z$rsnr[[1L]],
    rep = z$rep_idx[[1L]], RMSE = NA_real_, Mean_J = NA_real_, sec = NA_real_,
    error = as.character(message), stringsAsFactors = FALSE
  )
}

run_job_batches(
  jobs,
  FUN = function(z) fit_one(z$fn_id[[1L]], z$rsnr[[1L]], z$rep_idx[[1L]]),
  out_csv = RAW_CSV,
  n_cores = N_CORES,
  label = "Friedman",
  error_row = friedman_error_row,
  progress = function(z, completed, total, elapsed) {
    cat(sprintf(
      "  [%3d/%d] fn=%d RSNR=%g rep=%3d RMSE=%s sec=%.1f\n",
      completed, total, z$function_id[[1L]], z$RSNR[[1L]], z$rep[[1L]],
      ifelse(is.na(z$RMSE[[1L]]), "NA", sprintf("%.4f", z$RMSE[[1L]])),
      z$sec[[1L]]
    ))
  }
)

raw <- read_results(RAW_CSV)
if (is.null(raw) || nrow(raw) == 0L) stop("No Friedman MLABS results were produced.")
raw_key <- paste(raw$function_id, raw$RSNR, raw$rep, sep = "|")
raw <- latest_by_key(raw, raw_key)
raw <- raw[!is.na(raw$RMSE), , drop = FALSE]
if (nrow(raw) == 0L) stop("All Friedman MLABS fits failed; no summary can be produced.")

summary_df <- do.call(rbind, lapply(
  split(raw, list(raw$function_id, raw$RSNR), drop = TRUE),
  function(z) data.frame(
    Method = "MLABS",
    function_id = z$function_id[[1L]],
    RSNR = z$RSNR[[1L]],
    Mean_RMSE = mean(z$RMSE),
    SD_RMSE = sd_or_na(z$RMSE),
    Mean_J = mean_or_na(z$Mean_J),
    Mean_sec = mean_or_na(z$sec),
    N = nrow(z)
  )
))
summary_df <- summary_df[order(summary_df$function_id, summary_df$RSNR), ]
rownames(summary_df) <- NULL
write.csv(summary_df, SUM_CSV, row.names = FALSE)
print(summary_df, row.names = FALSE)
cat("Saved:", SUM_CSV, "\n")
