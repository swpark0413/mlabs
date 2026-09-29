# =============================================================================
# Reproduce the MLABS results for the four two-dimensional surface experiments
# =============================================================================
# Paper surfaces
#   1. Radial
#   2. Mexican hat
#   3. Genz
#   4. Oblique
#
# Paper design:
#   * 30 x 30 regular training grid on [0,1]^2
#   * RSNR = 1 and 5
#   * 2,500 independent test points
#   * 100 independent replications per setting
#
# This public script fits MLABS only. Hyperparameter search,
# and competing-method code are excluded.
# Selected paper settings are read from parameters.R.
#
# Data are generated in memory by simulation_data.R with fixed per-job seeds.
# No surface CSV files are required. These are new random realizations of the
# documented design, not guaranteed to match the historical data realizations.
# =============================================================================

THIS_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
})
source(file.path(THIS_DIR, "common.R"))
source_parameters(THIS_DIR)
source(file.path(THIS_DIR, "simulation_data.R"))

RSNR_VEC   <- c(1, 5)
NREP       <- 100L
RESULT_DIR <- file.path(dirname(THIS_DIR), "results", SIM_DATA_VERSION)
RAW_CSV    <- file.path(RESULT_DIR, "surface_mlabs_raw.csv")
SUM_CSV    <- file.path(RESULT_DIR, "surface_mlabs_summary.csv")
N_CORES    <- get_n_cores()

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# MLABS fitting
# -----------------------------------------------------------------------------
# Selected settings are validated against the mlabs 1.0.0 interface in common.R.

# Preserve the MCMC seed rules from the two supplied historical scripts.
# In source B the supplied grid contains one configuration, hence grid_id = 1.
surface_mcmc_seed <- function(spec, rsnr, rep_idx) {
  if (identical(spec$source, "A")) {
    return(800000L + spec$source_id * 10000L + as.integer(rsnr) * 1000L + rep_idx)
  }
  grid_id <- 1L
  820000L + spec$source_id * 100000L + as.integer(rsnr) * 50000L +
    grid_id * 200L + rep_idx
}

fit_one <- function(surface_key, rsnr, rep_idx) {
  spec <- SURFACE_SPECS[[surface_key]]
  dat <- load_surface_rep(surface_key, rsnr, rep_idx)
  params <- get_surface_params(surface_key, rsnr)
  validate_mlabs_spec(
    params,
    sprintf("Surface %s / RSNR %g", spec$label, rsnr)
  )

  set.seed(surface_mcmc_seed(spec, rsnr, rep_idx))
  t0 <- proc.time()[["elapsed"]]

  out <- tryCatch({
    fit <- fit_mlabs_spec(
      y = dat$Y_train,
      X = dat$X_train,
      spec = params,
      context = sprintf("Surface %s / RSNR %g", spec$label, rsnr)
    )

    pred <- predict_mlabs_mean(fit, dat$X_test)
    list(
      RMSE = rmse(pred, dat$mu_test),
      Mean_J = safe_mean_j(fit),
      error = NA_character_
    )
  }, error = function(e) {
    list(RMSE = NA_real_, Mean_J = NA_real_, error = conditionMessage(e))
  })

  data.frame(
    Method = "MLABS",
    surface_id = spec$paper_id,
    surface = surface_key,
    surface_label = spec$label,
    RSNR = rsnr,
    rep = rep_idx,
    RMSE = out$RMSE,
    Mean_J = out$Mean_J,
    sec = as.numeric(proc.time()[["elapsed"]] - t0),
    source_script = spec$source,
    source_surface_id = spec$source_id,
    error = out$error,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Resume-safe execution
# -----------------------------------------------------------------------------
existing <- read_results(RAW_CSV)
done <- if (is.null(existing)) {
  character(0)
} else {
  existing_key <- paste(existing$surface, existing$RSNR, existing$rep, sep = "|")
  completed_latest_keys(existing, existing_key)
}

jobs <- do.call(rbind, lapply(SURFACE_KEYS, function(surface_key) {
  do.call(rbind, lapply(RSNR_VEC, function(rsnr) {
    data.frame(
      surface_key = surface_key,
      rsnr = rsnr,
      rep_idx = seq_len(NREP),
      stringsAsFactors = FALSE
    )
  }))
}))
jobs$key <- paste(jobs$surface_key, jobs$rsnr, jobs$rep_idx, sep = "|")
jobs <- jobs[!(jobs$key %in% done), , drop = FALSE]

cat(sprintf(
  "MLABS surface reproduction: %d jobs remaining; workers=%d\n",
  nrow(jobs), N_CORES
))
cat("Surface mapping:\n")
for (surface_key in SURFACE_KEYS) {
  s <- SURFACE_SPECS[[surface_key]]
  cat(sprintf(
    "  %-12s -> paper #%d; source %s surf%d\n",
    s$label, s$paper_id, s$source, s$source_id
  ))
}

surface_error_row <- function(z, message) {
  key <- as.character(z$surface_key[[1L]])
  spec <- SURFACE_SPECS[[key]]
  data.frame(
    Method = "MLABS", surface_id = spec$paper_id, surface = key,
    surface_label = spec$label, RSNR = z$rsnr[[1L]], rep = z$rep_idx[[1L]],
    RMSE = NA_real_, Mean_J = NA_real_, sec = NA_real_,
    source_script = spec$source, source_surface_id = spec$source_id,
    error = as.character(message), stringsAsFactors = FALSE
  )
}

# Batch the expensive MCMC jobs so completed batches are saved immediately.
run_job_batches(
  jobs,
  FUN = function(z) fit_one(z$surface_key[[1L]], z$rsnr[[1L]], z$rep_idx[[1L]]),
  out_csv = RAW_CSV,
  n_cores = N_CORES,
  label = "Surface experiments",
  error_row = surface_error_row,
  progress = function(z, completed, total, elapsed) {
    cat(sprintf(
      "  %-12s RSNR=%g rep=%3d RMSE=%s Mean_J=%s sec=%s\n",
      z$surface_label[[1L]], z$RSNR[[1L]], z$rep[[1L]],
      ifelse(is.na(z$RMSE[[1L]]), "NA", sprintf("%.4f", z$RMSE[[1L]])),
      ifelse(is.na(z$Mean_J[[1L]]), "NA", sprintf("%.2f", z$Mean_J[[1L]])),
      ifelse(is.na(z$sec[[1L]]), "NA", sprintf("%.1f", z$sec[[1L]]))
    ))
  }
)

# -----------------------------------------------------------------------------
# Concise MLABS summary for reporting in the paper
# -----------------------------------------------------------------------------
raw <- read_results(RAW_CSV)
if (is.null(raw) || nrow(raw) == 0L) stop("No surface results were produced.")
raw_key <- paste(raw$surface, raw$RSNR, raw$rep, sep = "|")
raw <- latest_by_key(raw, raw_key)
raw <- raw[!is.na(raw$RMSE), , drop = FALSE]
if (nrow(raw) == 0L) stop("All surface MLABS fits failed; no summary can be produced.")

summary_df <- do.call(rbind, lapply(
  split(raw, list(raw$surface, raw$RSNR), drop = TRUE),
  function(z) data.frame(
    Method = "MLABS",
    surface_id = z$surface_id[[1L]],
    surface = z$surface[[1L]],
    surface_label = z$surface_label[[1L]],
    RSNR = z$RSNR[[1L]],
    Mean_RMSE = mean(z$RMSE),
    SD_RMSE = sd_or_na(z$RMSE),
    Mean_J = mean_or_na(z$Mean_J),
    Mean_sec = mean_or_na(z$sec),
    N = nrow(z),
    stringsAsFactors = FALSE
  )
))
summary_df <- summary_df[order(summary_df$surface_id, summary_df$RSNR), ]
rownames(summary_df) <- NULL

write.csv(summary_df, SUM_CSV, row.names = FALSE)
print(summary_df, row.names = FALSE)
cat("Saved:", SUM_CSV, "\n")
