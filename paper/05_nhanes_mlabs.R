# =============================================================================
# Reproduce the MLABS NHANES POP/LTL analysis
# =============================================================================
# This public script contains the MLABS analysis only. BKMR comparison code and
# hyperparameter search are intentionally excluded.
#
# Data:
#   NHANES 2001-2002 POP / leukocyte telomere length data
#   processed file: studypop.csv
#   pinned public commit for reproducibility
#
# Analysis inherited from the supplied final experiment code:
#   * n = 1,003 complete cases
#   * 18 lipid-adjusted POP exposures
#   * log-transform + standardize POPs and LTL
#   * linear adjustment for standardized age and squared standardized age
#   * one fixed five-fold CV split, saved at the subject (SEQN) level
#   * the selected MLABS settings from the supplied final NHANES code are fixed below
# =============================================================================

THIS_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
})
source(file.path(THIS_DIR, "common.R"))

K_FOLDS    <- 5L
BASE_SEED  <- 413L
RESULT_DIR <- file.path(dirname(THIS_DIR), "results")
RAW_CSV    <- file.path(RESULT_DIR, "nhanes_mlabs_folds.csv")
SUM_CSV    <- file.path(RESULT_DIR, "nhanes_mlabs_summary.csv")
FOLD_FILE  <- file.path(RESULT_DIR, "nhanes_cv_fold_assignment.csv")

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- MLABS partially linear implementation ---------------------------------
# The supplied experiment used mlabs_cov_collapsed() and predict_mlabs_cov().

load_mlabs_cov_if_needed <- function() {
  if (exists("mlabs_cov_collapsed", mode = "function", inherits = TRUE) &&
      exists("predict_mlabs_cov", mode = "function", inherits = TRUE)) {
    return(invisible(TRUE))
  }

  candidates <- c(
    Sys.getenv("MLABS_COV_SOURCE", ""),
    file.path(THIS_DIR, "mlabs_cov.R"),
    file.path(dirname(THIS_DIR), "R", "mlabs_cov.R"),
    file.path(dirname(THIS_DIR), "mlabs_cov.R")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0L) source(hit[[1L]])

  if (!exists("mlabs_cov_collapsed", mode = "function", inherits = TRUE) ||
      !exists("predict_mlabs_cov", mode = "function", inherits = TRUE)) {
    stop(
      "mlabs_cov_collapsed() and predict_mlabs_cov() are required. ",
      "The archive should contain reproducibility/mlabs_cov.R. ",
      "Alternatively set MLABS_COV_SOURCE=/path/to/mlabs_cov.R."
    )
  }
  invisible(TRUE)
}
load_mlabs_cov_if_needed()

# ---- data -------------------------------------------------------------------
# Pinned processed dataset used by the public mixture-analysis materials.

DATA_URL <- paste0(
  "https://raw.githubusercontent.com/lizzyagibson/SHARP.Mixtures.Workshop/",
  "1f2da3a14bb096d99b2c45a69d11053b0ef60088/Data/studypop.csv"
)

LOCAL_DATA <- Sys.getenv("MLABS_NHANES_DATA", "")
if (nzchar(LOCAL_DATA)) {
  if (!file.exists(LOCAL_DATA)) stop("MLABS_NHANES_DATA does not exist: ", LOCAL_DATA)
  dat0 <- read.csv(LOCAL_DATA, check.names = FALSE)
} else {
  dat0 <- read.csv(DATA_URL, check.names = FALSE)
}
if (nrow(dat0) != 1330L) {
  stop("Unexpected raw sample size: ", nrow(dat0), " (expected 1330).")
}

# Preserve the complete-case cohort used in the supplied analysis code.
dat <- na.omit(dat0)
if (nrow(dat) != 1003L) {
  stop("Unexpected complete-case sample size: ", nrow(dat), " (expected 1003).")
}

pops <- c(
  "LBX074LA", "LBX099LA", "LBX118LA", "LBX138LA", "LBX153LA", "LBX170LA",
  "LBX180LA", "LBX187LA", "LBX194LA", "LBXHXCLA", "LBXPCBLA",
  "LBXD03LA", "LBXD05LA", "LBXD07LA", "LBXF03LA", "LBXF04LA", "LBXF05LA",
  "LBXF08LA"
)

required_cols <- c("SEQN", "TELOMEAN", pops, "age_cent", "age_sq")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0L) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}
if (any(dat$TELOMEAN <= 0) || any(as.matrix(dat[, pops, drop = FALSE]) <= 0)) {
  stop("TELOMEAN and all POP concentrations must be positive before log transform.")
}

# The supplied final code standardizes these variables once on the analysis cohort.
Z <- scale(log(as.matrix(dat[, pops, drop = FALSE])))
colnames(Z) <- pops
y <- as.numeric(scale(log(dat$TELOMEAN)))

age_z <- as.numeric(scale(dat$age_cent))
X <- cbind(age = age_z, age2 = age_z^2)
Xcov <- cbind(Intercept = 1, X)

n <- length(y)
cat(sprintf(
  "NHANES MLABS reproduction: raw n=%d; analysis n=%d; exposures=%d; adjustment covariates=%d\n",
  nrow(dat0), n, ncol(Z), ncol(X)
))
cat(sprintf(
  "Data check: max |age_sq - age_cent^2| = %.3g\n",
  max(abs(dat$age_sq - dat$age_cent^2))
))

# ---- fixed five-fold split ---------------------------------------------------
# Save the subject-level fold assignment so every rerun uses exactly the same
# split rather than relying only on the RNG implementation/version.
if (file.exists(FOLD_FILE)) {
  fold_df <- read.csv(FOLD_FILE, stringsAsFactors = FALSE)
  if (!all(c("SEQN", "fold") %in% names(fold_df))) {
    stop("Invalid fold file: expected columns SEQN and fold: ", FOLD_FILE)
  }
  if (anyDuplicated(fold_df$SEQN)) stop("Duplicate SEQN values in fold file: ", FOLD_FILE)
  if (nrow(fold_df) != n) stop("Fold file row count does not match the analysis cohort.")
  idx <- match(dat$SEQN, fold_df$SEQN)
  if (anyNA(idx)) stop("Saved fold file does not cover all current SEQN values.")
  folds <- as.integer(fold_df$fold[idx])
  if (!all(folds %in% seq_len(K_FOLDS))) stop("Invalid fold labels in: ", FOLD_FILE)
} else {
  set.seed(BASE_SEED)
  folds <- sample(rep(seq_len(K_FOLDS), length.out = n))
  write.csv(
    data.frame(SEQN = dat$SEQN, fold = folds),
    FOLD_FILE, row.names = FALSE
  )
}
cat("Fold sizes:", paste(tabulate(folds, nbins = K_FOLDS), collapse = ", "), "\n")

# ---- selected NHANES parameters ---------------------------------------------
# These are the final settings from the supplied NHANES analysis

MLABS_CFG <- list(
  nmcmc = 200000L,
  nburn = 100000L,
  nthin = 10L,
  NB_max = 100L,
  maxInt = 2L,
  allowed_degrees = c(1L, 2L),
  aJ = 10,
  shrink_c = 1.0,
  exp_ratio = 1.0,
  normalize = TRUE,
  lambda_J = 0,
  k_weights = NULL,
  tau_beta2 = 100,
  control = list(
    bJ = 1,
    sigr = 0.01,
    sigR = 0.01,
    alpha_dir = 1.0,
    init_J = NULL,
    jump_prob = c(0.4, 0.4, 0.2),
    pi_warmup = NULL,
    ref_sigma = 1.0,
    knot_temperature = NULL,
    boost = 5.0,
    xi_step_size = 0.05,
    n_refresh = 5L,
    odd_anchor_uniform = 0.05
  )
)

fit_fold <- function(k) {
  tr <- folds != k
  te <- folds == k

  set.seed(BASE_SEED + 200L + k)
  t0 <- proc.time()[["elapsed"]]

  out <- tryCatch({
    fit <- mlabs_cov_collapsed(
      y = y[tr],
      Xcov = Xcov[tr, , drop = FALSE],
      Z = Z[tr, , drop = FALSE],
      nmcmc = MLABS_CFG$nmcmc,
      nburn = MLABS_CFG$nburn,
      nthin = MLABS_CFG$nthin,
      NB_max = MLABS_CFG$NB_max,
      maxInt = MLABS_CFG$maxInt,
      allowed_degrees = MLABS_CFG$allowed_degrees,
      aJ = MLABS_CFG$aJ,
      shrink_c = MLABS_CFG$shrink_c,
      exp_ratio = MLABS_CFG$exp_ratio,
      normalize = MLABS_CFG$normalize,
      lambda_J = MLABS_CFG$lambda_J,
      k_weights = MLABS_CFG$k_weights,
      tau_beta2 = MLABS_CFG$tau_beta2,
      verbose = FALSE,
      control = MLABS_CFG$control
    )
    pred <- as.numeric(predict_mlabs_cov(
      fit,
      Xcov[te, , drop = FALSE],
      Z[te, , drop = FALSE]
    ))
    list(
      RMSE = rmse(pred, y[te]),
      Mean_J = safe_mean_j(fit),
      error = NA_character_
    )
  }, error = function(e) {
    list(RMSE = NA_real_, Mean_J = NA_real_, error = conditionMessage(e))
  })

  data.frame(
    Method = "MLABS",
    fold = k,
    RMSE = out$RMSE,
    Mean_J = out$Mean_J,
    sec = as.numeric(proc.time()[["elapsed"]] - t0),
    error = out$error,
    stringsAsFactors = FALSE
  )
}

# Resume safely from already completed folds.
existing <- read_results(RAW_CSV)
done <- if (is.null(existing)) {
  integer(0)
} else {
  as.integer(completed_latest_keys(existing, as.character(existing$fold)))
}
todo <- setdiff(seq_len(K_FOLDS), done)

if (length(todo) > 0L) {
  cat("Fitting folds:", paste(todo, collapse = ", "), "\n")
  # Keep folds sequential: the supplied final code used one deterministic MCMC
  # seed per fold and the NHANES analysis has only five jobs.
  for (k in todo) {
    z <- fit_fold(k)
    append_csv(z, RAW_CSV)
    if (is.na(z$RMSE)) {
      message(sprintf("fold %d failed: %s", k, z$error))
    } else {
      message(sprintf("fold %d: RMSE=%.4f, time=%.1f sec", k, z$RMSE, z$sec))
    }
  }
}

raw <- read_results(RAW_CSV)
if (is.null(raw) || nrow(raw) == 0L) stop("No NHANES MLABS results were produced.")

# If a fold was rerun manually and appears more than once, keep the latest row.
raw <- latest_by_key(raw, raw$fold)
raw <- raw[order(raw$fold), , drop = FALSE]
valid <- raw[!is.na(raw$RMSE), , drop = FALSE]
if (nrow(valid) == 0L) stop("All NHANES MLABS folds failed; no summary can be produced.")

summary_df <- data.frame(
  Method = "MLABS",
  RMSE_mean = mean(valid$RMSE),
  RMSE_sd = sd_or_na(valid$RMSE),
  Mean_J = mean_or_na(valid$Mean_J),
  Time_mean_sec = mean_or_na(valid$sec),
  N_folds = nrow(valid),
  stringsAsFactors = FALSE
)
write.csv(summary_df, SUM_CSV, row.names = FALSE)

cat("\nPer-fold results\n")
print(raw[, c("fold", "RMSE", "Mean_J", "sec")], row.names = FALSE)
cat("\nFive-fold summary\n")
print(summary_df, row.names = FALSE)
cat("Saved:", SUM_CSV, "\n")
