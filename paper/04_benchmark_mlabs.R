# =============================================================================
# Reproduce the MLABS results for the six real-data benchmark datasets
# =============================================================================
# Paper design:
#   * Bodyfat, Boston housing, Concrete compressive strength,
#     Residential building, Tecator meat, and Chemical manufacturing process
#   * 20 repetitions of five-fold cross-validation
#   * The predefined fold indices stored in each .RData file are reused exactly.
#
# This public script fits MLABS only. Competing methods and tuning/search code are
# intentionally excluded. Final selected paper settings are read from parameters.R.
# =============================================================================

THIS_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
})
source(file.path(THIS_DIR, "common.R"))
source_parameters(THIS_DIR)

NREP      <- 20L
NFOLD     <- 5L
N_CORES   <- get_n_cores()
# Default: mlabs_repro/benchmark_data/, beside reproducibility/.
# To use another folder, set MLABS_BENCHMARK_DATA to its absolute path,
# or replace DATA_DIR below with e.g. "D:/MLABS/benchmark_data".
DATA_DIR <- path.expand(Sys.getenv(
  "MLABS_BENCHMARK_DATA",
  unset = file.path(dirname(THIS_DIR), "benchmark_data")
))
RESULT_DIR <- file.path(dirname(THIS_DIR), "results")
RAW_CSV   <- file.path(RESULT_DIR, "benchmark_mlabs_raw.csv")
SUM_CSV   <- file.path(RESULT_DIR, "benchmark_mlabs_summary.csv")

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

DATASETS <- list(
  bodyfat = list(
    label = "Bodyfat",
    file = "reg_testset_bodyfat.RData", obj = "bodyfat", y = "BodyFat"
  ),
  boston = list(
    label = "Boston housing",
    file = "reg_testset_boston.RData", obj = "boston", y = "medv"
  ),
  conc = list(
    label = "Concrete compressive strength",
    file = "reg_testset_conc.RData", obj = "conc", y = "y"
  ),
  resident = list(
    label = "Residential building",
    file = "reg_testset_resident.RData", obj = "resident", y = "y"
  ),
  tecator = list(
    label = "Tecator meat",
    file = "reg_testset_tecator.RData", obj = "tecator", y = "y"
  ),
  cmp = list(
    label = "Chemical manufacturing process",
    file = "reg_testset_cmp.RData", obj = "cmp", y = "Yield"
  )
)

load_dataset <- function(data_id) {
  spec <- DATASETS[[data_id]]
  path <- file.path(DATA_DIR, spec$file)
  if (!file.exists(path)) stop("Required benchmark data file not found: ", path)

  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  raw <- data.matrix(get(spec$obj, envir = e))
  cv_idx <- get(paste0(spec$obj, "_cv_idx"), envir = e)

  y_idx <- which(colnames(raw) == spec$y)
  if (length(y_idx) != 1L) {
    stop(sprintf("[%s] response column '%s' was not found uniquely.", data_id, spec$y))
  }

  y <- as.numeric(raw[, y_idx])
  X <- raw[, -y_idx, drop = FALSE]
  storage.mode(X) <- "double"
  n <- length(y)

  if (length(cv_idx) != NREP) {
    stop(sprintf("[%s] expected %d CV repetitions, found %d.",
                 data_id, NREP, length(cv_idx)))
  }
  for (r in seq_len(NREP)) {
    if (length(cv_idx[[r]]) != NFOLD) {
      stop(sprintf("[%s] CV repetition %d does not contain %d folds.",
                   data_id, r, NFOLD))
    }
    ids <- sort(as.integer(unlist(cv_idx[[r]], use.names = FALSE)))
    if (!identical(ids, seq_len(n))) {
      stop(sprintf("[%s] CV repetition %d does not partition 1:%d.", data_id, r, n))
    }
  }

  list(
    id = data_id, label = spec$label, X = X, y = y,
    cv_idx = cv_idx, n = n, p = ncol(X)
  )
}

# Report missing inputs before starting any of the long MCMC runs.
required_files <- file.path(DATA_DIR, vapply(DATASETS, `[[`, character(1), "file"))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Missing benchmark .RData files:\n  ", paste(missing_files, collapse = "\n  "),
    "\nPlace the six files in benchmark_data/ or set MLABS_BENCHMARK_DATA.",
    call. = FALSE
  )
}

get_fold <- function(dat, rep_idx, fold_idx) {
  te <- as.integer(dat$cv_idx[[rep_idx]][[fold_idx]])
  tr <- setdiff(seq_len(dat$n), te)
  list(
    X_train = dat$X[tr, , drop = FALSE],
    Y_train = dat$y[tr],
    X_test  = dat$X[te, , drop = FALSE],
    Y_test  = dat$y[te]
  )
}

# Selected settings are validated against the mlabs 1.0.0 interface in common.R.

fit_one <- function(dat, rep_idx, fold_idx) {
  params <- get_benchmark_params(dat$id)
  validate_mlabs_spec(params, paste("benchmark", dat$id))
  fold <- get_fold(dat, rep_idx, fold_idx)

  data_index <- match(dat$id, names(DATASETS))
  # The archived all-grid experiment used + grid_id * 200 in its MCMC seed.
  set.seed(820000L + data_index * 100000L + rep_idx * 7L + fold_idx)
  t0 <- proc.time()[["elapsed"]]

  out <- tryCatch({
    fit <- fit_mlabs_spec(
      y = fold$Y_train,
      X = fold$X_train,
      spec = params,
      context = paste("benchmark", dat$id)
    )
    pred <- predict_mlabs_mean(fit, fold$X_test)
    list(RMSE = rmse(pred, fold$Y_test), Mean_J = safe_mean_j(fit), error = NA_character_)
  }, error = function(e) {
    list(RMSE = NA_real_, Mean_J = NA_real_, error = conditionMessage(e))
  })

  data.frame(
    Method = "MLABS",
    data = dat$id,
    dataset = dat$label,
    rep = rep_idx,
    fold = fold_idx,
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
  existing_key <- paste(existing$data, existing$rep, existing$fold, sep = "|")
  completed_latest_keys(existing, existing_key)
}

all_jobs <- do.call(rbind, lapply(names(DATASETS), function(data_id) {
  expand.grid(
    data = data_id,
    rep = seq_len(NREP),
    fold = seq_len(NFOLD),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}))
all_jobs$key <- paste(all_jobs$data, all_jobs$rep, all_jobs$fold, sep = "|")
all_jobs <- all_jobs[!(all_jobs$key %in% done), , drop = FALSE]

cat(sprintf("MLABS benchmark reproduction: %d jobs remaining; workers=%d\n",
            nrow(all_jobs), N_CORES))

# Load each dataset once in the parent process when fitting is needed.
loaded <- NULL
if (nrow(all_jobs) > 0L) {
  loaded <- lapply(names(DATASETS), function(data_id) {
    tryCatch(load_dataset(data_id), error = identity)
  })
  names(loaded) <- names(DATASETS)
}

benchmark_error_row <- function(z, message) {
  id <- as.character(z$data[[1L]])
  data.frame(
    Method = "MLABS", data = id, dataset = DATASETS[[id]]$label,
    rep = z$rep[[1L]], fold = z$fold[[1L]], RMSE = NA_real_,
    Mean_J = NA_real_, sec = NA_real_, error = as.character(message),
    stringsAsFactors = FALSE
  )
}

run_job_batches(
  all_jobs,
  FUN = function(z) {
    dat <- loaded[[z$data[[1L]]]]
    if (inherits(dat, "error")) stop(conditionMessage(dat), call. = FALSE)
    fit_one(dat, z$rep[[1L]], z$fold[[1L]])
  },
  out_csv = RAW_CSV,
  n_cores = N_CORES,
  label = "Real-data benchmarks",
  error_row = benchmark_error_row,
  progress = function(z, completed, total, elapsed) {
    cat(sprintf(
      "  [%3d/%d] %-10s rep=%2d fold=%d RMSE=%s sec=%.1f\n",
      completed, total, z$data[[1L]], z$rep[[1L]], z$fold[[1L]],
      ifelse(is.na(z$RMSE[[1L]]), "NA", sprintf("%.4f", z$RMSE[[1L]])),
      z$sec[[1L]]
    ))
  }
)

raw <- read_results(RAW_CSV)
if (is.null(raw) || nrow(raw) == 0L) stop("No real-benchmark MLABS results were produced.")
raw_key <- paste(raw$data, raw$rep, raw$fold, sep = "|")
raw <- latest_by_key(raw, raw_key)
raw <- raw[!is.na(raw$RMSE), , drop = FALSE]
if (nrow(raw) == 0L) stop("All real-benchmark MLABS fits failed; no summary can be produced.")

summary_df <- do.call(rbind, lapply(split(raw, raw$data), function(z) {
  data.frame(
    Method = "MLABS",
    data = z$data[[1L]],
    dataset = z$dataset[[1L]],
    Mean_RMSE = mean(z$RMSE),
    SD_RMSE = sd_or_na(z$RMSE),
    Mean_J = mean_or_na(z$Mean_J),
    Mean_sec = mean_or_na(z$sec),
    N = nrow(z)
  )
}))
summary_df <- summary_df[match(intersect(names(DATASETS), summary_df$data), summary_df$data), ]
rownames(summary_df) <- NULL
write.csv(summary_df, SUM_CSV, row.names = FALSE)
print(summary_df, row.names = FALSE)
cat("Saved:", SUM_CSV, "\n")
