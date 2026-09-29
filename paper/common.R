# Common utilities for the MLABS reproducibility scripts.

if (!requireNamespace("mlabs", quietly = TRUE)) {
  stop(
    "The 'mlabs' package is required. Install mlabs version 1.0.0 before ",
    "running the reproducibility scripts."
  )
}

pkg_version <- as.character(utils::packageVersion("mlabs"))
if (pkg_version != "1.0.0") {
  warning(
    "These scripts were written for mlabs 1.0.0; installed version is ",
    pkg_version, ". Results may differ if the interface or sampler changed."
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))))
  }
  getwd()
}

source_parameters <- function(this_dir) {
  path <- file.path(this_dir, "parameters.R")
  if (!file.exists(path)) {
    stop(
      "Missing reproducibility/parameters.R. Copy parameters_TEMPLATE.R to ",
      "parameters.R and fill in the selected settings used in the paper."
    )
  }
  source(path, local = parent.frame())
  invisible(path)
}

rmse <- function(pred, truth) {
  pred <- as.numeric(pred)
  truth <- as.numeric(truth)
  if (length(pred) != length(truth)) {
    stop("Prediction and truth vectors have different lengths.")
  }
  sqrt(mean((pred - truth)^2))
}

safe_mean_j <- function(fit) {
  if (!is.null(fit$J)) return(mean(as.numeric(fit$J), na.rm = TRUE))
  NA_real_
}

mean_or_na <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

sd_or_na <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

get_n_cores <- function() {
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (length(detected) != 1L || is.na(detected) || detected < 1L) detected <- 1L
  detected <- max(1L, as.integer(detected) - 1L)

  requested <- suppressWarnings(as.integer(Sys.getenv("MLABS_N_CORES", detected)))
  if (length(requested) != 1L || is.na(requested) || requested < 1L) requested <- 1L
  min(detected, requested)
}

parallel_map <- function(X, FUN, n_cores = get_n_cores()) {
  if (length(X) == 0L) return(list())
  if (n_cores <= 1L || .Platform$OS.type == "windows") {
    if (.Platform$OS.type == "windows" && n_cores > 1L) {
      message("Windows detected: using sequential execution for portability.")
    }
    return(lapply(X, FUN))
  }
  parallel::mclapply(
    X, FUN,
    mc.cores = n_cores,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
}

append_csv <- function(x, path) {
  if (!is.data.frame(x) || nrow(x) < 1L || ncol(x) < 1L ||
      anyDuplicated(names(x)) || any(!nzchar(names(x)))) {
    stop("CSV results must be a non-empty data frame with unique named columns.")
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  first_write <- !file.exists(path) || file.info(path)$size == 0L
  if (first_write) {
    write.csv(x, path, row.names = FALSE)
  } else {
    expected <- names(read.csv(path, nrows = 0L, check.names = FALSE))
    if (!identical(names(x), expected)) {
      stop(
        "Refusing to append a result with a different CSV schema. Expected: ",
        paste(expected, collapse = ", "), "; received: ",
        paste(names(x), collapse = ", "), "."
      )
    }
    write.table(x, path, sep = ",", row.names = FALSE,
                col.names = FALSE, append = TRUE)
  }
  invisible(path)
}

run_job_batches <- function(jobs, FUN, out_csv, n_cores = get_n_cores(),
                            label = "MLABS", progress = NULL,
                            error_row) {
  if (is.null(jobs) || nrow(jobs) == 0L) return(invisible(NULL))
  if (!is.function(error_row)) stop("error_row must be a function.")

  make_error_row <- function(job, message) {
    row <- error_row(job, as.character(message)[[1L]])
    if (!is.data.frame(row) || nrow(row) != 1L) {
      stop("error_row must return a one-row data frame.")
    }
    row
  }
  template <- make_error_row(jobs[1L, , drop = FALSE], "result template")
  expected_names <- names(template)
  safe_fun <- function(job) {
    tryCatch(FUN(job), error = function(e) make_error_row(job, conditionMessage(e)))
  }
  normalize_result <- function(result, job) {
    if (inherits(result, "try-error")) {
      return(make_error_row(job, as.character(result)))
    }
    if (!is.data.frame(result) || nrow(result) != 1L ||
        !identical(names(result), expected_names)) {
      return(make_error_row(
        job,
        paste0("Worker returned an invalid result schema; expected: ",
               paste(expected_names, collapse = ", "))
      ))
    }
    result
  }

  batch_size <- max(1L, as.integer(n_cores))
  starts <- seq.int(1L, nrow(jobs), by = batch_size)
  completed <- 0L
  t_start <- proc.time()[["elapsed"]]

  for (st in starts) {
    en <- min(st + batch_size - 1L, nrow(jobs))
    batch <- split(jobs[st:en, , drop = FALSE], seq_len(en - st + 1L))
    results <- parallel_map(batch, safe_fun, n_cores = n_cores)
    if (length(results) != length(batch)) {
      stop("Parallel execution returned a different number of results than jobs.")
    }
    results <- Map(normalize_result, results, batch)

    for (z in results) {
      append_csv(z, out_csv)
      completed <- completed + 1L
      if (is.function(progress)) {
        progress(z, completed, nrow(jobs), proc.time()[["elapsed"]] - t_start)
      }
    }
  }

  message(sprintf("%s: completed %d job(s).", label, completed))
  invisible(NULL)
}

read_results <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0L) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

latest_by_key <- function(df, key) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  if (length(key) != nrow(df)) stop("Key length does not match number of rows.")
  df[!duplicated(key, fromLast = TRUE), , drop = FALSE]
}

completed_latest_keys <- function(df, key, success_col = "RMSE") {
  if (is.null(df) || nrow(df) == 0L) return(character(0))
  if (length(key) != nrow(df)) stop("Key length does not match number of rows.")
  if (!is.character(success_col) || length(success_col) != 1L ||
      !(success_col %in% names(df))) {
    stop("success_col must name a result column.")
  }
  keep <- !duplicated(key, fromLast = TRUE)
  latest <- df[keep, , drop = FALSE]
  latest_key <- as.character(key[keep])
  value <- suppressWarnings(as.numeric(latest[[success_col]]))
  latest_key[is.finite(value)]
}

# Validate one selected setting written in the same structure as mlabs 1.0.0:
#   top-level entries -> formal arguments of mlabs_collapsed()
#   control           -> named list accepted by mlabs_control()
validate_mlabs_spec <- function(spec, context = "MLABS") {
  if (!is.list(spec) || is.null(names(spec))) {
    stop(context, " settings must be a named list.")
  }

  required <- c(
    "nmcmc", "nburn", "nthin", "NB_max", "maxInt", "allowed_degrees",
    "aJ", "shrink_c", "exp_ratio", "normalize", "lambda_J", "k_weights",
    "control"
  )
  missing <- setdiff(required, names(spec))
  if (length(missing)) {
    stop(context, " settings are incomplete. Missing: ",
         paste(missing, collapse = ", "))
  }

  top_allowed <- c(
    "nmcmc", "nburn", "nthin", "NB_max", "maxInt", "allowed_degrees",
    "aJ", "shrink_c", "exp_ratio", "normalize", "lambda_J", "k_weights",
    "control"
  )
  unknown_top <- setdiff(names(spec), top_allowed)
  if (length(unknown_top)) {
    stop(context, " contains unknown top-level setting(s): ",
         paste(unknown_top, collapse = ", "))
  }

  if (!is.list(spec$control) || (length(spec$control) && is.null(names(spec$control)))) {
    stop(context, "$control must be a named list (possibly empty).")
  }
  control_allowed <- names(formals(mlabs::mlabs_control))
  unknown_control <- setdiff(names(spec$control), control_allowed)
  if (length(unknown_control)) {
    stop(context, "$control contains unknown setting(s): ",
         paste(unknown_control, collapse = ", "))
  }

  # Let the package perform detailed validation; these checks give earlier,
  # clearer messages for common mistakes in the paper settings file.
  if (!is.numeric(spec$nmcmc) || length(spec$nmcmc) != 1L ||
      !is.numeric(spec$nburn) || length(spec$nburn) != 1L ||
      spec$nmcmc <= 0 || spec$nburn < 0 || spec$nburn >= spec$nmcmc) {
    stop(context, ": require 0 <= nburn < nmcmc.")
  }
  if (!is.numeric(spec$nthin) || length(spec$nthin) != 1L || spec$nthin < 1) {
    stop(context, ": nthin must be >= 1.")
  }
  if (!is.logical(spec$normalize) || length(spec$normalize) != 1L || is.na(spec$normalize)) {
    stop(context, ": normalize must be TRUE or FALSE.")
  }
  if (!is.numeric(spec$allowed_degrees) || !length(spec$allowed_degrees) ||
      any(spec$allowed_degrees < 0) || any(spec$allowed_degrees > 3) ||
      any(spec$allowed_degrees != floor(spec$allowed_degrees))) {
    stop(context, ": allowed_degrees must contain integers from 0 through 3.")
  }

  # Important mlabs 1.0.0 semantic: control$sigmab is used only when
  # shrink_c = NULL. Warn rather than silently pretending both matter.
  if (!is.null(spec$shrink_c) && "sigmab" %in% names(spec$control)) {
    warning(context, ": control$sigmab is ignored by mlabs 1.0.0 when shrink_c is non-NULL.")
  }

  invisible(TRUE)
}

fit_mlabs_spec <- function(y, X, spec, context = "MLABS") {
  validate_mlabs_spec(spec, context)
  mlabs::mlabs_collapsed(
    y = y,
    X = X,
    nmcmc = spec$nmcmc,
    nburn = spec$nburn,
    nthin = spec$nthin,
    NB_max = spec$NB_max,
    maxInt = spec$maxInt,
    allowed_degrees = spec$allowed_degrees,
    aJ = spec$aJ,
    shrink_c = spec$shrink_c,
    exp_ratio = spec$exp_ratio,
    normalize = spec$normalize,
    lambda_J = spec$lambda_J,
    k_weights = spec$k_weights,
    verbose = FALSE,
    control = spec$control
  )
}

predict_mlabs_mean <- function(fit, newX) {
  as.numeric(stats::predict(fit, newdata = newX, type = "mean"))
}
