# ## Column names used inside ggplot2 aes() calls; declared so that R CMD check
# ## does not report them as undefined global variables.
# utils::globalVariables(c("count", "degree_label", "important", "label",
#                          "mean_count", "observed", "order_label", "proportion",
#                          "std_resid", "type", "value", "variable"))

.check_scalar <- function(value, name, lower = -Inf, upper = Inf,
                          integer = FALSE, strict = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value > upper || (if (strict) value <= lower else value < lower) ||
      (integer && value != floor(value)))
    stop(sprintf("Invalid value for %s.", name), call. = FALSE)
  invisible(value)
}

.validate_fit <- function(fit) {
  fields <- c("x", "y", "beta0", "beta", "xi", "C", "nu", "K", "J", "sigma")
  if (!inherits(fit, "mlabs") || !all(fields %in% names(fit)))
    stop("fit must be an mlabs object returned by mlabs_collapsed().")
  ns <- length(fit$J)
  if (!ns || any(!is.finite(fit$J)) || any(fit$J < 0) || any(fit$J != floor(fit$J)) ||
      any(vapply(fit[c("beta", "xi", "C", "nu", "K", "sigma")], length, integer(1)) != ns))
    stop("The stored posterior draws are empty or inconsistent.")
  invisible(fit)
}

.variable_names <- function(fit, X_train = fit$x, var_names) {
  .validate_fit(fit)
  if (is.null(dim(X_train)) || ncol(X_train) != ncol(fit$x))
    stop("X_train must have the same number of columns as fit$x.")
  if (is.null(var_names)) var_names <- colnames(X_train)
  if (is.null(var_names)) var_names <- paste0("x", seq_len(ncol(fit$x)))
  if (length(var_names) != ncol(fit$x) || anyNA(var_names) || anyDuplicated(var_names))
    stop("var_names must contain one distinct, non-missing name per predictor.")
  as.character(var_names)
}

.training_X <- function(fit) {
  if (is.null(fit$scale)) return(fit$x)
  sweep(sweep(fit$x, 2, fit$scale$x_max - fit$scale$x_min, "*"),
        2, fit$scale$x_min, "+")
}

.training_y <- function(fit) {
  if (is.null(fit$scale)) fit$y else fit$y * fit$scale$y_sd + fit$scale$y_mean
}

.fitted_values <- function(fit) .predict_mlabs(fit, .training_X(fit))

.check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly=TRUE))
    stop(sprintf("Package '%s' is required. install.packages('%s')", pkg, pkg),
         call.=FALSE)
}

.get_sigma2_chain <- function(fit) {
  .validate_fit(fit)
  multiplier <- if (is.null(fit$scale)) 1 else fit$scale$y_sd
  as.numeric(fit$sigma * multiplier)^2
}

.get_loglik_chain <- function(fit) {
  ll <- fit$loglik
  if (is.null(ll) && !is.null(fit$mse)) {
    n <- nrow(fit$x)
    ll <- -n / 2 * (log(2 * pi * fit$sigma^2) + fit$mse / fit$sigma^2)
  }
  if (!is.null(ll) && !is.null(fit$scale)) ll <- ll - nrow(fit$x) * log(fit$scale$y_sd)
  if (is.null(ll)) NULL else as.numeric(ll)
}

.get_mse_chain <- function(fit) {
  if (is.null(fit$mse)) return(NULL)
  multiplier <- if (is.null(fit$scale)) 1 else fit$scale$y_sd^2
  as.numeric(fit$mse) * multiplier
}

.get_function_draws <- function(fit, observation_index) {
  n_save <- length(fit$J)
  out <- matrix(NA_real_, n_save, length(observation_index))
  X_eval <- fit$x[observation_index, , drop = FALSE]
  for (s in seq_len(n_save)) {
    value <- rep(fit$beta0, length(observation_index))
    if (fit$J[s] > 0L) {
      B <- build_design_matrix(X_eval, fit$xi[[s]], fit$C[[s]], fit$nu[[s]],
                               as.integer(fit$K[[s]]), fit$J[s])
      value <- value + as.numeric(B %*% fit$beta[[s]])
    }
    if (!is.null(fit$scale))
      value <- value * fit$scale$y_sd + fit$scale$y_mean
    out[s, ] <- value
  }
  colnames(out) <- paste0("f[", observation_index, "]")
  out
}

summary.mlabs <- function(object,
                           probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                           include_pi = TRUE,
                           include_fitted = TRUE,
                           fitted_index = NULL, ...) {
  .validate_fit(object)
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1) || anyDuplicated(probs))
    stop("probs must contain distinct finite probabilities between zero and one.")
  probs <- sort(as.numeric(probs))
  if (!is.logical(include_pi) || length(include_pi) != 1L || is.na(include_pi))
    stop("include_pi must be TRUE or FALSE.")
  if (!is.logical(include_fitted) || length(include_fitted) != 1L || is.na(include_fitted))
    stop("include_fitted must be TRUE or FALSE.")
  if (is.null(fitted_index)) fitted_index <- seq_len(nrow(object$x))
  if (!is.numeric(fitted_index) || any(!is.finite(fitted_index)) ||
      any(fitted_index != floor(fitted_index)) ||
      any(fitted_index < 1L | fitted_index > nrow(object$x)) ||
      anyDuplicated(fitted_index))
    stop("fitted_index must contain distinct valid observation indices.")
  fitted_index <- as.integer(fitted_index)

  n_save <- length(object$J)
  chains <- list(J = as.numeric(object$J))
  if (!is.null(object$M)) chains$M <- as.numeric(object$M)
  chains[["sigma^2"]] <- .get_sigma2_chain(object)
  mse <- .get_mse_chain(object)
  if (!is.null(mse)) chains$MSE <- mse
  loglik <- .get_loglik_chain(object)
  if (!is.null(loglik)) chains$logLik <- loglik

  if (include_pi && !is.null(object$pi_var)) {
    pi_var <- as.matrix(object$pi_var)
    if (!is.numeric(pi_var) || nrow(pi_var) != n_save ||
        ncol(pi_var) != ncol(object$x) || any(!is.finite(pi_var)))
      stop("The stored pi_var draws are inconsistent with the fitted model.")
    variable_names <- .variable_names(object, object$x, NULL)
    for (j in seq_len(ncol(pi_var)))
      chains[[paste0("pi[", variable_names[j], "]")]] <- pi_var[, j]
  }

  lengths <- vapply(chains, length, integer(1))
  if (any(lengths != n_save) || any(!vapply(chains, function(z) all(is.finite(z)), logical(1))))
    stop("The stored posterior sample lengths or values are inconsistent.")
  draws <- do.call(cbind, chains)
  fitted_names <- character(0)
  if (include_fitted && length(fitted_index)) {
    fitted_draws <- .get_function_draws(object, fitted_index)
    fitted_names <- colnames(fitted_draws)
    draws <- cbind(draws, fitted_draws)
  }
  statistics <- cbind(
    Mean = colMeans(draws),
    Variance = apply(draws, 2, stats::var)
  )
  quantiles <- t(apply(draws, 2, stats::quantile, probs = probs, names = FALSE))
  colnames(quantiles) <- paste0(format(100 * probs, trim = TRUE, scientific = FALSE), "%")

  structure(
    list(statistics = statistics,
         quantiles = quantiles,
         n_draws = n_save,
         saved_iterations = object$saved_iterations %||% seq_len(n_save),
         include_pi = include_pi,
         include_fitted = include_fitted,
         fitted_index = if (include_fitted) fitted_index else integer(0),
         fitted_names = fitted_names,
         probs = probs),
    class = "summary.mlabs"
  )
}

print.summary.mlabs <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 max_fitted = 20L, ...) {
  .check_scalar(digits, "digits", 1, 15, integer = TRUE)
  .check_scalar(max_fitted, "max_fitted", 0, integer = TRUE)
  global_names <- setdiff(rownames(x$statistics), x$fitted_names)
  cat("\nPosterior sample summary\n")
  cat("Saved draws:", x$n_draws, "\n\n")
  cat("1. Empirical mean and variance\n\n")
  print(x$statistics[global_names, , drop = FALSE], digits = digits)
  cat("\n2. Quantiles\n\n")
  print(x$quantiles[global_names, , drop = FALSE], digits = digits)
  if (length(x$fitted_names) && max_fitted > 0L) {
    shown <- head(x$fitted_names, max_fitted)
    fitted_table <- cbind(x$statistics[shown, , drop = FALSE],
                          x$quantiles[shown, , drop = FALSE])
    cat("\n3. Regression-function values at observed predictors\n\n")
    print(fitted_table, digits = digits)
    if (length(shown) < length(x$fitted_names))
      cat("\nDisplayed", length(shown), "of", length(x$fitted_names),
          "function values. Use print(x, max_fitted = n) to display more.\n")
  }
  cat("\nThe trans-dimensional basis coefficients, degrees, and knots are not\n",
      "included because their identities can change across saved draws.\n", sep = "")
  invisible(x)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

marginal_var_prob <- function(fit, X_train = fit$x,
                                    var_names = NULL,
                                    threshold  = 0.5) {
  p <- ncol(X_train)
  var_names <- .variable_names(fit, X_train, var_names)

  parsed  <- .parse_basis_samples(fit)

  if (is.null(parsed))
    stop("Stored basis information is unavailable in fit.")

  n_save <- parsed$n_save

  t1_included <- matrix(0L, nrow=n_save, ncol=p)

  t2_basis_count <- matrix(0L, nrow=n_save, ncol=p)
  t2_total       <- integer(n_save)

  t3_basis_count <- matrix(0L, nrow=n_save, ncol=p)
  t3_total       <- integer(n_save)

  .iter_bases(parsed, function(t, xi, deg) {
    K_j <- length(xi)
    t2_total[t] <<- t2_total[t] + 1L

    valid_xi <- xi[xi >= 1L & xi <= p]
    if (length(valid_xi) == 0) return()

    t1_included[t, valid_xi] <<- 1L

    t2_basis_count[t, valid_xi] <<- t2_basis_count[t, valid_xi] + 1L

    if (K_j >= 2L) {
      t3_total[t] <<- t3_total[t] + 1L
      t3_basis_count[t, valid_xi] <<- t3_basis_count[t, valid_xi] + 1L
    }
  })

  mip <- colMeans(t1_included)

  abr_per_iter <- t2_basis_count / pmax(t2_total, 1L)
  abr <- colMeans(abr_per_iter)

  icr_per_iter <- t3_basis_count / pmax(t3_total, 1L)
  icr_per_iter[t3_total == 0, ] <- NA_real_
  icr <- colMeans(icr_per_iter, na.rm=TRUE)
  icr[is.nan(icr)] <- NA_real_

  result_df <- data.frame(
    variable  = var_names,

    MIP       = round(mip, 4),

    ABR       = round(abr, 4),

    ICR       = round(icr, 4),

    important = mip >= threshold,
    stringsAsFactors = FALSE
  )
  result_df <- result_df[order(-result_df$MIP), ]

  make_post_summary <- function(mat) {

    do.call(rbind, lapply(seq_len(p), function(j) {
      v <- mat[, j]
      v <- v[!is.na(v)]
      if (length(v) == 0)
        return(data.frame(variable=var_names[j],
                          mean=NA, sd=NA, q025=NA, median=NA, q975=NA))
      data.frame(
        variable = var_names[j],
        mean     = round(mean(v), 4),
        sd       = round(sd(v),   4),
        q025     = round(quantile(v, 0.025), 4),
        median   = round(quantile(v, 0.500), 4),
        q975     = round(quantile(v, 0.975), 4),
        stringsAsFactors = FALSE
      )
    }))
  }

  post_MIP <- make_post_summary(t1_included)
  post_ABR <- make_post_summary(abr_per_iter)
  post_ICR <- make_post_summary(icr_per_iter)

  structure(
    list(
      summary   = result_df,
      post_MIP  = post_MIP,
      post_ABR  = post_ABR,
      post_ICR  = post_ICR,

      samples_MIP = t1_included,
      samples_ABR = abr_per_iter,
      samples_ICR = icr_per_iter,

      var_names = var_names,
      threshold = threshold,
      n_save    = n_save,
      p         = p
    ),
    class = "marginal_prob"
  )
}

print.marginal_prob <- function(x, ...) {
  cat("==============================================================\n")
  cat(" MLABS Marginal Variable Inclusion Probability\n")
  cat(sprintf(" p = %d  |  Saved MCMC draws = %d  |  threshold = %.2f\n",
              x$p, x$n_save, x$threshold))
  cat("==============================================================\n")
  cat(sprintf(" %-12s  %6s  %6s  %6s  %s\n",
              "Variable", "MIP", "ABR", "ICR", "Above threshold"))
  cat(strrep("-", 58), "\n")
  df <- x$summary
  for (i in seq_len(nrow(df))) {
    cat(sprintf(" %-12s  %6.4f  %6.4f  %6s  %s\n",
                df$variable[i],
                df$MIP[i],
                df$ABR[i],
                ifelse(is.na(df$ICR[i]), "  NA  ", sprintf("%6.4f", df$ICR[i])),
                ifelse(df$important[i], "*", "")))
  }
  cat(strrep("-", 58), "\n")
  cat(" MIP: Marginal Inclusion Probability  P(xj in model | data)\n")
  cat(" ABR: Average Basis Ratio             E[bases containing j / J]\n")
  cat(" ICR: Interaction-cond. Ratio         E[bases containing j / interaction bases | interaction bases exist]\n\n")
  invisible(x)
}

plot_marginal_var_prob <- function(mp,
                                   threshold  = NULL,
                                   title      = "Marginal Variable Inclusion Probability",
                                   save_path  = NULL,
                                   width      = 1600,
                                   height     = 1200) {
  .check_pkg("ggplot2"); .check_pkg("patchwork")

  thr   <- threshold %||% mp$threshold
  BLUE  <- "#2166ac"; ORANGE <- "#d95f02"; GREEN <- "#1b7837"
  GRAY  <- "#636363"

  df    <- mp$summary

  var_lvl <- df$variable[order(df$MIP)]

  make_bar <- function(col, fill_col, x_lab, sub_title) {
    d <- data.frame(variable = df$variable,
                    value    = df[[col]],
                    important= df$important)
    d$label <- ifelse(is.na(d$value), "NA", sprintf("%.3f", d$value))
    d$value[is.na(d$value)] <- 0
    d$variable <- factor(d$variable, levels=var_lvl)
    ggplot2::ggplot(d, ggplot2::aes(x=variable, y=value, fill=important)) +
      ggplot2::geom_col(width=0.7, alpha=0.85) +
      ggplot2::geom_text(ggplot2::aes(
        label=label),
        hjust=-0.15, size=3.0) +
      ggplot2::geom_hline(yintercept=thr, linetype="dashed",
                          color="red", linewidth=0.7) +
      ggplot2::scale_fill_manual(
        values=c(`TRUE`=fill_col, `FALSE`="#d0d0d0"), guide="none") +
      ggplot2::scale_y_continuous(limits=c(0, 1.05)) +
      ggplot2::coord_flip() +
      ggplot2::labs(title=sub_title, x=NULL, y=x_lab) +
      ggplot2::theme_bw(base_size=10) +
      ggplot2::theme(axis.text.y=ggplot2::element_text(size=8))
  }

  make_box <- function(samples_mat, fill_col, x_lab, sub_title) {

    long <- do.call(rbind, lapply(seq_len(mp$p), function(j) {
      v <- samples_mat[, j]
      v <- v[!is.na(v)]
      if (length(v) == 0) return(NULL)
      data.frame(variable=mp$var_names[j], value=v,
                 stringsAsFactors=FALSE)
    }))
    if (is.null(long)) return(ggplot2::ggplot() + ggplot2::theme_void())
    long$variable <- factor(long$variable, levels=var_lvl)

    ggplot2::ggplot(long, ggplot2::aes(x=variable, y=value)) +
      ggplot2::geom_boxplot(fill=fill_col, alpha=0.6,
                            outlier.size=0.5, outlier.alpha=0.3,
                            width=0.6) +
      ggplot2::geom_hline(yintercept=thr, linetype="dashed",
                          color="red", linewidth=0.7) +
      ggplot2::scale_y_continuous(limits=c(0, 1)) +
      ggplot2::coord_flip() +
      ggplot2::labs(title=sub_title, x=NULL, y=x_lab) +
      ggplot2::theme_bw(base_size=10) +
      ggplot2::theme(axis.text.y=ggplot2::element_text(size=8))
  }

  p_mip_bar <- make_bar("MIP", BLUE,
                        "P(xj in model | data)",
                        "MIP - Posterior mean")
  p_mip_box <- make_box(mp$samples_MIP, BLUE,
                        "Inclusion per draw (0/1)",
                        "MIP - Distribution across draws")

  p_abr_bar <- make_bar("ABR", ORANGE,
                        "E[bases containing j / J]",
                        "ABR - Posterior mean")
  p_abr_box <- make_box(mp$samples_ABR, ORANGE,
                        "Basis fraction per draw",
                        "ABR - Distribution across draws")

  p_icr_bar <- make_bar("ICR", GREEN,
                        "E[bases containing j / interaction bases | interaction bases exist]",
                        "ICR - Posterior mean")
  p_icr_box <- make_box(mp$samples_ICR, GREEN,
                        "ICR per eligible draw",
                        "ICR - Distribution across draws")

  final <- (p_mip_bar | p_mip_box) /
    (p_abr_bar | p_abr_box) /
    (p_icr_bar | p_icr_box) +
    patchwork::plot_annotation(
      title    = title,
      subtitle = sprintf(
        "threshold = %.2f  |  n_save = %d  |  p = %d\n%s",
        thr, mp$n_save, mp$p,
        paste(
          "MIP: P(xj in model | data)",
          "ABR: E[basis fraction for j]",
          "ICR: E[basis fraction for j | interaction bases exist]",
          sep="   "
        )
      ),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size=13, face="bold"),
        plot.subtitle = ggplot2::element_text(size=9,  color=GRAY)
      )
    )

  if (!is.null(save_path)) {
    png(save_path, width=width, height=height, res=130)
    print(final)
    dev.off()
    cat(sprintf("[Saved] %s\n", save_path))
  } else {
    print(final)
  }
  invisible(mp)
}

.parse_basis_samples <- function(fit) {
  .validate_fit(fit)
  list(data = fit, n_save = length(fit$J))
}

.iter_bases <- function(parsed, callback) {
  fit <- parsed$data
  for (t in seq_len(parsed$n_save)) {
    for (j in seq_len(fit$J[t])) {
      callback(t, as.integer(fit$nu[[t]][[j]]) + 1L,
               as.integer(fit$C[[t]][[j]]))
    }
  }
  invisible(NULL)
}

basis_analysis <- function(fit, X_train = fit$x, var_names=NULL, min_Kj=2L) {
  .check_scalar(min_Kj, "min_Kj", 1, integer = TRUE)
  p <- ncol(X_train)
  var_names <- .variable_names(fit, X_train, var_names)

  parsed <- .parse_basis_samples(fit)
  if (is.null(parsed))
    stop("Stored basis information is unavailable in fit.\n",
         "Check the stored basis fields.")

  n_save <- parsed$n_save

  pi_j_samples <- numeric(n_save)

  main_count   <- integer(p)

  inter_count  <- matrix(0L, p, p,
                         dimnames=list(var_names, var_names))

  degree_tally <- list()

  order_tally <- list()

  inter_Kj_count <- matrix(0L, p, p,
                           dimnames=list(var_names, var_names))

  iter_Kj_count  <- integer(n_save)
  iter_total     <- integer(n_save)

  .iter_bases(parsed, function(t, xi, deg) {
    K_j <- length(xi)
    iter_total[t] <<- iter_total[t] + 1L

    if (K_j >= min_Kj) {
      iter_Kj_count[t] <<- iter_Kj_count[t] + 1L

      if (K_j >= 2L) {
        pairs <- combn(xi, 2L)
        for (z in seq_len(ncol(pairs))) {
          i1 <- pairs[1L, z]; i2 <- pairs[2L, z]
          inter_Kj_count[i1, i2] <<- inter_Kj_count[i1, i2] + 1L
          inter_Kj_count[i2, i1] <<- inter_Kj_count[i2, i1] + 1L
        }
      }
    }

    if (K_j == 1L && xi >= 1L && xi <= p) {
      main_count[xi] <<- main_count[xi] + 1L
    }

    if (K_j >= 2L) {
      pairs <- combn(xi, 2L)
      for (z in seq_len(ncol(pairs))) {
        i1 <- pairs[1L, z]; i2 <- pairs[2L, z]
        inter_count[i1, i2] <<- inter_count[i1, i2] + 1L
        inter_count[i2, i1] <<- inter_count[i2, i1] + 1L
      }
    }

    for (dk in as.character(deg))
      degree_tally[[dk]] <<- (degree_tally[[dk]] %||% 0L) + 1L

    ok <- as.character(K_j)
    order_tally[[ok]] <<- (order_tally[[ok]] %||% 0L) + 1L
  })

  for (t in seq_len(n_save))
    pi_j_samples[t] <- if (iter_total[t] > 0)
      iter_Kj_count[t] / iter_total[t] else 0

  pi_summary <- list(
    samples = pi_j_samples,
    mean    = mean(pi_j_samples),
    sd      = sd(pi_j_samples),
    median  = median(pi_j_samples),
    q025    = quantile(pi_j_samples, 0.025),
    q975    = quantile(pi_j_samples, 0.975)
  )

  main_df <- data.frame(
    variable   = var_names,
    mean_count = main_count / n_save,
    total_count= main_count,
    stringsAsFactors = FALSE
  )
  main_df <- main_df[order(-main_df$mean_count), ]

  pairs <- which(upper.tri(inter_count), arr.ind = TRUE)
  inter_df <- data.frame(
    var1 = var_names[pairs[, 1]], var2 = var_names[pairs[, 2]],
    label = paste(var_names[pairs[, 1]], var_names[pairs[, 2]], sep = " x "),
    mean_count_all = inter_count[pairs] / n_save,
    mean_count_Kj = inter_Kj_count[pairs] / n_save,
    total_count_all = inter_count[pairs], total_count_Kj = inter_Kj_count[pairs])
  inter_df <- inter_df[order(-inter_df$mean_count_all), , drop = FALSE]

  degree_df <- data.frame(
    degree      = as.integer(names(degree_tally)),
    total_count = as.integer(unlist(degree_tally)),
    stringsAsFactors = FALSE
  )
  degree_df$proportion <- degree_df$total_count / sum(degree_df$total_count)
  degree_df <- degree_df[order(degree_df$degree), ]

  order_df <- data.frame(
    order       = as.integer(names(order_tally)),
    total_count = as.integer(unlist(order_tally)),
    stringsAsFactors = FALSE
  )
  order_df$proportion <- order_df$total_count / sum(order_df$total_count)
  order_df <- order_df[order(order_df$order), ]

  Kj_count_df <- data.frame(
    iter        = seq_len(n_save),
    Kj_count    = iter_Kj_count,
    total_bases = iter_total,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      p          = p,
      var_names  = var_names,
      min_Kj     = min_Kj,
      n_save     = n_save,
      pi_j       = pi_summary,
      main_effect= main_df,
      interaction= inter_df,
      degree_dist= degree_df,
      order_dist = order_df,
      Kj_count   = Kj_count_df,
      inter_mat_all = inter_count / n_save,
      inter_mat_Kj  = inter_Kj_count / n_save
    ),
    class = "basis_analysis"
  )
}

print.basis_analysis <- function(x, top_n=10, digits=4, ...) {
  cat("======================================================\n")
  cat(sprintf(" MLABS Basis Analysis  (K_j >= %d basis threshold)\n", x$min_Kj))
  cat(sprintf(" p = %d  |  Saved MCMC draws = %d\n", x$p, x$n_save))
  cat("======================================================\n")

  pi <- x$pi_j
  cat(sprintf("\n-- basis_ratio (K_j >= %d fraction of bases) Posterior distribution --------------\n",
              x$min_Kj))
  cat(sprintf("  Posterior mean   = %.4f\n", pi$mean))
  cat(sprintf("  Posterior median = %.4f\n", pi$median))
  cat(sprintf("  Posterior SD     = %.4f\n", pi$sd))
  cat(sprintf("  95%% CI      = [%.4f, %.4f]\n", pi$q025, pi$q975))

  cat("\n-- Main-effect counts (mean bases per draw) ------------------\n")
  me <- head(x$main_effect, top_n)
  for (i in seq_len(nrow(me)))
    cat(sprintf("  %-12s  %.3f  %s\n", me$variable[i], me$mean_count[i],
                paste(rep("*", min(20, round(me$mean_count[i]*5))),
                      collapse="")))

  if (!is.null(x$interaction) && nrow(x$interaction) > 0) {
    cat("\n-- Top pair co-occurrence counts (mean bases per draw) -------------\n")
    top_inter <- head(x$interaction[order(-x$interaction$mean_count_all), ],
                      min(top_n, nrow(x$interaction)))
    for (i in seq_len(nrow(top_inter)))
      cat(sprintf("  %-20s  all=%.3f  selected=%.3f\n",
                  top_inter$label[i],
                  top_inter$mean_count_all[i],
                  top_inter$mean_count_Kj[i]))
  }

  cat("\n-- Factor degree distribution -------------------------------\n")
  for (i in seq_len(nrow(x$degree_dist)))
    cat(sprintf("  degree %d : %5.1f%%  (total %d)\n",
                x$degree_dist$degree[i],
                x$degree_dist$proportion[i]*100,
                x$degree_dist$total_count[i]))
  cat("\n")
  invisible(x)
}

plot_basis_analysis <- function(ba, top_n=15,
                                title="MLABS Basis Analysis",
                                save_path=NULL,
                                width=1600, height=1200) {
  .check_pkg("ggplot2"); .check_pkg("patchwork"); .check_pkg("scales")

  BLUE  <- "#2166ac"; ORANGE <- "#d95f02"; GRAY <- "#636363"

  od <- ba$order_dist
  od$order_label <- if (nrow(od)) paste0("K = ", od$order) else character(0)
  p1 <- ggplot2::ggplot(od, ggplot2::aes(x=order_label, y=proportion)) +
    ggplot2::geom_col(width=0.6, fill=BLUE, alpha=0.85) +
    ggplot2::geom_text(ggplot2::aes(label=sprintf("%.1f%%", proportion*100)),
                       vjust=-0.4, size=3.5) +
    ggplot2::scale_y_continuous(labels=function(z) paste0(round(100*z), "%"),
                                limits=c(0, max(od$proportion)*1.18)) +
    ggplot2::labs(title="Interaction-order distribution",
                  subtitle="Fraction of all bases at each interaction order",
                  x=expression(K[j]), y="Fraction") +
    ggplot2::theme_bw(base_size=11)

  me <- ba$main_effect
  me$variable <- factor(me$variable,
                        levels=me$variable[order(me$mean_count)])
  p2 <- ggplot2::ggplot(me,
                        ggplot2::aes(x=variable, y=mean_count)) +
    ggplot2::geom_col(fill=BLUE, alpha=0.85, width=0.7) +
    ggplot2::geom_text(ggplot2::aes(label=sprintf("%.2f", mean_count)),
                       hjust=-0.15, size=3.2) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title="Main-effect basis counts",
      subtitle="Mean number of univariate bases per draw",
      x="Predictor", y="Mean bases per draw") +
    ggplot2::theme_bw(base_size=11)

  if (!is.null(ba$interaction) && nrow(ba$interaction) > 0) {
    inter_top <- head(ba$interaction[order(-ba$interaction$mean_count_all), ],
                      top_n)
    inter_long <- rbind(
      data.frame(label=inter_top$label,
                 count=inter_top$mean_count_Kj,
                 type="Selected order", stringsAsFactors=FALSE),
      data.frame(label=inter_top$label,
                 count=inter_top$mean_count_all - inter_top$mean_count_Kj,
                 type="Other orders", stringsAsFactors=FALSE)
    )
    inter_long$label <- factor(inter_long$label,
                               levels=inter_top$label[
                                 order(inter_top$mean_count_all)])
    p3 <- ggplot2::ggplot(inter_long,
                          ggplot2::aes(x=label, y=count, fill=type)) +
      ggplot2::geom_col(width=0.7, alpha=0.85) +
      ggplot2::scale_fill_manual(
        values=c("Selected order"=BLUE, "Other orders"="#92c5de"),
        name="") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title=sprintf("Top %d pair co-occurrence counts", top_n),
        subtitle="Mean bases containing each pair per draw",
        x="Predictor pair", y="Mean bases per draw") +
      ggplot2::theme_bw(base_size=11) +
      ggplot2::theme(legend.position="bottom")
  } else {
    p3 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x=0.5, y=0.5,
                        label="No interaction bases\n(check the selected interaction order)",
                        size=5) +
      ggplot2::theme_void()
  }

  dd <- ba$degree_dist
  dd$degree_label <- if (nrow(dd)) paste0("degree ", dd$degree) else character(0)
  p4a <- ggplot2::ggplot(dd,
                         ggplot2::aes(x=degree_label, y=proportion,
                                      fill=degree_label)) +
    ggplot2::geom_col(width=0.6, alpha=0.85, show.legend=FALSE) +
    ggplot2::geom_text(ggplot2::aes(label=sprintf("%.1f%%",proportion*100)),
                       vjust=-0.4, size=3.5) +
    ggplot2::scale_fill_brewer(palette="Set2") +
    ggplot2::scale_y_continuous(labels=function(x) paste0(round(100*x), "%"),
                                limits=c(0, max(dd$proportion)*1.18)) +
    ggplot2::labs(title="Factor degree distribution",
                  x="Degree", y="Fraction") +
    ggplot2::theme_bw(base_size=11)

  final <- (p1 | p2) / (p3 | p4a) +
    patchwork::plot_annotation(
      title   = title,
      subtitle= sprintf("n_save=%d  |  p=%d  |  K_j>=%d basis analysis",
                        ba$n_save, ba$p, ba$min_Kj),
      theme   = ggplot2::theme(
        plot.title   = ggplot2::element_text(size=14, face="bold"),
        plot.subtitle= ggplot2::element_text(size=10, color=GRAY)
      )
    )

  if (!is.null(save_path)) {
    png(save_path, width=width, height=height, res=130)
    print(final)
    dev.off()
    cat(sprintf("[Saved] %s\n", save_path))
  } else {
    print(final)
  }
  invisible(ba)
}

## Chains that plot_acf() can display, assembled in the order given by pars.
.acf_chains <- function(fit, pars, which_pi, fitted_index) {
  chains <- list()
  for (par in pars) {
    if (par == "sigma2") {
      chains[["sigma^2"]] <- .get_sigma2_chain(fit)

    } else if (par == "M") {
      chains[["M"]] <- as.numeric(fit$M)

    } else if (par == "pi") {
      pi_var <- as.matrix(fit$pi_var)
      p <- ncol(pi_var)
      variable_names <- .variable_names(fit, fit$x, NULL)
      if (is.null(which_pi)) {
        which_pi <- order(colMeans(pi_var), decreasing = TRUE)[seq_len(min(p, 3L))]
      } else if (is.character(which_pi)) {
        which_pi <- match(which_pi, variable_names)
        if (anyNA(which_pi))
          stop("which_pi must name predictors of the fitted design.", call. = FALSE)
      } else {
        which_pi <- as.integer(which_pi)
        if (any(!is.finite(which_pi)) || any(which_pi < 1L) || any(which_pi > p))
          stop("which_pi must contain valid predictor indices.", call. = FALSE)
      }
      for (j in which_pi)
        chains[[paste0("pi[", variable_names[j], "]")]] <- pi_var[, j]

    } else if (par == "f") {
      n <- nrow(fit$x)
      if (is.null(fitted_index)) {
        y_train <- .training_y(fit)
        fitted_index <- vapply(stats::quantile(y_train, c(0.25, 0.5, 0.75)),
                               function(q) which.min(abs(y_train - q)), integer(1))
        fitted_index <- unique(as.integer(fitted_index))
      } else {
        fitted_index <- as.integer(fitted_index)
        if (any(!is.finite(fitted_index)) || any(fitted_index < 1L) ||
            any(fitted_index > n) || anyDuplicated(fitted_index))
          stop("fitted_index must contain distinct valid observation indices.",
               call. = FALSE)
      }
      f_draws <- .get_function_draws(fit, fitted_index)
      for (k in seq_along(fitted_index))
        chains[[colnames(f_draws)[k]]] <- f_draws[, k]

    } else if (par == "J") {
      chains[["J"]] <- as.numeric(fit$J)

    } else if (par == "logLik") {
      chains[["log-likelihood"]] <- .get_loglik_chain(fit)

    } else if (par == "mse") {
      chains[["MSE"]] <- .get_mse_chain(fit)
    }
  }
  chains
}

plot_acf <- function(fit, pars = c("sigma2", "M", "pi"),
                     which_pi = NULL, fitted_index = NULL,
                     lag_max = 50L, title = "MLABS Autocorrelation Diagnostics") {
  .validate_fit(fit)
  .check_pkg("ggplot2")
  .check_pkg("patchwork")
  .check_scalar(lag_max, "lag_max", 1, integer = TRUE)

  allowed <- c("sigma2", "M", "pi", "f", "J", "logLik", "mse")
  if (!is.character(pars) || !length(pars) || anyDuplicated(pars) ||
      !all(pars %in% allowed))
    stop("pars must be distinct values from: ",
         paste(allowed, collapse = ", "), ".", call. = FALSE)

  chains <- .acf_chains(fit, pars, which_pi, fitted_index)
  if (!length(chains))
    stop("No chain is available for the requested pars.", call. = FALSE)

  make_panel <- function(values, label) {
    if (is.null(values) || length(values) < 2L ||
        any(!is.finite(values)) || stats::sd(values) == 0) {
      return(
        ggplot2::ggplot() +
          ggplot2::annotate(
            "text", x = 0, y = 0,
            label = "ACF is unavailable for this chain."
          ) +
          ggplot2::labs(title = label) +
          ggplot2::theme_void()
      )
    }

    maximum_lag <- min(as.integer(lag_max), length(values) - 1L)
    acf_result <- stats::acf(values, lag.max = maximum_lag, plot = FALSE)
    data <- data.frame(
      lag = as.numeric(acf_result$lag[, , 1]),
      acf = as.numeric(acf_result$acf[, , 1])
    )
    bound <- stats::qnorm(0.975) / sqrt(length(values))

    ggplot2::ggplot(data, ggplot2::aes(x = lag, y = acf)) +
      ggplot2::geom_col(fill = "#2166ac", width = 0.6, alpha = 0.85) +
      ggplot2::geom_hline(
        yintercept = c(-bound, bound), linetype = "dashed",
        color = "red", linewidth = 0.6
      ) +
      ggplot2::labs(title = label, x = "Lag", y = "ACF") +
      ggplot2::theme_bw(base_size = 11)
  }

  plots <- Map(make_panel, chains, names(chains))
  patchwork::wrap_plots(plots, ncol = min(3L, length(plots))) +
    patchwork::plot_annotation(title = title)
}

diagnosis <- function(fit, title = "MLABS MCMC Trace Diagnostics",
                      smooth = TRUE, span = 0.25) {
  .validate_fit(fit)
  .check_pkg("ggplot2")
  .check_pkg("patchwork")
  if (!is.logical(smooth) || length(smooth) != 1L || is.na(smooth))
    stop("smooth must be TRUE or FALSE.")
  .check_scalar(span, "span", 0, 1, strict = TRUE)
  n_save <- length(fit$J)
  iteration <- fit$saved_iterations %||% seq_len(n_save)
  if (length(iteration) != n_save || any(!is.finite(iteration)))
    stop("The stored saved_iterations are inconsistent with the posterior draws.")
  loglik <- .get_loglik_chain(fit)
  if (is.null(loglik) || length(loglik) != n_save || any(!is.finite(loglik)))
    stop("Finite log-likelihood draws are required for diagnosis().")
  sigma2 <- .get_sigma2_chain(fit)

  make_trace <- function(value, y_label, color) {
    data <- data.frame(iteration = iteration, value = value)
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = iteration, y = value)) +
      ggplot2::geom_line(color = color, linewidth = 0.4, alpha = 0.4) +
      ggplot2::labs(x = "MCMC iteration", y = y_label) +
      ggplot2::theme_bw(base_size = 11)
    if (smooth && n_save >= 3L)
      plot <- plot + ggplot2::geom_smooth(
        method = "loess", formula = y ~ x, se = FALSE,
        span = max(span, min(1, 5 / n_save)),
        method.args = list(degree = 1),
        color = color, linewidth = 1.0
      )
    plot
  }

  plots <- list(
    make_trace(as.numeric(fit$J), "Number of basis terms (J)", "#2166ac"),
    make_trace(loglik, "Log-likelihood", "#1b7837"),
    make_trace(sigma2, expression(sigma^2), "#d95f02")
  )
  patchwork::wrap_plots(plots, ncol = 1) +
    patchwork::plot_annotation(title = title)
}

plot_fitted <- function(fit, Y_train = .training_y(fit), title="Fitted vs Observed") {
  .check_pkg("ggplot2")

  fitted <- .fitted_values(fit)
  if (length(Y_train) != length(fitted) || any(!is.finite(Y_train)))
    stop("Y_train must contain one finite response per training observation.")
  if (is.null(fitted))
    stop("Training predictions are unavailable in fit.")

  df <- data.frame(observed=Y_train, fitted=as.numeric(fitted))
  lim <- range(c(df$observed, df$fitted), na.rm=TRUE)

  ggplot2::ggplot(df, ggplot2::aes(x=observed, y=fitted)) +
    ggplot2::geom_point(alpha=0.5, size=1.8, color="#2166ac") +
    ggplot2::geom_abline(slope=1, intercept=0,
                         linetype="dashed", color="red", linewidth=0.8) +
    ggplot2::coord_fixed(xlim=lim, ylim=lim) +
    ggplot2::labs(title=title, x="Observed", y="Fitted") +
    ggplot2::theme_bw(base_size=12)
}

plot_residuals <- function(fit, Y_train = .training_y(fit), title="Residual Diagnostics") {
  .check_pkg("ggplot2"); .check_pkg("patchwork")

  fitted <- .fitted_values(fit)
  if (length(Y_train) != length(fitted) || any(!is.finite(Y_train)))
    stop("Y_train must contain one finite response per training observation.")
  if (is.null(fitted))
    stop("Training predictions are unavailable in fit.")

  resid <- Y_train - fitted
  df    <- data.frame(fitted=fitted, resid=resid,
                      std_resid=if (sd(resid) > 0) (resid - mean(resid)) / sd(resid) else rep(0, length(resid)))

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x=fitted, y=resid)) +
    ggplot2::geom_point(alpha=0.5, size=1.8, color="#2166ac") +
    ggplot2::geom_hline(yintercept=0, linetype="dashed",
                        color="red", linewidth=0.7) +
    ggplot2::geom_smooth(method="loess", se=FALSE, color="orange",
                         linewidth=0.8, formula=y~x) +
    ggplot2::labs(title="Residual vs Fitted", x="Fitted", y="Residual") +
    ggplot2::theme_bw(base_size=11)

  p2 <- ggplot2::ggplot(df, ggplot2::aes(sample=std_resid)) +
    ggplot2::stat_qq(alpha=0.5, size=1.5, color="#2166ac") +
    ggplot2::stat_qq_line(color="red", linetype="dashed", linewidth=0.7) +
    ggplot2::labs(title="Normal Q-Q", x="Theoretical quantile", y="Standardized residual") +
    ggplot2::theme_bw(base_size=11)

  p3 <- ggplot2::ggplot(df, ggplot2::aes(x=resid)) +
    ggplot2::geom_histogram(fill="#2166ac", alpha=0.7, bins=30,
                            color="white") +
    ggplot2::labs(title="Residual Histogram", x="Residual", y="Frequency") +
    ggplot2::theme_bw(base_size=11)

  (p1 | p2 | p3) +
    patchwork::plot_annotation(title=title)
}
