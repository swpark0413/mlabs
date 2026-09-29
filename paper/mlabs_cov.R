# =============================================================================
# Covariate-integrated MLABS
# =============================================================================
# Partially linear model
#   y = Xcov %*% beta_cov + f(Z) + eps
#
# This helper adapts the covariate-adjusted sampler used in the NHANES
# analysis to the internal RJMCMC functions
# =============================================================================

if (!requireNamespace("mlabs", quietly = TRUE)) {
  stop("The 'mlabs' package is required by mlabs_cov.R.")
}

ns <- asNamespace("mlabs")

mlabs_cov_collapsed <- function(y, Xcov, Z,
                                nmcmc = 6000, nburn = 3000, nthin = 5,
                                NB_max = 100, maxInt = 2,
                                allowed_degrees = c(0, 1, 2),
                                aJ = 5,
                                shrink_c = 3.0,
                                exp_ratio = 0.1,
                                normalize = TRUE,
                                lambda_J = 0.0,
                                k_weights = NULL,
                                tau_beta2 = 100,
                                verbose = FALSE,
                                control = list()) {
  y <- as.numeric(y)
  Z <- as.matrix(Z)
  Xcov <- as.matrix(Xcov)

  n <- length(y)
  if (nrow(Z) != n || nrow(Xcov) != n) {
    stop("y, Xcov, and Z must contain the same number of observations.")
  }
  if (n == 0L || ncol(Z) == 0L || ncol(Xcov) == 0L) {
    stop("y, Xcov, and Z must be non-empty.")
  }
  if (any(!is.finite(y)) || any(!is.finite(Z)) || any(!is.finite(Xcov))) {
    stop("y, Xcov, and Z must contain finite values only.")
  }

  nmcmc <- as.integer(nmcmc)
  nburn <- as.integer(nburn)
  nthin <- as.integer(nthin)
  NB_max <- as.integer(NB_max)
  maxInt <- as.integer(maxInt)
  if (nmcmc < 1L || nburn < 0L || nburn >= nmcmc || nthin < 1L) {
    stop("Require nmcmc >= 1, 0 <= nburn < nmcmc, and nthin >= 1.")
  }
  if ((nmcmc - nburn) %/% nthin < 1L) {
    stop("No posterior draw would be stored; adjust nmcmc, nburn, or nthin.")
  }
  if (NB_max < 2L || maxInt < 1L || maxInt > ncol(Z)) {
    stop("Require NB_max >= 2 and 1 <= maxInt <= ncol(Z).")
  }

  if (!is.numeric(allowed_degrees) || !length(allowed_degrees) ||
      any(!is.finite(allowed_degrees)) || any(allowed_degrees < 0) ||
      any(allowed_degrees > 3) || any(allowed_degrees != floor(allowed_degrees))) {
    stop("allowed_degrees must contain integers between 0 and 3.")
  }
  allowed_degrees <- as.integer(sort(unique(allowed_degrees)))

  if (!is.numeric(aJ) || length(aJ) != 1L || !is.finite(aJ) || aJ <= 0) {
    stop("aJ must be a positive number.")
  }
  if (!is.null(shrink_c) &&
      (!is.numeric(shrink_c) || length(shrink_c) != 1L ||
       !is.finite(shrink_c) || shrink_c <= 0)) {
    stop("shrink_c must be NULL or a positive number.")
  }
  if (!is.numeric(exp_ratio) || length(exp_ratio) != 1L ||
      !is.finite(exp_ratio) || exp_ratio <= 0) {
    stop("exp_ratio must be a finite, strictly positive number.")
  }
  if (!is.numeric(lambda_J) || length(lambda_J) != 1L ||
      !is.finite(lambda_J) || lambda_J < 0) {
    stop("lambda_J must be a non-negative number.")
  }
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) {
    stop("normalize must be TRUE or FALSE.")
  }
  if (!is.numeric(tau_beta2) || length(tau_beta2) != 1L ||
      !is.finite(tau_beta2) || tau_beta2 <= 0) {
    stop("tau_beta2 must be a positive number.")
  }

  if (!is.null(k_weights)) {
    k_weights <- as.numeric(k_weights)
    if (length(k_weights) != maxInt || any(!is.finite(k_weights)) ||
        any(k_weights < 0) || sum(k_weights) <= 0) {
      stop("k_weights must be a non-negative vector of length maxInt with positive sum.")
    }
    k_weights <- k_weights / sum(k_weights)
  }

  # Resolve sampler controls exactly as mlabs 1.0.0 does.
  ctrl <- ns$.resolve_control(control)
  if (!is.null(ctrl$beta0)) {
    stop("control$beta0 is not used here: the MLABS component intercept is fixed at 0.")
  }

  bJ <- ctrl$bJ
  sigr <- ctrl$sigr
  sigR <- ctrl$sigR
  alpha_dir <- ctrl$alpha_dir
  sigmab <- ctrl$sigmab
  init_J <- ctrl$init_J
  jump_prob <- ctrl$jump_prob
  pi_warmup <- ctrl$pi_warmup
  ref_sigma <- ctrl$ref_sigma
  knot_temperature <- ctrl$knot_temperature
  boost <- ctrl$boost
  xi_step_size <- ctrl$xi_step_size
  n_refresh <- ctrl$n_refresh
  odd_anchor_uniform <- ctrl$odd_anchor_uniform

  p <- ncol(Z)
  pc <- ncol(Xcov)
  if (qr(Xcov)$rank < pc) {
    stop("Xcov must have full column rank for the conjugate Gaussian update.")
  }

  # The intercept/linear covariates remain on their supplied scale.  Z and y are
  # normalized exactly once inside this sampler when normalize=TRUE.
  scale_info <- NULL
  z_mins <- apply(Z, 2, min)
  z_maxs <- apply(Z, 2, max)
  z_rng <- z_maxs - z_mins
  if (any(!is.finite(z_rng) | z_rng <= 0)) {
    stop("Every column of Z must have a finite positive range.")
  }
  if (normalize) {
    Z <- sweep(sweep(Z, 2, z_mins, "-"), 2, z_rng, "/")

    y_mean <- mean(y)
    y_sd <- stats::sd(y)
    if (!is.finite(y_sd) || y_sd <= 0) {
      stop("y must have a positive standard deviation when normalize=TRUE.")
    }
    y <- as.numeric((y - y_mean) / y_sd)
    scale_info <- list(
      z_min = z_mins,
      z_max = z_maxs,
      z_rng = z_rng,
      y_mean = y_mean,
      y_sd = y_sd
    )
  }

  beta0 <- 0.0
  XtX <- crossprod(Xcov)

  if (is.null(init_J)) {
    init_J <- min(NB_max, max(10L, min(NB_max %/% 2L, 2L * p %/% 5L + 10L)))
  } else if (init_J > NB_max) {
    stop(sprintf("control$init_J (%d) must not exceed NB_max (%d).",
                 init_J, NB_max), call. = FALSE)
  }

  # OLS initialization for the linear component; the spline component starts on
  # the corresponding partial residual.
  beta_cov <- as.numeric(qr.solve(Xcov, y))
  y_work <- as.numeric(y - Xcov %*% beta_cov)
  f_scale <- stats::sd(y_work)
  if (!is.finite(f_scale) || f_scale <= 0) {
    stop("The initial partial residual must have positive standard deviation.")
  }

  if (is.null(sigmab)) sigmab <- f_scale
  if (is.null(knot_temperature)) {
    knot_temperature <- max(1.0, f_scale / ref_sigma)
  }
  if (is.null(pi_warmup)) pi_warmup <- nburn %/% 2L
  if (is.null(jump_prob)) jump_prob <- c(0.4, 0.4, 0.2)

  kw <- if (is.null(k_weights)) numeric(0) else as.numeric(k_weights)
  boundary <- ns$compute_boundary(Z, exp_ratio)
  lower <- boundary[1, ] - boundary[3, ]
  upper <- boundary[2, ] + boundary[3, ]
  if (any(!is.finite(lower) | !is.finite(upper) |
          !is.finite(upper - lower) |
          lower >= boundary[1, ] | upper <= boundary[2, ])) {
    stop("exp_ratio must give finite knot bounds strictly outside the predictor range.")
  }

  degree_wt <- rep(1 / length(allowed_degrees), length(allowed_degrees))
  M <- stats::rgamma(1, 1, 1)
  sigmay <- f_scale
  pi_var <- rep(1 / p, p)
  J <- init_J
  K_vec <- if (is.null(k_weights)) {
    sample(seq_len(maxInt), J, replace = TRUE)
  } else {
    sample(seq_len(maxInt), J, replace = TRUE, prob = k_weights)
  }

  C_list <- vector("list", J)
  nu_list <- vector("list", J)
  xi_list <- vector("list", J)
  for (j in seq_len(J)) {
    Kj <- K_vec[j]
    nu_j <- sample(0:(p - 1L), Kj, replace = FALSE)
    C_j <- allowed_degrees[sample.int(length(allowed_degrees), Kj,
                                      replace = TRUE, prob = degree_wt)]
    xi_j <- vector("list", Kj)
    for (k in seq_len(Kj)) {
      x_k <- Z[, nu_j[k] + 1L]
      bd_k <- boundary[, nu_j[k] + 1L]
      xi_j[[k]] <- ns$generate_knots_weighted(
        x_k, bd_k, rep(0, n), as.integer(C_j[k]), FALSE,
        odd_anchor_uniform = odd_anchor_uniform
      )
    }
    C_list[[j]] <- as.integer(C_j)
    nu_list[[j]] <- as.integer(nu_j)
    xi_list[[j]] <- xi_j
  }

  Bmat <- ns$build_design_matrix(
    Z, xi_list, C_list, nu_list, as.integer(K_vec), J
  )

  lambda_eff <- exp(-lambda_J)
  tau2 <- if (is.null(shrink_c)) NA_real_ else (shrink_c * f_scale)^2
  if (is.null(shrink_c)) {
    beta <- stats::rnorm(J, 0, sigmab)
  } else {
    phiM_init <- sqrt(tau2 / max(M * lambda_eff, .Machine$double.eps))
    beta <- stats::rnorm(J, 0, phiM_init)
  }

  npost <- (nmcmc - nburn) %/% nthin
  saved_iterations <- nburn + seq_len(npost) * nthin
  saveIter <- 1L
  betaL <- xiL <- CL <- nuL <- KL <- vector("list", npost)
  betacovL <- matrix(NA_real_, npost, pc)
  sigmaL <- numeric(npost)
  JL <- integer(npost)
  ML <- numeric(npost)
  piL <- matrix(NA_real_, npost, p)
  mseL <- numeric(npost)

  rho_cache <- ns$compute_relevance(Z, y_work)
  b_acc <- d_acc <- b_att <- d_att <- 0L

  if (verbose) {
    cat("====== COVARIATE-INTEGRATED MLABS START:", as.character(Sys.time()), "======\n")
  }

  for (mc in seq_len(nmcmc)) {
    rho_interval <- if (mc <= nburn %/% 2L) 1L else 10L
    y_work <- as.numeric(y - Xcov %*% beta_cov)

    # Match the mlabs 1.0.0 coefficient-scale hierarchy.
    sigmab_t <- if (is.null(shrink_c)) {
      sigmab
    } else {
      sqrt(tau2 / max(M * lambda_eff, .Machine$double.eps))
    }

    bf <- ns$backfitting_update(
      y_work, Z, beta, Bmat, C_list, xi_list, nu_list,
      as.integer(K_vec), J, beta0, sigmay, sigmab_t, boundary,
      n, allowed_degrees,
      xi_step_size = xi_step_size,
      knot_temperature = knot_temperature,
      n_refresh = as.integer(n_refresh),
      odd_anchor_uniform = odd_anchor_uniform
    )
    beta <- bf$beta
    Bmat <- bf$Bmat
    C_list <- bf$C
    xi_list <- bf$xi

    if (J >= NB_max) {
      move <- "death"
    } else if (J <= 1L) {
      move <- "birth"
    } else {
      u <- stats::runif(1)
      move <- if (u < jump_prob[1]) {
        "birth"
      } else if (u < jump_prob[1] + jump_prob[2]) {
        "death"
      } else {
        "none"
      }
    }

    if (move == "birth") {
      b_att <- b_att + 1L
      br <- ns$birth_step_collapsed(
        y_work, Z, beta, Bmat, C_list, xi_list, nu_list,
        as.integer(K_vec), J, M, beta0, sigmay, sigmab_t,
        boundary, maxInt, n, p, jump_prob, pi_var, allowed_degrees,
        knot_temperature, rho_cache, as.integer(mc), as.integer(rho_interval),
        boost, ref_sigma, as.integer(NB_max), lambda_J, kw,
        odd_anchor_uniform = odd_anchor_uniform
      )
      beta <- br$beta
      Bmat <- br$Bmat
      C_list <- br$C
      xi_list <- br$xi
      nu_list <- br$nu
      K_vec <- br$K
      J <- br$J
      rho_cache <- br$rho
      if (isTRUE(br$accepted)) b_acc <- b_acc + 1L

    } else if (move == "death") {
      d_att <- d_att + 1L
      dr <- ns$death_step_collapsed(
        y_work, Z, beta, Bmat, C_list, xi_list, nu_list,
        as.integer(K_vec), J, M, beta0, sigmay, sigmab_t,
        boundary, n, p, jump_prob, pi_var, allowed_degrees,
        knot_temperature, rho_cache, as.integer(mc), as.integer(rho_interval),
        boost, ref_sigma, as.integer(NB_max), lambda_J,
        odd_anchor_uniform = odd_anchor_uniform
      )
      beta <- dr$beta
      Bmat <- dr$Bmat
      C_list <- dr$C
      xi_list <- dr$xi
      nu_list <- dr$nu
      K_vec <- dr$K
      J <- dr$J
      rho_cache <- dr$rho
      if (isTRUE(dr$accepted)) d_acc <- d_acc + 1L
    }

    # Conjugate Gaussian update for the linear covariate coefficients.
    if (J < 1L || J > NB_max) {
      stop("Internal error: basis count is outside [1, NB_max].", call. = FALSE)
    }

    f_hat <- as.numeric(Bmat %*% beta)
    precision <- XtX / sigmay^2 + diag(1 / tau_beta2, pc)
    V <- chol2inv(chol(precision))
    m <- V %*% (crossprod(Xcov, y - f_hat) / sigmay^2)
    beta_cov <- as.numeric(m + t(chol(V)) %*% stats::rnorm(pc))

    # Shared residual scale and M update.
    y_work <- as.numeric(y - Xcov %*% beta_cov)
    sigmay <- ns$update_sigmay(y_work, beta, Bmat, beta0, sigr, sigR)
    if (is.null(shrink_c)) {
      M <- stats::rgamma(1, shape = aJ + J, rate = bJ + lambda_eff)
    } else {
      beta_ss <- sum(beta^2)
      M <- stats::rgamma(
        1,
        shape = aJ + 1.5 * J,
        rate = bJ + lambda_eff * (1 + beta_ss / (2 * tau2))
      )
    }

    if (mc > pi_warmup) {
      pi_var <- ns$update_pi_mh_wor(
        pi_var, nu_list, as.integer(K_vec), J, p, alpha_dir
      )
    }

    if (mc > nburn && ((mc - nburn) %% nthin) == 0L) {
      betaL[[saveIter]] <- as.numeric(beta)
      xiL[[saveIter]] <- lapply(
        xi_list, function(xj) lapply(xj, function(knots) knots + 0)
      )
      CL[[saveIter]] <- lapply(C_list, function(degrees) degrees + 0L)
      nuL[[saveIter]] <- lapply(nu_list, function(vars) vars + 0L)
      KL[[saveIter]] <- K_vec + 0L
      betacovL[saveIter, ] <- beta_cov
      sigmaL[saveIter] <- sigmay
      JL[saveIter] <- J
      ML[saveIter] <- M
      piL[saveIter, ] <- pi_var
      mseL[saveIter] <- mean((y_work - as.numeric(Bmat %*% beta))^2)
      saveIter <- saveIter + 1L
    }

    if (verbose && mc %% 1000L == 0L) {
      cat(sprintf(
        "  it %d  J=%d sigma=%.3f b_acc=%.2f d_acc=%.2f\n",
        mc, J, sigmay, b_acc / max(1L, b_att), d_acc / max(1L, d_att)
      ))
    }
  }

  if (verbose) {
    cat("====== COVARIATE-INTEGRATED MLABS FINISH:", as.character(Sys.time()), "======\n")
  }

  stopifnot(saveIter - 1L == npost)

  structure(
    list(
      y = y,
      Z = Z,
      Xcov = Xcov,
      beta = betaL,
      xi = xiL,
      C = CL,
      nu = nuL,
      K = KL,
      beta_cov = betacovL,
      sigma = sigmaL,
      J = JL,
      M = ML,
      pi_var = piL,
      mse = mseL,
      beta0 = beta0,
      scale = scale_info,
      prior = list(
        aJ = aJ,
        bJ = bJ,
        shrink_c = shrink_c,
        exp_ratio = exp_ratio,
        lambda_J = lambda_J,
        k_weights = k_weights,
        tau_beta2 = tau_beta2,
        allowed_degrees = allowed_degrees,
        maxInt = maxInt,
        NB_max = NB_max,
        odd_anchor_uniform = odd_anchor_uniform
      ),
      control = ctrl,
      saved_iterations = saved_iterations,
      mcmc = list(
        nmcmc = nmcmc,
        nburn = nburn,
        nthin = nthin,
        n_save = npost
      ),
      acceptance = list(
        birth = b_acc / max(1L, b_att),
        death = d_acc / max(1L, d_att)
      )
    ),
    class = "mlabs_cov"
  )
}

predict_mlabs_cov <- function(fit, Xcov_new, Z_new) {
  if (!inherits(fit, "mlabs_cov")) {
    stop("fit must be an object returned by mlabs_cov_collapsed().")
  }

  Z_new <- as.matrix(Z_new)
  Xcov_new <- as.matrix(Xcov_new)
  if (nrow(Z_new) != nrow(Xcov_new)) {
    stop("Xcov_new and Z_new must contain the same number of rows.")
  }
  if (ncol(Xcov_new) != ncol(fit$beta_cov)) {
    stop("Xcov_new has a different number of columns from the fitted model.")
  }

  if (!is.null(fit$scale)) {
    if (ncol(Z_new) != length(fit$scale$z_min)) {
      stop("Z_new has a different number of columns from the fitted model.")
    }
    Z_new <- sweep(sweep(Z_new, 2, fit$scale$z_min, "-"),
                   2, fit$scale$z_rng, "/")
  }

  S <- length(fit$beta)
  if (S == 0L) stop("No posterior draws are stored in fit.")

  fmat <- matrix(0, S, nrow(Z_new))
  for (s in seq_len(S)) {
    Js <- fit$J[s]
    if (Js == 0L) next
    Bs <- ns$build_design_matrix(
      Z_new, fit$xi[[s]], fit$C[[s]], fit$nu[[s]],
      as.integer(fit$K[[s]]), Js
    )
    Bs[!is.finite(Bs)] <- 0
    fmat[s, ] <- as.numeric(Bs %*% fit$beta[[s]])
  }

  pred <- as.numeric(
    Xcov_new %*% colMeans(fit$beta_cov) + colMeans(fmat)
  )
  if (!is.null(fit$scale)) {
    pred <- pred * fit$scale$y_sd + fit$scale$y_mean
  }
  pred
}
