## Draws of the regression function at newX, on the original response scale.
## Rows are retained posterior draws, columns are rows of newX.
.function_draws <- function(fit, newX) {
  newX <- as.matrix(newX)
  if (ncol(newX) != ncol(fit$x))
    stop("newdata must have the same number of columns as the fitted design.",
         call. = FALSE)

  if (!is.null(fit$scale)) {
    newX <- sweep(newX, 2, fit$scale$x_min, "-")
    newX <- sweep(newX, 2, fit$scale$x_max - fit$scale$x_min, "/")
  }

  nsamples <- length(fit$beta)
  n_new    <- nrow(newX)
  preds    <- matrix(0, nsamples, n_new)

  for (s in seq_len(nsamples)) {
    J_s <- fit$J[s]
    if (J_s == 0) {
      preds[s, ] <- fit$beta0
      next
    }
    Bmat_s <- build_design_matrix(
      newX,
      fit$xi[[s]], fit$C[[s]], fit$nu[[s]],
      as.integer(fit$K[[s]]), J_s)
    Bmat_s[is.nan(Bmat_s)] <- 0
    preds[s, ] <- fit$beta0 + as.vector(Bmat_s %*% fit$beta[[s]])
  }

  if (!is.null(fit$scale) && !is.null(fit$scale$y_mean))
    preds <- preds * fit$scale$y_sd + fit$scale$y_mean

  preds
}

## Posterior draws of the noise standard deviation on the original scale.
.sigma_draws <- function(fit) {
  multiplier <- if (is.null(fit$scale)) 1 else fit$scale$y_sd
  as.numeric(fit$sigma) * multiplier
}

.summarise_draws <- function(draws, level) {
  alpha <- (1 - level) / 2
  quantiles <- apply(draws, 2, stats::quantile, probs = c(alpha, 1 - alpha),
                     names = FALSE)
  out <- cbind(mean  = colMeans(draws),
               sd    = apply(draws, 2, stats::sd),
               lower = quantiles[1, ],
               upper = quantiles[2, ])
  rownames(out) <- NULL
  out
}

.predict_mlabs <- function(fit, newX,
                           type = c("mean", "function", "predictive"),
                           draws = FALSE, level = 0.95) {
  .validate_fit(fit)
  type <- match.arg(type)
  if (!is.logical(draws) || length(draws) != 1L || is.na(draws))
    stop("draws must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1)
    stop("level must be a single number in (0, 1).", call. = FALSE)

  f_draws <- .function_draws(fit, newX)
  if (type == "mean") return(colMeans(f_draws))

  if (type == "predictive") {
    sigma <- .sigma_draws(fit)
    if (length(sigma) != nrow(f_draws) || any(!is.finite(sigma)))
      stop("The stored noise standard deviations are inconsistent with the draws.",
           call. = FALSE)
    f_draws <- f_draws +
      matrix(stats::rnorm(length(f_draws), 0, rep(sigma, ncol(f_draws))),
             nrow(f_draws), ncol(f_draws))
  }

  if (draws) f_draws else .summarise_draws(f_draws, level)
}

predict.mlabs <- function(object, newdata,
                          type = c("mean", "function", "predictive"),
                          draws = FALSE, level = 0.95, ...) {
  .predict_mlabs(object, newdata, type = match.arg(type),
                 draws = draws, level = level)
}

mlabs_collapsed <- function(y, X,
                            nmcmc = 200000, nburn = 100000, nthin = 10,
                            NB_max = 100,
                            maxInt = 2,
                            allowed_degrees = c(0, 1, 2),
                            aJ = 5,
                            shrink_c = 3.0,
                            exp_ratio = 0.1,
                            normalize = TRUE,
                            lambda_J = 0.0,
                            k_weights = NULL,
                            verbose = TRUE,
                            control = list()) {

  ## ---- main arguments ---------------------------------------------------
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(X)
  if (!n || nrow(X) != n)
    stop("y and X must contain the same number of observations.")
  if (any(!is.finite(y)) || any(!is.finite(X)))
    stop("y and X must contain finite values only.")
  if (!p)
    stop("X must contain at least one predictor column.")

  nmcmc  <- .check_scalar(nmcmc, "nmcmc", 1, integer = TRUE)
  nburn  <- .check_scalar(nburn, "nburn", 0, integer = TRUE)
  nthin  <- .check_scalar(nthin, "nthin", 1, integer = TRUE)
  NB_max <- .check_scalar(NB_max, "NB_max", 1, integer = TRUE)
  maxInt <- .check_scalar(maxInt, "maxInt", 1, integer = TRUE)
  if (nburn >= nmcmc)
    stop("nburn must be smaller than nmcmc.")
  if ((nmcmc - nburn) %/% nthin < 1L)
    stop("No draw would be stored: increase nmcmc - nburn or decrease nthin.")
  if (maxInt > p)
    stop("maxInt cannot exceed the number of predictors.")

  if (!is.numeric(allowed_degrees) || !length(allowed_degrees) ||
      any(!is.finite(allowed_degrees)) || any(allowed_degrees < 0) ||
      any(allowed_degrees > 3) ||
      any(allowed_degrees != floor(allowed_degrees)))
    stop("allowed_degrees must contain integers between 0 and 3.")
  allowed_degrees <- as.integer(sort(unique(allowed_degrees)))

  aJ        <- .check_scalar(aJ, "aJ", 0, strict = TRUE)
  exp_ratio <- .check_scalar(exp_ratio, "exp_ratio", 0)
  lambda_J  <- .check_scalar(lambda_J, "lambda_J", 0)
  if (!is.null(shrink_c)) shrink_c <- .check_scalar(shrink_c, "shrink_c", 0, strict = TRUE)
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize))
    stop("normalize must be TRUE or FALSE.")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("verbose must be TRUE or FALSE.")

  if (!is.null(k_weights)) {
    k_weights <- as.numeric(k_weights)
    if (length(k_weights) != maxInt || any(!is.finite(k_weights)) ||
        any(k_weights < 0) || sum(k_weights) <= 0)
      stop("k_weights must be a non-negative vector of length maxInt with a positive sum.")
    k_weights <- k_weights / sum(k_weights)
  }

  ## ---- tuning parameters ------------------------------------------------
  ctrl <- .resolve_control(control)

  bJ                 <- ctrl$bJ
  sigr               <- ctrl$sigr
  sigR               <- ctrl$sigR
  alpha_dir          <- ctrl$alpha_dir
  beta0              <- ctrl$beta0
  sigmab             <- ctrl$sigmab
  init_J             <- ctrl$init_J
  jump_prob          <- ctrl$jump_prob
  pi_warmup          <- ctrl$pi_warmup
  ref_sigma          <- ctrl$ref_sigma
  knot_temperature   <- ctrl$knot_temperature
  boost              <- ctrl$boost
  xi_step_size       <- ctrl$xi_step_size
  n_refresh          <- ctrl$n_refresh
  odd_anchor_uniform <- ctrl$odd_anchor_uniform

  ## ---- standardisation --------------------------------------------------
  scale_info <- NULL
  if (normalize) {
    x_mins <- apply(X, 2, min)
    x_maxs <- apply(X, 2, max)
    if (any(x_maxs - x_mins <= 0))
      stop("Every predictor must have a positive range when normalize = TRUE.")
    X <- sweep(X, 2, x_mins, "-")
    X <- sweep(X, 2, x_maxs - x_mins, "/")

    y_mean <- mean(y)
    y_sd   <- sd(y)
    if (!is.finite(y_sd) || y_sd <= 0)
      stop("y must have a positive standard deviation when normalize = TRUE.")
    y      <- (y - y_mean) / y_sd

    scale_info <- list(
      x_min  = x_mins,
      x_max  = x_maxs,
      y_mean = y_mean,
      y_sd   = y_sd
    )
  }

  ## ---- data-dependent defaults -----------------------------------------
  if (is.null(beta0))  beta0  <- mean(y)
  if (is.null(init_J)) {
    half_NB <- NB_max %/% 2L
    init_J <- max(10L, min(half_NB, 2L * p %/% 5L + 10L))
  } else if (init_J >= NB_max) {
    warning(sprintf(
      "control$init_J (%d) >= NB_max (%d): birth moves are blocked from the start. Increase NB_max or decrease init_J.",
      init_J, NB_max))
  }

  if (is.null(knot_temperature))
    knot_temperature <- max(1.0, sd(y) / ref_sigma)

  if (is.null(sigmab))    sigmab    <- sd(y)
  if (is.null(pi_warmup)) pi_warmup <- nburn %/% 2L

  is_adaptive_jp <- is.null(jump_prob)
  if (is_adaptive_jp) jump_prob <- c(0.4, 0.4, 0.2)

  kw <- if (is.null(k_weights)) numeric(0) else as.numeric(k_weights)

  if (verbose) {
    cat("Set beta0  =", beta0,  "\n")
    cat("Set sigmab =", sigmab,
        if (is.null(shrink_c)) "(used: shrink_c = NULL)\n"
        else "(unused: shrink_c is set)\n")
    cat("Allowed degrees S =", allowed_degrees, "\n")
    cat("Alpha (Dirichlet) =", alpha_dir, "\n")
    cat("init_J =", init_J, "(p =", p, ")\n")
    cat("NB_max =", NB_max, "\n")
    cat("nmcmc / nburn =", nmcmc, "/", nburn, "\n")
    cat("knot_temperature =", round(knot_temperature, 3), "\n")
    cat("boost =", boost, "  xi_step_size =", xi_step_size,
        "  ref_sigma =", ref_sigma, "\n")
    cat("pi_warmup =", pi_warmup, "\n")
    cat("n_refresh =", n_refresh, "\n")
  }

  ## ---- initial state ----------------------------------------------------
  boundary <- compute_boundary(X, exp_ratio)

  M <- rgamma(1, 1, 1)
  sigmay <- sd(y)

  pi_var <- rep(1 / p, p)

  J <- init_J
  K_vec <- if (is.null(k_weights)) sample(1:maxInt, J, TRUE)
  else sample(1:maxInt, J, TRUE, prob = k_weights)
  if (length(allowed_degrees) == 1) {
    allowed_degrees <- rep(allowed_degrees, maxInt)
  }

  C_list <- vector("list", J)
  nu_list <- vector("list", J)
  xi_list <- vector("list", J)

  for (j in 1:J) {
    Kj <- K_vec[j]
    nu_j <- sample(0:(p - 1), Kj, replace = FALSE)
    C_j <- sample(allowed_degrees, Kj, replace = TRUE)
    xi_j <- vector("list", Kj)
    for (k in 1:Kj) {
      x_k <- X[, nu_j[k] + 1]
      bd_k <- boundary[, nu_j[k] + 1]
      xi_j[[k]] <- generate_knots_weighted(x_k, bd_k, rep(0, n),
                                           as.integer(C_j[k]), FALSE,
                                           odd_anchor_uniform = odd_anchor_uniform)
    }
    C_list[[j]] <- as.integer(C_j)
    nu_list[[j]] <- as.integer(nu_j)
    xi_list[[j]] <- xi_j
  }

  Bmat <- build_design_matrix(X, xi_list, C_list, nu_list,
                              as.integer(K_vec), J)
  if (is.null(shrink_c)) {
    beta <- rnorm(J, 0, sigmab)
  } else {
    tau2_init <- (shrink_c * sd(y))^2
    phiM_init <- sqrt(tau2_init / max(M * exp(-lambda_J), .Machine$double.eps))
    beta <- rnorm(J, 0, phiM_init)
  }

  npost <- (nmcmc - nburn) %/% nthin
  saveIter <- 1

  betaL <- vector("list", npost)
  xiL   <- vector("list", npost)
  CL    <- vector("list", npost)
  nuL   <- vector("list", npost)
  KL    <- vector("list", npost)
  sigmaL <- numeric(npost)
  ML     <- numeric(npost)
  JL     <- integer(npost)
  MSE    <- numeric(npost)
  piL    <- matrix(0, npost, p)

  if (verbose) {
    cat("====== MLABS START:", as.character(Sys.time()), "======\n")
  }

  ## ---- reversible-jump sampler -----------------------------------------
  birth_accepts <- 0
  death_accepts <- 0
  birth_attempts <- 0
  death_attempts <- 0

  rho_cache      <- compute_relevance(X, y - beta0)
  rho_upd_burnin <- 1L
  rho_upd_main   <- 10L
  tau2 <- if (is.null(shrink_c)) NA_real_ else (shrink_c * sd(y))^2
  for (mc in 1:nmcmc) {
    rho_interval <- if (mc <= nburn %/% 2L) rho_upd_burnin else rho_upd_main

    if (is.null(shrink_c)) {
      sigmab_t <- sigmab
    } else {
      sigmab_t <- sqrt(tau2 / (M * exp(-lambda_J)))
    }
    bf_res <- backfitting_update(y, X, beta, Bmat,
                                 C_list, xi_list, nu_list,
                                 as.integer(K_vec), J,
                                 beta0, sigmay, sigmab_t,
                                 boundary, n, allowed_degrees,
                                 xi_step_size = xi_step_size,
                                 knot_temperature = knot_temperature,
                                 n_refresh = as.integer(n_refresh),
                                 odd_anchor_uniform = odd_anchor_uniform)
    beta    <- bf_res$beta
    Bmat    <- bf_res$Bmat
    C_list  <- bf_res$C
    xi_list <- bf_res$xi

    pb <- jump_prob[1]
    pd <- jump_prob[2]

    if (J >= NB_max) {
      move <- "death"
    } else if (J <= 1) {
      move <- "birth"
    } else {
      u <- runif(1)
      if (u < pb) {
        move <- "birth"
      } else if (u < pb + pd) {
        move <- "death"
      } else {
        move <- "none"
      }
    }

    if (move == "birth") {
      birth_attempts <- birth_attempts + 1
      br <- birth_step_collapsed(y, X, beta, Bmat,
                                 C_list, xi_list, nu_list,
                                 as.integer(K_vec), J, M,
                                 beta0, sigmay, sigmab_t,
                                 boundary, maxInt, n, p,
                                 jump_prob, pi_var, allowed_degrees,
                                 knot_temperature,
                                 rho_cache, as.integer(mc),
                                 as.integer(rho_interval),
                                 boost, ref_sigma, as.integer(NB_max),
                                 lambda_J,
                                 kw, odd_anchor_uniform = odd_anchor_uniform)
      beta      <- br$beta
      Bmat      <- br$Bmat
      C_list    <- br$C
      xi_list   <- br$xi
      nu_list   <- br$nu
      K_vec     <- br$K
      J         <- br$J
      rho_cache <- br$rho
      if (br$accepted) birth_accepts <- birth_accepts + 1

    } else if (move == "death") {
      death_attempts <- death_attempts + 1
      dr <- death_step_collapsed(y, X, beta, Bmat,
                                 C_list, xi_list, nu_list,
                                 as.integer(K_vec), J, M,
                                 beta0, sigmay, sigmab_t,
                                 boundary, n, p,
                                 jump_prob, pi_var, allowed_degrees,
                                 knot_temperature,
                                 rho_cache, as.integer(mc),
                                 as.integer(rho_interval),
                                 boost, ref_sigma, as.integer(NB_max),
                                 lambda_J, odd_anchor_uniform = odd_anchor_uniform)
      beta      <- dr$beta
      Bmat      <- dr$Bmat
      C_list    <- dr$C
      xi_list   <- dr$xi
      nu_list   <- dr$nu
      K_vec     <- dr$K
      J         <- dr$J
      rho_cache <- dr$rho
      if (dr$accepted) death_accepts <- death_accepts + 1
    }

    sigmay <- update_sigmay(y, beta, Bmat, beta0, sigr, sigR)
    if (is.null(shrink_c)) {
      M <- rgamma(1, shape = aJ + J, rate = bJ + exp(-lambda_J))
    } else {
      lambda_eff <- exp(-lambda_J)
      beta_ss    <- sum(beta^2)
      M <- rgamma(1,
                  shape = aJ + 1.5 * J,
                  rate  = bJ + lambda_eff * (1 + beta_ss / (2 * tau2)))
    }

    if (mc > pi_warmup) {
      pi_var <- update_pi_mh_wor(pi_var, nu_list, as.integer(K_vec), J, p, alpha_dir)
    }

    if (mc > nburn && (mc %% nthin) == 0) {
      for (jj in seq_len(J)) {
        C_jj <- C_list[[jj]]
        xi_jj <- xi_list[[jj]]
        K_jj <- K_vec[jj]
        for (kk in seq_len(K_jj)) {
          expected_nk <- C_jj[kk] + 2L
          actual_nk   <- length(xi_jj[[kk]])
          if (actual_nk != expected_nk) {
            bd_kk <- boundary[, nu_list[[jj]][kk] + 1L]
            lo <- bd_kk[1] - bd_kk[3]; hi <- bd_kk[2] + bd_kk[3]
            xi_list[[jj]][[kk]] <- sort(runif(expected_nk, lo, hi))
          }
        }
      }
      betaL[[saveIter]]  <- as.numeric(beta)
      xiL[[saveIter]]    <- lapply(
        xi_list, function(xj) lapply(xj, function(knots) knots + 0)
      )
      CL[[saveIter]]     <- lapply(C_list, function(degrees) degrees + 0L)
      nuL[[saveIter]]    <- lapply(nu_list, function(vars) vars + 0L)
      KL[[saveIter]]     <- K_vec + 0L
      sigmaL[saveIter]   <- sigmay
      ML[saveIter]       <- M
      JL[saveIter]       <- J
      piL[saveIter, ]    <- pi_var
      fit <- beta0 + Bmat %*% beta
      MSE[saveIter] <- mean((y - fit)^2)
      saveIter <- saveIter + 1
    }

    if (verbose && mc %% 1000 == 0) {
      cat("  Iter:", mc, "J:", J, "sigma:", round(sigmay, 4),
          "birth_acc:", round(birth_accepts / max(1, birth_attempts), 3),
          "death_acc:", round(death_accepts / max(1, death_attempts), 3), "\n")
    }
  }

  if (verbose) {
    cat("====== MLABS FINISH:", as.character(Sys.time()), "======\n")
  }

  ## ---- output -----------------------------------------------------------
  prior <- list(odd_anchor_uniform = odd_anchor_uniform,
                aJ = aJ, bJ = bJ, sigr = sigr, sigR = sigR,
                jump_prob = jump_prob, adaptive_jump_prob = is_adaptive_jp,
                exp_ratio = exp_ratio,
                boundary = boundary, allowed_degrees = allowed_degrees,
                sigmab = sigmab, shrink_c = shrink_c, lambda_J = lambda_J,
                NB_max = NB_max, alpha_dir = alpha_dir,
                maxInt = maxInt, k_weights = k_weights,
                ref_sigma = ref_sigma, knot_temperature = knot_temperature)

  res <- list(x = X, y = y, prior = prior, control = ctrl, beta0 = beta0,
              beta = betaL, xi = xiL, C = CL, nu = nuL, K = KL,
              sigma = sigmaL, J = JL, M = ML, pi_var = piL,
              mse = MSE, scale = scale_info, type = "regression",
              mcmc = list(nmcmc = nmcmc, nburn = nburn, nthin = nthin,
                          n_save = npost),
              acceptance = list(
                birth = birth_accepts / max(1, birth_attempts),
                death = death_accepts / max(1, death_attempts)))
  class(res) <- "mlabs"
  res
}
