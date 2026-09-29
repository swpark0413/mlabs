# =============================================================================
# MLABS settings for reproducibility
# =============================================================================
# Configurations for the Surface, Friedman, and real-data benchmark experiments.
# Sampler controls are passed through control=.
# ref_sigma is specified on the standardized response scale when normalize=TRUE.
# Bodyfat uses normalize=FALSE.
# When shrink_c is non-NULL, control$sigmab is ignored.
# odd_anchor_uniform is fixed at 0.05.
# =============================================================================

.make_spec <- function(nmcmc = 200000L,
                       nburn = 100000L,
                       nthin = 10L,
                       NB_max = 100L,
                       maxInt = 2L,
                       allowed_degrees,
                       aJ,
                       shrink_c,
                       exp_ratio,
                       normalize = TRUE,
                       lambda_J = 0,
                       k_weights = NULL,
                       bJ = 1,
                       sigr = 0.01,
                       sigR = 0.01,
                       alpha_dir = 5,
                       init_J = 15L,
                       jump_prob = c(0.4, 0.45, 0.15),
                       ref_sigma = 1,
                       boost = 5,
                       xi_step_size = 0.03,
                       n_refresh = 5L,
                       odd_anchor_uniform = 0.05) {
  list(
    nmcmc = as.integer(nmcmc),
    nburn = as.integer(nburn),
    nthin = as.integer(nthin),
    NB_max = as.integer(NB_max),
    maxInt = as.integer(maxInt),
    allowed_degrees = as.integer(allowed_degrees),
    aJ = aJ,
    shrink_c = shrink_c,
    exp_ratio = exp_ratio,
    normalize = normalize,
    lambda_J = lambda_J,
    k_weights = if (is.null(k_weights)) NULL else as.numeric(k_weights),
    control = list(
      bJ = bJ,
      sigr = sigr,
      sigR = sigR,
      alpha_dir = alpha_dir,
      init_J = as.integer(init_J),
      jump_prob = as.numeric(jump_prob),
      ref_sigma = ref_sigma,
      boost = boost,
      xi_step_size = xi_step_size,
      n_refresh = as.integer(n_refresh),
      odd_anchor_uniform = odd_anchor_uniform
    )
  )
}

# ----------------------------------------------------------------------------
# Four two-dimensional surfaces
# ----------------------------------------------------------------------------
get_surface_params <- function(surface, rsnr) {
  surface <- as.character(surface)
  rsnr <- as.numeric(rsnr)

  # Radial: archived Surface 1
  if (identical(surface, "radial") && rsnr == 1) {
    return(.make_spec(
      NB_max = 100L, maxInt = 2L, allowed_degrees = 2L,
      aJ = 5, shrink_c = 5, exp_ratio = 0.1,
      k_weights = c(0.9, 0.1)
    ))
  }
  if (identical(surface, "radial") && rsnr == 5) {
    return(.make_spec(
      NB_max = 100L, maxInt = 2L, allowed_degrees = 2L,
      aJ = 10, shrink_c = 3, exp_ratio = 0.1,
      k_weights = c(0.5, 0.5)
    ))
  }

  # Mexican hat: archived Surface 2
  if (identical(surface, "mexican_hat") && rsnr %in% c(1, 5)) {
    return(.make_spec(
      NB_max = 100L, maxInt = 2L, allowed_degrees = c(0L, 1L, 2L),
      aJ = 5, shrink_c = 3, exp_ratio = 0.1,
      k_weights = c(0, 1)
    ))
  }

  # Genz discontinuous: archived Surface 3 from source B
  if (identical(surface, "genz") && rsnr %in% c(1, 5)) {
    return(.make_spec(
      NB_max = 100L, maxInt = 2L, allowed_degrees = c(0L, 1L, 2L),
      aJ = 5, shrink_c = 3, exp_ratio = 0.1,
      k_weights = c(0, 1)
    ))
  }

  # Oblique discontinuous: archived Surface 3 from source A
  if (identical(surface, "oblique") && rsnr == 1) {
    return(.make_spec(
      NB_max = 100L, maxInt = 2L, allowed_degrees = c(0L, 1L),
      aJ = 10, shrink_c = 3, exp_ratio = 0.1,
      k_weights = c(0.3, 0.7)
    ))
  }
  if (identical(surface, "oblique") && rsnr == 5) {
    # User-specified override: degrees=(0,2), NB_max=50, aJ=5.
    # Remaining entries follow the archived oblique RSNR=5 report.
    return(.make_spec(
      NB_max = 50L, maxInt = 2L, allowed_degrees = c(0L, 2L),
      aJ = 5, shrink_c = 5, exp_ratio = 0.1,
      k_weights = c(0.3, 0.7)
    ))
  }

  stop("No selected surface setting for surface=", surface, ", RSNR=", rsnr)
}

# ----------------------------------------------------------------------------
# Friedman benchmark functions
# ----------------------------------------------------------------------------
get_friedman_params <- function(fn_id, rsnr) {
  fn_id <- as.integer(fn_id)
  rsnr <- as.numeric(rsnr)

  # Friedman 1
  if (fn_id == 1L && rsnr %in% c(1, 5)) {
    return(.make_spec(
      NB_max = 100L, maxInt = 2L, allowed_degrees = 2L,
      aJ = 10, shrink_c = 5, exp_ratio = 1,
      k_weights = NULL
    ))
  }

  # Friedman 2
  if (fn_id == 2L && rsnr == 1) {
    return(.make_spec(
      NB_max = 100L, maxInt = 2L, allowed_degrees = 2L,
      aJ = 10, shrink_c = 5, exp_ratio = 2,
      k_weights = NULL
    ))
  }
  if (fn_id == 2L && rsnr == 5) {
    return(.make_spec(
      NB_max = 100L, maxInt = 3L, allowed_degrees = 2L,
      aJ = 10, shrink_c = 5, exp_ratio = 2,
      k_weights = NULL
    ))
  }

  # Friedman 3
  if (fn_id == 3L && rsnr %in% c(1, 5)) {
    return(.make_spec(
      NB_max = 100L, maxInt = 3L, allowed_degrees = 2L,
      aJ = 10, shrink_c = 3, exp_ratio = 1,
      k_weights = NULL
    ))
  }

  stop("No selected Friedman setting for function=", fn_id, ", RSNR=", rsnr)
}

# ----------------------------------------------------------------------------
# Six real-data benchmark datasets
# ----------------------------------------------------------------------------
# Use birth/death moves with jump_prob=c(0.50, 0.50, 0).
# Set control$knot_temperature=1.
# control$sigmab is unused when shrink_c is non-NULL.
# ----------------------------------------------------------------------------

get_benchmark_params <- function(data_id) {
  data_id <- as.character(data_id)

  common_control <- list(
    bJ = 1, sigr = 0.01, sigR = 0.01, alpha_dir = 5,
    init_J = 15L, jump_prob = c(0.50, 0.50, 0),
    ref_sigma = 1, knot_temperature = 1,
    boost = 5, xi_step_size = 0.03, n_refresh = 5L,
    odd_anchor_uniform = 0.05
  )

  make_benchmark <- function(maxInt, allowed_degrees, aJ, shrink_c,
                             exp_ratio, normalize) {
    z <- .make_spec(
      nmcmc = 200000L, nburn = 100000L, nthin = 10L,
      NB_max = 100L, maxInt = maxInt,
      allowed_degrees = allowed_degrees,
      aJ = aJ, shrink_c = shrink_c, exp_ratio = exp_ratio,
      normalize = normalize, lambda_J = 0, k_weights = NULL,
      bJ = common_control$bJ, sigr = common_control$sigr,
      sigR = common_control$sigR, alpha_dir = common_control$alpha_dir,
      init_J = common_control$init_J, jump_prob = common_control$jump_prob,
      ref_sigma = common_control$ref_sigma, boost = common_control$boost,
      xi_step_size = common_control$xi_step_size,
      n_refresh = common_control$n_refresh,
      odd_anchor_uniform = common_control$odd_anchor_uniform
    )
    z$control$knot_temperature <- common_control$knot_temperature
    z
  }

  # Bodyfat: archived sigmab = 2 * sd(Y_train); normalize = FALSE.
  if (data_id == "bodyfat") {
    return(make_benchmark(
      maxInt = 2L, allowed_degrees = 2L, aJ = 10,
      shrink_c = 5, exp_ratio = 3, normalize = FALSE
    ))
  }

  # Boston housing: use the newer real_boston.pdf result, not the older
  # Boston configuration embedded in bodyfat.pdf. Archived sigmab multiplier=2.
  if (data_id == "boston") {
    return(make_benchmark(
      maxInt = 3L, allowed_degrees = 1L, aJ = 30,
      shrink_c = 3, exp_ratio = 1, normalize = TRUE
    ))
  }

  # Concrete compressive strength (conc): archived sigmab multiplier=1.
  if (data_id == "conc") {
    return(make_benchmark(
      maxInt = 3L, allowed_degrees = 1L, aJ = 50,
      shrink_c = 3, exp_ratio = 1, normalize = TRUE
    ))
  }

  # Residential building: archived sigmab multiplier=2.
  if (data_id == "resident") {
    return(make_benchmark(
      maxInt = 3L, allowed_degrees = c(2L, 3L), aJ = 50,
      shrink_c = 3, exp_ratio = 2, normalize = TRUE
    ))
  }

  # Tecator: archived sigmab multiplier=2.
  if (data_id == "tecator") {
    return(make_benchmark(
      maxInt = 2L, allowed_degrees = 3L, aJ = 5,
      shrink_c = 20, exp_ratio = 1, normalize = TRUE
    ))
  }

  # Chemical manufacturing process (cmp): archived sigmab multiplier=1.
  if (data_id == "cmp") {
    return(make_benchmark(
      maxInt = 3L, allowed_degrees = 1L, aJ = 30,
      shrink_c = 5, exp_ratio = 1, normalize = TRUE
    ))
  }

  stop("No selected benchmark setting for data_id=", data_id)
}

