mlabs_control <- function(bJ = 1.0,
                          sigr = 0.01,
                          sigR = 0.01,
                          alpha_dir = 5.0,
                          beta0 = NULL,
                          sigmab = NULL,
                          init_J = NULL,
                          jump_prob = NULL,
                          pi_warmup = NULL,
                          ref_sigma = 1.0,
                          knot_temperature = NULL,
                          boost = 5.0,
                          xi_step_size = 0.03,
                          n_refresh = 5L,
                          odd_anchor_uniform = 0.05) {

  .pos <- function(value, name, strict = TRUE) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        (strict && value <= 0) || (!strict && value < 0))
      stop(sprintf("control$%s must be a single %s number.", name,
                   if (strict) "positive" else "non-negative"), call. = FALSE)
    as.numeric(value)
  }
  .cnt <- function(value, name, lower = 0L) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != floor(value) || value < lower)
      stop(sprintf("control$%s must be a single integer >= %d.", name, lower),
           call. = FALSE)
    as.integer(value)
  }

  bJ        <- .pos(bJ, "bJ")
  sigr      <- .pos(sigr, "sigr")
  sigR      <- .pos(sigR, "sigR")
  alpha_dir <- .pos(alpha_dir, "alpha_dir")
  ref_sigma <- .pos(ref_sigma, "ref_sigma")
  boost     <- .pos(boost, "boost", strict = FALSE)
  n_refresh <- .cnt(n_refresh, "n_refresh", 0L)

  if (!is.numeric(xi_step_size) || length(xi_step_size) != 1L ||
      !is.finite(xi_step_size) || xi_step_size <= 0 || xi_step_size >= 1)
    stop("control$xi_step_size must be a single number in (0, 1).", call. = FALSE)

  if (!is.numeric(odd_anchor_uniform) || length(odd_anchor_uniform) != 1L ||
      !is.finite(odd_anchor_uniform) ||
      odd_anchor_uniform <= 0 || odd_anchor_uniform > 1)
    stop("control$odd_anchor_uniform must be a single number in (0, 1].",
         call. = FALSE)

  if (!is.null(beta0)) {
    if (!is.numeric(beta0) || length(beta0) != 1L || !is.finite(beta0))
      stop("control$beta0 must be NULL or a single finite number.", call. = FALSE)
    beta0 <- as.numeric(beta0)
  }
  if (!is.null(sigmab))           sigmab           <- .pos(sigmab, "sigmab")
  if (!is.null(knot_temperature)) {
    knot_temperature <- .pos(knot_temperature, "knot_temperature")
    if (knot_temperature < 1)
      stop("control$knot_temperature must be >= 1.", call. = FALSE)
  }
  if (!is.null(init_J))    init_J    <- .cnt(init_J, "init_J", 1L)
  if (!is.null(pi_warmup)) pi_warmup <- .cnt(pi_warmup, "pi_warmup", 0L)

  if (!is.null(jump_prob)) {
    jump_prob <- as.numeric(jump_prob)
    if (length(jump_prob) != 3L || any(!is.finite(jump_prob)) ||
        any(jump_prob < 0) || abs(sum(jump_prob) - 1) > 1e-8)
      stop("control$jump_prob must be a non-negative vector of length three summing to one.",
           call. = FALSE)
  }

  list(bJ = bJ,
       sigr = sigr,
       sigR = sigR,
       alpha_dir = alpha_dir,
       beta0 = beta0,
       sigmab = sigmab,
       init_J = init_J,
       jump_prob = jump_prob,
       pi_warmup = pi_warmup,
       ref_sigma = ref_sigma,
       knot_temperature = knot_temperature,
       boost = boost,
       xi_step_size = xi_step_size,
       n_refresh = n_refresh,
       odd_anchor_uniform = odd_anchor_uniform)
}

.resolve_control <- function(control) {
  if (is.null(control)) control <- list()
  if (!is.list(control))
    stop("control must be a (possibly empty) named list; see ?mlabs_control.",
         call. = FALSE)
  if (!length(control)) return(mlabs_control())

  nms <- names(control)
  if (is.null(nms) || any(!nzchar(nms)))
    stop("Every element of control must be named; see ?mlabs_control.",
         call. = FALSE)
  if (anyDuplicated(nms))
    stop("control contains duplicated names: ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "), call. = FALSE)

  known <- names(formals(mlabs_control))
  unknown <- setdiff(nms, known)
  if (length(unknown))
    stop("Unknown control argument(s): ", paste(unknown, collapse = ", "),
         ".\n  Allowed: ", paste(known, collapse = ", "), call. = FALSE)

  do.call(mlabs_control, control)
}
