rtrunc_gamma <- function(n, lower = -Inf, upper = Inf, Tg_shape, Tg_rate) {
  pfun <- function(x) pgamma(x, shape = Tg_shape, rate = Tg_rate)
  qfun <- function(p) qgamma(p, shape = Tg_shape, rate = Tg_rate)
  a <- if (is.finite(lower)) pfun(lower) else 0
  b <- if (is.finite(upper)) pfun(upper) else 1
  if (b <= a) {
    ## No mass in the requested interval — return boundary value to avoid errors
    return(rep(if (is.finite(lower)) lower else upper, n))
  }
  u <- runif(n, a, b)
  qfun(u)
}

## Build a sampling grid for the upfront sanity check on time-varying
## probability parameters. For any inputs that came from make_time_varying(),
## we include the breakpoints (and midpoints) so the grid covers every
## changepoint. Falls back to an evenly spaced grid on [0, tf] (or 0–365 if
## tf is infinite) when no time-varying inputs are supplied.
build_sanity_grid <- function(params, tf = NULL) {
  breakpoints <- numeric(0)
  for (p in params) {
    if (inherits(p, "time_varying_fn")) {
      bp <- attr(p, "times", exact = TRUE)
      if (is.null(bp)) {
        ## Fallback for time_varying_fns built before make_time_varying
        ## exposed breakpoints as attributes.
        bp <- tryCatch(environment(p)$x, error = function(e) NULL)
      }
      if (is.numeric(bp)) {
        breakpoints <- c(breakpoints, bp)
      }
    }
  }
  breakpoints <- sort(unique(breakpoints))

  if (length(breakpoints) >= 2L) {
    mids <- (breakpoints[-length(breakpoints)] + breakpoints[-1]) / 2
    sort(unique(c(breakpoints, mids)))
  } else if (length(breakpoints) == 1L) {
    bp <- breakpoints[1]
    sort(unique(c(0, bp - 1, bp, bp + 1, max(bp + 10, 100))))
  } else {
    upper <- if (!is.null(tf) && is.finite(tf)) tf else 365
    seq(0, upper, length.out = 50)
  }
}

## Sanity-check a probability parameter (scalar or function) by sampling it
## on `grid` and confirming every resolved value is in [0, 1]. NULL inputs
## are skipped so this can be called uniformly for optional parameters.
check_probability_on_grid <- function(param, grid, param_name) {
  if (is.null(param)) return(invisible(NULL))
  values <- resolve_time_varying(param, grid, param_name)
  if (any(values < 0 | values > 1)) {
    bad <- which(values < 0 | values > 1)
    show <- bad[seq_len(min(3L, length(bad)))]
    stop(sprintf(
      "`%s` resolves to values outside [0, 1] at t = %s (got %s).",
      param_name,
      paste(round(grid[show], 3), collapse = ", "),
      paste(round(values[show], 4), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(NULL)
}

## Sanity-check a strictly positive parameter (scalar or function) by
## sampling it on `grid` and confirming every resolved value is > 0.
check_positive_on_grid <- function(param, grid, param_name) {
  if (is.null(param)) return(invisible(NULL))
  values <- resolve_time_varying(param, grid, param_name)
  if (any(values <= 0)) {
    bad <- which(values <= 0)
    show <- bad[seq_len(min(3L, length(bad)))]
    stop(sprintf(
      "`%s` must resolve to positive value(s); failed at t = %s (got %s).",
      param_name,
      paste(round(grid[show], 3), collapse = ", "),
      paste(round(values[show], 4), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(NULL)
}
