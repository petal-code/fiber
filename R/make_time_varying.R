#' Create a Time-Varying Parameter Function
#'
#' Constructs a function of time from a set of breakpoints, using either
#' piecewise-linear interpolation or a step (constant) function.
#' This is a convenience wrapper around [stats::approxfun()] designed for
#' use with time-varying simulation parameters such as
#' `prob_hospitalised_genPop` and `hospitalisation_delay_factor`.
#'
#' @param times Numeric vector of time points (non-decreasing, length >= 2).
#' @param values Numeric vector of corresponding values (same length as `times`).
#' @param method Character string: `"linear"` (default) for piecewise-linear
#'   interpolation, or `"constant"` for a left-continuous step function
#'   (each value holds until the next breakpoint).
#'
#' @return A function `f(t)` that accepts a single numeric time and returns the
#'   interpolated or stepped value. Values are clamped to `values[1]` for
#'   `t < min(times)` and `values[length(values)]` for `t > max(times)`.
#'   The returned function has class `c("time_varying_fn", "function")`.
#'
#' @examples
#' # Piecewise-linear: probability rises from 0.05 to 0.5 over 60 days
#' prob_hosp <- make_time_varying(
#'   times  = c(0, 14, 30, 60),
#'   values = c(0.05, 0.1, 0.3, 0.5)
#' )
#' prob_hosp(20)
#'
#' # Step function: constant between changepoints
#' prob_hosp_step <- make_time_varying(
#'   times  = c(0, 14, 30, 60),
#'   values = c(0.05, 0.1, 0.3, 0.5),
#'   method = "constant"
#' )
#' prob_hosp_step(20)
#'
#' @export
make_time_varying <- function(times, values, method = "linear") {

  ## --- Input validation ---
  if (!is.numeric(times)) {
    stop("`times` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(values)) {
    stop("`values` must be a numeric vector.", call. = FALSE)
  }
  if (length(times) != length(values)) {
    stop("`times` and `values` must have the same length.", call. = FALSE)
  }
  if (length(times) < 2L) {
    stop("`times` and `values` must have length >= 2.", call. = FALSE)
  }
  if (is.unsorted(times, strictly = FALSE)) {
    stop("`times` must be non-decreasing.", call. = FALSE)
  }
  method <- match.arg(method, choices = c("linear", "constant"))


  ## --- Build interpolation function ---
  if (method == "linear") {
    fn <- stats::approxfun(times, values, rule = 2)
  } else {
    ## f = 0 gives left-continuous step: value changes at the breakpoint
    fn <- stats::approxfun(times, values, method = "constant", rule = 2, f = 0)
  }

  class(fn) <- c("time_varying_fn", "function")
  return(fn)
}
