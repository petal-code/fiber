## Internal utility: resolve a parameter that is either a scalar or a function(t)
## Not exported.
resolve_time_varying <- function(param, t, param_name = "parameter") {
  if (is.function(param)) {
    val <- param(t)
  } else if (is.numeric(param) && length(param) == 1L) {
    val <- param
  } else {
    stop(sprintf("`%s` must be a single numeric value or a function(t).", param_name),
         call. = FALSE)
  }
  if (!is.numeric(val) || length(val) != 1L || is.na(val)) {
    stop(sprintf("`%s` evaluated to an invalid value at t = %s.", param_name, t),
         call. = FALSE)
  }
  return(val)
}
