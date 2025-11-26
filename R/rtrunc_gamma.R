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
