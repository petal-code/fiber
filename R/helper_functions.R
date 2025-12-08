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

prob_hosp_given_symptoms <- function(prob_hosp, prob_symptomatic) {
  if (prob_hosp > prob_symptomatic) {
    stop("P(hospitalisation) cannot exceed P(symptomatic).", call. = FALSE)
  }
  if (prob_symptomatic == 0 & prob_hosp > 0) {
    stop("P(symptomatic) = 0 but P(hospitalisation) > 0: impossible.", call. = FALSE)
  }
  return(prob_hosp / prob_symptomatic)
}

prob_death_given_symptoms <- function(prob_death, prob_symptomatic) {
  if (prob_death > prob_symptomatic) {
    stop("P(death) cannot exceed P(symptomatic).", call. = FALSE)
  }
  if (prob_symptomatic == 0 & prob_death > 0) {
    stop("P(symptomatic) = 0 but P(death) > 0: impossible.", call. = FALSE)
  }
  return(prob_death / prob_symptomatic)
}
