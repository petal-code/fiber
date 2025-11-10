calculate_implied_hospital_transmission_share <- function(
    p_hosp,
    shape_h, rate_h,          # time to hospitalisation ~ Gamma(shape_h, rate_h)
    shape_out, rate_out,      # time to outcome        ~ Gamma(shape_out, rate_out)
    gt_shape, gt_rate,        # generation time        ~ Gamma(gt_shape, gt_rate)
    kappa = 1,                # hospital keep-prob (can be scalar or vector in [0,1])
    n = 20000,                # Monte Carlo draws
    seed = NULL               # optional reproducibility
) {
  if (!is.null(seed)) set.seed(seed)

  # simulate parent paths (independent times)
  H      <- rbinom(n, 1, p_hosp) == 1
  t_out  <- rgamma(n, shape = shape_out, rate = rate_out)
  t_hosp <- ifelse(H, rgamma(n, shape = shape_h, rate = rate_h), Inf)
  # outcome on/before admission ⇒ no hospital window
  t_hosp[is.finite(t_hosp) & (t_out <= t_hosp)] <- Inf

  # generation-time mass components
  F_out <- pgamma(t_out, shape = gt_shape, rate = gt_rate)
  # qc = mass before min(t_hosp, t_out) relative to mass before t_out
  qc <- ifelse(!is.finite(t_hosp) | F_out <= 0, 1,
               pmin(pmax(pgamma(t_hosp, shape = gt_shape, rate = gt_rate) / F_out, 0), 1))
  qh <- 1 - qc

  # helper: expected hospital share for a given k
  exp_share <- function(k) {
    mean(ifelse(qh > 0, (k * qh) / (qc + k * qh), 0))
  }

  # allow vector kappa
  kappa <- as.numeric(kappa)
  pi_hat <- vapply(kappa, exp_share, numeric(1))

  # Monte Carlo SE (use delta on the plug-in estimator)
  g_k <- function(k) ifelse(qh > 0, (k * qh) / (qc + k * qh), 0)
  se <- vapply(kappa, function(k) sd(g_k(k)) / sqrt(n), numeric(1))

  # maximum achievable share (when kappa=1)
  max_pi <- mean(qh)

  out <- data.frame(
    kappa = kappa,
    pi_hat = pi_hat,
    mc_se = se,
    max_pi = rep(max_pi, length(kappa))
  )
  rownames(out) <- NULL
  out
}

## assume implied_hospital_share() is already defined

res <- implied_hospital_share(
  p_hosp   = 1,     # Pr(hospitalised)
  shape_h  = 500,     # time to hospitalisation ~ Gamma(3, 0.6)
  rate_h   = 50,
  shape_out= 1000,       # time to outcome ~ Gamma(4, 0.5)
  rate_out = 50,
  gt_shape = 500,     # generation time ~ Gamma(2.5, 0.5)
  gt_rate  = 50,
  kappa    = seq(0, 1, 0.02),  # hospital keep probabilities
  n        = 20000,   # Monte Carlo draws
  seed     = 42
)

print(res)

## (Optional) quick plot of implied hospital share vs kappa
plot(res$kappa, res$pi_hat, type = "b", xlab = "kappa (hospital keep prob)",
     ylab = "Implied hospital share (pi_hat)")
abline(h = unique(res$max_pi), lty = 2)
