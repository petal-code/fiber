#' Implied hospital transmission share under random hospitalisation/outcome times
#'
#' Computes the **expected fraction** of transmissions that occur in hospital
#' (among retained/kept events) when community keep-probability is fixed at 1
#' and hospital keep-probability is \code{kappa}. Times to hospitalisation and
#' outcome are assumed **independent** and Gamma-distributed; the generation
#' time (GT) is also Gamma-distributed. The expectation is evaluated by Monte
#' Carlo over the parent population.
#'
#' For a given parent with GT mass fractions \eqn{q_c} (pre-admission) and
#' \eqn{q_h} (post-admission), the hospital share among kept events is
#' \deqn{f(\kappa) = \frac{\kappa q_h}{q_c + \kappa q_h}.}
#' This function estimates \eqn{\Pi(\kappa) = E[f(\kappa)]} across parents.
#'
#' @param p_hosp Numeric in \code{[0,1]}. Probability a parent is hospitalised.
#' @param shape_h,rate_h Positive numerics. Gamma(shape = \code{shape_h}, rate = \code{rate_h})
#'   for the time to hospitalisation (when hospitalised).
#' @param shape_out,rate_out Positive numerics. Gamma(shape = \code{shape_out}, rate = \code{rate_out})
#'   for the time to outcome (death/recovery).
#' @param gt_shape,gt_rate Positive numerics. Gamma parameters for the generation time.
#' @param kappa Numeric scalar or vector in \code{[0,1]}. Hospital keep-probability
#'   used in thinning (community keep is fixed at 1).
#' @param n Integer \eqn{\ge} 1. Number of Monte Carlo parent draws.
#' @param seed Optional integer for reproducibility passed to \code{set.seed()}.
#'
#' @return A \code{data.frame} with columns:
#' \itemize{
#'   \item \code{kappa}: input hospital keep-probability.
#'   \item \code{pi_hat}: estimated expected hospital share \eqn{\Pi(\kappa)}.
#'   \item \code{mc_se}: Monte Carlo standard error of \code{pi_hat}.
#'   \item \code{max_pi}: feasibility ceiling \eqn{E[q_h]} (achieved at \code{kappa = 1}).
#' }
#'
#' @details
#' For each parent draw:
#' \enumerate{
#'   \item Sample hospitalisation indicator \eqn{H \sim \mathrm{Bernoulli}(p_\mathrm{hosp})}.
#'   \item Sample \eqn{t_\mathrm{out} \sim \mathrm{Gamma}(\text{shape\_out}, \text{rate\_out})}.
#'   \item If \eqn{H=1}, sample \eqn{t_\mathrm{hosp} \sim \mathrm{Gamma}(\text{shape\_h}, \text{rate\_h})};
#'         else set \eqn{t_\mathrm{hosp} = \infty}.
#'   \item If \eqn{t_\mathrm{out} \le t_\mathrm{hosp}}, set \eqn{t_\mathrm{hosp} = \infty}
#'         (no hospital transmission window).
#'   \item Compute GT CDF mass \eqn{F_\mathrm{GT}(\cdot)} at those times and set
#'         \eqn{q_c = F_\mathrm{GT}(\min(t_\mathrm{hosp}, t_\mathrm{out}))/F_\mathrm{GT}(t_\mathrm{out})},
#'         \eqn{q_h = 1 - q_c}.
#' }
#'
#' @examples
#' # Example: sweep kappa from 0 to 1 in steps of 0.02
#' res <- calculate_implied_hospital_transmission_share(
#'   p_hosp    = 1.0,       # everyone hospitalised in expectation
#'   shape_h   = 500, rate_h = 50,    # time to hospitalisation ~ Gamma(500, 50)
#'   shape_out = 1000, rate_out = 50, # time to outcome        ~ Gamma(1000, 50)
#'   gt_shape  = 500, gt_rate = 50,   # generation time        ~ Gamma(500, 50)
#'   kappa     = seq(0, 1, by = 0.02),
#'   n         = 20000,
#'   seed      = 42
#' )
#' print(head(res))
#' plot(res$kappa, res$pi_hat, type = "b",
#'      xlab = "kappa (hospital keep probability)",
#'      ylab = "Implied hospital share (pi_hat)")
#' abline(h = unique(res$max_pi), lty = 2)  # feasibility ceiling at kappa = 1
#'
#' @export
calculate_implied_hospital_transmission_share <- function(
    p_hosp,
    shape_h, rate_h,          # time to hospitalisation ~ Gamma(shape_h, rate_h)
    shape_out, rate_out,      # time to outcome        ~ Gamma(shape_out, rate_out)
    gt_shape, gt_rate,        # generation time        ~ Gamma(gt_shape, gt_rate)
    kappa = 1,                # hospital keep-prob (scalar or vector in [0,1])
    n = 20000,                # Monte Carlo draws
    seed = NULL               # optional reproducibility
) {
  ## ---------------------------
  ## Basic argument validation
  ## ---------------------------
  if (!is.numeric(p_hosp) || length(p_hosp) != 1L || is.na(p_hosp) || p_hosp < 0 || p_hosp > 1)
    stop("`p_hosp` must be a single number in [0, 1].", call. = FALSE)
  for (nm in c("shape_h","rate_h","shape_out","rate_out","gt_shape","gt_rate")) {
    val <- get(nm, inherits = FALSE)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0)
      stop(sprintf("`%s` must be a single positive number.", nm), call. = FALSE)
  }
  if (!is.numeric(kappa) || any(is.na(kappa)) || any(kappa < 0 | kappa > 1))
    stop("`kappa` must be a number (or vector) in [0, 1].", call. = FALSE)
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1)
    stop("`n` must be a single integer >= 1.", call. = FALSE)

  ## -----------------------------------------
  ## Reproducibility (optional, user-provided)
  ## -----------------------------------------
  if (!is.null(seed)) set.seed(seed)

  ## ------------------------------------------------
  ## 1) Simulate parent paths (independent times)
  ## ------------------------------------------------
  # Hospitalisation indicator: H = TRUE if hospitalised
  H      <- rbinom(n, size = 1, prob = p_hosp) == 1

  # Time to outcome (always drawn)
  t_out  <- rgamma(n, shape = shape_out, rate = rate_out)

  # Time to hospitalisation when hospitalised, else +Inf (no admission)
  t_hosp <- ifelse(H, rgamma(n, shape = shape_h, rate = rate_h), Inf)

  # If outcome occurs on/before admission, there is no hospital transmission window
  t_hosp[is.finite(t_hosp) & (t_out <= t_hosp)] <- Inf

  ## -------------------------------------------------------------------
  ## 2) Compute GT mass fractions qc (community) and qh (hospital)
  ## -------------------------------------------------------------------
  # F_out = GT CDF at outcome time; denom for normalising mass on [0, t_out]
  F_out <- pgamma(t_out, shape = gt_shape, rate = gt_rate)

  # qc = F_GT(min(t_hosp, t_out)) / F_GT(t_out); if no hospital window or denom ~ 0 ⇒ qc = 1
  qc <- ifelse(
    !is.finite(t_hosp) | F_out <= 0,
    1,
    pmin(pmax(pgamma(t_hosp, shape = gt_shape, rate = gt_rate) / F_out, 0), 1)
  )
  qh <- 1 - qc  # post-admission mass fraction

  ## -------------------------------------------------------------------
  ## 3) Expected hospital share for each kappa (community keep = 1)
  ## -------------------------------------------------------------------
  # For each parent, hospital share among kept = (kappa * qh) / (qc + kappa * qh)
  exp_share <- function(k) {
    mean(ifelse(qh > 0, (k * qh) / (qc + k * qh), 0))
  }

  # Vectorised over kappa values
  kappa <- as.numeric(kappa)
  pi_hat <- vapply(kappa, exp_share, numeric(1))

  # Monte Carlo SE via sample SD of g_k(t) / sqrt(n)
  g_k <- function(k) ifelse(qh > 0, (k * qh) / (qc + k * qh), 0)
  mc_se <- vapply(kappa, function(k) sd(g_k(k)) / sqrt(n), numeric(1))

  # Feasibility ceiling: expected hospital mass fraction if kappa = 1
  max_pi <- mean(qh)

  ## -------------------
  ## 4) Return results
  ## -------------------
  out <- data.frame(
    kappa  = kappa,
    pi_hat = pi_hat,
    mc_se  = mc_se,
    max_pi = rep(max_pi, length(kappa))
  )
  rownames(out) <- NULL
  out
}
