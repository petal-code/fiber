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
#' @export
calculate_implied_hospital_transmission_share_genPop <- function(
    p_hosp,                   # probability of an infection going to hospital
    shape_h, rate_h,          # time to hospitalisation ~ Gamma(shape_h, rate_h)
    shape_out, rate_out,      # time to outcome        ~ Gamma(shape_out, rate_out)
    gt_shape, gt_rate,        # generation time        ~ Gamma(gt_shape, gt_rate)
    kappa = 1,                # hospital keep-prob (scalar or vector in [0,1]) - to reflect reduced transmission in hospitals (quarantine)
    n = 20000,                # Monte Carlo draws to assess fraction hospital transmission
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
    # pgamma call calculates fraction of Tg probability mass that lies before admission relative to all mass before outcome
    pmin(pmax(pgamma(t_hosp, shape = gt_shape, rate = gt_rate) / F_out, 0), 1) # p_min and p_max to ensure [0, 1]
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


#' Implied hospital transmission share for HCW parents with distinct pre- and post-admission thinning
#'
#' Computes the **expected fraction** of transmissions that occur in hospital
#' (among retained/kept events) for **healthcare-worker (HCW)** parents when
#' community keep-probability is fixed at 1, but **hospital keep-probabilities
#' differ pre- vs post-admission**. Times to hospitalisation and outcome are
#' assumed **independent** and Gamma-distributed; the generation time (GT) is
#' also Gamma-distributed. The expectation is evaluated by Monte Carlo over the
#' parent population.
#'
#' For a given parent, let \eqn{q_{\mathrm{pre}} = \Pr(T \le t_{\mathrm{hosp}} \mid T \le t_{\mathrm{out}})}
#' and \eqn{q_{\mathrm{post}} = 1 - q_{\mathrm{pre}}}. Pre-admission candidate mass splits by the HCW
#' workplace probability \eqn{p_{\mathrm{preHosp}}}:
#' \eqn{q_{\mathrm{pre,h}} = p_{\mathrm{preHosp}} q_{\mathrm{pre}}}, \eqn{q_{\mathrm{pre,c}} = (1-p_{\mathrm{preHosp}}) q_{\mathrm{pre}}}.
#' Post-admission candidates are always hospital: \eqn{q_{\mathrm{post,h}} = q_{\mathrm{post}}}.
#' With community keep = 1, pre-admission hospital keep \eqn{k_{\mathrm{pre}} = 1 - \mathrm{ppe\_efficacy\_hcw}},
#' and post-admission hospital keep \eqn{k_{\mathrm{post}} = 1 - \mathrm{hospital\_quarantine\_efficacy}}, the
#' hospital share among kept events for that parent is
#' \deqn{f(k_{\mathrm{pre}}, k_{\mathrm{post}}) =
#' \frac{k_{\mathrm{pre}}\,q_{\mathrm{pre,h}} + k_{\mathrm{post}}\,q_{\mathrm{post,h}}}
#'      {q_{\mathrm{pre,c}} + k_{\mathrm{pre}}\,q_{\mathrm{pre,h}} + k_{\mathrm{post}}\,q_{\mathrm{post,h}}}.}
#' This function estimates \eqn{\Pi = E[f(k_{\mathrm{pre}}, k_{\mathrm{post}})]} across parents.
#'
#' @param p_hosp Numeric in \code{[0,1]}. Probability a parent is hospitalised.
#' @param shape_h,rate_h Positive numerics. Gamma(shape = \code{shape_h}, rate = \code{rate_h})
#'   for the time to hospitalisation (when hospitalised).
#' @param shape_out,rate_out Positive numerics. Gamma(shape = \code{shape_out}, rate = \code{rate_out})
#'   for the time to outcome (death/recovery).
#' @param gt_shape,gt_rate Positive numerics. Gamma parameters for the generation time.
#' @param prob_hospital_cond_hcw_preAdm Numeric in \code{[0,1]}. Probability that a **pre-admission**
#'   candidate infection occurs in the **hospital** setting (while the HCW is working).
#' @param ppe_efficacy_hcw Numeric in \code{[0,1]}. Efficacy of PPE/IPC at reducing **pre-admission
#'   hospital** transmission; pre-admission hospital keep-probability is \code{1 - ppe_efficacy_hcw}.
#'   May be a scalar or a numeric vector.
#' @param hospital_quarantine_efficacy Numeric in \code{[0,1]}. Efficacy of **post-admission hospital**
#'   quarantine at reducing transmission; post-admission hospital keep-probability is
#'   \code{1 - hospital_quarantine_efficacy}. May be a scalar or a numeric vector.
#' @param n Integer \eqn{\ge} 1. Number of Monte Carlo parent draws.
#' @param seed Optional integer for reproducibility passed to \code{set.seed()}.
#' @export
calculate_implied_hospital_transmission_share_hcw_separate <- function(
    p_hosp,                          # Pr(parent hospitalised)
    shape_h, rate_h,                 # time to hospitalisation ~ Gamma(shape_h, rate_h)
    shape_out, rate_out,             # time to outcome        ~ Gamma(shape_out, rate_out)
    gt_shape, gt_rate,               # generation time        ~ Gamma(gt_shape, gt_rate)
    prob_hospital_cond_hcw_preAdm,   # Pr(pre-admission candidate is hospital while working)
    ppe_efficacy_hcw,                # in [0,1]; allows vector
    hospital_quarantine_efficacy,    # in [0,1]; allows vector
    n = 20000,                       # Monte Carlo parents
    seed = NULL                      # optional reproducibility
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
  if (!is.numeric(prob_hospital_cond_hcw_preAdm) || length(prob_hospital_cond_hcw_preAdm) != 1L ||
      is.na(prob_hospital_cond_hcw_preAdm) || prob_hospital_cond_hcw_preAdm < 0 || prob_hospital_cond_hcw_preAdm > 1)
    stop("`prob_hospital_cond_hcw_preAdm` must be a single number in [0, 1].", call. = FALSE)
  if (!is.numeric(ppe_efficacy_hcw) || any(is.na(ppe_efficacy_hcw)) || any(ppe_efficacy_hcw < 0 | ppe_efficacy_hcw > 1))
    stop("`ppe_efficacy_hcw` must be a number (or vector) in [0, 1].", call. = FALSE)
  if (!is.numeric(hospital_quarantine_efficacy) || any(is.na(hospital_quarantine_efficacy)) ||
      any(hospital_quarantine_efficacy < 0 | hospital_quarantine_efficacy > 1))
    stop("`hospital_quarantine_efficacy` must be a number (or vector) in [0, 1].", call. = FALSE)
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1)
    stop("`n` must be a single integer >= 1.", call. = FALSE)

  ## -----------------------------------------
  ## Reproducibility (optional)
  ## -----------------------------------------
  if (!is.null(seed)) set.seed(seed)

  ## ------------------------------------------------
  ## 1) Simulate parent paths (independent times)
  ## ------------------------------------------------
  H      <- rbinom(n, size = 1, prob = p_hosp) == 1            # hospitalised?
  t_out  <- rgamma(n, shape = shape_out, rate = rate_out)       # outcome time
  t_hosp <- ifelse(H, rgamma(n, shape = shape_h, rate = rate_h), Inf)  # admit time or Inf

  # No post-admission window if outcome <= admission
  t_hosp[is.finite(t_hosp) & (t_out <= t_hosp)] <- Inf

  ## -------------------------------------------------------------------
  ## 2) Compute GT mass fractions: q_pre (<= t_hosp | <= t_out), q_post
  ## -------------------------------------------------------------------
  F_out <- pgamma(t_out, shape = gt_shape, rate = gt_rate)      # mass on [0, t_out]
  q_pre <- ifelse(
    !is.finite(t_hosp) | F_out <= 0,
    1,
    pmin(pmax(pgamma(t_hosp, shape = gt_shape, rate = gt_rate) / F_out, 0), 1)
  )
  q_post <- 1 - q_pre

  # Split pre-admission mass by HCW workplace probability
  p_preHosp <- prob_hospital_cond_hcw_preAdm
  q_pre_h  <- p_preHosp * q_pre         # pre-admission hospital candidates
  q_pre_c  <- (1 - p_preHosp) * q_pre   # pre-admission community candidates
  q_post_h <- q_post                    # post-admission hospital candidates

  ## -------------------------------------------------------------------
  ## 3) Evaluate Pi for each (k_pre, k_post) pair
  ## -------------------------------------------------------------------
  k_pre_vals  <- as.numeric(1 - ppe_efficacy_hcw)
  k_post_vals <- as.numeric(1 - hospital_quarantine_efficacy)

  grid <- expand.grid(k_pre = k_pre_vals, k_post = k_post_vals, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  # parent-level share function
  g_pair <- function(k_pre, k_post) {
    num <- k_pre * q_pre_h + k_post * q_post_h
    den <- q_pre_c + num
    ifelse(den > 0, num / den, 0)
  }

  # compute pi_hat and MC SE for each pair
  pi_hat <- mc_se <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g <- g_pair(grid$k_pre[i], grid$k_post[i])
    pi_hat[i] <- mean(g)
    mc_se[i]  <- sd(g) / sqrt(n)
  }

  # Feasibility ceiling at k_pre = k_post = 1
  max_pi <- mean(q_pre_h + q_post_h)

  ## -------------------
  ## 4) Return results
  ## -------------------
  out <- data.frame(
    ppe_efficacy_hcw = rep(1 - grid$k_pre, 1),
    hospital_quarantine_efficacy = rep(1 - grid$k_post, 1),
    k_pre  = grid$k_pre,
    k_post = grid$k_post,
    pi_hat = pi_hat,
    mc_se  = mc_se,
    max_pi = rep(max_pi, nrow(grid)),
    row.names = NULL,
    check.names = FALSE
  )
  out
}
