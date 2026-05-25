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


#' Obeldesivir PEP efficacy as a function of days post challenge/exposure
#'
#' This helper implements the current NHP-derived working curve used for the
#' first OBV/ODV branch. The defaults are the maximum-likelihood estimates from
#' `ODV_PEP_plot1_only_with_ribbon_GOF.R` using the scaled logistic model with
#' `dpc_zero = 15` and `k = 1`. For this model branch, efficacy is set to zero
#' after `max_dpc = 10`, so the fitted line is intentionally cut at 10 DPC.
#'
#' @param dpc Numeric vector. Days post challenge/exposure at which OBV is first
#'   received.
#' @param E0 Fitted efficacy at DPC 0 on the hazard scale.
#' @param d50 Fitted logistic midpoint.
#' @param k Fixed logistic steepness.
#' @param dpc_zero DPC at which efficacy is forced to zero in the underlying fit.
#' @param max_dpc Maximum DPC retained in this model branch; later dosing is
#'   assigned zero efficacy.
#'
#' @return Numeric vector of efficacy values between 0 and 1.
#' @export
obv_pep_efficacy_from_dpc <- function(dpc,
                                      E0 = 0.82342697,
                                      d50 = 5.59722823,
                                      k = 1,
                                      dpc_zero = 15,
                                      max_dpc = 10) {
  if (!is.numeric(dpc)) {
    stop("`dpc` must be numeric.", call. = FALSE)
  }
  if (length(dpc) == 0L) {
    return(numeric(0))
  }
  if (any(!is.finite(dpc))) {
    stop("`dpc` must contain only finite values.", call. = FALSE)
  }
  if (any(dpc < 0)) {
    stop("`dpc` must be non-negative.", call. = FALSE)
  }
  for (nm in c("E0", "d50", "k", "dpc_zero", "max_dpc")) {
    val <- get(nm, inherits = FALSE)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || !is.finite(val)) {
      stop(sprintf("`%s` must be a single finite numeric value.", nm), call. = FALSE)
    }
  }
  if (E0 < 0 || E0 > 1) {
    stop("`E0` must be in [0, 1].", call. = FALSE)
  }
  if (dpc_zero <= 0 || max_dpc < 0) {
    stop("`dpc_zero` must be positive and `max_dpc` must be non-negative.", call. = FALSE)
  }

  g_logistic <- function(d) 1 / (1 + exp(k * (d - d50)))

  g0 <- g_logistic(0)
  gz <- g_logistic(dpc_zero)
  gd <- g_logistic(dpc)
  denom <- g0 - gz
  if (!is.finite(denom) || abs(denom) < 1e-10) {
    stop("Invalid OBV efficacy curve: denominator is effectively zero.", call. = FALSE)
  }

  eff <- E0 * (gd - gz) / denom
  eff[dpc >= dpc_zero] <- 0
  eff[dpc > max_dpc] <- 0
  pmin(1, pmax(0, eff))
}

empty_obv_pep_counts <- function() {
  list(
    eligible  = 0L,
    received  = 0L,
    adherent  = 0L,
    prevented = 0L
  )
}

empty_offspring_dataframe <- function() {
  out <- data.frame(
    infection_location      = character(0),
    time_infection_relative = numeric(0),
    class                   = character(0),
    obv_pep_eligible        = logical(0),
    obv_pep_received        = logical(0),
    obv_pep_adherent        = logical(0),
    obv_pep_dpc             = numeric(0),
    obv_pep_efficacy        = numeric(0),
    stringsAsFactors = FALSE
  )
  attr(out, "obv_pep_counts") <- empty_obv_pep_counts()
  out
}

resolve_obv_efficacy <- function(obv_pep_efficacy, dpc) {
  if (length(dpc) == 0L) {
    return(numeric(0))
  }

  if (is.null(obv_pep_efficacy)) {
    value <- obv_pep_efficacy_from_dpc(dpc)
  } else if (is.function(obv_pep_efficacy)) {
    value <- obv_pep_efficacy(dpc)
  } else {
    value <- obv_pep_efficacy
  }

  if (!is.numeric(value)) {
    stop("`obv_pep_efficacy` must be NULL, a function(dpc), or a numeric value in [0, 1].", call. = FALSE)
  }
  if (length(value) == 1L && length(dpc) > 1L) {
    value <- rep(value, length(dpc))
  }
  if (length(value) != length(dpc)) {
    stop("`obv_pep_efficacy` must return one value per DPC, or a single value to recycle.", call. = FALSE)
  }
  if (any(!is.finite(value)) || any(value < 0 | value > 1)) {
    stop("`obv_pep_efficacy` must resolve to value(s) in [0, 1].", call. = FALSE)
  }

  as.numeric(value)
}

#' Apply the obeldesivir PEP gate to candidate offspring
#'
#' This is an infection-prevention gate. It is applied after existing setting
#' gates such as IPC/PPE and ETU/quarantine, and after offspring class has been
#' assigned, because the first OBV branch targets exposed HCWs. A candidate
#' exposure is eligible when the candidate offspring is an HCW and the exposure
#' occurred in one of `obv_pep_target_locations` (hospital by default). Among
#' eligible exposures, `obv_pep_coverage_hcw` decides whether OBV is received,
#' `obv_pep_adherence` decides whether the course is effectively adhered to, and
#' `obv_pep_efficacy(dpc)` thins infections among adherent recipients.
#'
#' @param infection_location Character vector of candidate exposure settings.
#' @param offspring_class Character vector of candidate offspring classes.
#' @param infection_time_absolute Numeric vector of candidate exposure clock times.
#' @param obv_pep_enabled Logical scalar.
#' @param obv_pep_coverage_hcw Numeric between 0 and 1 or function(t). Probability an
#'   eligible HCW exposure receives OBV.
#' @param obv_pep_adherence Numeric between 0 and 1 or function(t). Probability an OBV
#'   recipient adheres sufficiently for efficacy to apply.
#' @param obv_pep_dpc Non-negative numeric or function(t). Days post challenge /
#'   exposure at which first dose is received. For the simplest working scenario,
#'   set this to 1.
#' @param obv_pep_efficacy NULL, numeric between 0 and 1, or function(dpc). If NULL, the
#'   NHP-derived `obv_pep_efficacy_from_dpc()` helper is used.
#' @param obv_pep_target_locations Character vector of exposure locations eligible
#'   for OBV PEP. Defaults should usually be "hospital" for HCW PEP.
#'
#' @return A list containing `keep`, row-level `metadata`, and aggregate `counts`.
#' @noRd
apply_obv_pep_gate <- function(infection_location,
                               offspring_class,
                               infection_time_absolute,
                               obv_pep_enabled = FALSE,
                               obv_pep_coverage_hcw = 0,
                               obv_pep_adherence = 1,
                               obv_pep_dpc = 1,
                               obv_pep_efficacy = NULL,
                               obv_pep_target_locations = "hospital") {

  n <- length(infection_location)
  if (length(offspring_class) != n || length(infection_time_absolute) != n) {
    stop("OBV PEP gate inputs must have the same length.", call. = FALSE)
  }

  metadata <- data.frame(
    obv_pep_eligible = rep(FALSE, n),
    obv_pep_received = rep(FALSE, n),
    obv_pep_adherent = rep(FALSE, n),
    obv_pep_dpc      = rep(NA_real_, n),
    obv_pep_efficacy = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )
  keep <- rep(TRUE, n)
  counts <- empty_obv_pep_counts()

  if (n == 0L || !isTRUE(obv_pep_enabled)) {
    return(list(keep = keep, metadata = metadata, counts = counts))
  }

  if (!is.character(obv_pep_target_locations) || length(obv_pep_target_locations) < 1L) {
    stop("`obv_pep_target_locations` must be a non-empty character vector.", call. = FALSE)
  }

  validate_probability_or_time_varying <- function(param, param_name) {
    if (is.function(param)) {
      return(invisible(NULL))
    }
    if (!is.numeric(param) || length(param) != 1L || is.na(param) || param < 0 || param > 1) {
      stop(sprintf("`%s` must be a function(t) or a single numeric in [0, 1].", param_name),
           call. = FALSE)
    }
    invisible(NULL)
  }

  validate_nonnegative_or_time_varying <- function(param, param_name) {
    if (is.function(param)) {
      return(invisible(NULL))
    }
    if (!is.numeric(param) || length(param) != 1L || is.na(param) || param < 0) {
      stop(sprintf("`%s` must be a function(t) or a single non-negative numeric.", param_name),
           call. = FALSE)
    }
    invisible(NULL)
  }

  resolve_probability <- function(param, t, param_name) {
    value <- resolve_time_varying(param = param, t = t, param_name = param_name)
    if (any(value < 0 | value > 1)) {
      stop(sprintf("`%s` must resolve to value(s) in [0, 1].", param_name), call. = FALSE)
    }
    value
  }

  resolve_nonnegative <- function(param, t, param_name) {
    value <- resolve_time_varying(param = param, t = t, param_name = param_name)
    if (any(!is.finite(value)) || any(value < 0)) {
      stop(sprintf("`%s` must resolve to finite, non-negative value(s).", param_name), call. = FALSE)
    }
    value
  }

  validate_probability_or_time_varying(obv_pep_coverage_hcw, "obv_pep_coverage_hcw")
  validate_probability_or_time_varying(obv_pep_adherence, "obv_pep_adherence")
  validate_nonnegative_or_time_varying(obv_pep_dpc, "obv_pep_dpc")

  eligible <- offspring_class == "HCW" & infection_location %in% obv_pep_target_locations
  metadata$obv_pep_eligible <- eligible
  counts$eligible <- sum(eligible)

  if (!any(eligible)) {
    return(list(keep = keep, metadata = metadata, counts = counts))
  }

  eligible_idx <- which(eligible)
  coverage_t <- resolve_probability(
    obv_pep_coverage_hcw,
    infection_time_absolute[eligible_idx],
    "obv_pep_coverage_hcw"
  )
  received <- as.logical(rbinom(n = length(eligible_idx), size = 1, prob = coverage_t))
  metadata$obv_pep_received[eligible_idx] <- received
  counts$received <- sum(received)

  if (!any(received)) {
    return(list(keep = keep, metadata = metadata, counts = counts))
  }

  received_idx <- eligible_idx[received]
  adherence_t <- resolve_probability(
    obv_pep_adherence,
    infection_time_absolute[received_idx],
    "obv_pep_adherence"
  )
  adherent <- as.logical(rbinom(n = length(received_idx), size = 1, prob = adherence_t))
  metadata$obv_pep_adherent[received_idx] <- adherent
  counts$adherent <- sum(adherent)

  dpc <- resolve_nonnegative(
    obv_pep_dpc,
    infection_time_absolute[received_idx],
    "obv_pep_dpc"
  )
  metadata$obv_pep_dpc[received_idx] <- dpc

  efficacy <- rep(0, length(received_idx))
  prevented <- rep(FALSE, length(received_idx))
  if (any(adherent)) {
    efficacy[adherent] <- resolve_obv_efficacy(obv_pep_efficacy, dpc[adherent])
    prevented[adherent] <- as.logical(rbinom(n = sum(adherent), size = 1,
                                             prob = efficacy[adherent]))
  }
  metadata$obv_pep_efficacy[received_idx] <- efficacy
  keep[received_idx[prevented]] <- FALSE
  counts$prevented <- sum(prevented)

  list(keep = keep, metadata = metadata, counts = counts)
}
