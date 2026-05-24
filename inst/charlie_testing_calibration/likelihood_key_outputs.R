#' Likelihood of observed summary statistics given a model trajectory
#'
#' For a single model run, compute the (log-)likelihood of one or more
#' observed summary statistics under user-chosen distributions.
#'
#' Counts (`n_cases_total`, `n_deaths_total`, `n_deaths_HCW`) use either
#' Poisson or Negative Binomial with the modelled count as the mean
#' (NB requires a user-supplied `size` overdispersion parameter).
#' Duration uses a Normal with the modelled value as the mean and
#' user-supplied `sd_duration`. Only one of the three duration metrics
#' may be active at a time.
#'
#' @param tdf            Transmission data frame from `branching_process_main()`.
#' @param observed       Named list with entries matching the names returned
#'                       by `key_outputs()` for every metric being fitted.
#' @param fit_n_cases_total,fit_n_deaths_total,fit_n_deaths_HCW Logical
#'                       fit flags for the three count metrics.
#' @param fit_duration_first,fit_duration_10th,fit_duration_25th Logical
#'                       fit flags for the three duration metrics
#'                       (at most one may be TRUE).
#' @param dist_n_cases_total,dist_n_deaths_total,dist_n_deaths_HCW Either
#'                       `"poisson"` or `"nbinom"`.
#' @param size_n_cases_total,size_n_deaths_total,size_n_deaths_HCW Positive
#'                       scalar NB `size` (overdispersion); required when the
#'                       corresponding distribution is `"nbinom"`.
#' @param sd_duration    Positive scalar SD for the Normal likelihood on the
#'                       active duration metric.
#' @param log            Return log-likelihood (default TRUE) or likelihood.
#' @param details        If FALSE (default) return the total scalar. If TRUE
#'                       return a list `list(loglik, contributions)`, where
#'                       `contributions` is a length-6 named numeric vector of
#'                       per-metric (log-)likelihoods (NA for inactive metrics).
#' @param subset         Passed through to `key_outputs()`.
#'
#' @return Either a scalar (default) or a list with the total plus a per-metric
#'   contributions vector. Total is `-Inf` (or `0` when `log = FALSE`) if any
#'   active contribution is `-Inf`, e.g. NA modelled duration.
#'
#' @export
likelihood_key_outputs <- function(
    tdf,
    observed,
    fit_n_cases_total   = FALSE,
    fit_n_deaths_total  = FALSE,
    fit_n_deaths_HCW    = FALSE,
    fit_duration_first  = FALSE,
    fit_duration_10th   = FALSE,
    fit_duration_25th   = FALSE,
    dist_n_cases_total  = "poisson",
    dist_n_deaths_total = "poisson",
    dist_n_deaths_HCW   = "poisson",
    size_n_cases_total  = NULL,
    size_n_deaths_total = NULL,
    size_n_deaths_HCW   = NULL,
    sd_duration         = NULL,
    log                 = TRUE,
    details             = FALSE,
    subset              = "realised_subset"
) {

  ## ---- Duration flags: at most one ---------------------------------------
  duration_flags <- c(
    duration_first_death_to_last_outcome = fit_duration_first,
    duration_10th_death_to_last_outcome  = fit_duration_10th,
    duration_25th_death_to_last_outcome  = fit_duration_25th
  )
  if (sum(duration_flags) > 1L) {
    stop("At most one of fit_duration_first, fit_duration_10th, ",
         "fit_duration_25th may be TRUE.", call. = FALSE)
  }

  count_flags <- c(
    n_cases_total  = fit_n_cases_total,
    n_deaths_total = fit_n_deaths_total,
    n_deaths_HCW   = fit_n_deaths_HCW
  )

  if (!any(c(count_flags, duration_flags))) {
    stop("At least one fit_* flag must be TRUE.", call. = FALSE)
  }

  ## ---- Validate count distribution + overdispersion combos ---------------
  validate_count <- function(fit_flag, dist, size, name) {
    if (!fit_flag) return(invisible(NULL))
    if (!(length(dist) == 1L && dist %in% c("poisson", "nbinom"))) {
      stop(sprintf('dist_%s must be "poisson" or "nbinom".', name),
           call. = FALSE)
    }
    if (dist == "nbinom") {
      if (is.null(size) || !is.numeric(size) ||
          length(size) != 1L || !is.finite(size) || size <= 0) {
        stop(sprintf('size_%s must be a positive scalar when dist_%s = "nbinom".',
                     name, name), call. = FALSE)
      }
    }
  }
  validate_count(fit_n_cases_total,  dist_n_cases_total,  size_n_cases_total,  "n_cases_total")
  validate_count(fit_n_deaths_total, dist_n_deaths_total, size_n_deaths_total, "n_deaths_total")
  validate_count(fit_n_deaths_HCW,   dist_n_deaths_HCW,   size_n_deaths_HCW,   "n_deaths_HCW")

  if (any(duration_flags)) {
    if (is.null(sd_duration) || !is.numeric(sd_duration) ||
        length(sd_duration) != 1L || !is.finite(sd_duration) ||
        sd_duration <= 0) {
      stop("sd_duration must be a positive scalar when a duration metric ",
           "is being fitted.", call. = FALSE)
    }
  }

  ## ---- Validate observed list -------------------------------------------
  needed <- c(names(count_flags)[count_flags],
              names(duration_flags)[duration_flags])
  if (is.null(observed) || !is.list(observed) || is.null(names(observed))) {
    stop("`observed` must be a named list.", call. = FALSE)
  }
  missing_obs <- setdiff(needed, names(observed))
  if (length(missing_obs)) {
    stop("`observed` is missing entries for: ",
         paste(missing_obs, collapse = ", "), call. = FALSE)
  }
  for (nm in needed) {
    v <- observed[[nm]]
    if (is.null(v) || !is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      stop(sprintf("observed$%s must be a single finite numeric value.", nm),
           call. = FALSE)
    }
  }

  ## ---- Compute the model summary ----------------------------------------
  out <- key_outputs(tdf, subset = subset)

  ## ---- Per-metric log-likelihood contributions --------------------------
  contributions <- c(
    n_cases_total                        = NA_real_,
    n_deaths_total                       = NA_real_,
    n_deaths_HCW                         = NA_real_,
    duration_first_death_to_last_outcome = NA_real_,
    duration_10th_death_to_last_outcome  = NA_real_,
    duration_25th_death_to_last_outcome  = NA_real_
  )

  ll_count <- function(modelled, obs_val, dist, size) {
    if (is.na(modelled) || modelled < 0) return(-Inf)
    if (dist == "poisson") {
      dpois(obs_val, lambda = modelled, log = TRUE)
    } else {
      dnbinom(obs_val, mu = modelled, size = size, log = TRUE)
    }
  }

  if (fit_n_cases_total) {
    contributions["n_cases_total"] <- ll_count(
      out$n_cases_total, observed$n_cases_total,
      dist_n_cases_total, size_n_cases_total)
  }
  if (fit_n_deaths_total) {
    contributions["n_deaths_total"] <- ll_count(
      out$n_deaths_total, observed$n_deaths_total,
      dist_n_deaths_total, size_n_deaths_total)
  }
  if (fit_n_deaths_HCW) {
    contributions["n_deaths_HCW"] <- ll_count(
      out$n_deaths_HCW, observed$n_deaths_HCW,
      dist_n_deaths_HCW, size_n_deaths_HCW)
  }

  if (any(duration_flags)) {
    dur_name     <- names(duration_flags)[duration_flags]
    modelled_dur <- out[[dur_name]]
    obs_dur      <- observed[[dur_name]]
    contributions[dur_name] <- if (is.na(modelled_dur)) {
      -Inf
    } else {
      dnorm(obs_dur, mean = modelled_dur, sd = sd_duration, log = TRUE)
    }
  }

  ## ---- Combine: -Inf wins, otherwise sum of active contributions --------
  ## "Active" is a structural property (the user asked to fit this metric),
  ## not "the value is non-NA" - the two are not the same and conflating them
  ## could silently drop an unexpected NA contribution.
  active_mask <- c(
    n_cases_total                        = fit_n_cases_total,
    n_deaths_total                       = fit_n_deaths_total,
    n_deaths_HCW                         = fit_n_deaths_HCW,
    duration_first_death_to_last_outcome = fit_duration_first,
    duration_10th_death_to_last_outcome  = fit_duration_10th,
    duration_25th_death_to_last_outcome  = fit_duration_25th
  )
  active_v <- contributions[active_mask]

  ## Active contributions must be finite or -Inf. NA / NaN / +Inf means a bug
  ## or unhandled edge case in the per-metric likelihood - fail loudly.
  bad <- is.na(active_v) | (is.infinite(active_v) & active_v > 0)
  if (any(bad)) {
    stop("Unexpected non-finite log-likelihood contribution(s) for active ",
         "metric(s): ", paste(names(active_v)[bad], collapse = ", "),
         ". Each active contribution must be finite or -Inf.",
         call. = FALSE)
  }

  loglik <- if (any(is.infinite(active_v) & active_v < 0)) -Inf else sum(active_v)

  if (!log) {
    loglik        <- exp(loglik)
    contributions <- exp(contributions)
  }

  if (details) list(loglik = loglik, contributions = contributions) else loglik
}
