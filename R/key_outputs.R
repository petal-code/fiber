#' Generate six key outputs from a model run
#'
#' Summarises a transmission data frame (`tdf`) returned by
#' `branching_process_main()` into six headline outputs: total cases,
#' total deaths, total healthcare-worker deaths, and three measures of
#' epidemic duration based on the time between the Nth death (N = 1, 10, 25)
#' and the final outcome of the last case.
#'
#' @param tdf Transmission data frame from a model run. Must contain
#'   `class`, `outcome`, `time_outcome_absolute` and (when
#'   `subset = "realised_subset"`) `offspring_generated`.
#' @param subset Either `"realised_subset"` (only cases that generated
#'   offspring, default) or `"total_tdf"` (all rows).
#'
#' @return A named list with six elements:
#'   \itemize{
#'     \item `n_cases_total`
#'     \item `n_deaths_total`
#'     \item `n_deaths_HCW`
#'     \item `duration_first_death_to_last_outcome`
#'     \item `duration_10th_death_to_last_outcome`
#'     \item `duration_25th_death_to_last_outcome`
#'   }
#'   Duration entries are `NA_real_` when there are fewer deaths than the
#'   required rank.
#'
#' @export
key_outputs <- function(tdf, subset = "realised_subset") {

  if (!(subset %in% c("total_tdf", "realised_subset"))) {
    stop('subset must be either "total_tdf" or "realised_subset"')
  }

  if (subset == "realised_subset") {
    if (is.null(tdf$offspring_generated)) {
      stop('tdf must contain column "offspring_generated" for realised_subset')
    }
    subset_vector <- tdf$offspring_generated == TRUE
  } else {
    subset_vector <- rep(TRUE, nrow(tdf))
  }

  if (!any(subset_vector)) {
    stop("No rows selected by subset; nothing to summarise.")
  }

  ## Case and death counts
  n_cases_total  <- sum(subset_vector)
  n_deaths_total <- sum(tdf$outcome & subset_vector, na.rm = TRUE)
  n_deaths_HCW   <- sum(tdf$outcome & tdf$class == "HCW" &
                          subset_vector, na.rm = TRUE)

  ## Times of all deaths in the subset, sorted ascending
  death_subset <- tdf$outcome & subset_vector & !is.na(tdf$time_outcome_absolute)
  death_times  <- sort(tdf$time_outcome_absolute[death_subset])

  ## Final outcome time across all included cases (deaths and recoveries)
  outcome_times    <- tdf$time_outcome_absolute[subset_vector]
  last_outcome_time <- suppressWarnings(max(outcome_times, na.rm = TRUE))
  if (!is.finite(last_outcome_time)) last_outcome_time <- NA_real_

  duration_from_nth_death <- function(n) {
    if (length(death_times) < n || is.na(last_outcome_time)) return(NA_real_)
    last_outcome_time - death_times[n]
  }

  list(
    n_cases_total                        = n_cases_total,
    n_deaths_total                       = n_deaths_total,
    n_deaths_HCW                         = n_deaths_HCW,
    duration_first_death_to_last_outcome = duration_from_nth_death(1),
    duration_10th_death_to_last_outcome  = duration_from_nth_death(10),
    duration_25th_death_to_last_outcome  = duration_from_nth_death(25)
  )
}
