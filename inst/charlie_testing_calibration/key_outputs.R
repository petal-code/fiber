#' Generate six key outputs from a model run
#'
#' Summarises a transmission data frame (`tdf`) returned by
#' `branching_process_main()` into six headline outputs: total cases,
#' total deaths, total healthcare-worker deaths, and three measures of
#' epidemic duration based on the time between the Nth death (N = 1, 10, 25)
#' and the final outcome of the last case.
#'
#' Counts every infection that was added to `tdf` (every row with a non-NA
#' `time_infection_absolute`), including leaf cases whose own onward
#' transmission was not simulated. If any such unprocessed cases are present,
#' that almost always indicates the simulation cap (`check_final_size`) was
#' hit before the outbreak completed; a warning is emitted with the count.
#'
#' @param tdf Transmission data frame from a model run. Must contain
#'   `class`, `outcome`, `time_infection_absolute`, `time_outcome_absolute`,
#'   and `offspring_generated`.
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
key_outputs <- function(tdf) {

  ## "Real" cases = rows actually populated by the simulation. Preallocated
  ## but unused rows have time_infection_absolute = NA.
  is_real <- !is.na(tdf$time_infection_absolute)

  if (!any(is_real)) {
    stop("No real cases found in tdf; nothing to summarise.", call. = FALSE)
  }

  ## If any real case did not have its offspring simulated, the cap was
  ## probably hit before the outbreak finished. Warn with the count.
  n_real        <- sum(is_real)
  n_unprocessed <- sum(is_real & !tdf$offspring_generated)
  if (n_unprocessed > 0L) {
    warning(sprintf(
      paste0("Simulation cap likely hit: %d of %d infections did not have ",
             "their onward transmission simulated (offspring_generated = FALSE). ",
             "The outbreak may not have run to completion - consider raising ",
             "`check_final_size`."),
      n_unprocessed, n_real
    ), call. = FALSE)
  }

  ## Case and death counts (over all real cases, processed or not)
  n_cases_total  <- n_real
  n_deaths_total <- sum(is_real & tdf$outcome, na.rm = TRUE)
  n_deaths_HCW   <- sum(is_real & tdf$outcome & tdf$class == "HCW",
                        na.rm = TRUE)

  ## Times of all deaths, sorted ascending
  death_subset <- is_real & tdf$outcome & !is.na(tdf$time_outcome_absolute)
  death_times  <- sort(tdf$time_outcome_absolute[death_subset])

  ## Final outcome time across all real cases (deaths and recoveries)
  outcome_times     <- tdf$time_outcome_absolute[is_real]
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
