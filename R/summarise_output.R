#' Summarise a model run into descriptive statistics
#'
#' Computes counts, CFRs, transmission-setting breakdown, attack rates, and
#' outbreak timing/duration metrics from a transmission data frame returned
#' by `branching_process_main()`.
#'
#' Counts every infection in `tdf` (every row with a non-NA
#' `time_infection_absolute`), including leaf cases whose own onward
#' transmission was not simulated. If any such unprocessed cases are present,
#' that almost always indicates the simulation cap (`check_final_size`) was
#' hit before the outbreak completed; a warning is emitted with the count.
#'
#' Includes the three Nth-death-to-last-outcome durations previously computed
#' by `key_outputs()`, so a single call returns everything needed by both
#' descriptive workflows and the likelihood machinery.
#'
#' @param tdf      Transmission data frame from `branching_process_main()`.
#' @param sim_info Optional list with `population` and `hcw_total` (the
#'   `sim_info` element returned by `branching_process_main()`). Required
#'   for attack-rate calculations; otherwise those are returned as `NA_real_`.
#'
#' @return Named list of summary statistics. Headline fields:
#'   \itemize{
#'     \item Counts: `n_cases_total`, `n_cases_genPop`, `n_cases_HCW`,
#'                   `n_deaths_total`, `n_deaths_genPop`, `n_deaths_HCW`.
#'     \item Timing: `outbreak_start_time`, `outbreak_end_time`,
#'                   `outbreak_duration_cont`, `outbreak_duration_days`,
#'                   `time_first_death`, `time_10th_death`, `time_25th_death`,
#'                   `duration_first_death_to_last_outcome`,
#'                   `duration_10th_death_to_last_outcome`,
#'                   `duration_25th_death_to_last_outcome`.
#'     \item Setting: `n_comm`, `n_hosp`, `n_funeral`, plus proportions.
#'     \item CFRs: `cfr_overall`, `cfr_genPop`, `cfr_HCW`.
#'     \item Attack: `attack_rate_overall`, `attack_rate_genPop`,
#'                   `hcw_attack_rate`, `deaths_per_1000_pop`.
#'   }
#'
#' @export
summarise_output <- function(tdf, sim_info = NULL) {

  ##--------------------------------------------------------------
  ## 0. "Real" cases + cap-hit diagnostic
  ##--------------------------------------------------------------
  is_real <- !is.na(tdf$time_infection_absolute)

  if (!any(is_real)) {
    stop("No real cases found in tdf; nothing to summarise.", call. = FALSE)
  }

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

  ##--------------------------------------------------------------
  ## 1. Outbreak timing: first infection, Nth deaths, last outcome
  ##--------------------------------------------------------------
  first_infection_time <- min(tdf$time_infection_absolute[is_real],
                              na.rm = TRUE)

  outcome_times     <- tdf$time_outcome_absolute[is_real]
  last_outcome_time <- suppressWarnings(max(outcome_times, na.rm = TRUE))
  if (!is.finite(last_outcome_time)) last_outcome_time <- NA_real_

  death_subset <- is_real & tdf$outcome & !is.na(tdf$time_outcome_absolute)
  death_times  <- sort(tdf$time_outcome_absolute[death_subset])

  nth_death_time <- function(n) {
    if (length(death_times) < n) return(NA_real_)
    death_times[n]
  }
  time_first_death <- nth_death_time(1)
  time_10th_death  <- nth_death_time(10)
  time_25th_death  <- nth_death_time(25)

  duration_from_nth_death <- function(n) {
    t_n <- nth_death_time(n)
    if (is.na(t_n) || is.na(last_outcome_time)) return(NA_real_)
    last_outcome_time - t_n
  }
  duration_first_death_to_last_outcome <- duration_from_nth_death(1)
  duration_10th_death_to_last_outcome  <- duration_from_nth_death(10)
  duration_25th_death_to_last_outcome  <- duration_from_nth_death(25)

  if (is.finite(first_infection_time) && is.finite(last_outcome_time)) {
    outbreak_duration_cont <- last_outcome_time - first_infection_time
    outbreak_duration_days <- round(outbreak_duration_cont)
  } else {
    outbreak_duration_cont <- NA_real_
    outbreak_duration_days <- NA_real_
  }

  ##--------------------------------------------------------------
  ## 2. Counts by class (genPop vs HCW)
  ##--------------------------------------------------------------
  class_vec   <- tdf$class
  outcome_vec <- tdf$outcome

  n_cases_total  <- n_real
  n_cases_genPop <- sum(is_real & class_vec == "genPop", na.rm = TRUE)
  n_cases_HCW    <- sum(is_real & class_vec == "HCW",    na.rm = TRUE)

  n_deaths_total  <- sum(is_real & outcome_vec, na.rm = TRUE)
  n_deaths_genPop <- sum(is_real & outcome_vec & class_vec == "genPop",
                         na.rm = TRUE)
  n_deaths_HCW    <- sum(is_real & outcome_vec & class_vec == "HCW",
                         na.rm = TRUE)

  ## CFRs
  cfr_overall <- if (n_cases_total  > 0) n_deaths_total  / n_cases_total  else NA_real_
  cfr_genPop  <- if (n_cases_genPop > 0) n_deaths_genPop / n_cases_genPop else NA_real_
  cfr_HCW     <- if (n_cases_HCW    > 0) n_deaths_HCW    / n_cases_HCW    else NA_real_

  ##--------------------------------------------------------------
  ## 3. Attack rates
  ##--------------------------------------------------------------
  population <- NA_real_
  hcw_total  <- NA_real_

  if (!is.null(sim_info)) {
    if (!is.null(sim_info$population)) population <- sim_info$population
    if (!is.null(sim_info$hcw_total))  hcw_total  <- sim_info$hcw_total
  }

  attack_rate_overall <- if (is.finite(population) && population > 0) {
    n_cases_total / population
  } else NA_real_

  genpop_pop <- if (is.finite(population) && is.finite(hcw_total)) {
    population - hcw_total
  } else NA_real_

  attack_rate_genPop <- if (is.finite(genpop_pop) && genpop_pop > 0) {
    n_cases_genPop / genpop_pop
  } else NA_real_

  hcw_attack_rate <- if (is.finite(hcw_total) && hcw_total > 0) {
    n_cases_HCW / hcw_total
  } else NA_real_

  deaths_per_1000_pop <- if (is.finite(population) && population > 0) {
    n_deaths_total / population * 1000
  } else NA_real_

  ##--------------------------------------------------------------
  ## 4. Transmission setting breakdown
  ##--------------------------------------------------------------
  if (!is.null(tdf$infection_location)) {
    setting_vec <- tdf$infection_location

    n_comm    <- sum(is_real & setting_vec == "community", na.rm = TRUE)
    n_hosp    <- sum(is_real & setting_vec == "hospital",  na.rm = TRUE)
    n_funeral <- sum(is_real & setting_vec == "funeral",   na.rm = TRUE)

    n_with_setting <- n_comm + n_hosp + n_funeral

    prop_comm    <- if (n_with_setting > 0) n_comm    / n_with_setting else NA_real_
    prop_hosp    <- if (n_with_setting > 0) n_hosp    / n_with_setting else NA_real_
    prop_funeral <- if (n_with_setting > 0) n_funeral / n_with_setting else NA_real_
  } else {
    n_comm <- n_hosp <- n_funeral <- NA_real_
    prop_comm <- prop_hosp <- prop_funeral <- NA_real_
  }

  ##--------------------------------------------------------------
  ## 5. Return
  ##--------------------------------------------------------------
  list(
    ## Outbreak timing
    outbreak_start_time      = first_infection_time,
    outbreak_end_time        = last_outcome_time,
    outbreak_duration_cont   = outbreak_duration_cont,
    outbreak_duration_days   = outbreak_duration_days,

    ## Times of key deaths
    time_first_death         = time_first_death,
    time_10th_death          = time_10th_death,
    time_25th_death          = time_25th_death,

    ## Durations from Nth death to last outcome
    duration_first_death_to_last_outcome = duration_first_death_to_last_outcome,
    duration_10th_death_to_last_outcome  = duration_10th_death_to_last_outcome,
    duration_25th_death_to_last_outcome  = duration_25th_death_to_last_outcome,

    ## Case + death counts by class
    n_cases_total            = n_cases_total,
    n_cases_genPop           = n_cases_genPop,
    n_cases_HCW              = n_cases_HCW,
    n_deaths_total           = n_deaths_total,
    n_deaths_genPop          = n_deaths_genPop,
    n_deaths_HCW             = n_deaths_HCW,

    ## Setting breakdown
    n_comm                   = n_comm,
    n_hosp                   = n_hosp,
    n_funeral                = n_funeral,
    prop_comm                = prop_comm,
    prop_hosp                = prop_hosp,
    prop_funeral             = prop_funeral,

    ## CFRs
    cfr_overall              = cfr_overall,
    cfr_genPop               = cfr_genPop,
    cfr_HCW                  = cfr_HCW,

    ## Population + attack rates
    population               = population,
    hcw_total                = hcw_total,
    attack_rate_overall      = attack_rate_overall,
    attack_rate_genPop       = attack_rate_genPop,
    hcw_attack_rate          = hcw_attack_rate,
    deaths_per_1000_pop      = deaths_per_1000_pop
  )
}
