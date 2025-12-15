summarise_output <- function(
    tdf,
    subset   = "realised_subset",
    sim_info = NULL   # optional list from branching_process_main
) {

  ##--------------------------------------------------------------
  ## 0. Subset choice and basic checks
  ##--------------------------------------------------------------
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

  ##--------------------------------------------------------------
  ## 1. Outbreak timing (start / end / duration)
  ##--------------------------------------------------------------
  ## Continuous times
  first_infection_time <- min(tdf$time_infection_absolute[subset_vector],
                              na.rm = TRUE)
  last_outcome_time <- max(tdf$time_outcome_absolute[subset_vector],
                           na.rm = TRUE)

  if (!is.finite(first_infection_time) || !is.finite(last_outcome_time)) {
    stop("Could not determine outbreak start/end times (all times NA?).")
  }

  outbreak_duration_cont <- last_outcome_time - first_infection_time
  outbreak_duration_days <- round(outbreak_duration_cont)

  ##--------------------------------------------------------------
  ## 2. Basic counts by class (genPop vs HCW)
  ##--------------------------------------------------------------
  class_vec <- tdf$class

  n_cases_total <- sum(subset_vector)
  n_cases_genPop <- sum(class_vec == "genPop" & subset_vector, na.rm = TRUE)
  n_cases_HCW   <- sum(class_vec == "HCW"    & subset_vector, na.rm = TRUE)

  ## outcome == TRUE corresponds to death in the current model
  outcome_vec <- tdf$outcome

  n_deaths_total  <- sum(outcome_vec & subset_vector, na.rm = TRUE)
  n_deaths_genPop <- sum(outcome_vec & class_vec == "genPop" &
                           subset_vector, na.rm = TRUE)
  n_deaths_HCW    <- sum(outcome_vec & class_vec == "HCW" &
                           subset_vector, na.rm = TRUE)

  ## CFRs (case fatality risks)
  cfr_overall <- if (n_cases_total > 0)  n_deaths_total  / n_cases_total  else NA_real_
  cfr_genPop  <- if (n_cases_genPop > 0) n_deaths_genPop / n_cases_genPop else NA_real_
  cfr_HCW     <- if (n_cases_HCW > 0)    n_deaths_HCW    / n_cases_HCW    else NA_real_

  ##--------------------------------------------------------------
  ## 3. Attack rates (require sim_info with population & hcw_total)
  ##--------------------------------------------------------------
  population <- NA_real_
  hcw_total  <- NA_real_

  if (!is.null(sim_info)) {
    if (!is.null(sim_info$population)) {
      population <- sim_info$population
    }
    if (!is.null(sim_info$hcw_total)) {
      hcw_total <- sim_info$hcw_total
    }
  }

  ## Overall attack rate
  attack_rate_overall <- if (is.finite(population) && population > 0) {
    n_cases_total / population
  } else NA_real_

  ## Approximate genPop population as N - HCWs if both are supplied
  genpop_pop <- if (is.finite(population) && is.finite(hcw_total)) {
    population - hcw_total
  } else NA_real_

  attack_rate_genPop <- if (is.finite(genpop_pop) && genpop_pop > 0) {
    n_cases_genPop / genpop_pop
  } else NA_real_

  ## HCW attack rate = proportion of HCW stock infected
  hcw_attack_rate <- if (is.finite(hcw_total) && hcw_total > 0) {
    n_cases_HCW / hcw_total
  } else NA_real_

  ## deaths per 1,000 population (overall)
  deaths_per_1000_pop <- if (is.finite(population) && population > 0) {
    n_deaths_total / population * 1000
  } else NA_real_

  ##--------------------------------------------------------------
  ## 4. Transmission setting breakdown (community / hospital / funeral)
  ##--------------------------------------------------------------
  ## Assume tdf$setting contains strings like "community", "hospital",
  ## "funeral". If not present, we just skip this bit.
  if (!is.null(tdf$setting)) {
    setting_vec <- tdf$setting

    n_comm    <- sum(setting_vec == "community" & subset_vector, na.rm = TRUE)
    n_hosp    <- sum(setting_vec == "hospital"  & subset_vector, na.rm = TRUE)
    n_funeral <- sum(setting_vec == "funeral"   & subset_vector, na.rm = TRUE)

    n_with_setting <- n_comm + n_hosp + n_funeral

    prop_comm    <- if (n_with_setting > 0) n_comm    / n_with_setting else NA_real_
    prop_hosp    <- if (n_with_setting > 0) n_hosp    / n_with_setting else NA_real_
    prop_funeral <- if (n_with_setting > 0) n_funeral / n_with_setting else NA_real_

  } else {

    n_comm <- n_hosp <- n_funeral <- NA_real_
    prop_comm <- prop_hosp <- prop_funeral <- NA_real_
  }

  ##--------------------------------------------------------------
  ## 5. Return a named list
  ##--------------------------------------------------------------
  out <- list(
    ## Outbreak timing
    outbreak_start_time      = first_infection_time,
    outbreak_end_time        = last_outcome_time,
    outbreak_duration_cont   = outbreak_duration_cont,
    outbreak_duration_days   = outbreak_duration_days,

    ##  case and death counts by class
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

  return(out)
}
