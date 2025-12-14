hcw_loss <- function(
    tdf, # takes output from branching_process_main
    subset = "realised_subset",
    outbreak_end_time = NULL,   # optional: if NULL, taken from tdf
    sim_info = NULL # takes summary of inputs to branching_process_main for given simulation
) {

  ## Check correct specification of argument "subset"
  if (!(subset %in% c("total_tdf", "realised_subset"))) {
    stop("subset must be either total_tdf or realised_subset")
  }
  if (subset == "realised_subset") {
    subset_vector <- tdf$offspring_generated == TRUE
  } else {
    subset_vector <- rep(TRUE, nrow(tdf))
  }

  ##################################################################
  ### Step 1: Define the "end of the outbreak"

  ## By default we take this as the latest observed outcome time in the
  ## transmission tree, within the chosen subset.
  ##################################################################
  max_poss_outbreak_end_time <- max(
    tdf$time_outcome_absolute[subset_vector],
    na.rm = TRUE
  )

  if (is.null(outbreak_end_time)) {

    outbreak_end_time <- max_poss_outbreak_end_time

  } else {

    ## User-specified continuous end time must not exceed the latest
    ## observed outcome time.
    if (!(outbreak_end_time <= max_poss_outbreak_end_time)) {
      stop("outbreak_end_time must be <= max(c(tdf$time_outcome_absolute[subset_vector]), na.rm = TRUE)")
    }
  }

  ## Something has gone wrong and we can't determine an end time
  if (!is.finite(outbreak_end_time)) {
    stop("hcw_loss: could not determine outbreak end time (all times are NA).")
  }

  ##################################################################
  ### Step 2: Identify HCWs in the transmission tree

  ## HCWs are defined as those with class == "HCW" and a known
  ## infection time, whose outcome occurs on/before the end of the
  ## outbreak, and who are included in the chosen subset.
  ##
  ## We count:
  ##   - how many HCWs were ever infected (n_hcw_infected)
  ##   - how many HCWs died (n_hcw_deaths; outcome == TRUE)
  ##################################################################
  outbreak_end_time_subset <- !is.na(tdf$time_outcome_absolute) &
    tdf$time_outcome_absolute <= outbreak_end_time


  is_hcw <- tdf$class == "HCW" &
    !is.na(tdf$time_infection_absolute) &
    outbreak_end_time_subset

  is_hcw_subset <- is_hcw & subset_vector

  n_hcw_infected <- sum(is_hcw_subset, na.rm = TRUE)

  ## outcome == TRUE corresponds to death for this model
  n_hcw_deaths <- sum(tdf$outcome[is_hcw_subset], na.rm = TRUE)

  ## Proportion of infected HCWs who die
  hcw_cfr <- if (n_hcw_infected > 0) {
    n_hcw_deaths / n_hcw_infected
  } else {
    NA_real_
  }

  ##################################################################
  ### Step 3: Compute HCW-days lost

  ## Definition (for each HCW i):
  ##   HCW-days-lost_i = outbreak_end_time - time_infection_absolute_i
  ##
  ## i.e. we assume HCWs cease to provide care from the time they are
  ## infected until the end of the outbreak.
  ##
  ## We return:
  ##   - hcw_days_lost           : per-HCW continuous days lost
  ##   - total_hcw_days_lost_cont: continuous total (sum of continuous)
  ##   - total_hcw_days_lost     : rounded continuous
  ##################################################################

  if (n_hcw_infected > 0) {

    hcw_days_lost <- pmax(
      outbreak_end_time - tdf$time_infection_absolute[is_hcw_subset],
      0
    )

  } else {

    ## No HCWs infected: return empty vector
    hcw_days_lost <- numeric(0)

  }

  total_hcw_days_lost_cont <- sum(hcw_days_lost)
  total_hcw_days_lost      <- round(total_hcw_days_lost_cont)



  #################################################################
  ### step 4: Scenario-specific measures


  ## Scenario-level derived measures (if sim_info provided)
  if (!is.null(sim_info) &&
      !is.null(sim_info$hcw_total) &&
      !is.null(sim_info$population) &&
      is.finite(sim_info$hcw_total) && sim_info$hcw_total > 0 &&
      is.finite(sim_info$population) && sim_info$population > 0) {


    ## total HCWs and potential HCW-days
    total_hcw      <- sim_info$hcw_total
    total_hcw_days <- sim_info$hcw_total * round(outbreak_end_time)

    ## proportion of total HCW stock that died
    prop_hcw_died_total <- n_hcw_deaths / total_hcw

    ### proportion hcw deaths relative to total population
    hcw_deaths_per_1000_pop <- n_hcw_deaths / sim_info$population * 1000

    ## proportion of potential HCW-days lost:
    prop_potential_hcw_days_lost <- if (total_hcw_days > 0) {
      total_hcw_days_lost_cont / total_hcw_days
    } else {
      NA_real_
    }

  } else {
    total_hcw                  <- NA_real_
    total_hcw_days             <- NA_real_
    prop_hcw_died_total        <- NA_real_
    hcw_deaths_per_1000_pop    <- NA_real_
    prop_potential_hcw_days_lost <- NA_real_


  }





  ##################################################################
  ### Step 5: Construct per-HCW summary dataframe

  ## For downstream diagnostics / plotting, we return a small
  ## dataframe with:
  ##   - id              : individual id in the transmission tree
  ##   - class           : should always be "HCW" here
  ##   - time_infection_absolute
  ##   - dead            : logical, TRUE if this HCW died (outcome == TRUE)
  ##   - hcw_days_lost   : continuous HCW-days lost for this HCW
  ##################################################################

  if (n_hcw_infected > 0) {

    hcw_details <- tdf[is_hcw_subset,
                       c("id", "class", "time_infection_absolute", "outcome")]

    names(hcw_details)[names(hcw_details) == "outcome"] <- "dead"
    hcw_details$hcw_days_lost <- hcw_days_lost

    ## Order by infection time then id for readability
    hcw_details <- hcw_details[order(hcw_details$time_infection_absolute,
                                     hcw_details$id), ]
    rownames(hcw_details) <- NULL

  } else {

    hcw_details <- tdf[0, c("id", "class", "time_infection_absolute"), drop = FALSE]
    hcw_details$dead          <- logical(0)
    hcw_details$hcw_days_lost <- numeric(0)

  }

  ##################################################################
  ### Step 6: Return summary as a list
  ##################################################################

  out <- list(

    ### HCW stock
    total_hcw                = total_hcw,
    total_hcw_days           = total_hcw_days,

    ### HCW counts
    n_hcw_infected           = n_hcw_infected,
    n_hcw_deaths             = n_hcw_deaths,
    hcw_cfr                  = hcw_cfr,

    ### proportions ###

    prop_hcw_died_total       = prop_hcw_died_total,      # proportion of total HCW who died
    hcw_deaths_per_1000_pop   = hcw_deaths_per_1000_pop,
    prop_potential_hcw_days_lost = prop_potential_hcw_days_lost, ## proportion of HCW_days_lost

    ## Overall outbreak timing
    outbreak_end_time        = round(outbreak_end_time),    # integer "final day"
    outbreak_end_time_cont   = outbreak_end_time,           # continuous final time

    ## HCW-days lost
    total_hcw_days_lost      = total_hcw_days_lost,         # integer total
    total_hcw_days_lost_cont = total_hcw_days_lost_cont,    # continuous total

    ## Per-HCW details
    hcw_details              = hcw_details
  )

  return(out)
}
