hcw_loss <- function(
    tdf,
    subset = "realised_subset",
    outbreak_end_time = NULL   # optional: if NULL, taken from tdf
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
  ##################################################################
  ## By default we take this as the latest observed time in the
  ## transmission tree (either infection or outcome time).
  ## We keep BOTH:
  ##   - a continuous time (outbreak_end_time_cont)
  ##   - an integer day (outbreak_end_time), using ceiling().
  ##################################################################
  max_poss_outbreak_end_time <- max(c(tdf$time_outcome_absolute[subset_vector]), na.rm = TRUE)

  if (is.null(outbreak_end_time)) {

    outbreak_end_time <- max_poss_outbreak_end_time

  } else {

    ## User-specified continuous end time
    if (!(outbreak_end_time <= max_poss_outbreak_end_time)) {
      stop("outbreak_end_time must be <= max(c(tdf$time_outcome_absolute[subset_vector]), na.rm = TRUE)")
    }
  }

  ## Something has gone wrong and we can't determine an end time
  if (!is.finite(outbreak_end_time)) {
    stop("summarise_hcw_impact: could not determine outbreak end time (all times are NA).")
  }

  ##################################################################
  ### Step 2: Identify HCWs in the transmission tree
  ##################################################################
  ## HCWs are defined as those with class == "HCW" and a known
  ## infection time. We count:
  ##   - how many HCWs were ever infected
  ##   - how many HCWs died (outcome == TRUE means death)
  ##################################################################
  outbreak_end_time_subset <- (tdf$time_outcome_absolute <= outbreak_end_time)
  is_hcw <- tdf$class == "HCW" & !is.na(tdf$time_infection_absolute) & outbreak_end_time_subset
  is_hcw_subset <- is_hcw & subset_vector

  n_hcw_infected <- sum(is_hcw_subset, na.rm = TRUE)

  ## outcome == TRUE corresponds to death for this model
  n_hcw_deaths <- sum(tdf$outcome[is_hcw_subset], na.rm = TRUE)


  ##################################################################
  ### Step 3: Compute HCW-days lost
  ##################################################################
  ## Definition (for each HCW i):
  ##   HCW-days-lost_i = outbreak_end_time_cont - time_infection_absolute_i
  ##
  ## i.e. we assume HCWs cease to provide care from the time they are
  ## infected until the end of the outbreak.
  ##
  ## We do this in two ways:
  ##   - continuous (hcw_days_lost_cont)
  ##   - integer days (hcw_days_lost_int), using round().
  ##
  ## If we want a more conservative definition ("any fraction of a
  ## day counts as a full day lost"), we can swap round() for
  ## ceiling() below.
  ##################################################################

  if (n_hcw_infected > 0) {

    hcw_days_lost <- pmax(
      outbreak_end_time - tdf$time_infection_absolute[is_hcw_subset],
      0
    )

  } else {

    ## No HCWs infected: return empty vectors
    hcw_days_lost_cont <- numeric(0)

  }

  total_hcw_days_lost <- sum(round(hcw_days_lost))


  ##################################################################
  ### Step 4: Construct per-HCW summary dataframe
  ##################################################################
  ## For downstream diagnostics / plotting, we return a small
  ## dataframe with:
  ##   - id          : individual id in the transmission tree
  ##   - class       : should always be "HCW" here
  ##   - time_infection_absolute
  ##   - dead        : logical, TRUE if this HCW died (outcome == TRUE)
  ##   - hcw_days_lost : integer HCW-days lost for this HCW
  ##################################################################

  if (n_hcw_infected > 0) {

    hcw_details <- tdf[is_hcw_subset, c("id", "class", "time_infection_absolute", "outcome")]
    names(hcw_details)[names(hcw_details) == "outcome"] <- "dead"
    hcw_details$hcw_days_lost <- hcw_days_lost

    ## Order by infection time then id for readability
    hcw_details <- hcw_details[order(hcw_details$time_infection_absolute,
                                     hcw_details$id), ]
    rownames(hcw_details) <- NULL

  } else {

    hcw_details <- tdf[0, c("id", "class", "time_infection_absolute"), drop = FALSE]
    hcw_details$dead           <- logical(0)
    hcw_details$hcw_days_lost  <- integer(0)

  }


  ##################################################################
  ### Step 5: Return summary as a list (for easy extraction)
  ##################################################################

  out <- list(
    n_hcw_infected           = n_hcw_infected,
    n_hcw_deaths             = n_hcw_deaths,
    outbreak_end_time        = round(outbreak_end_time),    # integer "final day"
    outbreak_end_time_cont   = outbreak_end_time,           # continuous final time
    total_hcw_days_lost      = round(total_hcw_days_lost),  # integer total
    total_hcw_days_lost_cont = total_hcw_days_lost,         # continuous total
    hcw_details              = hcw_details
  )

  return(out)
}
