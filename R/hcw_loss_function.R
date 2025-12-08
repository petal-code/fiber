hcw_loss <- function(
    tdf,
    outbreak_end_time = NULL   # optional: if NULL, taken from tdf
) {

  ##################################################################
  ### Step 1: Define the "end of the outbreak"
  ##################################################################
  ## By default we take this as the latest observed time in the
  ## transmission tree (either infection or outcome time).
  ## We keep BOTH:
  ##   - a continuous time (outbreak_end_time_cont)
  ##   - an integer day (outbreak_end_time), using ceiling().
  ##################################################################

  if (is.null(outbreak_end_time)) {

    outbreak_end_time_cont <- max(
      c(tdf$time_infection_absolute, tdf$time_outcome_absolute),
      na.rm = TRUE
    )

  } else {

    ## User-specified continuous end time
    outbreak_end_time_cont <- outbreak_end_time

  }

  ## Something has gone wrong and we can't determine an end time
  if (!is.finite(outbreak_end_time_cont)) {
    stop("summarise_hcw_impact: could not determine outbreak end time (all times are NA).")
  }

  ## Integer "final day of outbreak" (e.g. day 98 if last event is at 97.6 days)
  outbreak_end_time_int <- ceiling(outbreak_end_time_cont)


  ##################################################################
  ### Step 2: Identify HCWs in the transmission tree
  ##################################################################
  ## HCWs are defined as those with class == "HCW" and a known
  ## infection time. We count:
  ##   - how many HCWs were ever infected
  ##   - how many HCWs died (outcome == TRUE means death)
  ##################################################################

  is_hcw <- tdf$class == "HCW" & !is.na(tdf$time_infection_absolute)

  n_hcw_infected <- sum(is_hcw, na.rm = TRUE)

  ## outcome == TRUE corresponds to death for this model
  n_hcw_deaths <- sum(tdf$class == "HCW" & tdf$outcome, na.rm = TRUE)


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

    hcw_days_lost_cont <- pmax(
      outbreak_end_time_cont - tdf$time_infection_absolute[is_hcw],
      0
    )

    ## Integer days: here we use round() to avoid systematic
    ## upwards bias from rounding every individual up.
    ## Replace round() with ceiling() if prefer the
    ## conservative "partial day = full day" definition.
    hcw_days_lost_int <- round(hcw_days_lost_cont)

  } else {

    ## No HCWs infected: return empty vectors
    hcw_days_lost_cont <- numeric(0)
    hcw_days_lost_int  <- integer(0)

  }

  total_hcw_days_lost_cont <- sum(hcw_days_lost_cont)
  total_hcw_days_lost_int  <- sum(hcw_days_lost_int)


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

    hcw_details <- tdf[is_hcw, c("id", "class", "time_infection_absolute", "outcome")]
    names(hcw_details)[names(hcw_details) == "outcome"] <- "dead"

    hcw_details$hcw_days_lost <- hcw_days_lost_int

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
    outbreak_end_time        = outbreak_end_time_int,    # integer "final day"
    outbreak_end_time_cont   = outbreak_end_time_cont,   # continuous final time
    total_hcw_days_lost      = total_hcw_days_lost_int,  # integer total
    total_hcw_days_lost_cont = total_hcw_days_lost_cont, # continuous total
    hcw_details              = hcw_details
  )

  return(out)
}
