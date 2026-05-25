#' Generate incidence time series from a model trajectory
#'
#' Aggregates events from a transmission data frame into daily or weekly
#' incidence counts for one of four metrics.
#'
#' Counts every infection that was added to `tdf` (every row with a non-NA
#' `time_infection_absolute`), including leaf cases whose own onward
#' transmission was not simulated. This matches the convention used by
#' `key_outputs()`. To check whether the simulation finished, inspect
#' `all(tdf$offspring_generated)` separately.
#'
#' @param tdf     Transmission data frame from `branching_process_main()`.
#' @param metric  One of `"cases"`, `"deaths"`, `"cases_HCW"`, `"deaths_HCW"`.
#'                Cases are anchored at `time_infection_absolute`; deaths at
#'                `time_outcome_absolute` (filtered by `outcome == TRUE`).
#' @param period  `"daily"` (1-day bins) or `"weekly"` (7-day bins).
#' @param dense   If TRUE (default), the output spans every bin from day 0 to
#'                the last event, with zero counts for empty bins (useful for
#'                plotting). If FALSE, only bins with at least one event are
#'                returned.
#'
#' @return A data frame with two columns:
#'   \itemize{
#'     \item `time` - start of the bin in days from outbreak day 0.
#'     \item `incidence` - integer count of events falling in that bin.
#'   }
#'   If no events match, returns a zero-row data frame with the same schema.
#'
#' @examples
#' \dontrun{
#'   res     <- branching_process_main(...)
#'   daily_d <- generate_incidence(res$tdf, metric = "deaths", period = "daily")
#'   weekly  <- generate_incidence(res$tdf, metric = "cases_HCW",
#'                                 period = "weekly")
#' }
#'
#' @export
generate_incidence <- function(
    tdf,
    metric = c("cases", "deaths", "cases_HCW", "deaths_HCW"),
    period = c("daily", "weekly"),
    dense  = TRUE
) {

  metric <- match.arg(metric)
  period <- match.arg(period)

  ## "Real" cases = rows actually populated by the simulation. Preallocated
  ## but unused rows have time_infection_absolute = NA.
  is_real <- !is.na(tdf$time_infection_absolute)

  ## ---- Pick time anchor and event filter per metric ---------------------
  is_case <- metric %in% c("cases", "cases_HCW")
  is_hcw  <- metric %in% c("cases_HCW", "deaths_HCW")

  if (is_case) {
    event_time <- tdf$time_infection_absolute
    event_flag <- rep(TRUE, nrow(tdf))                # every row is a case
  } else {
    event_time <- tdf$time_outcome_absolute
    event_flag <- !is.na(tdf$outcome) & tdf$outcome   # deaths only, NA-safe
  }

  if (is_hcw) {
    event_flag <- event_flag & (tdf$class == "HCW")
  }

  event_subset <- is_real & event_flag & !is.na(event_time)
  event_times  <- event_time[event_subset]

  empty_out <- data.frame(time = integer(0), incidence = integer(0))
  if (length(event_times) == 0L) return(empty_out)

  ## ---- Bin --------------------------------------------------------------
  bin_width <- if (period == "daily") 1L else 7L
  bin_idx   <- as.integer(floor(event_times / bin_width))   # 0-based bin index

  if (any(bin_idx < 0L)) {
    stop("Found event time < 0; cannot bin against day 0 origin.",
         call. = FALSE)
  }

  if (dense) {
    bin_max <- max(bin_idx)
    counts  <- tabulate(bin_idx + 1L, nbins = bin_max + 1L)
    data.frame(
      time      = seq.int(0L, bin_max) * bin_width,
      incidence = counts
    )
  } else {
    tab <- table(bin_idx)
    data.frame(
      time      = as.integer(names(tab)) * bin_width,
      incidence = as.integer(tab)
    )
  }
}
