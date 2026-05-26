#' Compute reproduction number(s) from a transmission tree
#'
#' Given the transmission tree (\code{tdf}) returned by
#' \code{\link{branching_process_main}}, compute time-binned estimates of the
#' reproduction number. Two flavours are supported:
#'
#' \describe{
#'   \item{Case (cohort) reproduction number, \eqn{R_c(t)}:}{
#'     For each case infected at calendar time \eqn{t}, the realised offspring
#'     count (\code{n_offspring}) is averaged over the cohort. The time index
#'     is the INFECTOR's infection time, so \eqn{R_c} is forward-looking --
#'     it's the realised reproduction number for the cohort of cases infected
#'     at \eqn{t}.
#'   }
#'   \item{Instantaneous reproduction number, \eqn{R(t)}:}{
#'     A Cori-style renewal-equation estimator:
#'     \deqn{R(t) = \frac{I(t)}{\sum_{s \ge 1} I(t - s)\, w(s)}}
#'     where \eqn{I(t)} is incidence at time \eqn{t} (count of cases
#'     infected in bin \eqn{t}) and \eqn{w(s)} is the generation-interval
#'     probability mass function at lag \eqn{s} in units of \code{bin_width}.
#'     The time index is the TRANSMISSION event, so \eqn{R(t)} is
#'     backward-looking -- it's the \eqn{R} that explains incidence at \eqn{t}
#'     given prior incidence.
#'   }
#' }
#'
#' Cases with \code{offspring_generated == FALSE} are excluded with a warning
#' indicating the count removed. These are typically trailing cases whose
#' offspring weren't computed before the simulation hit its size cap
#' (\code{check_final_size}) -- including them would artificially deflate
#' \eqn{R_c} for the latest cohorts.
#'
#' When estimating the generation interval empirically (the default), all
#' parent-child time gaps in the (filtered) tdf are used, rounded up to the
#' nearest \code{bin_width} lag. Gaps of 0 are bumped to 1 lag because the
#' renewal equation sums over \eqn{s \ge 1}.
#'
#' @param tdf Transmission tree data frame from \code{branching_process_main()}.
#'   Must contain \code{time_infection_absolute}, \code{n_offspring},
#'   \code{offspring_generated}, \code{id}, and \code{parent}.
#' @param bin_width Numeric, positive. Bin width for both the output time
#'   grid and (when estimated empirically) the generation interval pmf.
#'   Defaults to 7 (weekly).
#' @param type One of \code{"both"} (default), \code{"case"}, or
#'   \code{"instantaneous"}. Which estimate(s) to compute.
#' @param generation_interval Optional numeric vector. The pmf of the
#'   generation interval at lags 1, 2, ..., \code{length(generation_interval)},
#'   in units of \code{bin_width}. Need not sum to 1 -- it's normalised
#'   internally. If \code{NULL} (default), estimated empirically from
#'   parent-child infection-time gaps in \code{tdf}.
#' @param min_cases_per_bin Integer, default 5. Bins with fewer cases are
#'   returned with NA for that bin's R estimate.
#'
#' @return A data frame with one row per bin (covering
#'   \code{min(time)..max(time)} on a regular grid of step \code{bin_width}):
#'   \itemize{
#'     \item \code{time}: left edge of the bin (calendar time, same units as
#'       \code{time_infection_absolute})
#'     \item \code{n_cases}: number of cases with \code{time_infection_absolute}
#'       falling in this bin
#'     \item \code{R_case}: case reproduction number for this bin (NA when
#'       \code{n_cases < min_cases_per_bin} or when type excludes \code{"case"})
#'     \item \code{R_instantaneous}: Cori-style instantaneous R for this bin
#'       (NA when the renewal-equation denominator is too small or when type
#'       excludes \code{"instantaneous"})
#'   }
#'
#' @examples
#' \dontrun{
#' out <- branching_process_main(...)
#' rn  <- compute_reproduction_number(out$tdf, bin_width = 7)
#' head(rn)
#' }
#'
#' @export
compute_reproduction_number <- function(
    tdf,
    bin_width = 7,
    type = c("both", "case", "instantaneous"),
    generation_interval = NULL,
    min_cases_per_bin = 5
) {
  type <- match.arg(type)

  ## --- Input validation --------------------------------------------------
  needed <- c("time_infection_absolute", "n_offspring",
              "offspring_generated", "id", "parent")
  missing_cols <- setdiff(needed, names(tdf))
  if (length(missing_cols) > 0L) {
    stop(sprintf("`tdf` is missing required column(s): %s",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (!is.numeric(bin_width) || length(bin_width) != 1L || bin_width <= 0) {
    stop("`bin_width` must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(min_cases_per_bin) || length(min_cases_per_bin) != 1L ||
      min_cases_per_bin < 1) {
    stop("`min_cases_per_bin` must be a single positive integer.", call. = FALSE)
  }

  ## --- Filter: keep only real, fully-processed cases --------------------
  is_real <- !is.na(tdf$time_infection_absolute)
  is_unprocessed <- is_real & (is.na(tdf$offspring_generated) |
                                  !tdf$offspring_generated)
  n_excluded <- sum(is_unprocessed)
  if (n_excluded > 0L) {
    warning(sprintf(
      "%d individual(s) had offspring_generated == FALSE; excluded from R(t) calculation.",
      n_excluded
    ), call. = FALSE)
  }
  tdf <- tdf[is_real & !is_unprocessed, , drop = FALSE]
  if (nrow(tdf) == 0L) {
    stop("No cases remain after filtering to offspring_generated == TRUE.",
         call. = FALSE)
  }

  ## --- Build the bin grid ------------------------------------------------
  bin_left <- floor(tdf$time_infection_absolute / bin_width) * bin_width
  grid_min <- min(bin_left)
  grid_max <- max(bin_left)
  bin_grid <- seq(grid_min, grid_max, by = bin_width)
  n_bins <- length(bin_grid)

  ## --- Per-bin case counts and case R -----------------------------------
  bin_idx <- match(bin_left, bin_grid)
  n_cases <- as.integer(tabulate(bin_idx, nbins = n_bins))

  R_case <- rep(NA_real_, n_bins)
  if (type %in% c("both", "case")) {
    for (i in seq_len(n_bins)) {
      sel <- bin_idx == i
      if (sum(sel) >= min_cases_per_bin) {
        R_case[i] <- mean(tdf$n_offspring[sel])
      }
    }
  }

  ## --- Instantaneous R via Cori renewal equation ------------------------
  R_inst <- rep(NA_real_, n_bins)
  if (type %in% c("both", "instantaneous")) {

    ## Get or estimate the generation-interval pmf at lags 1, 2, ...
    if (is.null(generation_interval)) {
      parent_pos <- match(tdf$parent, tdf$id)
      gaps <- tdf$time_infection_absolute - tdf$time_infection_absolute[parent_pos]
      gaps <- gaps[!is.na(gaps) & gaps > 0]
      if (length(gaps) < 5L) {
        warning(sprintf(
          "Only %d parent-child gap(s) available to estimate the generation interval; instantaneous R may be unreliable.",
          length(gaps)
        ), call. = FALSE)
      }
      if (length(gaps) == 0L) {
        warning("No valid parent-child gaps in tdf; instantaneous R will be NA throughout.",
                call. = FALSE)
        return(assemble_rn_output(bin_grid, n_cases, R_case, R_inst, type))
      }
      ## Convert gaps (continuous) to integer lags in units of bin_width.
      ## Gaps of 0 (rare, would mean simultaneous infection times) are
      ## bumped to lag 1.
      lags <- pmax(1L, ceiling(gaps / bin_width))
      max_lag <- max(lags)
      gi_pmf <- tabulate(lags, nbins = max_lag) / length(lags)
    } else {
      if (!is.numeric(generation_interval) ||
          any(!is.finite(generation_interval)) ||
          any(generation_interval < 0) ||
          sum(generation_interval) <= 0) {
        stop("`generation_interval` must be a finite non-negative numeric vector with positive sum.",
             call. = FALSE)
      }
      gi_pmf <- generation_interval / sum(generation_interval)
      max_lag <- length(gi_pmf)
    }

    ## Apply Cori renewal equation per bin
    for (t_idx in seq_len(n_bins)) {
      if (n_cases[t_idx] < min_cases_per_bin) next
      denom <- 0
      for (s in seq_len(max_lag)) {
        prev_idx <- t_idx - s
        if (prev_idx >= 1L) {
          denom <- denom + n_cases[prev_idx] * gi_pmf[s]
        }
      }
      if (denom > 0) {
        R_inst[t_idx] <- n_cases[t_idx] / denom
      }
    }
  }

  assemble_rn_output(bin_grid, n_cases, R_case, R_inst, type)
}

## Internal helper: pack outputs into a data frame, including only the
## R columns the caller asked for.
assemble_rn_output <- function(time, n_cases, R_case, R_inst, type) {
  out <- data.frame(time = time, n_cases = n_cases, stringsAsFactors = FALSE)
  if (type %in% c("both", "case")) out$R_case <- R_case
  if (type %in% c("both", "instantaneous")) out$R_instantaneous <- R_inst
  out
}
