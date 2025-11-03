#' Generate offspring infections via healthcare transmission
#'
#' This is a first sketch of a route-specific offspring function.
#' For a single index case, it draws the number of healthcare-acquired
#' infections they cause, and returns a data frame describing those offspring.
#'
#' At this stage, this function:
#'   * uses a negative binomial distribution for the number of offspring
#'   * samples which susceptible individuals become infected, uniformly at random
#'   * is written to be easy to extend later (e.g. adding a contact matrix)
#'
#' It is *not* yet wired into the branching_process_basic() function.
#'
#' @param index_id Integer ID of the index case (e.g. row/ID in the tree).
#' @param index_group Character label for the index case group/type
#'   (e.g. "HCW", "patient"). Currently only stored; later can control rates.
#' @param mn_offspring_hc Mean number of healthcare offspring per index case.
#' @param disp_offspring_hc Overdispersion for healthcare offspring
#'   (size parameter of the negative binomial; must be > 0).
#' @param susceptible_pool Data frame describing who can be infected in
#'   healthcare settings. For now we assume it has at least:
#'     - a column `id`   : unique identifier of each susceptible
#'     - a column `group`: group/type of each susceptible (e.g. "HCW", "patient")
#' @param hc_contact_matrix Optional contact matrix for healthcare transmission,
#'   e.g. rows = index group, cols = recipient group. This is not used yet
#'   in this sketch, but left here for future extension.
#'
#' @return A list with two elements:
#'   - n_offspring  : integer, number of healthcare infections generated
#'   - offspring_df : data frame with one row per healthcare infection and columns:
#'         * transmission_route
#'         * index_id
#'         * index_group
#'         * recipient_id
#'         * recipient_group
#'
#' @keywords internal
offspring_healthcare <- function(index_id,
                                 index_group,
                                 mn_offspring_hc,
                                 disp_offspring_hc,
                                 susceptible_pool,
                                 hc_contact_matrix = NULL) {

  ## -----------------------------
  ## 1. Basic input sanity checks
  ## -----------------------------

  # Check dispersion is positive; negative binomial requires size > 0
  if (!is.numeric(disp_offspring_hc) || length(disp_offspring_hc) != 1L || disp_offspring_hc <= 0) {
    stop("disp_offspring_hc must be a single positive number.")
  }

  # Check mean is non-negative
  if (!is.numeric(mn_offspring_hc) || length(mn_offspring_hc) != 1L || mn_offspring_hc < 0) {
    stop("mn_offspring_hc must be a single non-negative number.")
  }

  # Ensure susceptible_pool is a data.frame with required columns
  if (!is.data.frame(susceptible_pool)) {
    stop("susceptible_pool must be a data.frame.")
  }
  required_cols <- c("id", "group")
  missing_cols <- setdiff(required_cols, names(susceptible_pool))
  if (length(missing_cols) > 0) {
    stop("susceptible_pool is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # If there are no susceptibles in healthcare, we cannot generate infections
  if (nrow(susceptible_pool) == 0L) {
    return(list(
      n_offspring  = 0L,
      offspring_df = data.frame(
        transmission_route = character(0),
        index_id           = integer(0),
        index_group        = character(0),
        recipient_id       = integer(0),
        recipient_group    = character(0),
        stringsAsFactors   = FALSE
      )
    ))
  }

  ## --------------------------------------------
  ## 2. Draw number of healthcare offspring (NB)
  ## --------------------------------------------

  # Negative binomial draw for number of healthcare infections
  # - mu   = mn_offspring_hc (mean)
  # - size = disp_offspring_hc (overdispersion / "size" parameter)
  n_off <- stats::rnbinom(n = 1L, mu = mn_offspring_hc, size = disp_offspring_hc)

  # If the draw is zero, return an empty offspring data frame
  if (n_off == 0L) {
    offspring_df <- data.frame(
      transmission_route = character(0),
      index_id           = integer(0),
      index_group        = character(0),
      recipient_id       = integer(0),
      recipient_group    = character(0),
      stringsAsFactors   = FALSE
    )

    return(list(
      n_offspring  = 0L,
      offspring_df = offspring_df
    ))
  }

  ## --------------------------------------------------
  ## 3. Cap n_off by available susceptibles (safety)
  ## --------------------------------------------------

  # We cannot infect more people than exist in the susceptible pool
  if (n_off > nrow(susceptible_pool)) {
    n_off <- nrow(susceptible_pool)
  }

  ## ------------------------------------------------------------------
  ## 4. Sample which susceptible individuals become infected
  ## ------------------------------------------------------------------

  # For this sketch we keep it simple: we sample uniformly at random from
  # the susceptible_pool. In future, hc_contact_matrix can be used to
  # weight sampling by recipient group.
  sampled_idx <- sample(seq_len(nrow(susceptible_pool)),
                        size    = n_off,
                        replace = FALSE)

  sampled <- susceptible_pool[sampled_idx, , drop = FALSE]

  ## -----------------------------------------
  ## 5. Construct offspring characteristics
  ## -----------------------------------------

  # Build a df describing each healthcare infection.
  # Each row = one infected individual in healthcare.
  offspring_df <- data.frame(
    transmission_route = rep("healthcare", n_off),
    index_id           = rep(index_id,      n_off),
    index_group        = rep(index_group,   n_off),
    recipient_id       = sampled$id,
    recipient_group    = sampled$group,
    stringsAsFactors   = FALSE
  )

  ## ---------------------------
  ## 6. Return as a simple list
  ## ---------------------------

  # We return both the count and the detailed data frame.
  # This mirrors the pattern in the mpoxRingVax offspring function, but
  # scoped down to a single transmission route.
  list(
    n_offspring  = n_off,
    offspring_df = offspring_df
  )
}
