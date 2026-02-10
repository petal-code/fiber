#' Generate offspring for a deceased parent via an unsafe funeral event
#'
#' Simulates secondary infections (offspring) generated at an **unsafe funeral**
#' for a deceased parent. Funeral offspring are produced only if:
#'   (i) the parent died, AND
#'   (ii) an unsafe funeral occurs (Bernoulli with probability p_unsafe_comm or p_unsafe_hosp).
#'
#' The number of funeral offspring is drawn from a Negative Binomial distribution.
#' Infection times are generated as a **Gamma-distributed delay** after the parent’s
#' outcome time, with small variance to represent a concentrated exposure event
#' with modest jitter. All funeral infections occur in the \code{"funeral"} setting.
#'
#' Each realised offspring is then assigned a class (\code{"HCW"} or \code{"genPop"})
#' using \code{prob_hcw_cond_funeral}.
#'
#' @param parent_hospitalised Logical scalar. Whether the parent was hospitalised.
#' @param parent_time_to_hospitalisation Numeric scalar
#'   Time from infection to admission; NA if not hospitalised.
#' @param parent_time_to_outcome Numeric scalar.
#'   Time from infection to outcome (death/recovery), have this already.
#' @param parent_died Logical scalar. Whether the parent died (TRUE) or recovered (FALSE).
#'
#' @param parent_class Character scalar, either "genPop" or "HCW". Class of the parent.
#'
#' @param p_unsafe_funeral_comm_hcw Numeric in [0,1]. Probability that a **community** death of a HCW leads
#'   to an unsafe funeral.
#' @param p_unsafe_funeral_hosp_hcw Numeric in [0,1]. Probability that a **hospital** death of a HCW  leads
#'   to an unsafe funeral (may be small but non-zero for completeness).
#' @param p_unsafe_funeral_comm_genPop Numeric in [0,1]. Probability that a **community** death of a genPop leads
#'   to an unsafe funeral.
#' @param p_unsafe_funeral_hosp_genPop Numeric in [0,1]. Probability that a **hospital** death of a genPop leads
#'   to an unsafe funeral (may be small but non-zero for completeness).
#'
#' @param mn_offspring_funeral Positive numeric. Mean of NB distribution for number
#'   of offspring at an unsafe funeral.
#' @param overdisp_offspring_funeral Positive numeric. Negative Binomial size (overdispersion).
#'
#' @param Tg_shape_funeral Positive numeric. Shape of Gamma delay distribution
#'   for time from outcome to funeral infections.
#' @param Tg_rate_funeral Positive numeric. Rate of Gamma delay distribution. ## mean GT = shape/rate
#' @param safe_funeral_efficacy numeric in [0,1]. Efficacy in preventing transmission at a safe funeral where 1 = perfect efficacy/
#' no onward transmission; 0 = no efficacy (equivalent to an unsafe funeral). Acts as a thinning parameter,
#' reducing number of offspring generated.
#'
#' @param prob_hcw_cond_funeral_hcw Numeric in [0,1]. Probability a funeral infection is an HCW.
#' @param prob_hcw_cond_funeral_genPop Numeric in [0,1]. Probability a funeral infection is a GenPop
#'
#' @return A data.frame with one row per realised funeral offspring and columns:
#'   \code{id}, \code{parent_class}, \code{setting}, \code{time_infection}, \code{class}.
#'   Returns a 0-row data.frame if no unsafe funeral occurs or zero offspring.
#' @export
offspring_function_funeral <- function(

  ## *Parent* characteristics and properties
  parent_info,

  ## Offspring distribution for unsafe funeral event
  mn_offspring_funeral = NULL, # mean number of offspring at unsafe funeral
  overdisp_offspring_funeral = NULL, # overdispersion of the above number of offspring

  ## Timing of funeral infections (delay after outcome)
  Tg_shape_funeral = NULL, # gamma shape parameter for Tg distribution at funerals ### have high shape, high rate to get low variance ##
  Tg_rate_funeral = NULL,  #gamma rate parameter for Tg distribution at funerals

  ### efficacy of a safe funeral (thinning funeral offspring)
  safe_funeral_efficacy = NULL, ## efficacy of a safe burial in reducing transmission in a funeral setting

  ## HCW vs genPop at funeral
  prob_hcw_cond_funeral_hcw = NULL, ### probability that the unsafe funeral infector infects a HCW
  prob_hcw_cond_funeral_genPop = NULL,

  ## HCW availability constraint
  hcw_available = NULL                      # number of uninfected HCWs remaining
) {

  ### note: high probability that this is genPop to genPop (which is what we want) but should allow for parent to be HCW

  ## Step 0: Extract parent info
  parent_hospitalised = parent_info$hospitalisation                          # whether the parent (infector) is hospitalised or not
  parent_time_to_hospitalisation = parent_info$time_hospitalisation_relative # if parent is hospitalised, the time of hospitalisation (relative to infection)
  parent_time_to_outcome = parent_info$time_outcome_relative                 # the time when the parent dies/recovers (relative to time of infection)
  parent_died = parent_info$outcome                                          # was the outcome at parent_time_to_outcome death? (yes/no)
  parent_class = parent_info$class                                           # "genPop" or "HCW" this matters now because unlike other two functions (although improbable) parent could be either
  parent_funeral = parent_info$funeral_safety                                # whether the funeral was safe or not (for those dying)

  #########################################################################################
  ## Input checks
  #########################################################################################

  if (is.null(parent_died) || !is.logical(parent_died) || length(parent_died) != 1L || is.na(parent_died))
    stop("`parent_died` must be a single logical value.", call. = FALSE)

  if (is.null(parent_hospitalised) || !is.logical(parent_hospitalised) || length(parent_hospitalised) != 1L || is.na(parent_hospitalised))
    stop("`parent_hospitalised` must be TRUE or FALSE.", call. = FALSE)

  if (is.null(parent_class) ||
      !parent_class %in% c("genPop", "HCW") ||
      length(parent_class) != 1L) {
    stop("`parent_class` must be either 'genPop' or 'HCW'.", call. = FALSE)
  }

  if (is.null(parent_time_to_outcome) || !is.numeric(parent_time_to_outcome) ||
      length(parent_time_to_outcome) != 1L || is.na(parent_time_to_outcome) ||
      parent_time_to_outcome < 0)
    stop("`parent_time_to_outcome` must be a single non-negative numeric value.", call. = FALSE)

  if (isTRUE(parent_hospitalised)) {
    if (is.null(parent_time_to_hospitalisation) || !is.numeric(parent_time_to_hospitalisation) ||
        length(parent_time_to_hospitalisation) != 1L || is.na(parent_time_to_hospitalisation) ||
        parent_time_to_hospitalisation < 0)
      stop("If parent is hospitalised, `parent_time_to_hospitalisation` must be non-negative.", call. = FALSE)
  } else {
    if (!is.na(parent_time_to_hospitalisation))
      stop("If parent is NOT hospitalised, `parent_time_to_hospitalisation` must be NA.", call. = FALSE)
  }

  # Other positive parameters
  for (nm in c("mn_offspring_funeral", "overdisp_offspring_funeral",
               "Tg_shape_funeral", "Tg_rate_funeral")) {
    val <- get(nm, inherits = FALSE)
    if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0)
      stop(sprintf("`%s` must be a single positive numeric value.", nm), call. = FALSE)
  }

  if (is.null(prob_hcw_cond_funeral_hcw) || !is.numeric(prob_hcw_cond_funeral_hcw) ||
      length(prob_hcw_cond_funeral_hcw) != 1L || is.na(prob_hcw_cond_funeral_hcw) ||
      prob_hcw_cond_funeral_hcw < 0 || prob_hcw_cond_funeral_hcw > 1)
    stop("`prob_hcw_cond_funeral_hcw` must be in [0, 1].", call. = FALSE)

  if (is.null(prob_hcw_cond_funeral_genPop) || !is.numeric(prob_hcw_cond_funeral_genPop) ||
      length(prob_hcw_cond_funeral_genPop) != 1L || is.na(prob_hcw_cond_funeral_genPop) ||
      prob_hcw_cond_funeral_genPop < 0 || prob_hcw_cond_funeral_genPop > 1)
    stop("`prob_hcw_cond_funeral_genPop` must be in [0, 1].", call. = FALSE)

  ## Check hcw_available parameter
  if (is.null(hcw_available) || length(hcw_available) != 1L ||
      !is.numeric(hcw_available) || is.na(hcw_available) || hcw_available < 0) {
    stop("`hcw_available` must be a single non-negative numeric value.", call. = FALSE)
  }


  #########################################################################################
  ## Logic for whether a safe or unsafe funeral occurs
  #########################################################################################

  # Step 1: Ensure parent is dead. If parent survived, no funeral transmission
  if (!isTRUE(parent_died)) {
    return(data.frame(infection_location = character(0), time_infection_relative = numeric(0), class = character(0), stringsAsFactors=FALSE))
  }

  # Step 2: Information on whether the parent had an unsafe or safe funeral
  has_unsafe_funeral <- parent_funeral # as.logical(rbinom(n = 1, size = 1, prob = p_unsafe_funeral)) # prev: Bernoulli trial for determining whether an unsafe funeral occurs
  has_unsafe_funeral <- ifelse(has_unsafe_funeral == "unsafe", TRUE, FALSE)

  #########################################################################################
  ## Produce funeral offspring
  #########################################################################################

  # Step 3: Draw raw number of infections from NB at the funeral (assuming initially an unsafe one)
  num_offspring_raw <- rnbinom(
    n    = 1,
    mu   = mn_offspring_funeral,
    size = overdisp_offspring_funeral
  )

  # Step 4: Thin the number of infections if the funeral is a safe one
  if (!has_unsafe_funeral) {
    keep_infection <- as.logical(rbinom(n = num_offspring_raw, size = 1, prob = 1 - safe_funeral_efficacy))
    num_offspring <- sum(keep_infection)
  } else {
    num_offspring <- num_offspring_raw
  }

  if (num_offspring == 0L) {
    return(data.frame(infection_location = character(0), time_infection_relative = numeric(0), class = character(0), stringsAsFactors=FALSE))
  }

  # Step 5: Generate infection times = outcome time + Gamma distributed 'delay', typically with
  #         little variance to represent a singular point from which infections arose
  delay_funeral <- rgamma(n = num_offspring, shape = Tg_shape_funeral, rate  = Tg_rate_funeral)
  infection_times <- parent_time_to_outcome + delay_funeral

  # Step 6: Assign setting ("funeral") and class (HCW or genPop) - currently, we have separate probabilities
  #         where prob_hcw_cond_funeral depends on class of the parent
  infection_settings <- rep("funeral", num_offspring)
  offspring_class <- rep("genPop", num_offspring)
  if (parent_class == "genPop") {
    flip_hcw <- as.logical(rbinom(n = num_offspring, size = 1, prob = prob_hcw_cond_funeral_genPop))
    offspring_class[flip_hcw] <- "HCW"
  } else if (parent_class == "HCW") {
    flip_hcw <- as.logical(rbinom(n = num_offspring, size = 1, prob = prob_hcw_cond_funeral_hcw))
    offspring_class[flip_hcw] <- "HCW"
  } else {
    stop("Step 6 of funeral offspring function is broken")
  }

  # Step 7: Cap HCW infections based on hcw_available - if more HCWs were generated than are available, randomly convert excess back to genPop
  hcw_idx <- which(offspring_class == "HCW")
  n_hcw_generated <- length(hcw_idx)
  if (n_hcw_generated > hcw_available) {
    # Randomly select which HCWs to convert back to genPop
    # Note: use sample.int() with indexing to avoid R's sample() single-value gotcha
    # where sample(n, size=k) samples from 1:n instead of c(n) when length(n)==1
    n_excess <- n_hcw_generated - hcw_available
    convert_idx <- hcw_idx[sample.int(n_hcw_generated, size = n_excess)]
    offspring_class[convert_idx] <- "genPop"
  }

  offspring_df <- data.frame(infection_location = infection_settings,
                             time_infection_relative = infection_times,
                             class = offspring_class,
                             stringsAsFactors = FALSE)
  return(offspring_df)
}
