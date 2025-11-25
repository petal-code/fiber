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
#' @param p_unsafe_funeral_comm Numeric in [0,1]. Probability that a **community** death leads
#'   to an unsafe funeral.
#' @param p_unsafe_funeral_hosp Numeric in [0,1]. Probability that a **hospital** death leads
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
#' @param prob_hcw_cond_funeral Numeric in [0,1]. Probability a funeral infection
#'   is an HCW.
#'
#' @return A data.frame with one row per realised funeral offspring and columns:
#'   \code{id}, \code{parent_class}, \code{setting}, \code{time_infection}, \code{class}.
#'   Returns a 0-row data.frame if no unsafe funeral occurs or zero offspring.
#' @export

offspring_function_funeral <- function(

  ## *Parent* characteristics and properties
  parent_hospitalised = NULL, ### is the parent hospitalised (yes/no)
  parent_time_to_hospitalisation = NULL, ### when over the course of the infection was the parent hospitalised
  parent_time_to_outcome = NULL, ### when does the infection resolve (death or recovered)
  parent_died = NULL, ### was the outcome at parent_time_to_outcome death? (yes/no)
  parent_class = NULL,   # "genPop" or "HCW" this matters now because unlike other two functions (although improbable) parent could be either

  ## Funeral occurrence
  p_unsafe_funeral_comm = NULL, ## probability the parent had an unsafe funeral Cond. on a community death
  p_unsafe_funeral_hosp = NULL, ## probability parent had an unsafe funeral cond. on a hospital death

  ## Offspring distribution for unsafe funeral event
  mn_offspring_funeral = NULL, # mean number of offspring at unsafe funeral
  overdisp_offspring_funeral = NULL, # overdispersion of the above number of offspring

  ## Timing of funeral infections (delay after outcome)
  Tg_shape_funeral = NULL, # gamma shape parameter for Tg distribution at funerals ### have high shape, high rate to get low variance ##
  Tg_rate_funeral = NULL,  #gamma rate parameter for Tg distribution at funerals

  ### efficacy of a safe funeral (thinning funeral offspring)
  safe_funeral_efficacy = NULL, ## efficacy of a safe burial in reducing transmission in a funeral setting

  ## HCW vs genPop at funeral
  prob_hcw_cond_funeral = NULL ### probability that the unsafe funeral infector infects a HCW

) {

  ### note: high probability that this is genPop to genPop (which is what we want) but should allow for parent to be HCW

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

  # p_unsafe_funeral_comm/hosp
  for (nm in c("p_unsafe_funeral_comm", "p_unsafe_funeral_hosp")) {
    val <- get(nm, inherits = FALSE)
    if (is.null(val) || length(val) != 1L || !is.numeric(val) || is.na(val) || val < 0 || val > 1)
      stop(sprintf("`%s` must be a numeric scalar in [0, 1].", nm), call. = FALSE)
  }

 # safe funeral efficacy
  if (is.null(safe_funeral_efficacy) || !is.numeric(safe_funeral_efficacy) ||
      safe_funeral_efficacy < 0 || safe_funeral_efficacy > 1)
    stop("`safe_funeral_efficacy` must be between 0 and 1.", call. = FALSE)


  # Other positive parameters
  for (nm in c("mn_offspring_funeral", "overdisp_offspring_funeral",
               "Tg_shape_funeral", "Tg_rate_funeral")) {
    val <- get(nm, inherits = FALSE)
    if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0)
      stop(sprintf("`%s` must be a single positive numeric value.", nm), call. = FALSE)
  }

  if (is.null(prob_hcw_cond_funeral) || !is.numeric(prob_hcw_cond_funeral) ||
      length(prob_hcw_cond_funeral) != 1L || is.na(prob_hcw_cond_funeral) ||
      prob_hcw_cond_funeral < 0 || prob_hcw_cond_funeral > 1)
    stop("`prob_hcw_cond_funeral` must be in [0, 1].", call. = FALSE)


  #########################################################################################
  ## Logic for whether a safe or unsafe funeral occurs
  #########################################################################################
  # step 1, ensure parent is dead
  # If parent survived, no funeral transmission
  if (!isTRUE(parent_died)) {
    return(data.frame(id=integer(0), parent_class=character(0), setting=character(0),
                      time_infection=numeric(0), class=character(0), stringsAsFactors=FALSE))
  }
  # step 2
  # Determine whether death occurred BEFORE hospitalisation (community death)
  community_death <- !isTRUE(parent_hospitalised) ||
    (parent_time_to_outcome < parent_time_to_hospitalisation)


  # Assign the appropriate probability
  p_funeral <- if (isTRUE(community_death)) p_unsafe_funeral_comm else p_unsafe_funeral_hosp

  # step 3
  # Bernoulli trial: does an unsafe funeral occur?
  has_unsafe_funeral <- as.logical(rbinom(n = 1, size = 1, prob = p_funeral))


  if (!has_unsafe_funeral) {
    # Funeral was safe — thinning offspring mean according to efficacy
    mn_offspring_funeral <- mn_offspring_funeral * (1 - safe_funeral_efficacy)
  }



  #########################################################################################
  ## Produce funeral offspring
  #########################################################################################

  # Step 4:draw raw number of infections from NB at the unsafe funeral
  num_offspring_raw <- rnbinom(
    n    = 1,
    mu   = mn_offspring_funeral,
    size = overdisp_offspring_funeral
  )

  if (num_offspring_raw == 0L) {
    return(data.frame(id=integer(0), parent_class=character(0), setting=character(0),
                      time_infection=numeric(0), class=character(0), stringsAsFactors=FALSE))
  }

  # Step 5: infection times = outcome time + Gamma distributed 'delay' with little variance to represent a singular point from which infections arose
  delay_funeral <- rgamma(
    n     = num_offspring_raw,
    shape = Tg_shape_funeral,
    rate  = Tg_rate_funeral
  )
  infection_times <- parent_time_to_outcome + delay_funeral


  #########################################################################################
  # step 6 Assign setting ("funeral") and class (HCW or genPop)
  #########################################################################################

  infection_settings <- rep("funeral", num_offspring_raw)

  offspring_class <- rep("genPop", num_offspring_raw)
  flip_hcw <- as.logical(rbinom(n=num_offspring_raw, size=1, prob=prob_hcw_cond_funeral))
  offspring_class[flip_hcw] <- "HCW"


  #########################################################################################
  ## step 7 Output dataframe
  #########################################################################################

  offspring_df <- data.frame(
    id             = seq_len(num_offspring_raw),
    parent_class = rep(parent_class, num_offspring_raw), ## note this allows for both
    setting        = infection_settings,
    time_infection = infection_times,
    class          = offspring_class,
    stringsAsFactors = FALSE
  )

  return(offspring_df)
}
