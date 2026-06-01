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
#' @param parent_info One-row data.frame/list of parent attributes. Funeral transmission reads
#'   `hospitalisation`, `time_hospitalisation_relative`, `time_outcome_relative`, `outcome`
#'   (did the parent die?), `class` ("genPop"/"HCW"), `funeral_safety` ("safe"/"unsafe") and
#'   `time_infection_absolute`.
#'
#' @param mn_offspring_funeral Positive numeric or function(t). Mean of the NB distribution for the
#'   number of offspring at an unsafe funeral. May be supplied as a single positive scalar or as a
#'   function of absolute calendar time (e.g. built with \code{\link{make_time_varying}}); when
#'   time-varying it is resolved at the parent's absolute death (outcome) time, since the funeral
#'   occurs then rather than at the parent's infection time.
#' @param overdisp_offspring_funeral Positive numeric. Negative Binomial size (overdispersion).
#'
#' @param Tg_shape_funeral Positive numeric. Shape of Gamma delay distribution
#'   for time from outcome to funeral infections.
#' @param Tg_rate_funeral Positive numeric. Rate of Gamma delay distribution. ## mean GT = shape/rate
#' @param safe_funeral_efficacy numeric between 0 and 1. Efficacy in preventing transmission at a safe funeral where 1 = perfect efficacy/no onward transmission; 0 = no efficacy (equivalent to an unsafe funeral). This remains scalar; the time-varying probability that a funeral is safe or unsafe is resolved upstream.
#'
#' @param obv_pep_enabled Logical scalar. If TRUE, applies an OBV PEP infection-prevention
#'   gate around the safe-funeral thinning step. By default \code{obv_pep_target_locations
#'   = "hospital"}, so funeral exposures are not eligible; set
#'   \code{obv_pep_target_locations = c("hospital", "funeral", ...)} to broaden.
#' @param obv_pep_coverage Numeric in \code{[0,1]} or function(t). Probability a
#'   pre-thinning eligible candidate receives OBV at calendar time \code{t}.
#' @param obv_pep_adherence Numeric in \code{[0,1]} or function(t). Probability an OBV recipient
#'   adheres sufficiently for efficacy to apply.
#' @param obv_pep_dpc Non-negative numeric or function(t). Days post challenge/exposure to first
#'   dose. The simplest working assumption is \code{obv_pep_dpc = 1}.
#' @param obv_pep_efficacy NULL, numeric in \code{[0,1]}, or function(dpc). If NULL, uses
#'   \code{obv_pep_efficacy_from_dpc()}, which is cut to zero after 10 DPC.
#' @param obv_pep_target_class Character vector of offspring classes eligible for OBV PEP.
#'   Defaults to \code{"HCW"}.
#' @param obv_pep_target_locations Character vector of exposure settings eligible for OBV PEP.
#'   Defaults to \code{"hospital"} for HCW occupational exposures.
#' @param prob_hcw_cond_funeral_hcw Numeric between 0 and 1. Probability a funeral infection is an HCW.
#' @param prob_hcw_cond_funeral_genPop Numeric between 0 and 1. Probability a funeral infection is a GenPop
#'
#' @return A data.frame with one row per realised funeral offspring and columns:
#'   \code{id}, \code{parent_class}, \code{setting}, \code{time_infection}, \code{class}.
#'   Returns a 0-row data.frame if no unsafe funeral occurs or zero offspring.
#' @export
offspring_function_funeral <- function(

  ## *Parent* characteristics and properties
  parent_info,

  ## Offspring distribution for unsafe funeral event
  mn_offspring_funeral = NULL, # scalar or function(t): mean offspring at unsafe funeral (resolved at parent death time)
  overdisp_offspring_funeral = NULL, # overdispersion of the above number of offspring

  ## Timing of funeral infections (delay after outcome)
  Tg_shape_funeral = NULL, # gamma shape parameter for Tg distribution at funerals ### have high shape, high rate to get low variance ##
  Tg_rate_funeral = NULL,  #gamma rate parameter for Tg distribution at funerals

  ### efficacy of a safe funeral (thinning funeral offspring)
  safe_funeral_efficacy = NULL, ## efficacy of a safe burial in reducing transmission in a funeral setting

  ## Obeldesivir PEP for exposed HCWs
  obv_pep_enabled = FALSE,
  obv_pep_coverage = 0,
  obv_pep_adherence = 1,
  obv_pep_dpc = 1,
  obv_pep_efficacy = NULL,
  obv_pep_target_class = "HCW",
  obv_pep_target_locations = "hospital",

  ## HCW vs genPop at funeral
  prob_hcw_cond_funeral_hcw = NULL, ### probability that the unsafe funeral infector infects a HCW
  prob_hcw_cond_funeral_genPop = NULL
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

  # Other positive parameters (scalars). mn_offspring_funeral is handled
  # separately below because it may be a scalar OR a function(t).
  for (nm in c("overdisp_offspring_funeral",
               "Tg_shape_funeral", "Tg_rate_funeral")) {
    val <- get(nm, inherits = FALSE)
    if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0)
      stop(sprintf("`%s` must be a single positive numeric value.", nm), call. = FALSE)
  }

  # mn_offspring_funeral may be a single positive scalar or a function(t) of
  # absolute calendar time (resolved at the parent's death time, see Step 3).
  if (!is.function(mn_offspring_funeral)) {
    if (is.null(mn_offspring_funeral) || !is.numeric(mn_offspring_funeral) ||
        length(mn_offspring_funeral) != 1L || is.na(mn_offspring_funeral) ||
        mn_offspring_funeral <= 0)
      stop("`mn_offspring_funeral` must be a function(t) or a single positive numeric value.",
           call. = FALSE)
  }

  ## Keep safe funeral efficacy as a scalar here. The probability that a funeral
  ## is safe/unsafe is already resolved upstream in complete_offspring_info();
  ## this parameter only controls residual transmission if a funeral is labelled
  ## safe. Set safe_funeral_efficacy = 1 for fully transmission-blocking safe funerals.
  if (is.null(safe_funeral_efficacy) || !is.numeric(safe_funeral_efficacy) ||
      length(safe_funeral_efficacy) != 1L || is.na(safe_funeral_efficacy) ||
      safe_funeral_efficacy < 0 || safe_funeral_efficacy > 1) {
    stop("`safe_funeral_efficacy` must be a single numeric value in [0, 1].", call. = FALSE)
  }

  if (is.null(prob_hcw_cond_funeral_hcw) || !is.numeric(prob_hcw_cond_funeral_hcw) ||
      length(prob_hcw_cond_funeral_hcw) != 1L || is.na(prob_hcw_cond_funeral_hcw) ||
      prob_hcw_cond_funeral_hcw < 0 || prob_hcw_cond_funeral_hcw > 1)
    stop("`prob_hcw_cond_funeral_hcw` must be in [0, 1].", call. = FALSE)

  if (is.null(prob_hcw_cond_funeral_genPop) || !is.numeric(prob_hcw_cond_funeral_genPop) ||
      length(prob_hcw_cond_funeral_genPop) != 1L || is.na(prob_hcw_cond_funeral_genPop) ||
      prob_hcw_cond_funeral_genPop < 0 || prob_hcw_cond_funeral_genPop > 1)
    stop("`prob_hcw_cond_funeral_genPop` must be in [0, 1].", call. = FALSE)


  #########################################################################################
  ## Logic for whether a safe or unsafe funeral occurs
  #########################################################################################

  # Step 1: Ensure parent is dead. If parent survived, no funeral transmission
  if (!isTRUE(parent_died)) {
    return(empty_offspring_dataframe())
  }

  # Step 2: Information on whether the parent had an unsafe or safe funeral
  has_unsafe_funeral <- parent_funeral # as.logical(rbinom(n = 1, size = 1, prob = p_unsafe_funeral)) # prev: Bernoulli trial for determining whether an unsafe funeral occurs
  has_unsafe_funeral <- ifelse(has_unsafe_funeral == "unsafe", TRUE, FALSE)

  #########################################################################################
  ## Produce funeral offspring
  #########################################################################################

  # Step 3: Draw raw number of infections from NB at the funeral (assuming initially an unsafe one).
  #         mn_offspring_funeral may be time-varying. The funeral is a discrete event at the parent's
  #         death, so resolve transmissibility at the parent's absolute death (outcome) time rather
  #         than at their infection time.
  parent_death_time_absolute <- parent_info$time_infection_absolute + parent_time_to_outcome
  mn_offspring_funeral_t <- resolve_positive_time_varying(
    mn_offspring_funeral, parent_death_time_absolute, "mn_offspring_funeral")
  num_offspring_raw <- rnbinom(
    n    = 1,
    mu   = mn_offspring_funeral_t,
    size = overdisp_offspring_funeral
  )

  if (num_offspring_raw == 0L) {
    return(empty_offspring_dataframe())
  }

  # Step 4: Generate infection times = outcome time + Gamma distributed 'delay', typically with
  #         little variance to represent a singular point from which infections arose
  delay_funeral <- rgamma(n = num_offspring_raw, shape = Tg_shape_funeral, rate = Tg_rate_funeral)
  infection_times <- parent_time_to_outcome + delay_funeral

  # Step 5: Assign setting ("funeral") for all candidate offspring
  infection_settings <- rep("funeral", num_offspring_raw)

  # Step 6: Assign class (HCW or genPop) to each candidate offspring. Class is assigned BEFORE
  #         thinning to match the structure of offspring_function_genPop and offspring_function_hcw,
  #         even though funeral thinning is class-independent.
  offspring_class <- rep("genPop", num_offspring_raw)
  if (parent_class == "genPop") {
    flip_hcw <- as.logical(rbinom(n = num_offspring_raw, size = 1, prob = prob_hcw_cond_funeral_genPop))
  } else if (parent_class == "HCW") {
    flip_hcw <- as.logical(rbinom(n = num_offspring_raw, size = 1, prob = prob_hcw_cond_funeral_hcw))
  } else {
    stop("Step 6 of funeral offspring function is broken")
  }
  offspring_class[flip_hcw] <- "HCW"

  # Step 6b: Snapshot pre-thinning candidates for the two-phase OBV PEP gate.
  ## Funeral candidate exposure clock times are parent infection time + delay.
  pre_thinning <- list(
    infection_location      = infection_settings,
    offspring_class         = offspring_class,
    infection_time_absolute = parent_info$time_infection_absolute + infection_times
  )

  # Step 7: Thin if the funeral is a safe one (class-independent thinning by safe_funeral_efficacy).
  ##         OBV efficacy is applied separately, post-thinning, in Step 8.
  if (!has_unsafe_funeral) {
    keep_infection <- as.logical(rbinom(n = num_offspring_raw, size = 1, prob = 1 - safe_funeral_efficacy))
  } else {
    keep_infection <- rep(TRUE, num_offspring_raw)
  }

  # Step 8: OBV PEP gate. By default this does nothing for funeral exposures because
  #         obv_pep_target_locations = "hospital", but enabling e.g.
  #         obv_pep_target_locations = c("hospital", "funeral") turns it on here too.
  obv_gate <- apply_obv_pep_gate(
    pre_thinning             = pre_thinning,
    kept_indices             = which(keep_infection),
    obv_pep_enabled          = obv_pep_enabled,
    obv_pep_coverage     = obv_pep_coverage,
    obv_pep_adherence        = obv_pep_adherence,
    obv_pep_dpc              = obv_pep_dpc,
    obv_pep_efficacy         = obv_pep_efficacy,
    obv_pep_target_class     = obv_pep_target_class,
    obv_pep_target_locations = obv_pep_target_locations
  )

  ## Final realised set = kept AND not prevented by OBV.
  final_local        <- obv_gate$keep
  final_idx          <- which(keep_infection)[final_local]
  infection_times    <- infection_times[final_idx]
  infection_settings <- infection_settings[final_idx]
  offspring_class    <- offspring_class[final_idx]
  obv_metadata       <- obv_gate$metadata[final_local, , drop = FALSE]

  offspring_df <- data.frame(infection_location = infection_settings,
                             time_infection_relative = infection_times,
                             class = offspring_class,
                             obv_metadata,
                             stringsAsFactors = FALSE)
  attr(offspring_df, "obv_pep_num_treated") <- obv_gate$num_treated
  ## Stash the infections OBV prevented (no RNG drawn) so the caller can resolve
  ## their counterfactual would-be deaths after the simulation loop.
  attr(offspring_df, "obv_pep_prevented_info") <-
    extract_obv_prevented_info(pre_thinning, keep_infection, final_local)
  return(offspring_df)
}
