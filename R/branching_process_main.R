## Add in details on parameter arguments here
## Add in details about outputs etc here
## Basically all the required documentation for this to be
## functional inside a package

## Note: need to add descriptions of each of the function inputs, which are currently missing
## Note: currently the code has a mixture of distribution parameter inputs (e.g. generation time) and others
##       where there are actual parameters of the distribution inputs (e.g. Tg_shape_funeral and Tg_rate_funeral)
##       we should harmonise this at some point
## Note: we might make a function called "generate_seeding_case_attributes" or something that does everything that we do in step 2 currently
## Note: we should change the function so that 1s/0s of parent outcome are characters i.e. "death" / "recovery" explicitly
## Note: general note that we should be actively thinking about how to ensure we don't end up in weird edge cases where all
##       of our infections end up dying or recovering before they even need healthcare. that'll require us to be careful with how
##       we approach parameterising this (maybe we put a check in place??)
## Note: prob_hospitalised_hcw and prob_hospitalised_genPop - do they need to be made specific to the location of the infection as well?

#' Run a stochastic branching-process outbreak simulation
#'
#' Top-level entry point for fiber's filovirus branching-process model.
#' Iteratively generates offspring (community, hospital, and funeral
#' transmission) from active cases until the outbreak ends or the configured
#' final-size cap is reached. Time-varying scenario inputs (probabilities,
#' delay factors, IPC / ETU coverage, and mean-offspring transmissibility) can
#' be passed as scalars or functions of time produced by [make_time_varying()].
#'
#' @param mn_offspring_genPop Positive numeric or function(t). Mean of the Negative Binomial
#'   offspring distribution for general-population (genPop) parents; resolved at the parent's
#'   infection time. Scalar or a function of absolute calendar time (e.g. [make_time_varying()]).
#' @param overdisp_offspring_genPop Positive numeric. Negative Binomial size (overdispersion) for
#'   genPop parents.
#' @param Tg_shape_genPop Positive numeric. Shape of the Gamma generation-time distribution for
#'   genPop parents.
#' @param Tg_rate_genPop Positive numeric. Rate of the Gamma generation-time distribution for genPop
#'   parents (mean generation time = shape / rate).
#' @param mn_offspring_hcw Positive numeric or function(t). Mean Negative Binomial offspring
#'   distribution for healthcare-worker (HCW) parents; resolved at the parent's infection time.
#' @param overdisp_offspring_hcw Positive numeric. Negative Binomial size (overdispersion) for HCW
#'   parents.
#' @param Tg_shape_hcw Positive numeric. Shape of the Gamma generation-time distribution for HCW
#'   parents.
#' @param Tg_rate_hcw Positive numeric. Rate of the Gamma generation-time distribution for HCW
#'   parents.
#' @param mn_offspring_funeral Positive numeric or function(t). Mean Negative Binomial number of
#'   offspring at an unsafe funeral; resolved at the parent's death (outcome) time.
#' @param overdisp_offspring_funeral Positive numeric. Negative Binomial size (overdispersion) for
#'   funeral offspring.
#' @param Tg_shape_funeral Positive numeric. Shape of the Gamma outcome-to-funeral-infection delay
#'   distribution.
#' @param Tg_rate_funeral Positive numeric. Rate of the Gamma funeral-delay distribution (mean delay
#'   = shape / rate).
#' @param incubation_period Function(n). Random generator returning n incubation-period draws
#'   (infection to symptom onset).
#' @param onset_to_hospitalisation Function(n). Random generator returning n
#'   onset-to-hospitalisation delay draws.
#' @param hospitalisation_delay_factor Positive numeric or function(t). Multiplier applied to
#'   `onset_to_hospitalisation` draws; may vary with absolute calendar time. Defaults to 1.
#' @param onset_to_death Function(n). Random generator returning n onset-to-death delay draws.
#' @param onset_to_recovery Function(n). Random generator returning n onset-to-recovery delay draws.
#' @param hospitalisation_to_death Function(n). Random generator returning n
#'   hospitalisation-to-death delay draws.
#' @param hospitalisation_to_recovery Function(n). Random generator returning n
#'   hospitalisation-to-recovery delay draws.
#' @param prob_symptomatic Numeric in `[0, 1]`. Probability an infection is symptomatic.
#' @param prob_hospitalised_hcw Numeric in `[0, 1]` or function(t). P(hospitalised | symptomatic) for
#'   HCWs; resolved at symptom-onset time.
#' @param prob_hospitalised_genPop Numeric in `[0, 1]` or function(t). P(hospitalised | symptomatic)
#'   for genPop; resolved at symptom-onset time.
#' @param prob_death_comm Numeric in `[0, 1]`. P(death | symptomatic) for cases managed in the
#'   community.
#' @param prob_death_hosp Numeric in `[0, 1]`. P(death | symptomatic) for hospitalised cases; must be
#'   no greater than `prob_death_comm`.
#' @param prob_hcw_cond_genPop_comm Numeric in `[0, 1]`. Probability a community-located infection from
#'   a genPop parent is an HCW.
#' @param prob_hcw_cond_genPop_hospital Numeric in `[0, 1]`. Probability a hospital-located infection
#'   from a genPop parent is an HCW.
#' @param prob_hcw_cond_hcw_comm Numeric in `[0, 1]`. Probability a community-located infection from an
#'   HCW parent is an HCW.
#' @param prob_hcw_cond_hcw_hospital Numeric in `[0, 1]`. Probability a hospital-located infection from
#'   an HCW parent is an HCW.
#' @param prob_hospital_cond_hcw_preAdm Numeric in `[0, 1]`. Probability that an infection generated by
#'   an HCW parent before their own admission occurs in the hospital (while still working).
#' @param ppe_efficacy_hcw Numeric in `[0, 1]` or function(t). Per-layer efficacy of PPE/IPC at
#'   reducing hospital transmission; resolved at each candidate hospital transmission time.
#' @param hospital_quarantine_efficacy Optional numeric in `[0, 1]` or function(t). Direct
#'   post-admission hospital quarantine/ETU efficacy. If NULL, derived from `prop_etu`, `ipc_helper`
#'   and `etu_efficacy_baseline`.
#' @param prop_etu Numeric in `[0, 1]` or function(t). Proportion of hospitalised cases in ETU/ETC
#'   care; used when `hospital_quarantine_efficacy` is NULL.
#' @param ipc_helper Numeric in `[0, 1]` or function(t). IPC/response-maturity proxy; used when
#'   `hospital_quarantine_efficacy` is NULL.
#' @param etu_efficacy_baseline Numeric in `[0, 1]`. Baseline ETU/ETC transmission-blocking efficacy
#'   before IPC-maturity adjustment; used when `hospital_quarantine_efficacy` is NULL.
#' @param obv_pep_enabled Logical scalar. If TRUE, apply the obeldesivir (OBV) PEP
#'   infection-prevention gate.
#' @param obv_pep_coverage Numeric in `[0, 1]` or function(t). Probability an eligible candidate
#'   receives OBV PEP.
#' @param obv_pep_adherence Numeric in `[0, 1]` or function(t). Probability a received OBV course is
#'   adhered to.
#' @param obv_pep_dpc Non-negative numeric or function(t). Days post challenge/exposure to first
#'   dose.
#' @param obv_pep_efficacy NULL, numeric in `[0, 1]`, or function(dpc). OBV efficacy; NULL uses
#'   [obv_pep_efficacy_from_dpc()].
#' @param obv_pep_target_class Character vector. Offspring classes eligible for OBV PEP (default
#'   "HCW").
#' @param obv_pep_target_locations Character vector. Exposure settings eligible for OBV PEP (default
#'   "hospital").
#' @param p_unsafe_funeral_comm_hcw Numeric in `[0, 1]` or function(t). Probability of an unsafe
#'   funeral after a community death of an HCW.
#' @param p_unsafe_funeral_hosp_hcw Numeric in `[0, 1]` or function(t). Probability of an unsafe
#'   funeral after a hospital death of an HCW.
#' @param p_unsafe_funeral_comm_genPop Numeric in `[0, 1]` or function(t). Probability of an unsafe
#'   funeral after a community death of a genPop case.
#' @param p_unsafe_funeral_hosp_genPop Numeric in `[0, 1]` or function(t). Probability of an unsafe
#'   funeral after a hospital death of a genPop case.
#' @param safe_funeral_efficacy Numeric in `[0, 1]`. Efficacy of a safe burial at preventing funeral
#'   transmission (1 = fully blocking).
#' @param prob_hcw_cond_funeral_hcw Numeric in `[0, 1]`. Probability a funeral infection from an HCW
#'   parent is an HCW.
#' @param prob_hcw_cond_funeral_genPop Numeric in `[0, 1]`. Probability a funeral infection from a
#'   genPop parent is an HCW.
#' @param population Positive integer. Total population size.
#' @param hcw_per_capita Positive numeric. HCWs per capita; total HCWs is
#'   `round(hcw_per_capita * population)`.
#' @param check_final_size Positive integer. Final-size cap; the simulation stops once this many
#'   cases have been generated.
#' @param initial_immune Non-negative integer. Number initially immune (removed from the susceptible
#'   pool). Defaults to 0.
#' @param seeding_cases Positive integer. Number of initial genPop community seeding infections.
#' @param susceptible_deplete Logical scalar. Placeholder for susceptible-depletion behaviour (not
#'   yet implemented). Defaults to FALSE.
#' @param seed Optional integer. RNG seed for reproducibility.
#'
#' @return A list with the simulated transmission tree (`$tdf`) and
#'   auxiliary outputs (final-size flag, generation index, etc.).
#'
#' @export
branching_process_main <- function(

  ## Transmission
  mn_offspring_genPop = NULL,               # scalar or function(t): mean offspring distribution for genPop (resolved at parent infection time)
  overdisp_offspring_genPop = NULL,         # overdispersion of the offspring distribution for genPop
  Tg_shape_genPop = NULL,                   # gamma shape parameter for Tg distribution for general population
  Tg_rate_genPop = NULL,                    # gamma rate parameter for Tg distribution for general population
  mn_offspring_hcw = NULL,                  # scalar or function(t): mean offspring distribution for HCWs (resolved at parent infection time)
  overdisp_offspring_hcw = NULL,            # overdispersion of the offspring distribution for HCWs
  Tg_shape_hcw = NULL,                      # gamma shape parameter for Tg distribution for HCWs
  Tg_rate_hcw = NULL,                       # gamma rate parameter for Tg distribution for HCWs
  mn_offspring_funeral = NULL,              # scalar or function(t): mean offspring at unsafe funeral (resolved at parent death time)
  overdisp_offspring_funeral = NULL,        # overdispersion of the above number of offspring
  Tg_shape_funeral = NULL,                  # gamma shape parameter for Tg distribution at funerals ### have high shape, high rate to get low variance ##
  Tg_rate_funeral = NULL,                   # gamma rate parameter for Tg distribution at funerals

  ## Natural history
  incubation_period,              # DESCRIPTION HERE
  onset_to_hospitalisation,    # DESCRIPTION HERE
  hospitalisation_delay_factor = 1.0,   # scalar or function(t): multiplier on onset_to_hospitalisation draws
  onset_to_death,
  onset_to_recovery,
  hospitalisation_to_death,              # Note: Jacob to look up whether the time -> death is the same typically as time -> recovery (or do they need to be different)
  hospitalisation_to_recovery,           # Note: Jacob to look up whether the time -> death is the same typically as time -> recovery (or do they need to be different)

  # Disease severity and healthcare seeking. Hospitalisation and death
  # probabilities are CONDITIONAL on the offspring being symptomatic.
  prob_symptomatic = NULL,           # P(symptomatic | infected)
  prob_hospitalised_hcw = NULL,      # scalar or function(t): P(hospitalised | symptomatic) for HCWs
  prob_hospitalised_genPop = NULL,   # scalar or function(t): P(hospitalised | symptomatic) for genPop
  prob_death_comm = NULL,            # P(die | symptomatic, community)
  prob_death_hosp = NULL,            # P(die | symptomatic, hospitalised); must be <= prob_death_comm

  ## Probabilities for genPop infecting either genPop or HCWs, depending on the setting
  prob_hcw_cond_genPop_comm = NULL,         # prob that a community-located infection generated by genPop is a HCW
  prob_hcw_cond_genPop_hospital = NULL,     # prob that a hospital-located infection generated by genPop is a HCW

  ## Probabilities for HCW infecting either genPop or HCWs, depending on the setting
  prob_hcw_cond_hcw_comm = NULL,            # prob that a community-located infection generated by HCW is a HCW
  prob_hcw_cond_hcw_hospital = NULL,        # prob that a hospital-located infection generated by HCW is a HCW

  ## Setting model for HCWs
  prob_hospital_cond_hcw_preAdm = NULL,     # probability that an infection generated prior to parent hospitaliation occurs in the hospital (whilst HCW is working)
  ppe_coverage_hcw = NULL,                  # scalar or function(t): coverage/probability that a relevant HCW has PPE (time-varying coverage lever)
  ppe_efficacy = NULL,                      # scalar: efficacy of PPE at preventing infection conditional on having it
  prop_etu = NULL,                          # scalar/function(t): proportion of hospitalised cases in ETU/ETC care (time-varying coverage lever)
  etu_efficacy = NULL,                      # scalar: post-admission quarantine efficacy for ETU/ETC care
  general_hospital_quarantine_efficacy = NULL,  # scalar: post-admission quarantine efficacy for general (non-ETU) hospital care

  ## Obeldesivir PEP. The gate is applied around the Swiss-cheese thinning step
  ## in each offspring function: treatment status (received, adherent, DPC) is
  ## assigned to all pre-thinning eligible candidates, and efficacy is applied
  ## only to candidates that also survive PPE/quarantine. See apply_obv_pep_gate().
  ## Per-individual treatment-status decisions compose multiplicatively with
  ## the existing PPE source / PPE receiver / hospital quarantine layers in the
  ## offspring functions; the NB overdispersion of HCW-only offspring counts is
  ## not preserved exactly under this thinning.
  ##
  ## Reproducibility caveat: when obv_pep_enabled = TRUE, the gate consumes
  ## additional rbinom() draws (coverage, adherence, prevention) per offspring
  ## call. This means OBV-enabled and OBV-disabled simulations with the same
  ## `seed` will NOT produce identical genPop / hospital / funeral trajectories
  ## -- they diverge as soon as the first OBV draw fires. For paired
  ## counterfactual comparisons use many replicates with different seeds and
  ## compare distributions (see Scenario 2 in the OBV verification script).
  obv_pep_enabled = FALSE,                   # logical: apply OBV infection-prevention gate
  obv_pep_coverage = 0,                  # scalar/function(t): probability eligible candidate receives OBV
  obv_pep_adherence = 1,                     # scalar/function(t): probability received course is effectively adhered to
  obv_pep_dpc = 1,                           # scalar/function(t): days post challenge/exposure to first dose
  obv_pep_efficacy = NULL,                   # NULL/function(dpc)/scalar: efficacy; NULL uses obv_pep_efficacy_from_dpc()
  obv_pep_target_class = "HCW",              # character vector: offspring classes eligible for OBV PEP
  obv_pep_target_locations = "hospital",     # character vector: exposure settings eligible for OBV PEP

  ## Funeral occurrence
  p_unsafe_funeral_comm_hcw = NULL, ## scalar or function(t): probability of unsafe funeral after a community death, HCW
  p_unsafe_funeral_hosp_hcw = NULL, ## scalar or function(t): probability of unsafe funeral after a hospital death, HCW
  p_unsafe_funeral_comm_genPop = NULL, ## scalar or function(t): probability of unsafe funeral after a community death, genPop
  p_unsafe_funeral_hosp_genPop = NULL, ## scalar or function(t): probability of unsafe funeral after a hospital death, genPop
  safe_funeral_efficacy = NULL, ## efficacy of a safe burial in reducing transmission in a funeral setting

  ## HCW vs genPop at funeral
  prob_hcw_cond_funeral_hcw = NULL, ### probability that the unsafe funeral infector infects a HCW
  prob_hcw_cond_funeral_genPop = NULL, ## DESCRIPTION NEEDED HERE

  ## Misc
  population,
  hcw_per_capita = 10,
  check_final_size,
  initial_immune = 0,
  seeding_cases,
  susceptible_deplete = FALSE,  ## note - still need to add code around this as functionality
                                ## envisaging this will adapt mn_offspring to reflect susceptible depletion
  seed = NULL

) {

  ##################################################################
  ### Step 1: Set up everything we need for the simulation
  ##################################################################
  # Set seed for reproducibility
  set.seed(seed)

  ## Local helpers for resolving parameters that can be either scalars or
  ## functions of calendar time. These mirror the logic in complete_offspring_info
  ## because some time-varying parameters are needed before offspring are created.
  resolve_probability <- function(param, t, param_name) {
    value <- resolve_time_varying(param = param, t = t, param_name = param_name)
    if (any(value < 0 | value > 1)) {
      stop(sprintf("`%s` must resolve to value(s) in [0, 1].", param_name), call. = FALSE)
    }
    value
  }

  ## Post-admission hospital quarantine efficacy is always derived inside the
  ## offspring functions as a prop_etu(t) mixture of the ETU and general-hospital
  ## quarantine efficacies. Require the full set up front so missing inputs fail
  ## early with a clear message rather than failing downstream.
  if (is.null(prop_etu) ||
      is.null(etu_efficacy) ||
      is.null(general_hospital_quarantine_efficacy)) {
    stop(
      "Supply all of `prop_etu`, `etu_efficacy`, and `general_hospital_quarantine_efficacy`.",
      call. = FALSE
    )
  }

  ## Fixed scalar efficacies in [0, 1]. The time-varying response enters through the
  ## coverage levers (`ppe_coverage_hcw`, `prop_etu`), not the efficacies, which are
  ## deliberately scalar-only. The ETU and general-hospital efficacies are
  ## independently togglable; no ordering between them is enforced.
  validate_scalar_probability <- function(param, param_name) {
    if (!is.numeric(param) ||
        length(param) != 1L ||
        is.na(param) ||
        param < 0 ||
        param > 1) {
      stop(sprintf("`%s` must be a single numeric value between 0 and 1.", param_name),
           call. = FALSE)
    }
  }

  validate_scalar_probability(ppe_efficacy, "ppe_efficacy")
  validate_scalar_probability(etu_efficacy, "etu_efficacy")
  validate_scalar_probability(general_hospital_quarantine_efficacy,
                              "general_hospital_quarantine_efficacy")

  if (!is.logical(obv_pep_enabled) || length(obv_pep_enabled) != 1L || is.na(obv_pep_enabled)) {
    stop("`obv_pep_enabled` must be a single logical value.", call. = FALSE)
  }
  if (!is.character(obv_pep_target_class) || length(obv_pep_target_class) < 1L) {
    stop("`obv_pep_target_class` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(obv_pep_target_locations) || length(obv_pep_target_locations) < 1L) {
    stop("`obv_pep_target_locations` must be a non-empty character vector.", call. = FALSE)
  }
  if (isTRUE(obv_pep_enabled) &&
      is.numeric(obv_pep_coverage) && length(obv_pep_coverage) == 1L &&
      !is.na(obv_pep_coverage) && obv_pep_coverage == 0) {
    warning("`obv_pep_enabled = TRUE` but `obv_pep_coverage = 0`; the OBV PEP gate will be a no-op.",
            call. = FALSE)
  }

  ## OBV PEP per-call accumulator: 7 gate counters, see empty_obv_pep_num_treated().
  obv_num_treated <- empty_obv_pep_num_treated()
  ## Collected per-call snapshots of infections OBV prevented (no RNG drawn in
  ## the loop). Their counterfactual would-be deaths are resolved once, after the
  ## loop, to populate obv_num_treated$prevented_deaths without perturbing the
  ## simulated trajectory's RNG stream.
  obv_prevented_info_list <- list()
  ##################################################################
  ### Step 1b: Upfront sanity checks on time-varying parameters
  ###
  ### Sample each probability parameter on a grid built from the
  ### make_time_varying breakpoints (where supplied) plus midpoints, or
  ### a default 0..365 grid otherwise. This catches curves that resolve
  ### outside [0, 1] somewhere in the simulation horizon before the
  ### simulation starts, rather than mid-run with a cryptic error.
  ###
  ### Also enforces the cross-parameter constraint that
  ###     prob_death_hosp <= prob_death_comm
  ### so the second-chance ratio prob_death_hosp / prob_death_comm
  ### stays in [0, 1] (see complete_offspring_info Step 3.2).
  ##################################################################
  sanity_params <- list(
    prob_hospitalised_hcw         = prob_hospitalised_hcw,
    prob_hospitalised_genPop      = prob_hospitalised_genPop,
    p_unsafe_funeral_comm_hcw     = p_unsafe_funeral_comm_hcw,
    p_unsafe_funeral_hosp_hcw     = p_unsafe_funeral_hosp_hcw,
    p_unsafe_funeral_comm_genPop  = p_unsafe_funeral_comm_genPop,
    p_unsafe_funeral_hosp_genPop  = p_unsafe_funeral_hosp_genPop,
    ppe_coverage_hcw              = ppe_coverage_hcw,
    prop_etu                      = prop_etu,
    obv_pep_coverage          = obv_pep_coverage,
    obv_pep_adherence             = obv_pep_adherence
  )
  ## Build the sampling grid from ALL time-varying inputs -- the probabilities
  ## above plus the positive-valued curves below -- so the upfront check lands on
  ## every changepoint, including those of hospitalisation_delay_factor,
  ## obv_pep_dpc, and the time-varying mean-offspring parameters.
  grid_inputs <- c(
    sanity_params,
    list(
      hospitalisation_delay_factor = hospitalisation_delay_factor,
      obv_pep_dpc                  = obv_pep_dpc,
      mn_offspring_genPop          = mn_offspring_genPop,
      mn_offspring_hcw             = mn_offspring_hcw,
      mn_offspring_funeral         = mn_offspring_funeral
    )
  )
  sanity_grid <- build_sanity_grid(grid_inputs)

  for (nm in names(sanity_params)) {
    check_probability_on_grid(sanity_params[[nm]], sanity_grid, nm)
  }

  ## hospitalisation_delay_factor is strictly positive (a multiplier), not a probability.
  check_positive_on_grid(hospitalisation_delay_factor, sanity_grid,
                         "hospitalisation_delay_factor")

  ## mn_offspring_* are strictly positive NB means and may be scalars or functions
  ## of absolute calendar time. They are resolved inside the offspring functions
  ## (genPop/HCW at the parent's infection time, funeral at the parent's death
  ## time); here we only sanity-check positivity across the simulation horizon.
  check_positive_on_grid(mn_offspring_genPop,  sanity_grid, "mn_offspring_genPop")
  check_positive_on_grid(mn_offspring_hcw,     sanity_grid, "mn_offspring_hcw")
  check_positive_on_grid(mn_offspring_funeral, sanity_grid, "mn_offspring_funeral")

  ## obv_pep_dpc is non-negative (0 = same-day treatment is a meaningful boundary value).
  check_nonneg_on_grid(obv_pep_dpc, sanity_grid, "obv_pep_dpc")

  ## prob_death_hosp must not exceed prob_death_comm (so second_chance_death_prob <= 1).
  ## Currently both are scalars; if they become time-varying in future, this still
  ## works because resolve_time_varying recycles scalars across the grid.
  if (!is.null(prob_death_comm) && !is.null(prob_death_hosp)) {
    pdc <- resolve_time_varying(prob_death_comm, sanity_grid, "prob_death_comm")
    pdh <- resolve_time_varying(prob_death_hosp, sanity_grid, "prob_death_hosp")
    if (any(pdh > pdc)) {
      bad <- which(pdh > pdc)
      show <- bad[seq_len(min(3L, length(bad)))]
      stop(sprintf(
        "`prob_death_hosp` must be <= `prob_death_comm` at all times (so the second-chance ratio is a valid probability); violation at t = %s (prob_death_hosp = %s, prob_death_comm = %s).",
        paste(round(sanity_grid[show], 3), collapse = ", "),
        paste(round(pdh[show], 4), collapse = ", "),
        paste(round(pdc[show], 4), collapse = ", ")
      ), call. = FALSE)
    }
  }


  ## Initialise the susceptible population
  susc <- population - initial_immune

  ## Initialise the HCW population
  hcw_total <- round(hcw_per_capita * population)
  if (hcw_total <= 0) {
    stop("number of hcws is <= 0 as currently specified by hcw_per_capita and population")
  }
  hcw_available <- hcw_total

  ## Preallocate data frame -
  max_cases <- check_final_size
  tdf <- data.frame(
    id                             = integer(max_cases),   # id of the infected individual
    class                          = NA_character_,
    infection_location             = NA_character_,
    parent                         = integer(max_cases),   # ancestor of the infected individual i.e. the parent
    generation                     = integer(max_cases),   # generation of the infected individual i.e. how many infections precede them in the transmission chain
    time_infection_relative        = NA_real_,             # time of the infection relative to the parent
    time_infection_absolute        = NA_real_,             # time of the infection in absolute calendar time (i.e. since start of outbreak)
    incubation_period              = NA_real_,
    symptomatic                    = NA,                   # are they symptomatic?
    time_symptom_onset_relative    = NA_real_,             # time of symptom onset relative to the parent
    time_symptom_onset_absolute    = NA_real_,             # time of symptom onset in absolute calendar time (i.e. since start of outbreak)
    hospitalisation                = FALSE,
    time_hospitalisation_relative  = NA_real_,
    time_hospitalisation_absolute  = NA_real_,
    outcome                        = FALSE,         # what the outcome is for that individual
    outcome_location               = NA_character_,
    time_outcome_relative          = NA_real_,
    time_outcome_absolute          = NA_real_,
    funeral_safety                 = NA_character_,
    obv_pep_eligible               = rep(FALSE, max_cases),
    obv_pep_received               = rep(FALSE, max_cases),
    obv_pep_adherent               = rep(FALSE, max_cases),
    obv_pep_dpc                    = rep(NA_real_, max_cases),
    n_offspring                    = integer(max_cases),
    offspring_generated            = FALSE,
    stringsAsFactors = FALSE
  )

  #########################################################################
  ### Step 2: Initialise conditions and features of the seeding cases
  #########################################################################

  ## Deciding whether the seeding cases are symptomatic, and if so, when they develop symptoms
  seeding_cases_time_infection <- seq(from = 0, to = 0.01, length.out = seeding_cases)
  seeding_cases_incubation <- incubation_period(n = seeding_cases)
  seeding_cases_symptomatic <- as.logical(rbinom(n = seeding_cases, size = 1, prob = prob_symptomatic))
  seeding_cases_symptom_onset <- rep(NA_real_, seeding_cases)
  seeding_cases_symptom_onset[seeding_cases_symptomatic] <- seeding_cases_incubation[seeding_cases_symptomatic]
  seeding_cases_symptom_onset_absolute <- seeding_cases_time_infection + seeding_cases_symptom_onset

  ## Deciding on the outcome for the seeding cases, and if so, when that outcome occurs
  seeding_cases_outcome <- rep(FALSE, seeding_cases)
  seeding_cases_outcome[seeding_cases_symptomatic] <- as.logical(rbinom(n = sum(seeding_cases_symptomatic), size = 1, prob = prob_death_comm))
  seeding_cases_outcome_time <- rep(NA_real_, seeding_cases)
  seeding_cases_outcome_time[seeding_cases_outcome] <- seeding_cases_incubation[seeding_cases_outcome] + onset_to_death(n = sum(seeding_cases_outcome))
  seeding_cases_outcome_time[!seeding_cases_outcome] <- seeding_cases_incubation[!seeding_cases_outcome] + onset_to_recovery(n = sum(!seeding_cases_outcome))
  seeding_cases_outcome_time_absolute <- seeding_cases_time_infection + seeding_cases_outcome_time

  ## Deciding funeral safety for seed cases who die. Seed cases are genPop and
  ## their deaths occur in the community, so use p_unsafe_funeral_comm_genPop at
  ## each seed case's death time rather than forcing all seed funerals unsafe.
  seeding_cases_funeral_safety <- rep(NA_character_, seeding_cases)
  if (any(seeding_cases_outcome)) {
    seeding_cases_p_unsafe_funeral <- resolve_probability(p_unsafe_funeral_comm_genPop,
                                                          seeding_cases_outcome_time_absolute[seeding_cases_outcome],
                                                          "p_unsafe_funeral_comm_genPop")
    seeding_cases_funeral_safety[seeding_cases_outcome] <- ifelse(
      rbinom(n = sum(seeding_cases_outcome), size = 1, prob = seeding_cases_p_unsafe_funeral) == 1,
      "unsafe", "safe"
    )
  }

  ## Initialising the dataframe with the seed cases and their attributes
  tdf[1:seeding_cases, ] <- data.frame(
    id                             = seq_len(seeding_cases),
    class                          = rep("genPop", seeding_cases),
    infection_location             = rep("community", seeding_cases),
    parent                         = NA_character_,
    generation                     = 1,
    time_infection_relative        = seeding_cases_time_infection,
    time_infection_absolute        = seeding_cases_time_infection,
    incubation_period              = seeding_cases_incubation,
    symptomatic                    = seeding_cases_symptomatic,
    time_symptom_onset_relative    = seeding_cases_symptom_onset,
    time_symptom_onset_absolute    = seeding_cases_symptom_onset_absolute,
    hospitalisation                = rep(FALSE, seeding_cases),
    time_hospitalisation_relative  = NA_real_,
    time_hospitalisation_absolute  = NA_real_,
    outcome                        = seeding_cases_outcome,
    outcome_location               = rep("community", seeding_cases),
    time_outcome_relative          = seeding_cases_outcome_time,
    time_outcome_absolute          = seeding_cases_outcome_time_absolute,
    funeral_safety                 = seeding_cases_funeral_safety,
    obv_pep_eligible               = rep(FALSE, seeding_cases),
    obv_pep_received               = rep(FALSE, seeding_cases),
    obv_pep_adherent               = rep(FALSE, seeding_cases),
    obv_pep_dpc                    = rep(NA_real_, seeding_cases),
    n_offspring                    = NA_integer_,
    offspring_generated            = FALSE,
    stringsAsFactors = FALSE
  )

  #################################################################################
  ### Step 3: Loop through infections and generate offspring for each of them
  #################################################################################
  ## While we haven't hit the simulation cap size (check_final_size) and any infections exist where we have not yet generated the requisite offspring,
  ## continue to generate infections
  while (any(is.na(tdf$n_offspring)) && susc > 0 && nrow(tdf) <= check_final_size) {

    #############################################################################################
    ## Step 1: Get earliest infection not yet expanded to act as a parent, and their attributes
    #############################################################################################
    parent_time_infection <- min(tdf$time_infection_absolute[!tdf$offspring_generated & !is.na(tdf$time_infection_absolute)])
    idx <- which(tdf$time_infection_absolute == parent_time_infection & !tdf$offspring_generated)[1]
    parent_info <- tdf[idx, ]
    if (!(parent_info$class %in% c("genPop", "HCW"))) {
      stop("error with parent class")
    }

    ###################################################################################################################
    ### Step 2: Generate offspring associated with community and (if hospitalised) healthcare associated transmission
    ###################################################################################################################
    ## Pass scalar-or-time-varying response parameters into the offspring functions
    ## directly. Those functions know the candidate transmission times, so they can
    ## resolve PPE coverage and post-admission hospital quarantine/ETU efficacy at the
    ## actual absolute calendar time of each candidate hospital exposure. PPE thinning
    ## is ppe_coverage_hcw(t) * ppe_efficacy; hospital quarantine efficacy is
    ## calculated inside the offspring functions as a prop_etu(t)-weighted mixture
    ## of etu_efficacy and general_hospital_quarantine_efficacy.

    if (parent_info$class == "genPop") {
      offspring_community_healthcare_df <- offspring_function_genPop(parent_info = parent_info,
                                                                     mn_offspring_genPop = mn_offspring_genPop,
                                                                     overdisp_offspring_genPop = overdisp_offspring_genPop,
                                                                     Tg_shape_genPop = Tg_shape_genPop,
                                                                     Tg_rate_genPop = Tg_rate_genPop,
                                                                     prop_etu = prop_etu,
                                                                     etu_efficacy = etu_efficacy,
                                                                     general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
                                                                     obv_pep_enabled = obv_pep_enabled,
                                                                     obv_pep_coverage = obv_pep_coverage,
                                                                     obv_pep_adherence = obv_pep_adherence,
                                                                     obv_pep_dpc = obv_pep_dpc,
                                                                     obv_pep_efficacy = obv_pep_efficacy,
                                                                     obv_pep_target_class = obv_pep_target_class,
                                                                     obv_pep_target_locations = obv_pep_target_locations,
                                                                     ppe_coverage_hcw = ppe_coverage_hcw,
                                                                     ppe_efficacy = ppe_efficacy,
                                                                     prob_hcw_cond_genPop_comm = prob_hcw_cond_genPop_comm,
                                                                     prob_hcw_cond_genPop_hospital = prob_hcw_cond_genPop_hospital)
    } else if (parent_info$class == "HCW") {
      ## PPE/IPC efficacy is passed through unresolved. The HCW offspring
      ## function resolves it separately for each candidate pre-admission
      ## hospital transmission event.

      offspring_community_healthcare_df <- offspring_function_hcw(parent_info = parent_info,
                                                                  mn_offspring_hcw = mn_offspring_hcw,
                                                                  overdisp_offspring_hcw = overdisp_offspring_hcw,
                                                                  Tg_shape_hcw = Tg_shape_hcw,
                                                                  Tg_rate_hcw = Tg_rate_hcw,
                                                                  prob_hospital_cond_hcw_preAdm = prob_hospital_cond_hcw_preAdm,
                                                                  ppe_coverage_hcw = ppe_coverage_hcw,
                                                                  ppe_efficacy = ppe_efficacy,
                                                                  prop_etu = prop_etu,
                                                                  etu_efficacy = etu_efficacy,
                                                                  general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
                                                                  obv_pep_enabled = obv_pep_enabled,
                                                                  obv_pep_coverage = obv_pep_coverage,
                                                                  obv_pep_adherence = obv_pep_adherence,
                                                                  obv_pep_dpc = obv_pep_dpc,
                                                                  obv_pep_efficacy = obv_pep_efficacy,
                                                                  obv_pep_target_class = obv_pep_target_class,
                                                                  obv_pep_target_locations = obv_pep_target_locations,
                                                                  prob_hcw_cond_hcw_comm = prob_hcw_cond_hcw_comm,
                                                                  prob_hcw_cond_hcw_hospital = prob_hcw_cond_hcw_hospital)
    }

    ## Accumulate OBV PEP counters from this offspring call (community/healthcare).
    step <- attr(offspring_community_healthcare_df, "obv_pep_num_treated", exact = TRUE)
    obv_num_treated$pre_eligible  <- obv_num_treated$pre_eligible  + step$pre_eligible
    obv_num_treated$pre_treated   <- obv_num_treated$pre_treated   + step$pre_treated
    obv_num_treated$pre_adherent  <- obv_num_treated$pre_adherent  + step$pre_adherent
    obv_num_treated$post_eligible <- obv_num_treated$post_eligible + step$post_eligible
    obv_num_treated$post_treated  <- obv_num_treated$post_treated  + step$post_treated
    obv_num_treated$post_adherent <- obv_num_treated$post_adherent + step$post_adherent
    obv_num_treated$prevented     <- obv_num_treated$prevented     + step$prevented
    pinfo <- attr(offspring_community_healthcare_df, "obv_pep_prevented_info", exact = TRUE)
    if (!is.null(pinfo) && nrow(pinfo) > 0) {
      obv_prevented_info_list[[length(obv_prevented_info_list) + 1L]] <- pinfo
    }

    ## Count HCWs generated and subtract from available pool
    n_hcw_community_healthcare <- sum(offspring_community_healthcare_df$class == "HCW")
    hcw_available <- hcw_available - n_hcw_community_healthcare

    #############################################################################################
    ### Step 3: Generate offspring associated with funeral transmission
    #############################################################################################
    offspring_funeral_df <- offspring_function_funeral(parent_info = parent_info,
                                                       mn_offspring_funeral = mn_offspring_funeral,
                                                       overdisp_offspring_funeral = overdisp_offspring_funeral,
                                                       Tg_shape_funeral = Tg_shape_funeral,
                                                       Tg_rate_funeral = Tg_rate_funeral,
                                                       safe_funeral_efficacy = safe_funeral_efficacy,
                                                       obv_pep_enabled = obv_pep_enabled,
                                                       obv_pep_coverage = obv_pep_coverage,
                                                       obv_pep_adherence = obv_pep_adherence,
                                                       obv_pep_dpc = obv_pep_dpc,
                                                       obv_pep_efficacy = obv_pep_efficacy,
                                                       obv_pep_target_class = obv_pep_target_class,
                                                       obv_pep_target_locations = obv_pep_target_locations,
                                                       prob_hcw_cond_funeral_hcw = prob_hcw_cond_funeral_hcw,
                                                       prob_hcw_cond_funeral_genPop = prob_hcw_cond_funeral_genPop)

    ## Accumulate OBV PEP counters from this offspring call (funeral).
    step <- attr(offspring_funeral_df, "obv_pep_num_treated", exact = TRUE)
    obv_num_treated$pre_eligible  <- obv_num_treated$pre_eligible  + step$pre_eligible
    obv_num_treated$pre_treated   <- obv_num_treated$pre_treated   + step$pre_treated
    obv_num_treated$pre_adherent  <- obv_num_treated$pre_adherent  + step$pre_adherent
    obv_num_treated$post_eligible <- obv_num_treated$post_eligible + step$post_eligible
    obv_num_treated$post_treated  <- obv_num_treated$post_treated  + step$post_treated
    obv_num_treated$post_adherent <- obv_num_treated$post_adherent + step$post_adherent
    obv_num_treated$prevented     <- obv_num_treated$prevented     + step$prevented
    pinfo <- attr(offspring_funeral_df, "obv_pep_prevented_info", exact = TRUE)
    if (!is.null(pinfo) && nrow(pinfo) > 0) {
      obv_prevented_info_list[[length(obv_prevented_info_list) + 1L]] <- pinfo
    }

    ## Update hcw_available after funeral transmission
    n_hcw_funeral <- sum(offspring_funeral_df$class == "HCW")
    hcw_available <- hcw_available - n_hcw_funeral

    #################################################################################################################
    ### Step 4: Complete offspring information based on parent attributes and timings; and update parent information
    ##          (e.g. num_offspring, offspring_generated == TRUE etc)
    #################################################################################################################
    ## Completing offspring information if there are any
    if (nrow(rbind(offspring_community_healthcare_df, offspring_funeral_df)) > 0) {
      complete_offspring_df <- complete_offspring_info(parent_info = parent_info,
                                                       offspring_dataframe = rbind(offspring_community_healthcare_df, offspring_funeral_df),
                                                       prob_symptomatic = prob_symptomatic,
                                                       prob_hospitalised_hcw = prob_hospitalised_hcw,
                                                       prob_hospitalised_genPop = prob_hospitalised_genPop,
                                                       prob_death_comm = prob_death_comm,
                                                       prob_death_hosp = prob_death_hosp,
                                                       p_unsafe_funeral_comm_hcw = p_unsafe_funeral_comm_hcw,
                                                       p_unsafe_funeral_hosp_hcw = p_unsafe_funeral_hosp_hcw,
                                                       p_unsafe_funeral_comm_genPop = p_unsafe_funeral_comm_genPop,
                                                       p_unsafe_funeral_hosp_genPop = p_unsafe_funeral_hosp_genPop,
                                                       incubation_period = incubation_period,
                                                       onset_to_hospitalisation = onset_to_hospitalisation,
                                                       hospitalisation_delay_factor = hospitalisation_delay_factor,
                                                       hospitalisation_to_death = hospitalisation_to_death,
                                                       hospitalisation_to_recovery = hospitalisation_to_recovery,
                                                       onset_to_death = onset_to_death,
                                                       onset_to_recovery = onset_to_recovery)
      tdf$n_offspring[idx] <- nrow(complete_offspring_df)
    } else {
      complete_offspring_df <- tdf[0, , drop = FALSE]
      tdf$n_offspring[idx] <- 0
    }
    tdf$offspring_generated[idx] <- TRUE

    #################################################################################################################
    ### Step 5: Adding the complete offspring dataframe (complete_offspring_df) to the main dataframe (tdf)
    #################################################################################################################
    ## If offspring exist, append them
    if (nrow(complete_offspring_df) > 0) {
      current_max_row <- max(which(!is.na(tdf$time_infection_absolute)))
      current_max_id <- max(tdf$id[which(!is.na(tdf$time_infection_absolute))])
      complete_offspring_df$id <- (current_max_id + 1):(current_max_id + nrow(complete_offspring_df))
      tdf[(current_max_row + 1):(current_max_row + nrow(complete_offspring_df)), ] <- complete_offspring_df[, names(tdf), drop = FALSE]
    }
    ## Deplete susceptibles
    susc <- susc - tdf$n_offspring[idx]
  }

  ############################################################################################
  ### Final tidy of dataframe and then outputting it
  ############################################################################################
  tdf <- tdf[order(tdf$time_infection_absolute, tdf$id), ]
  rownames(tdf) <- NULL

  #########################################################################################
  ### Deferred OBV PEP "prevented deaths" counterfactual
  ###
  ### For each infection the OBV gate prevented, decide whether it WOULD have died had
  ### it occurred, by replaying it through the SAME outcome model as realised cases
  ### (complete_offspring_info: symptomatic -> potential hospitalisation -> community CFR
  ### -> hospital second-chance). This is a "direct" count -- the would-be death of each
  ### prevented index infection only -- mirroring `prevented`, which likewise excludes the
  ### averted onward transmission chains.
  ###
  ### Why after the loop: these draws consume RNG, so doing them inline would shift every
  ### subsequent draw and change the simulated trajectory. Run once the tree is finalised,
  ### nothing downstream in this call depends on them, and each branching_process_main()
  ### call re-seeds up front -- so the trajectory and every existing output are byte-for-byte
  ### identical to a run without this counter. Skipped entirely when nothing was prevented,
  ### so zero-prevention (incl. obv-disabled) runs draw nothing extra.
  #########################################################################################
  if (length(obv_prevented_info_list) > 0) {
    obv_prevented_info <- do.call(rbind, obv_prevented_info_list)
    if (nrow(obv_prevented_info) > 0) {
      ## Carry each prevented infection's absolute infection time as a "relative" time
      ## against a zero-time dummy parent, so complete_offspring_info's clock
      ## (parent_abs + relative) resolves time-varying parameters at the true time.
      prevented_offspring <- data.frame(
        infection_location      = obv_prevented_info$infection_location,
        time_infection_relative = obv_prevented_info$time_infection_absolute,
        class                   = obv_prevented_info$class,
        stringsAsFactors        = FALSE
      )
      dummy_parent <- data.frame(
        id                      = NA_integer_,
        generation              = NA_integer_,
        time_infection_absolute = 0,
        stringsAsFactors        = FALSE
      )
      prevented_completed <- complete_offspring_info(
        parent_info                  = dummy_parent,
        offspring_dataframe          = prevented_offspring,
        prob_symptomatic             = prob_symptomatic,
        prob_hospitalised_hcw        = prob_hospitalised_hcw,
        prob_hospitalised_genPop     = prob_hospitalised_genPop,
        prob_death_comm              = prob_death_comm,
        prob_death_hosp              = prob_death_hosp,
        p_unsafe_funeral_comm_hcw    = p_unsafe_funeral_comm_hcw,
        p_unsafe_funeral_hosp_hcw    = p_unsafe_funeral_hosp_hcw,
        p_unsafe_funeral_comm_genPop = p_unsafe_funeral_comm_genPop,
        p_unsafe_funeral_hosp_genPop = p_unsafe_funeral_hosp_genPop,
        incubation_period            = incubation_period,
        onset_to_hospitalisation     = onset_to_hospitalisation,
        hospitalisation_delay_factor = hospitalisation_delay_factor,
        hospitalisation_to_death     = hospitalisation_to_death,
        hospitalisation_to_recovery  = hospitalisation_to_recovery,
        onset_to_death               = onset_to_death,
        onset_to_recovery            = onset_to_recovery
      )
      obv_num_treated$prevented_deaths <- sum(prevented_completed$outcome, na.rm = TRUE)
    }
  }

  attr(tdf, "hcw_total") <- hcw_total
  attr(tdf, "hcw_infected") <- hcw_total - hcw_available
  attr(tdf, "hcw_remaining") <- hcw_available
  attr(tdf, "obv_pep_num_treated") <- obv_num_treated

  out <- list(
    tdf = tdf,
    sim_info = list(
      population          = population,
      hcw_per_capita      = hcw_per_capita,
      hcw_total           = hcw_total,
      seed                = seed,
      obv_pep_enabled     = obv_pep_enabled,
      obv_pep_num_treated = obv_num_treated
    )
  )

  return(out)
}

