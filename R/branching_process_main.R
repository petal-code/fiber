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

mn_offspring_genPop = 3
overdisp_offspring_genPop = 3
Tg_shape_genPop = 20
Tg_rate_genPop = 2
mn_offspring_hcw = 3
overdisp_offspring_hcw = 3
Tg_shape_hcw = 20
Tg_rate_hcw = 2
mn_offspring_funeral = 5
overdisp_offspring_funeral = 1
Tg_shape_funeral = 1000
Tg_rate_funeral = 100
incubation_period = function(n) { rgamma(n = n, shape = 5, rate = 2)}
onset_to_hospitalisation = function(n) { rgamma(n = n, shape = 20, rate = 2) }
onset_to_death = function(n) { rgamma(n = n, shape = 40, rate = 2)}
onset_to_recovery = function(n) { rgamma(n = n, shape = 40, rate = 2)}
hospitalisation_to_death = function(n) { rgamma(n = n, shape = 20, rate = 2)}
hospitalisation_to_recovery = function(n) { rgamma(n = n, shape = 20, rate = 2)}
prob_symptomatic = 0.8
prob_hospitalised_hcw = 0.7
prob_hospitalised_genPop = 0.7
prob_death_comm = 0.8
prob_death_hosp = 0.2
prob_hcw_cond_genPop_comm = 0.5
prob_hcw_cond_genPop_hospital = 0.5
prob_hcw_cond_hcw_comm = 0.5
prob_hcw_cond_hcw_hospital = 0.5
prob_hospital_cond_hcw_preAdm = 0.5
ppe_efficacy_hcw = 0.5
hospital_quarantine_efficacy = 0.5
p_unsafe_funeral_comm_hcw = 0.75
p_unsafe_funeral_hosp_hcw = 0.75
p_unsafe_funeral_comm_genPop = 0.75
p_unsafe_funeral_hosp_genPop = 0.75
safe_funeral_efficacy = 0.9
prob_hcw_cond_funeral_hcw = 0.5
prob_hcw_cond_funeral_genPop = 0.5
tf = Inf
population = 1e6
check_final_size = 10000
initial_immune = 0
seeding_cases = 10
susceptible_deplete = FALSE
seed = 400
source("R/complete_offspring_info.R")
source("R/offspring_function_funeral.R")
source("R/offspring_function_genPop.R")
source("R/offspring_function_hcw.R")
source("R/helper_functions.R")

branching_process_main <- function(

  ## Transmission
  mn_offspring_genPop = NULL,               # mean of the offspring distribution for genPop
  overdisp_offspring_genPop = NULL,         # overdispersion of the offspring distribution for genPop
  Tg_shape_genPop = NULL,                   # gamma shape parameter for Tg distribution for general population
  Tg_rate_genPop = NULL,                    # gamma rate parameter for Tg distribution for general population
  mn_offspring_hcw = NULL,                  # mean of the offspring distribution for HCWs
  overdisp_offspring_hcw = NULL,            # overdispersion of the offspring distribution for HCWs
  Tg_shape_hcw = NULL,                      # gamma shape parameter for Tg distribution for HCWs
  Tg_rate_hcw = NULL,                       # gamma rate parameter for Tg distribution for HCWs
  mn_offspring_funeral = NULL,              # mean number of offspring at unsafe funeral
  overdisp_offspring_funeral = NULL,        # overdispersion of the above number of offspring
  Tg_shape_funeral = NULL,                  # gamma shape parameter for Tg distribution at funerals ### have high shape, high rate to get low variance ##
  Tg_rate_funeral = NULL,                   # gamma rate parameter for Tg distribution at funerals

  ## Natural history
  incubation_period,              # DESCRIPTION HERE
  onset_to_hospitalisation,    # DESCRIPTION HERE
  onset_to_death,
  onset_to_recovery,
  hospitalisation_to_death,              # Note: Jacob to look up whether the time -> death is the same typically as time -> recovery (or do they need to be different)
  hospitalisation_to_recovery,           # Note: Jacob to look up whether the time -> death is the same typically as time -> recovery (or do they need to be different)

  # Disease severity and healthcare seeking
  prob_symptomatic = NULL,
  prob_hospitalised_hcw = NULL,      # Do these need to be made both individual (hcw vs genPop) and context of where the infection originates (community vs healthcare) specific??
  prob_hospitalised_genPop = NULL,   # Do these need to be made both individual (hcw vs genPop) and context of where the infection originates (community vs healthcare) specific??
  prob_death_comm = NULL,
  prob_death_hosp = NULL,

  ## Probabilities for genPop infecting either genPop or HCWs, depending on the setting
  prob_hcw_cond_genPop_comm = NULL,         # prob that a community-located infection generated by genPop is a HCW
  prob_hcw_cond_genPop_hospital = NULL,     # prob that a hospital-located infection generated by genPop is a HCW

  ## Probabilities for HCW infecting either genPop or HCWs, depending on the setting
  prob_hcw_cond_hcw_comm = NULL,            # prob that a community-located infection generated by HCW is a HCW
  prob_hcw_cond_hcw_hospital = NULL,        # prob that a hospital-located infection generated by HCW is a HCW

  ## Setting model for HCWs
  prob_hospital_cond_hcw_preAdm = NULL,     # probability that an infection generated prior to parent hospitaliation occurs in the hospital (whilst HCW is working)
  ppe_efficacy_hcw = NULL,                  # efficacy of PPE/IPC measures at reducing transmission (i.e. pre hospitalisation)
  hospital_quarantine_efficacy = NULL,      # efficacy of quarantine at reducing transmission (i.e. post hospitalisation)

  ## Funeral occurrence
  p_unsafe_funeral_comm_hcw = NULL, ## probability the parent had an unsafe funeral Cond. on a community death and parent class being HCW
  p_unsafe_funeral_hosp_hcw = NULL, ## probability parent had an unsafe funeral cond. on a hospital death and parent class being HCW
  p_unsafe_funeral_comm_genPop = NULL, ## probability the parent had an unsafe funeral Cond. on a community death and parent class being genPop
  p_unsafe_funeral_hosp_genPop = NULL, ## probability parent had an unsafe funeral cond. on a hospital death and parent class being genPop
  safe_funeral_efficacy = NULL, ## efficacy of a safe burial in reducing transmission in a funeral setting

  ## HCW vs genPop at funeral
  prob_hcw_cond_funeral_hcw = NULL, ### probability that the unsafe funeral infector infects a HCW
  prob_hcw_cond_funeral_genPop = NULL, ## DESCRIPTION NEEDED HERE

  ## Misc
  tf = Inf,
  population,
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

  ## Initialise the susceptible population
  susc <- population - initial_immune

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
    n_offspring                    = integer(max_cases),
    offspring_generated            = FALSE,
    stringsAsFactors = FALSE
  )

  #########################################################################
  ### Step 2: Initialise conditions and features of the seeding cases
  #########################################################################

  ## Deciding whether the seeding cases are symptomatic, and if so, when they develop symptoms
  seeding_cases_incubation <- incubation_period(n = seeding_cases)
  seeding_cases_symptomatic <- as.logical(rbinom(n = seeding_cases, size = 1, prob = prob_symptomatic))
  seeding_cases_symptom_onset <- rep(NA_real_, seeding_cases)
  seeding_cases_symptom_onset[seeding_cases_symptomatic] <- seeding_cases_incubation[seeding_cases_symptomatic]

  ## Deciding on the outcome for the seeding cases, and if so, when that outcome occurs
  seeding_cases_outcome <- rep(FALSE, seeding_cases)
  seeding_cases_outcome[seeding_cases_symptomatic] <- as.logical(rbinom(n = sum(seeding_cases_symptomatic), size = 1, prob = prob_death_comm))
  seeding_cases_outcome_time <- rep(NA_real_, seeding_cases)
  seeding_cases_outcome_time[seeding_cases_outcome] <- seeding_cases_incubation[seeding_cases_outcome] + onset_to_death(n = sum(seeding_cases_outcome))
  seeding_cases_outcome_time[!seeding_cases_outcome] <- seeding_cases_incubation[!seeding_cases_outcome] + onset_to_recovery(n = sum(!seeding_cases_outcome))

  ## Deciding on the times of funerals for dead individuals (all of whom are assumed to have unsafe funerals)
  seeding_cases_funeral_safety <- rep("unsafe", seeding_cases)

  ## Initialising the dataframe with the seed cases and their attributes
  tdf[1:seeding_cases, ] <- data.frame(
    id                             = seq_len(seeding_cases),
    class                          = rep("genPop", seeding_cases),
    infection_location             = rep("community", seeding_cases),
    parent                         = NA_character_,
    generation                     = 1,
    time_infection_relative        = seq(from = 0, to = 0.01, length.out = seeding_cases),
    time_infection_absolute        = seq(from = 0, to = 0.01, length.out = seeding_cases),
    incubation_period              = seeding_cases_incubation,
    symptomatic                    = seeding_cases_symptomatic,
    time_symptom_onset_relative    = seeding_cases_symptom_onset,
    time_symptom_onset_absolute    = seeding_cases_symptom_onset,
    hospitalisation                = rep(FALSE, seeding_cases),
    time_hospitalisation_relative  = NA_real_,
    time_hospitalisation_absolute  = NA_real_,
    outcome                        = seeding_cases_outcome,
    outcome_location               = rep("community", seeding_cases),
    time_outcome_relative          = seeding_cases_outcome_time,
    time_outcome_absolute          = seeding_cases_outcome_time,
    funeral_safety                 = seeding_cases_funeral_safety,
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
    if (parent_info$class == "genPop") {
      offspring_community_healthcare_df <- offspring_function_genPop(parent_info = parent_info,
                                                                     mn_offspring_genPop = mn_offspring_genPop,
                                                                     overdisp_offspring_genPop = overdisp_offspring_genPop,
                                                                     Tg_shape_genPop = Tg_shape_genPop,
                                                                     Tg_rate_genPop = Tg_rate_genPop,
                                                                     hospital_quarantine_efficacy = hospital_quarantine_efficacy,
                                                                     prob_hcw_cond_genPop_comm = prob_hcw_cond_genPop_comm,
                                                                     prob_hcw_cond_genPop_hospital = prob_hcw_cond_genPop_hospital)
    } else if (parent_info$class == "HCW") {
      offspring_community_healthcare_df <- offspring_function_hcw(parent_info = parent_info,
                                                                  mn_offspring_hcw = mn_offspring_hcw,
                                                                  overdisp_offspring_hcw = overdisp_offspring_hcw,
                                                                  Tg_shape_hcw = Tg_shape_hcw,
                                                                  Tg_rate_hcw = Tg_rate_hcw,
                                                                  prob_hospital_cond_hcw_preAdm = prob_hospital_cond_hcw_preAdm,
                                                                  ppe_efficacy_hcw = ppe_efficacy_hcw,
                                                                  hospital_quarantine_efficacy = hospital_quarantine_efficacy,
                                                                  prob_hcw_cond_hcw_comm = prob_hcw_cond_hcw_comm,
                                                                  prob_hcw_cond_hcw_hospital = prob_hcw_cond_hcw_hospital)
    }

    #############################################################################################
    ### Step 3: Generate offspring associated with funeral transmission
    #############################################################################################
    offspring_funeral_df <- offspring_function_funeral(parent_info = parent_info,
                                                       mn_offspring_funeral = mn_offspring_funeral,
                                                       overdisp_offspring_funeral = overdisp_offspring_funeral,
                                                       Tg_shape_funeral = Tg_shape_funeral,
                                                       Tg_rate_funeral = Tg_rate_funeral,
                                                       safe_funeral_efficacy = safe_funeral_efficacy,
                                                       prob_hcw_cond_funeral_hcw = prob_hcw_cond_funeral_hcw,
                                                       prob_hcw_cond_funeral_genPop = prob_hcw_cond_funeral_genPop)

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
                                                       hospitalisation_to_death = hospitalisation_to_death,
                                                       hospitalisation_to_recovery = hospitalisation_to_recovery,
                                                       onset_to_death = onset_to_death,
                                                       onset_to_recovery = onset_to_recovery)
      tdf$n_offspring[idx] <- nrow(complete_offspring_df)
    } else {
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
  return(tdf)
}

