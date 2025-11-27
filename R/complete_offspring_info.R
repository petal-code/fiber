# consider turning individual parent information into a list at some point

offspring_dataframe <- offspring_df
parent_id <- 2
parent_generation <- 2
parent_infection_time = 10
prob_symptomatic = 0.75
infection_to_onset <- function(n) { rgamma(n = n, shape = 4, rate = 2)}
onset_to_hospitalisation <- function(n) { rgamma(n = n, shape = 10, rate = 2)}
infection_to_death <- function(n) { rgamma(n = n, shape = 20, rate = 2)}
infection_to_recovery <- function(n) { rgamma(n = n, shape = 30, rate = 2)}
hospitalisation_to_death <- function(n) {rgamma(n = n, shape = 50, rate = 1)}
hospitalisation_to_recovery <- function(n) {rgamma(n = n, shape = 100, rate = 1)}
prob_symptomatic = 0.5
prob_hospitalised_hcw = 0.5
prob_hospitalised_genPop = 0.5
prob_death_comm = 0.5
prob_death_hosp = 0.25

complete_offspring_info <- function(
    offspring_dataframe = NULL,
    parent_id = NULL,
    parent_generation = NULL,
    parent_infection_time = NULL,
    prob_symptomatic = NULL,
    incubation_period,
    onset_to_hospitalisation,
    hospitalisation_to_death,
    hospitalisation_to_recovery,
    onset_to_death,
    onset_to_recovery,
    prob_symptomatic = NULL,
    prob_hospitalised_hcw = NULL,
    prob_hospitalised_genPop = NULL,
    prob_death_comm = NULL,
    prob_death_hosp = NULL)
{

  ## Step 1: Assigning symptom status, hospitalisation and outcome status
  ## Deciding whether the offspring cases are symptomatic, and if so, when they develop symptoms
  offspring_cases_incubation <- incubation_period(n = num_offspring)
  offspring_cases_symptomatic <- as.logical(rbinom(n = num_offspring, size = 1, prob = prob_symptomatic))
  offspring_cases_symptom_onset <- rep(NA_real_, num_offspring)
  offspring_cases_symptom_onset[offspring_cases_symptomatic] <- offspring_cases_incubation[offspring_cases_symptomatic]

  ## Step 2: Deciding on the hospitalisation status for the (symptomatic) offspring cases, and if so, when that outcome occurs
  prob_hosp_given_symptoms_hcw <- prob_hosp_given_symptoms(prob_hosp = prob_hospitalised_hcw, prob_symptomatic = prob_symptomatic)
  prob_hosp_given_symptoms_genPop <- prob_hosp_given_symptoms(prob_hosp = prob_hospitalised_genPop, prob_symptomatic = prob_symptomatic)
  prob_hosp <- ifelse(offspring_class[offspring_cases_symptomatic] == "hcw", prob_hosp_given_symptoms_hcw, prob_hosp_given_symptoms_genPop)
  offspring_cases_potentially_hospitalised <- rep(FALSE, num_offspring)
  offspring_cases_potentially_hospitalised[offspring_cases_symptomatic] <- as.logical(rbinom(n = sum(offspring_cases_symptomatic), size = 1, prob = prob_hosp))
  offspring_cases_potentially_hospitalised_time <- rep(NA_real_, num_offspring)
  offspring_cases_potentially_hospitalised_time[offspring_cases_potentially_hospitalised] <- onset_to_hospitalisation(n = sum(offspring_cases_potentially_hospitalised))

  ## Step 3: Deciding on the outcome for the offspring cases, and if so, when that outcome occurs
  ## Note: at a later date, we need to ensure we're capturing the fact that only symptomatic people will die (i.e. correlated outcomes)

  ## 3.1 Start by assuming community and calculate which event occurs and when it occurs
  prob_death_given_symptoms_comm <- prob_death_given_symptoms(prob_death = prob_death_comm, prob_symptomatic = prob_symptomatic)
  temp_comm_outcome_death <- rep(FALSE, num_offspring)
  temp_comm_outcome_death[offspring_cases_symptomatic] <-  as.logical(rbinom(n = sum(offspring_cases_symptomatic), size = 1, prob = prob_death_given_symptoms_comm))
  temp_comm_outcome_time <- rep(NA_real_, num_offspring)
  temp_comm_outcome_time[temp_comm_outcome_death] <- onset_to_death(n = sum(temp_comm_outcome_death))
  temp_comm_outcome_time[!temp_comm_outcome_death] <- onset_to_recovery(n = sum(!temp_comm_outcome_death))

  ## 3.2 If community event earlier than hosp, keep the community event and associated times; if hospitalisation precedes
  ##     the community event, hospitalisation is successful. If community event is recovery, keep as recovery in hospital.
  ##     If community event was death, give another "chance" to be saved
  hosp_success <- offspring_cases_potentially_hospitalised_time[offspring_cases_potentially_hospitalised] < temp_comm_outcome_time[offspring_cases_potentially_hospitalised]
  second_chance_death_prob <- prob_death_hosp / prob_death_comm  ## which should this be  ## check which one of these this should be - unclear to me currently, but overarching idea remains
  # second_chance_death_prob <- prob_death_hosp / prob_death_given_symptoms_comm ## which should this be
  hosp_outcome <- temp_comm_outcome_death[offspring_cases_potentially_hospitalised][hosp_success]
  hosp_outcome[which(hosp_outcome)] <- as.logical(rbinom(n = sum(hosp_outcome), size = 1, prob = second_chance_death_prob))

  # Generating the dataframes for i) potentially hospitalised; ii) actually hospitalised; and iii) outcome
  offspring_cases_actually_hospitalised <- offspring_cases_potentially_hospitalised
  offspring_cases_actually_hospitalised[offspring_cases_potentially_hospitalised] <- hosp_success
  offspring_cases_actually_hospitalised_time <- rep(NA_real_, num_offspring)
  offspring_cases_actually_hospitalised_time[offspring_cases_actually_hospitalised] <- offspring_cases_potentially_hospitalised_time[offspring_cases_actually_hospitalised]

  offspring_cases_outcome <- temp_comm_outcome_death
  offspring_cases_outcome[offspring_cases_actually_hospitalised] <- hosp_outcome
  offspring_cases_outcome_time <- rep(NA_real_, num_offspring)
  offspring_cases_outcome_time[!offspring_cases_actually_hospitalised] <- temp_comm_outcome_time
  offspring_cases_outcome_time[offspring_cases_outcome & offspring_cases_actually_hospitalised] <- offspring_cases_potentially_hospitalised_time[offspring_cases_outcome & offspring_cases_actually_hospitalised] + hospitalisation_to_death(n = sum(offspring_cases_outcome & offspring_cases_actually_hospitalised))
  offspring_cases_outcome_time[!offspring_cases_outcome & offspring_cases_actually_hospitalised] <- offspring_cases_potentially_hospitalised_time[!offspring_cases_outcome & offspring_cases_actually_hospitalised] + hospitalisation_to_recovery(n = sum(!offspring_cases_outcome & offspring_cases_actually_hospitalised))

  offspring_cases_outcome_location <- rep("community", num_offspring)
  offspring_cases_outcome_location[offspring_cases_actually_hospitalised] <- "hospital"

  ## Step 4: Update and output offspring dataframe
  offspring_dataframe$parent <- parent_id
  offspring_dataframe$generation <- parent_generation + 1
  offspring_dataframe$time_infection_absolute <- parent_infection_time + offspring_dataframe$time_infection_relative
  offspring_dataframe$incubation_period <- offspring_cases_incubation
  offspring_dataframe$symptomatic <- offspring_cases_symptomatic
  offspring_dataframe$time_symptom_onset_relative <- offspring_cases_symptom_onset
  offspring_dataframe$time_symptom_onset_absolute <- offspring_dataframe$time_infection_absolute + offspring_cases_symptom_onset
  offspring_dataframe$potentially_hospitalised <- offspring_cases_potentially_hospitalised
  offspring_dataframe$hospitalisation <- offspring_cases_actually_hospitalised
  offspring_dataframe$time_hospitalisation_relative <- offspring_cases_incubation + offspring_cases_actually_hospitalised_time
  offspring_dataframe$time_hospitalisation_absolute <- offspring_dataframe$time_infection_absolute + offspring_dataframe$time_hospitalisation_relative
  offspring_dataframe$outcome <- offspring_cases_outcome
  offspring_dataframe$outcome_location <- offspring_cases_outcome_location
  offspring_dataframe$time_outcome_relative <- offspring_cases_incubation + offspring_cases_outcome_time
  offspring_dataframe$time_outcome_absolute <- offspring_dataframe$time_infection_absolute + offspring_dataframe$time_outcome_relative

  ## Step 5: Return completed dataframe
  return(offspring_dataframe)

}
