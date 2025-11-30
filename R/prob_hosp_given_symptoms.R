prob_hosp_given_symptoms <- function(prob_hosp, prob_symptomatic) {
  if (prob_hosp > prob_symptomatic) {
    stop("P(hospitalisation) cannot exceed P(symptomatic).", call. = FALSE)
  }
  if (prob_symptomatic == 0 & prob_hosp > 0) {
    stop("P(symptomatic) = 0 but P(hospitalisation) > 0: impossible.", call. = FALSE)
  }
  return(prob_hosp / prob_symptomatic)
}

prob_death_given_symptoms <- function(prob_death, prob_symptomatic) {
  if (prob_death > prob_symptomatic) {
    stop("P(death) cannot exceed P(symptomatic).", call. = FALSE)
  }
  if (prob_symptomatic == 0 & prob_death > 0) {
    stop("P(symptomatic) = 0 but P(death) > 0: impossible.", call. = FALSE)
  }
  return(prob_death / prob_symptomatic)
}
