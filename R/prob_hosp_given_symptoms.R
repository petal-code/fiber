prob_hosp_given_symptoms <- function(prob_hosp, prob_symptomatic) {
  if (prob_hosp_all > prob_symptomatic) {
    stop("P(hospitalisation) cannot exceed P(symptomatic).", call. = FALSE)
  }
  if (prob_symptomatic == 0) {
    if (prob_hosp_all > 0) {
      stop("P(symptomatic) = 0 but P(hospitalisation) > 0: impossible.", call. = FALSE)
    }
    return(0)
  }
  return(prob_hosp / prob_symptomatic)
}
