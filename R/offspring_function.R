## need to account for people who die before they go to hospital (but who would otherwise do so)

offspring_function <- function(

  ## Characteristics and properties of the parent (who we are generating the offspring for)
  parent_infection_time = NULL,             # time of infection
  parent_hospitalised = NULL,               # whether or not the parent (infector) is hospitalised or not
  parent_time_to_hospitalisation = NULL,    # if hospitalised, the time (relative to infection) when the parent is hospitalised
  parent_outcome = NULL,                    # the disease outcome for the parent (died or recovered)
  parent_time_to_outcome = NULL,            # the time (relative to infection) when the parent dies/recovers
  parent_class = NULL,                      # whether the parent is a healthcare worker (HCW) or member of the general population (genPop)

  ## Parameters of the offspring distributions for healthcare workers (hcw) and general population (genPop)
  mn_offspring_hcw = NULL,
  overdisp_offspring_hcw = NULL,
  mn_offspring_genPop = NULL,
  overdisp_offspring_genPop = NULL,

  ##
  RR_community = NULL,
  Tg_rate = NULL,
  Tg_shape = NULL,

) {

  ## Checks to make sure inputs are in the required format
  if (is.null(parent_class) || length(parent_class) != 1L || is.na(parent_class) || !parent_class %in% c("genPop", "HCW")) {
    stop("`parent_class` must be a single string equal to 'genPop' or 'hcw'.", call. = FALSE)
  }
  if (is.null(parent_hospitalised) || length(parent_hospitalised) != 1L || !is.logical(parent_hospitalised) || is.na(parent_hospitalised)) {
    stop("`parent_hospitalised` must be a single logical value: TRUE or FALSE.", call. = FALSE)
  }

  ## Step 1: Draw overall number of offspring from offspring distribution - these distributions are
  ##         different depending on whether parent is HCW or genPop
  if (parent_class == "HCW") {
    num_offspring <- rnbinom(n = 1, mu = mn_offspring_hcw, size = overdisp_offspring_hcw)
  } else if (parent_class == "genPop") {
    num_offspring <- rnbinom(n = 1, mu = mn_offspring_genPop, size = overdisp_offspring_genPop)
  }

  ## Step 2: Apportion these infections across community vs healthcare settings and allocate their timings
  if (parent_class == "HCW") {

    ## HOW SHOULD WE HANDLE THESE??
    #### STILL TO DO

  } else if (parent_class == "genPop") {
    ## If parent is hospitalised, distribute infections across the community and healthcare setting according to time and comparative riskiness
    if (parent_hospitalised == TRUE) {
      time_in_hospital <- parent_time_to_outcome - parent_time_to_hospitalisation
      prob_offspring_in_community <- (parent_time_to_hospitalisation * RR_community) / ((parent_time_to_hospitalisation * RR_community) + time_in_hospital)
      num_offspring_community <- rbinom(n = num_offspring, prob = prob_offspring_in_community)
      num_offspring_healthcare <- num_offspring - num_offspring_community

      time_infection_offspring_community <- if (num_offspring_community > 0) {
        rtrunc_via_cdf(n = num_offspring_community,
                       lower = 0, upper = t_hosp,
                       tg_shape = Tg_shape, tg_rate = Tg_rate)
      } else numeric(0)

      times_infection_offspring_hosp <- if (num_offspring_healthcare > 0) {
        rtrunc_via_cdf(n = num_offspring_healthcare,
                       lower = t_hosp, upper = Inf,
                       tg_shape = Tg_shape, tg_rate = Tg_rate)
      } else numeric(0)

    ## If parent is not hospitalised, all infections necessarily occur in the community
    } else {
      num_offspring_community <- num_offspring
      num_offspring_healthcare <- 0
      time_infection_offspring_community <- rgamma(n = num_offspring_community,
                                                   shape = Tg_shape,
                                                   rate = Tg_rate)
    }
  }

  ## Step 3: Assigning the times of these infections



  ## Step 4: Assigning classes to these new infections (i.e. are they genPop or HCWs).
  ##         This depends on both the identity of the parent (genPop or HCW) AND
  ##         where the infection takes place (community vs healthcare setting)









  return(offspring_df)
}
