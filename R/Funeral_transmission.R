#' Generate offspring infections from funeral transmission
#'
#' @description
#' Models secondary infections arising from funerals following the death
#' of an infectious individual. Funerals are treated as discrete, high-risk
#' exposure events with overdispersed transmission potential.
#'
#' @param index_id               Numeric ID of the deceased infector.
#' @param deceased               Logical indicating if the infector has died.
#' @param mn_offspring_funeral   Mean number of funeral infections.
#' @param disp_offspring_funeral Dispersion parameter for Negative Binomial.
#' @param generation_time_fun    Function sampling infection generation times.
#' @param number_susceptible     Number of susceptible individuals remaining.
#' @param population_size        Total population size.
#' @param current_max_id         Current maximum infection ID.
#' @param funeral_contact_matrix Optional contact matrix representing
#'                               differential exposure across groups
#'                               (e.g., age, occupation). Default = NULL.
#'
#' @returns
#' A list containing:
#' * `n_offspring` – number of new infections generated
#' * `offspring_df` – data frame of new infections
#'
#' @export
offspring_funeral <- function(index_id,
                              deceased,
                              mn_offspring_funeral,
                              disp_offspring_funeral,
                              generation_time_fun,
                              number_susceptible,
                              population_size,
                              current_max_id,
                              funeral_contact_matrix = NULL) {

  ## --------------------------------------------------------------------
  ## 1. Check if funeral transmission should occur
  ## --------------------------------------------------------------------
  if (!isTRUE(deceased)) {
    return(list(
      n_offspring = 0,
      offspring_df = data.frame(
        id = integer(0),
        ancestor = integer(0),
        generation_time = numeric(0),
        transmission_route = character(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  ## --------------------------------------------------------------------
  ## 2. Adjust mean for susceptible depletion
  ## --------------------------------------------------------------------
  adjusted_mean <- mn_offspring_funeral * (number_susceptible / population_size)

  ## --------------------------------------------------------------------
  ## 3. Draw number of offspring using Negative Binomial
  ## --------------------------------------------------------------------
  n_offspring <- rnbinom(1, mu = adjusted_mean, size = disp_offspring_funeral)

  ## --------------------------------------------------------------------
  ## 4. If no secondary infections, return empty structure
  ## --------------------------------------------------------------------
  if (n_offspring == 0) {
    return(list(
      n_offspring = 0,
      offspring_df = data.frame(
        id = integer(0),
        ancestor = integer(0),
        generation_time = numeric(0),
        transmission_route = character(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  ## --------------------------------------------------------------------
  ## 5. Sample generation times (typically short for funeral transmission)
  ## --------------------------------------------------------------------
  child_gen_times <- generation_time_fun(n_offspring)

  ## --------------------------------------------------------------------
  ## 6. Placeholder for applying a contact matrix
  ## --------------------------------------------------------------------
  # If a contact matrix is supplied, weight infection probabilities
  # by demographic group (e.g., attendees aged 18+ more likely to be infected)
  # This is currently a placeholder; integration with demographic attributes
  # will depend on the model structure used in FIBER.
  if (!is.null(funeral_contact_matrix)) {
    # Example future use:
    # contact_weights <- funeral_contact_matrix[index_group, ]
    # probabilities <- contact_weights / sum(contact_weights)
    # infected_groups <- sample(seq_len(ncol(funeral_contact_matrix)),
    #                           size = n_offspring,
    #                           prob = probabilities,
    #                           replace = TRUE)
    # For now, this placeholder does not alter infection assignment.
    infected_groups <- rep(NA, n_offspring)
  } else {
    infected_groups <- rep(NA, n_offspring)
  }

  ## --------------------------------------------------------------------
  ## 7. Assign new IDs and compile offspring data frame
  ## --------------------------------------------------------------------
  new_ids <- seq(from = current_max_id + 1, length.out = n_offspring)

  offspring_df <- data.frame(
    id = new_ids,
    ancestor = rep(index_id, n_offspring),
    generation_time = child_gen_times,
    transmission_route = rep("funeral", n_offspring),
    contact_group = infected_groups,
    stringsAsFactors = FALSE
  )

  ## --------------------------------------------------------------------
  ## 8. Return list structure consistent with other offspring functions
  ## --------------------------------------------------------------------
  return(list(
    n_offspring = n_offspring,
    offspring_df = offspring_df
  ))
}
