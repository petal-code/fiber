#' Generate community offspring infections
#'
#' @description
#' Creates new infections arising through *community transmission*.
#' Intended as a simple, homogeneous mixing process suitable for integration into a branching process.
#'
#' @param index_id             Numeric ID of the infector.
#' @param mn_offspring_comm    Mean of the community offspring distribution (λ for Poisson).
#' @param generation_time_fun  Function that samples generation times (passed from main model).
#' @param number_susceptible   Number of susceptible individuals remaining in the population.
#' @param population_size      Total population size (for susceptible depletion scaling).
#' @param current_max_id       Current maximum infection ID (to increment new infections sequentially).
#'
#' @returns
#' A list with:
#' - `n_offspring` = number of community infections generated.
#' - `offspring_df` = data.frame of infection details (IDs, ancestry, timing, route).
#'
#' @export
offspring_community <- function(index_id,
                                mn_offspring_comm,
                                generation_time_fun,
                                number_susceptible,
                                population_size,
                                current_max_id) {

  ## --------------------------------------------------------------------
  ## 1. Adjust mean for susceptible depletion
  ## --------------------------------------------------------------------
  # Transmission rate scales with the proportion of susceptibles left.
  adjusted_mean <- mn_offspring_comm * (number_susceptible / population_size)

  ## --------------------------------------------------------------------
  ## 2. Draw number of community offspring (Poisson)
  ## --------------------------------------------------------------------
  n_offspring <- rpois(1, lambda = adjusted_mean)

  ## --------------------------------------------------------------------
  ## 3. If no offspring, return empty structure
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
  ## 4. If offspring exist, sample their infection times
  ## --------------------------------------------------------------------
  # Generation times are drawn from the supplied generation_time_fun.
  child_gen_times <- generation_time_fun(n_offspring)

  ## --------------------------------------------------------------------
  ## 5. Assign new IDs and build offspring data frame
  ## --------------------------------------------------------------------
  new_ids <- seq(from = current_max_id + 1, length.out = n_offspring)

  offspring_df <- data.frame(
    id = new_ids,
    ancestor = rep(index_id, n_offspring),
    generation_time = child_gen_times,
    transmission_route = rep("community", n_offspring),
    stringsAsFactors = FALSE
  )

  ## --------------------------------------------------------------------
  ## 6. Return as list (same interface as healthcare version)
  ## --------------------------------------------------------------------
  return(list(
    n_offspring = n_offspring,
    offspring_df = offspring_df
  ))
}
