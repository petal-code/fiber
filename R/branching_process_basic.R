

branching_process_basic <- function(

  ## Transmission
  offspring = c("pois", "nbinom"),
  mn_offspring,
  disp_offspring,            # only used if "nbinom"

  ## Natural history
  generation_time,
  infection_to_onset,

  ## quarantine ##

  quarantine_time = 0,
  quarantine_efficacy = 0,


  ## Misc
  t0 = 0,
  tf = Inf,
  population,
  check_final_size,
  initial_immune = 0,
  seeding_cases,
  seed
) {

  ## Set seed
  set.seed(seed)

  ## Susceptibles
  susc <- population - initial_immune

  ## Offspring distribution
  dist <- match.arg(offspring)
  if (dist == "pois") {
    offspring_fun <- function(n, susc) {
      rpois(n, lambda = mn_offspring * susc / population)
    }
  } else if (dist == "nbinom") {
    if (disp_offspring <= 1) {
      stop("For 'nbinom', require disp_offspring > 1 (use 'pois' otherwise).")
    }
    offspring_fun <- function(n, susc) {
      new_mn <- mn_offspring * susc / population
      size <- new_mn / (disp_offspring - 1)
      truncdist::rtrunc(n, spec = "nbinom", b = susc, mu = new_mn, size = size)
    }
  } else stop("offspring specification is wrong")

  ## Preallocate (minimal columns only)
  max_cases <- check_final_size
  tdf <- data.frame(
    id              = integer(max_cases),
    ancestor        = integer(max_cases),
    generation      = integer(max_cases),
    time_infection  = NA_real_,
    time_onset      = NA_real_,
    n_offspring     = integer(max_cases),
    offspring_generated = FALSE,
    stringsAsFactors = FALSE
  )

  ## Seed cases
  tdf[1:seeding_cases, ] <- data.frame(
    id              = seq_len(seeding_cases),
    ancestor        = NA_integer_,
    generation      = 1L,
    time_infection  = t0 + seq(from = 0, to = 0.01, length.out = seeding_cases),
    time_onset      = NA_real_,
    n_offspring     = NA_integer_,
    offspring_generated = FALSE,
    stringsAsFactors = FALSE
  )

  ## Main expansion loop
  while ( any(!tdf$offspring_generated & !is.na(tdf$time_infection)) ) {


    ## Get earliest infection not yet expanded
    time_infection_index <- min(tdf$time_infection[!tdf$offspring_generated & !is.na(tdf$time_infection)])

    # If no further infections to expand, stop
    if (is.infinite(time_infection_index)) break

    idx <- which(tdf$time_infection == time_infection_index & !tdf$offspring_generated)[1]

    id_index   <- tdf$id[idx]
    t_index    <- tdf$time_infection[idx]
    gen_index  <- tdf$generation[idx]
    current_max <- max(tdf$id, na.rm = TRUE)

    ## Natural history timing

    onset_rel <- infection_to_onset(1)
    tdf$time_onset[idx] <- t_index + onset_rel

    ## Draw offspring count
    n_off <- offspring_fun(1, susc)
    tdf$n_offspring[idx] <- n_off
    tdf$offspring_generated[idx] <- TRUE



    ### Quarantine comes in here ###

    if (n_off > 0 && quarantine_efficacy > 0) {

      # relative infection times for children
      child_rel_times <- generation_time(n_off)

      qres <- implement_quarantine(
        symptom_onset_time        = onset_rel,
        quarantine_time           = quarantine_time,
        n_offspring               = n_off,
        offspring_infection_times = child_rel_times,
        offspring_function_draw   = NULL,
        quarantine_efficacy       = quarantine_efficacy
      )

      # update n_off and timing after quarantine
      n_off <- qres$updated_n_offspring
      child_rel_times <- qres$updated_infection_times
    }

    # If quarantine is off, generate child infection times now
    if (n_off > 0 && (!exists("child_rel_times") || length(child_rel_times) == 0)) {
      child_rel_times <- generation_time(n_off)
    }

    ## If children exist, append them
    if (n_off > 0) {
      # stop if we’d exceed preallocated size
      if ((current_max + n_off) > max_cases)
        n_off <- max(0L, max_cases - current_max)
      child_rel_times <- child_rel_times[seq_len(n_off)]

      if (n_off > 0) {
        new_ids   <- current_max + seq_len(n_off)
        new_times <- t_index + child_rel_times


        rows <- new_ids
        tdf[rows, "id"]             <- new_ids
        tdf[rows, "ancestor"]       <- id_index
        tdf[rows, "generation"]     <- gen_index + 1L
        tdf[rows, "time_infection"] <- new_times
        tdf[rows, "time_onset"] <- t_index + infection_to_onset(n_off)
        tdf[rows, "n_offspring"]    <- NA_integer_
        tdf[rows, "offspring_generated"] <- FALSE
      }
    }

    ## Deplete susceptibles
    susc <- susc - tdf$n_offspring[idx]

    ## Optional: hard stop if we filled the table
    if (max(tdf$id, na.rm = TRUE) >= max_cases) break
  }

  ## Final tidy
  tdf <- tdf[!is.na(tdf$time_infection) & tdf$time_infection <= tf, ]
  tdf <- tdf[order(tdf$time_infection, tdf$id), ]
  tdf$offspring_generated <- NULL
  rownames(tdf) <- NULL
  tdf
}

