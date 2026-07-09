## estimate_leaky_onward(): post-hoc "leaky OBV" sensitivity estimator.
##
## Motivation. branching_process_main() models obeldesivir (OBV) PEP as a
## *sterilizing* gate: an engaged (adherent, efficacy-positive) exposure is
## removed from the tree entirely -- never infected, never a case, never
## transmitting. A reviewer's fair concern is that this is optimistic if the drug
## in reality only *reduces* infectiousness (and/or severity) without fully
## preventing infection. This function estimates what the averted infections
## would have contributed under a **leaky** reinterpretation of the SAME efficacy,
## without rebuilding the simulation framework: it runs entirely post-hoc on a
## finished run's `prevented_completed` output.
##
## Reinterpretation (see the design notes in the manuscript response):
##   * The existing efficacy(dpc) is read as P(the drug ENGAGES), i.e. who lands
##     in `prevented_completed` / is engaged at any generation below.
##   * `residual_transmissibility` (r in [0, 1]) is the fraction of normal
##     transmission an ENGAGED individual still produces. r = 0 recovers the
##     current transmission behaviour (engaged individuals transmit nothing) but,
##     unlike the sterilizing model, they still exist as cases and can die; r = 1
##     means OBV has no transmission effect at all.
##
## Validity. This post-hoc forward simulation is unbiased *contingent on the
## current no-depletion regime*: fiber does not feed susceptible or HCW depletion
## back into transmission (susceptible_deplete is an unimplemented placeholder and
## hcw_available is tracking-only), so the subtree rooted at each averted infection
## is an independent branching process and can be simulated separately and summed.
## IF depletion is ever wired into transmission, this decoupling breaks and the
## estimator must be revisited.

## Scale a mean-offspring parameter (scalar or function(t)) by a factor. Used to
## thin an engaged (leaky) parent's transmission to `residual_transmissibility`.
## A factor of 0 is handled by the caller (engaged parents are skipped) so the
## scaled value is never passed to the offspring functions' positivity checks.
.scale_mn <- function(mn, factor) {
  if (is.function(mn)) {
    force(mn); force(factor)
    function(t) mn(t) * factor
  } else {
    mn * factor
  }
}

## Pull the 3-column offspring core (infection_location, time_infection_relative,
## class) plus an `engaged` flag out of an offspring-function result. The realised
## rows (survived thinning AND not OBV-prevented) are engaged = FALSE; the rows the
## internal sterilizing gate *prevented* -- reconstructed from the
## `obv_pep_prevented_info` attribute -- are the engaged = TRUE leaky cases. Their
## union is exactly the "survived PPE/quarantine thinning" set, split by whether
## OBV engaged. Reuses apply_obv_pep_gate() as-is (no leaky variant needed).
.leaky_offspring_core <- function(off_df, parent_time_infection_absolute) {
  realised <- data.frame(
    infection_location      = off_df$infection_location,
    time_infection_relative = off_df$time_infection_relative,
    class                   = off_df$class,
    engaged                 = rep(FALSE, nrow(off_df)),
    stringsAsFactors        = FALSE
  )
  pinfo <- attr(off_df, "obv_pep_prevented_info", exact = TRUE)
  if (!is.null(pinfo) && nrow(pinfo) > 0) {
    prevented <- data.frame(
      infection_location      = pinfo$infection_location,
      ## complete_offspring_info() reads relative times; prevented_info carries
      ## absolute infection times, so convert back against this parent's clock.
      time_infection_relative = pinfo$time_infection_absolute - parent_time_infection_absolute,
      class                   = pinfo$class,
      engaged                 = rep(TRUE, nrow(pinfo)),
      stringsAsFactors        = FALSE
    )
    return(rbind(realised, prevented))
  }
  realised
}

#' Estimate the leaky onward burden of OBV-averted infections (post-hoc)
#'
#' A post-hoc sensitivity estimator for the transmission impact of obeldesivir
#' (OBV) PEP under a \emph{leaky} (reduced-infectiousness) reinterpretation of its
#' efficacy, addressing the concern that treating OBV as fully sterilizing is
#' optimistic. It runs \strong{outside} [branching_process_main()], consuming that
#' run's `prevented_completed` frame and parameter bundle, and never perturbs the
#' original trajectory or its RNG stream.
#'
#' \strong{What it does.} The infections OBV removed in the main run
#' (`prevented_completed`) are treated as the first cohort of \emph{engaged}
#' cases: they \emph{did} become infected, but transmit at only
#' `residual_transmissibility` (r) of the normal rate. Each such case is
#' forward-simulated through the same offspring (community / hospital / funeral)
#' and [complete_offspring_info()] machinery as the main model. Descendants are
#' ordinary infections and transmit fully -- but any descendant on whom the OBV
#' gate \emph{engages} (the same coverage / adherence / efficacy draws) becomes
#' leaky in turn, transmitting at r, recursively at every generation. Because
#' fiber currently feeds no susceptible / HCW depletion back into transmission,
#' each averted infection's subtree is an independent branching process, so this
#' post-hoc replay is an unbiased estimate of the onward burden.
#'
#' \strong{Interpreting r.} r = 0 reproduces the main model's transmission
#' behaviour (engaged individuals transmit nothing) yet, unlike the sterilizing
#' model, still counts them as cases that can die -- so at r = 0 the added
#' infections equal the number of averted index infections and the added deaths
#' equal `prevented_deaths`, with no onward chains. r = 1 means OBV has no
#' transmission effect, approaching the no-OBV world. Sweeping r brackets the
#' reviewer's "partial protection / reduced infectiousness" concern.
#'
#' \strong{Accounting (decision: index cases are added back).} The averted index
#' infections are counted as realised leaky cases and their (already-resolved,
#' fixed) deaths are added back, alongside their stochastic onward chains. r
#' scales an engaged individual's transmission on \emph{all} routes -- community,
#' hospital and funeral -- since it represents that individual's overall residual
#' infectiousness.
#'
#' @param prevented_completed The `prevented_completed` data frame returned by
#'   [branching_process_main()] (the averted index infections, each already
#'   carrying its counterfactual natural history). `NULL` or a 0-row frame is
#'   allowed and yields an all-zero result.
#' @param params Named list of the transmission / natural-history / OBV
#'   parameters the run used. The convenient source is `out$sim_info$params` from
#'   the same [branching_process_main()] call; it can also be hand-built with the
#'   same names as that function's arguments.
#' @param residual_transmissibility Single numeric in \code{[0, 1]}. Fraction of
#'   normal transmission an engaged (leaky) individual still produces.
#' @param n_replicates Positive integer. Number of independent forward-simulation
#'   replicates; the onward burden is stochastic, so the result reports its
#'   distribution rather than a single value. Defaults to 1.
#' @param max_tree_size Positive integer. Per-replicate cap on the total number of
#'   simulated infections (seeds + descendants); expansion stops if it is reached
#'   and a warning is issued (so a supercritical r near 1 cannot run away
#'   silently). Defaults to 1e5.
#' @param seed Optional integer. RNG seed for reproducibility of the replicates.
#' @param return_trees Logical. If TRUE, the full simulated infection frame for
#'   each replicate is returned in `$trees` (seeds first, `engaged` column
#'   flagging leaky cases). Defaults to FALSE.
#'
#' @return A named list:
#'   \describe{
#'     \item{`residual_transmissibility`, `n_replicates`}{Echoed inputs.}
#'     \item{`n_seeds`}{Number of averted index infections replayed (rows of
#'       `prevented_completed`).}
#'     \item{`seed_deaths`}{Fixed number of would-be deaths among those index
#'       infections (equals `prevented_deaths`); added back on every replicate.}
#'     \item{`per_replicate`}{Data frame, one row per replicate, with
#'       `onward_infections` / `onward_deaths` (descendants only),
#'       `added_infections` / `added_deaths` (index cases + descendants), and
#'       `capped` (did the replicate hit `max_tree_size`).}
#'     \item{`summary`}{Means and 2.5/50/97.5\% quantiles of the four count
#'       columns across replicates.}
#'     \item{`trees`}{Present only when `return_trees = TRUE`: a list of the
#'       per-replicate infection frames.}
#'   }
#' @export
estimate_leaky_onward <- function(prevented_completed,
                                  params,
                                  residual_transmissibility,
                                  n_replicates  = 1L,
                                  max_tree_size = 1e5,
                                  seed          = NULL,
                                  return_trees  = FALSE) {

  ## ---- Validate inputs ----------------------------------------------------
  if (!is.list(params) || is.data.frame(params)) {
    stop("`params` must be a named list (e.g. out$sim_info$params).", call. = FALSE)
  }
  if (!is.numeric(residual_transmissibility) || length(residual_transmissibility) != 1L ||
      is.na(residual_transmissibility) || residual_transmissibility < 0 ||
      residual_transmissibility > 1) {
    stop("`residual_transmissibility` must be a single numeric in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(n_replicates) || length(n_replicates) != 1L || is.na(n_replicates) ||
      n_replicates < 1 || n_replicates != round(n_replicates)) {
    stop("`n_replicates` must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(max_tree_size) || length(max_tree_size) != 1L || is.na(max_tree_size) ||
      max_tree_size < 1) {
    stop("`max_tree_size` must be a single positive integer.", call. = FALSE)
  }
  n_replicates  <- as.integer(n_replicates)
  max_tree_size <- as.integer(max_tree_size)
  r <- residual_transmissibility

  if (!is.null(seed)) {
    set.seed(seed)
  }

  ## Empty-input fast path: nothing was averted, so there is nothing to replay.
  n_seeds <- if (is.null(prevented_completed)) 0L else nrow(prevented_completed)
  if (n_seeds == 0L) {
    per_replicate <- data.frame(
      replicate         = seq_len(n_replicates),
      onward_infections = 0L,
      onward_deaths     = 0L,
      added_infections  = 0L,
      added_deaths      = 0L,
      capped            = FALSE,
      stringsAsFactors  = FALSE
    )
    zero_summ <- function() list(mean = 0, median = 0, q2.5 = 0, q97.5 = 0)
    out <- list(
      residual_transmissibility = r,
      n_replicates              = n_replicates,
      n_seeds                   = 0L,
      seed_deaths               = 0L,
      per_replicate             = per_replicate,
      summary = list(onward_infections = zero_summ(), onward_deaths = zero_summ(),
                     added_infections  = zero_summ(), added_deaths  = zero_summ())
    )
    if (isTRUE(return_trees)) out$trees <- rep(list(prevented_completed), n_replicates)
    return(out)
  }

  if (!is.data.frame(prevented_completed) ||
      !all(c("class", "infection_location", "time_infection_absolute",
             "time_infection_relative", "hospitalisation",
             "time_hospitalisation_relative", "time_outcome_relative",
             "outcome", "funeral_safety") %in% names(prevented_completed))) {
    stop("`prevented_completed` must be the data frame returned by branching_process_main() (its counterfactual natural-history columns are required).",
         call. = FALSE)
  }
  if (!all(prevented_completed$class %in% c("genPop", "HCW"))) {
    stop("`prevented_completed$class` must contain only 'genPop' or 'HCW'.", call. = FALSE)
  }

  ## Seed deaths are fixed (the index infections' outcomes were already resolved
  ## in the main run's replay); only their onward chains are redrawn per replicate.
  seed_deaths <- sum(prevented_completed$outcome, na.rm = TRUE)

  ## Build the offspring-function argument list common to every parent, differing
  ## only in the (r-scaled) mean-offspring parameters, which are set per-parent.
  hosp_delay_factor <- if (is.null(params$hospitalisation_delay_factor)) 1 else params$hospitalisation_delay_factor

  complete_one_brood <- function(parent_row, offspring_core) {
    complete_offspring_info(
      parent_info                  = parent_row,
      offspring_dataframe          = offspring_core,
      prob_symptomatic             = params$prob_symptomatic,
      prob_hospitalised_hcw        = params$prob_hospitalised_hcw,
      prob_hospitalised_genPop     = params$prob_hospitalised_genPop,
      prob_death_comm              = params$prob_death_comm,
      prob_death_hosp              = params$prob_death_hosp,
      p_unsafe_funeral_comm_hcw    = params$p_unsafe_funeral_comm_hcw,
      p_unsafe_funeral_hosp_hcw    = params$p_unsafe_funeral_hosp_hcw,
      p_unsafe_funeral_comm_genPop = params$p_unsafe_funeral_comm_genPop,
      p_unsafe_funeral_hosp_genPop = params$p_unsafe_funeral_hosp_genPop,
      incubation_period            = params$incubation_period,
      onset_to_hospitalisation     = params$onset_to_hospitalisation,
      hospitalisation_delay_factor = hosp_delay_factor,
      hospitalisation_to_death     = params$hospitalisation_to_death,
      hospitalisation_to_recovery  = params$hospitalisation_to_recovery,
      onset_to_death               = params$onset_to_death,
      onset_to_recovery            = params$onset_to_recovery
    )
  }

  ## Generate + complete the offspring of a single parent under the leaky rule.
  ## Returns a completed offspring frame (with an `engaged` column) or NULL when
  ## the parent produces nothing (including an engaged parent at r = 0, who is
  ## infectious at rate 0 and so generates no offspring on any route).
  expand_parent <- function(parent_row) {
    scale <- if (isTRUE(parent_row$engaged)) r else 1
    if (scale == 0) {
      return(NULL)
    }
    parent_abs <- parent_row$time_infection_absolute

    ## Community / hospital transmission (dispatch on parent class).
    if (parent_row$class == "genPop") {
      off_ch <- offspring_function_genPop(
        parent_info                          = parent_row,
        mn_offspring_genPop                  = .scale_mn(params$mn_offspring_genPop, scale),
        overdisp_offspring_genPop            = params$overdisp_offspring_genPop,
        Tg_shape_genPop                      = params$Tg_shape_genPop,
        Tg_rate_genPop                       = params$Tg_rate_genPop,
        prop_etu                             = params$prop_etu,
        etu_efficacy                         = params$etu_efficacy,
        general_hospital_quarantine_efficacy = params$general_hospital_quarantine_efficacy,
        obv_pep_enabled                      = params$obv_pep_enabled,
        obv_pep_coverage                     = params$obv_pep_coverage,
        obv_pep_adherence                    = params$obv_pep_adherence,
        obv_pep_dpc                          = params$obv_pep_dpc,
        obv_pep_dpc_sd                       = params$obv_pep_dpc_sd,
        obv_pep_efficacy                     = params$obv_pep_efficacy,
        obv_pep_efficacy_args                = params$obv_pep_efficacy_args,
        obv_pep_target_class                 = params$obv_pep_target_class,
        obv_pep_target_locations             = params$obv_pep_target_locations,
        ppe_coverage_hcw                     = params$ppe_coverage_hcw,
        ppe_efficacy                         = params$ppe_efficacy,
        prob_hcw_cond_genPop_comm            = params$prob_hcw_cond_genPop_comm,
        prob_hcw_cond_genPop_hospital        = params$prob_hcw_cond_genPop_hospital
      )
    } else {
      off_ch <- offspring_function_hcw(
        parent_info                          = parent_row,
        mn_offspring_hcw                     = .scale_mn(params$mn_offspring_hcw, scale),
        overdisp_offspring_hcw               = params$overdisp_offspring_hcw,
        Tg_shape_hcw                         = params$Tg_shape_hcw,
        Tg_rate_hcw                          = params$Tg_rate_hcw,
        prob_hospital_cond_hcw_preAdm        = params$prob_hospital_cond_hcw_preAdm,
        ppe_coverage_hcw                     = params$ppe_coverage_hcw,
        ppe_efficacy                         = params$ppe_efficacy,
        prop_etu                             = params$prop_etu,
        etu_efficacy                         = params$etu_efficacy,
        general_hospital_quarantine_efficacy = params$general_hospital_quarantine_efficacy,
        obv_pep_enabled                      = params$obv_pep_enabled,
        obv_pep_coverage                     = params$obv_pep_coverage,
        obv_pep_adherence                    = params$obv_pep_adherence,
        obv_pep_dpc                          = params$obv_pep_dpc,
        obv_pep_dpc_sd                       = params$obv_pep_dpc_sd,
        obv_pep_efficacy                     = params$obv_pep_efficacy,
        obv_pep_efficacy_args                = params$obv_pep_efficacy_args,
        obv_pep_target_class                 = params$obv_pep_target_class,
        obv_pep_target_locations             = params$obv_pep_target_locations,
        prob_hcw_cond_hcw_comm               = params$prob_hcw_cond_hcw_comm,
        prob_hcw_cond_hcw_hospital           = params$prob_hcw_cond_hcw_hospital
      )
    }

    ## Funeral transmission (a no-op unless the parent died and had an unsafe
    ## funeral). r scales this route too, since it is the same engaged individual.
    off_fun <- offspring_function_funeral(
      parent_info                  = parent_row,
      mn_offspring_funeral         = .scale_mn(params$mn_offspring_funeral, scale),
      overdisp_offspring_funeral   = params$overdisp_offspring_funeral,
      Tg_shape_funeral             = params$Tg_shape_funeral,
      Tg_rate_funeral              = params$Tg_rate_funeral,
      safe_funeral_efficacy        = params$safe_funeral_efficacy,
      obv_pep_enabled              = params$obv_pep_enabled,
      obv_pep_coverage             = params$obv_pep_coverage,
      obv_pep_adherence            = params$obv_pep_adherence,
      obv_pep_dpc                  = params$obv_pep_dpc,
      obv_pep_dpc_sd               = params$obv_pep_dpc_sd,
      obv_pep_efficacy             = params$obv_pep_efficacy,
      obv_pep_efficacy_args        = params$obv_pep_efficacy_args,
      obv_pep_target_class         = params$obv_pep_target_class,
      obv_pep_target_locations     = params$obv_pep_target_locations,
      prob_hcw_cond_funeral_hcw    = params$prob_hcw_cond_funeral_hcw,
      prob_hcw_cond_funeral_genPop = params$prob_hcw_cond_funeral_genPop
    )

    core <- rbind(
      .leaky_offspring_core(off_ch,  parent_abs),
      .leaky_offspring_core(off_fun, parent_abs)
    )
    if (nrow(core) == 0L) {
      return(NULL)
    }
    complete_one_brood(parent_row, core)
  }

  ## ---- Run replicates -----------------------------------------------------
  onward_infections <- integer(n_replicates)
  onward_deaths     <- integer(n_replicates)
  capped            <- logical(n_replicates)
  trees             <- if (isTRUE(return_trees)) vector("list", n_replicates) else NULL

  ## Column template the seeds and every completed brood share, so the tree can be
  ## assembled by rbind and split into parent rows. Seeds get engaged = TRUE.
  seed_frame <- prevented_completed
  seed_frame$engaged    <- TRUE
  seed_frame$id         <- seq_len(n_seeds)
  seed_frame$generation <- 1L

  for (rep_i in seq_len(n_replicates)) {
    ## Queue of parent rows to expand (FIFO); pure branching, so order is
    ## irrelevant to the counts. Seeds are the initial parents.
    queue        <- split(seed_frame, seq_len(nrow(seed_frame)))
    brood_frames <- list()
    total_size   <- n_seeds
    next_id      <- n_seeds + 1L
    was_capped   <- FALSE

    qi <- 1L
    while (qi <= length(queue)) {
      if (total_size >= max_tree_size) {
        was_capped <- TRUE
        break
      }
      parent_row <- queue[[qi]]
      qi <- qi + 1L

      completed <- expand_parent(parent_row)
      if (is.null(completed) || nrow(completed) == 0L) {
        next
      }
      n_new <- nrow(completed)
      completed$id <- next_id:(next_id + n_new - 1L)
      next_id <- next_id + n_new
      total_size <- total_size + n_new

      brood_frames[[length(brood_frames) + 1L]] <- completed
      queue <- c(queue, split(completed, seq_len(n_new)))
    }

    if (length(brood_frames) > 0L) {
      descendants <- do.call(rbind, brood_frames)
      onward_infections[rep_i] <- nrow(descendants)
      onward_deaths[rep_i]     <- sum(descendants$outcome, na.rm = TRUE)
    } else {
      descendants <- NULL
      onward_infections[rep_i] <- 0L
      onward_deaths[rep_i]     <- 0L
    }
    capped[rep_i] <- was_capped

    if (isTRUE(return_trees)) {
      ## Seeds first, then descendants; align columns by intersecting names so the
      ## seed frame's extra natural-history columns and the brood columns match.
      if (is.null(descendants)) {
        trees[[rep_i]] <- seed_frame
      } else {
        shared <- intersect(names(seed_frame), names(descendants))
        trees[[rep_i]] <- rbind(seed_frame[, shared, drop = FALSE],
                                descendants[, shared, drop = FALSE])
      }
    }
  }

  if (any(capped)) {
    warning(sprintf(
      "estimate_leaky_onward(): %d of %d replicate(s) hit max_tree_size = %d; onward burden is undercounted for those. Increase max_tree_size or lower residual_transmissibility.",
      sum(capped), n_replicates, max_tree_size), call. = FALSE)
  }

  added_infections <- n_seeds     + onward_infections
  added_deaths     <- seed_deaths + onward_deaths

  per_replicate <- data.frame(
    replicate         = seq_len(n_replicates),
    onward_infections = onward_infections,
    onward_deaths     = onward_deaths,
    added_infections  = added_infections,
    added_deaths      = added_deaths,
    capped            = capped,
    stringsAsFactors  = FALSE
  )

  summ <- function(x) {
    qs <- stats::quantile(x, probs = c(0.025, 0.5, 0.975), names = FALSE, type = 7)
    list(mean = mean(x), median = qs[2], q2.5 = qs[1], q97.5 = qs[3])
  }

  out <- list(
    residual_transmissibility = r,
    n_replicates              = n_replicates,
    n_seeds                   = n_seeds,
    seed_deaths               = seed_deaths,
    per_replicate             = per_replicate,
    summary = list(
      onward_infections = summ(onward_infections),
      onward_deaths     = summ(onward_deaths),
      added_infections  = summ(added_infections),
      added_deaths      = summ(added_deaths)
    )
  )
  if (isTRUE(return_trees)) out$trees <- trees
  out
}
