# 00_common_time_varying_scenario_setup.R
# -------------------------------------------------------------------------
# Common setup for initial time-varying FIBER branching-process runs.
#
# This script:
#   1. sources the updated time-varying model functions;
#   2. reads the scenario matrix from final_four_scenario_values.csv;
#   3. maps dense scenario-matrix columns onto branching_process_main() args;
#   4. supplies reasonable scalar defaults for parameters not in the matrix;
#   5. runs replicate simulations and writes summary outputs.
#
# No scenario-specific values are hard-coded into the package functions; the
# scenario differences enter only through make_time_varying() curves below.
# -------------------------------------------------------------------------

SCRIPT_DIR <- paste0(getwd(), "/R/")

source_model_functions <- function(function_dir = file.path(SCRIPT_DIR)) {
  files <- c(
    "helper_functions.R",
    "resolve_time_varying.R",
    "make_time_varying.R",
    "offspring_function_genPop.R",
    "offspring_function_hcw.R",
    "offspring_function_funeral.R",
    "complete_offspring_info.R",
    "branching_process_main.R",
    "summarise_output.R",
    "hcw_loss_function.R"
  )

  for (f in files) {
    path <- file.path(function_dir, f)
    if (!file.exists(path)) {
      stop("Could not find required model function file: ", path, call. = FALSE)
    }
    source(path)
  }

  invisible(TRUE)
}

source_model_functions()

# Check that the model functions in `model_functions/` are the updated versions
# that support deriving hospital_quarantine_efficacy(t) internally from
# prop_etu(t), ipc_helper(t), and etu_efficacy_baseline. This prevents the run
# folder silently using an older copy of branching_process_main().
required_branching_args <- c("prop_etu", "ipc_helper", "etu_efficacy_baseline")
missing_branching_args <- setdiff(required_branching_args, names(formals(branching_process_main)))
if (length(missing_branching_args) > 0L) {
  stop(
    "The model_functions folder appears to contain an older branching_process_main(). ",
    "Missing argument(s): ", paste(missing_branching_args, collapse = ", "), ". ",
    "Please update model_functions/ from the current time-varying FIBER branch.",
    call. = FALSE
  )
}

# -------------------------------------------------------------------------
# Small fallback helpers.
#
# In the package context these may already exist. They are defined here only
# if missing, so this run folder can be used as a standalone smoke-test setup.
# -------------------------------------------------------------------------

if (!exists("rtrunc_gamma", mode = "function")) {
  rtrunc_gamma <- function(n, lower, upper, Tg_shape, Tg_rate) {
    if (n <= 0L) return(numeric(0))
    if (!is.numeric(lower) || length(lower) != 1L || is.na(lower)) {
      stop("`lower` must be a single numeric value.", call. = FALSE)
    }
    if (!is.numeric(upper) || length(upper) != 1L || is.na(upper)) {
      stop("`upper` must be a single numeric value.", call. = FALSE)
    }
    if (upper <= lower) {
      return(rep(lower, n))
    }

    p_lower <- pgamma(lower, shape = Tg_shape, rate = Tg_rate)
    p_upper <- pgamma(upper, shape = Tg_shape, rate = Tg_rate)
    if (p_upper <= p_lower) {
      return(rep(lower, n))
    }

    qgamma(runif(n, min = p_lower, max = p_upper),
           shape = Tg_shape,
           rate = Tg_rate)
  }
}

if (!exists("prob_hosp_given_symptoms", mode = "function")) {
  prob_hosp_given_symptoms <- function(prob_hosp, prob_symptomatic) {
    value <- prob_hosp / prob_symptomatic
    pmin(pmax(value, 0), 1)
  }
}

if (!exists("prob_death_given_symptoms", mode = "function")) {
  prob_death_given_symptoms <- function(prob_death, prob_symptomatic) {
    value <- prob_death / prob_symptomatic
    pmin(pmax(value, 0), 1)
  }
}

# -------------------------------------------------------------------------
# Distribution helpers
# -------------------------------------------------------------------------

gamma_shape_rate_from_mean_sd <- function(mean, sd) {
  if (mean <= 0 || sd <= 0) {
    stop("Gamma mean and sd must both be positive.", call. = FALSE)
  }
  list(shape = (mean / sd)^2, rate = mean / (sd^2))
}

make_gamma_sampler <- function(mean, sd) {
  pars <- gamma_shape_rate_from_mean_sd(mean = mean, sd = sd)
  function(n) {
    if (n <= 0L) return(numeric(0))
    rgamma(n = n, shape = pars$shape, rate = pars$rate)
  }
}

# -------------------------------------------------------------------------
# Scenario matrix handling
# -------------------------------------------------------------------------

read_scenario_matrix <- function(matrix_path = "C:/Users/cwhittaker/Documents/Research Projects/fiber/inst/charlie_abc_attempt") {
  if (!file.exists(matrix_path)) {
    stop("Cannot find scenario matrix CSV: ", matrix_path, call. = FALSE)
  }

  x <- read.csv(matrix_path, stringsAsFactors = FALSE, check.names = FALSE)

  numeric_cols <- setdiff(names(x), c("scenario", "scenario_label"))
  for (nm in numeric_cols) {
    x[[nm]] <- as.numeric(x[[nm]])
  }

  required_cols <- c(
    "scenario", "scenario_label", "relative_day", "prob_hosp", "delay_hosp",
    "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp",
    "prob_unsafe_funeral_etu", "prop_etu", "ipc_helper"
  )

  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0L) {
    stop("Scenario matrix is missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  x
}

clip01 <- function(x) pmin(pmax(x, 0), 1)

make_curve <- function(times, values, method = "linear") {
  make_time_varying(times = times, values = values, method = method)
}

build_time_varying_args <- function(
    scenario_id,
    matrix = read_scenario_matrix(),
    curve_method = "linear",
    etu_efficacy_baseline = DEFAULT_SCALAR_INPUTS$etu_efficacy_baseline
) {
  scenario_matrix <- matrix[matrix$scenario == scenario_id, ]
  if (nrow(scenario_matrix) == 0L) {
    stop("No rows found for scenario_id = ", scenario_id, call. = FALSE)
  }

  scenario_matrix <- scenario_matrix[order(scenario_matrix$relative_day), ]
  times <- scenario_matrix$relative_day

  # Hospital deaths can occur in ordinary hospital care or ETUs. The matrix
  # includes an ETU-specific unsafe funeral probability, so use p_ETU as a
  # mixture weight. In the supplied matrices prob_unsafe_funeral_etu is usually
  # zero, which encodes highly controlled ETU burial pathways.
  p_unsafe_funeral_hosp_values <- clip01(
    (1 - scenario_matrix$prop_etu) * scenario_matrix$prob_unsafe_funeral_hosp +
      scenario_matrix$prop_etu * scenario_matrix$prob_unsafe_funeral_etu
  )

  # IMPORTANT: hospital_quarantine_efficacy(t) is no longer pre-computed here.
  # The core model now derives it inside the hospital-transmission offspring
  # functions as:
  #
  #   etu_efficacy(t) <- etu_efficacy_baseline +
  #     (1 - etu_efficacy_baseline) * ipc_helper(t)
  #
  #   hospital_quarantine_efficacy(t) <-
  #     prop_etu(t) * etu_efficacy(t) +
  #     (1 - prop_etu(t)) * ipc_helper(t)
  #
  # The analysis layer therefore supplies prop_etu(t), ipc_helper(t), and the
  # scalar etu_efficacy_baseline. This keeps the fitted time-varying curves and
  # literature-informed scalar assumptions traceable to the parameter table.

  list(
    scenario_label = unique(scenario_matrix$scenario_label)[1L],
    scenario_matrix = scenario_matrix,

    # Hospitalisation probability: same trajectory for GenPop and HCWs for now,
    # but passed separately to preserve package structure.
    prob_hospitalised_genPop = make_curve(times, clip01(scenario_matrix$prob_hosp), curve_method),
    prob_hospitalised_hcw    = make_curve(times, clip01(scenario_matrix$prob_hosp), curve_method),

    # delay_hosp is interpreted as an absolute mean delay in days. We therefore
    # use a unit-mean raw onset_to_hospitalisation sampler in make_base_args(),
    # and let this time-varying multiplier set the calendar-day mean.
    hospitalisation_delay_factor = make_curve(times, pmax(scenario_matrix$delay_hosp, 0.01), curve_method),

    # Unsafe funeral probabilities. GenPop and HCW are deliberately kept as
    # separate arguments even though they receive the same current trajectory.
    p_unsafe_funeral_comm_genPop = make_curve(times, clip01(scenario_matrix$prob_unsafe_funeral_comm), curve_method),
    p_unsafe_funeral_comm_hcw    = make_curve(times, clip01(scenario_matrix$prob_unsafe_funeral_comm), curve_method),
    p_unsafe_funeral_hosp_genPop = make_curve(times, p_unsafe_funeral_hosp_values, curve_method),
    p_unsafe_funeral_hosp_hcw    = make_curve(times, p_unsafe_funeral_hosp_values, curve_method),

    # IPC/PPE and ETU coverage inputs. ppe_efficacy_hcw uses the same fitted
    # IPC trajectory for pre-admission HCW exposures. prop_etu and ipc_helper
    # are passed into the core hospital-transmission logic to derive
    # hospital_quarantine_efficacy(t) at candidate hospital transmission times.
    ppe_efficacy_hcw      = make_curve(times, clip01(scenario_matrix$ipc_helper), curve_method),
    prop_etu              = make_curve(times, clip01(scenario_matrix$prop_etu), curve_method),
    ipc_helper            = make_curve(times, clip01(scenario_matrix$ipc_helper), curve_method),
    etu_efficacy_baseline = etu_efficacy_baseline
  )
}

# -------------------------------------------------------------------------
# Scalar defaults for parameters not present in the time-varying matrix.
#
# These are the central values from:
#   filovirus_model_parameter_table_4_scenarios_with_etu_baseline.xlsx
#
# The parameter table itself is the audit/reference document. The simulation
# uses the fitted dense response curves from final_four_scenario_values.csv and
# the scalar defaults below. This avoids trying to reconstruct fitted curves
# from min/max ranges or orange-point anchors in the publication table.
# -------------------------------------------------------------------------

DEFAULT_SCALAR_INPUTS <- list(
  # Hospital / ETU efficacy.
  etu_efficacy_baseline = 0.90,

  # Transmission means and dispersion.
  mn_offspring_genPop = 1.25,
  overdisp_offspring_genPop = 0.18,
  genpop_generation_mean = 15.4,
  genpop_generation_shape = 2.5,

  mn_offspring_hcw = 0.20,
  overdisp_offspring_hcw = 0.18,
  hcw_generation_mean = 15.4,
  hcw_generation_shape = 2.5,

  mn_offspring_funeral = 0.25,
  overdisp_offspring_funeral = 0.30,
  funeral_generation_shape = 20,
  funeral_generation_rate = 10,

  # Natural history.
  incubation_mean = 8.5,
  incubation_sd = 4.5,
  raw_onset_to_hospitalisation_mean = 1.0,
  raw_onset_to_hospitalisation_sd = 0.35,
  onset_to_death_mean = 9.3,
  onset_to_death_sd = 3.0,
  onset_to_recovery_mean = 13.0,
  onset_to_recovery_sd = 4.0,
  hospitalisation_to_death_mean = 4.5,
  hospitalisation_to_death_sd = 2.0,
  hospitalisation_to_recovery_mean = 8.0,
  hospitalisation_to_recovery_sd = 2.5,

  # Disease severity.
  prob_symptomatic = 1.0,
  prob_death_comm = 0.70,
  prob_death_hosp = 0.50,

  # Conditional class / setting assignment.
  prob_hcw_cond_genPop_comm = 0.005,
  prob_hcw_cond_genPop_hospital = 0.12,
  prob_hcw_cond_hcw_comm = 0.02,
  prob_hcw_cond_hcw_hospital = 0.20,
  prob_hospital_cond_hcw_preAdm = 0.50,

  # Funeral control / assignment.
  safe_funeral_efficacy = 1.0,
  prob_hcw_cond_funeral_hcw = 0.02,
  prob_hcw_cond_funeral_genPop = 0.005,

  # Population / simulation control.
  population = 1000000,
  hcw_per_capita = 0.005,
  seeding_cases = 5,
  check_final_size = 30000,
  initial_immune = 0,
  susceptible_deplete = FALSE,
  tf = Inf
)

make_base_args <- function(
    population = DEFAULT_SCALAR_INPUTS$population,
    seeding_cases = DEFAULT_SCALAR_INPUTS$seeding_cases,
    check_final_size = DEFAULT_SCALAR_INPUTS$check_final_size,
    hcw_per_capita = DEFAULT_SCALAR_INPUTS$hcw_per_capita,
    scalar_inputs = DEFAULT_SCALAR_INPUTS
) {
  list(
    # Transmission. Values are central estimates from the parameter table.
    mn_offspring_genPop = scalar_inputs$mn_offspring_genPop,
    overdisp_offspring_genPop = scalar_inputs$overdisp_offspring_genPop,
    Tg_shape_genPop = scalar_inputs$genpop_generation_shape,
    Tg_rate_genPop = scalar_inputs$genpop_generation_shape / scalar_inputs$genpop_generation_mean,

    mn_offspring_hcw = scalar_inputs$mn_offspring_hcw,
    overdisp_offspring_hcw = scalar_inputs$overdisp_offspring_hcw,
    Tg_shape_hcw = scalar_inputs$hcw_generation_shape,
    Tg_rate_hcw = scalar_inputs$hcw_generation_shape / scalar_inputs$hcw_generation_mean,

    mn_offspring_funeral = scalar_inputs$mn_offspring_funeral,
    overdisp_offspring_funeral = scalar_inputs$overdisp_offspring_funeral,
    Tg_shape_funeral = scalar_inputs$funeral_generation_shape,
    Tg_rate_funeral = scalar_inputs$funeral_generation_rate,

    # Natural history.
    incubation_period = make_gamma_sampler(
      mean = scalar_inputs$incubation_mean,
      sd = scalar_inputs$incubation_sd
    ),

    # Raw hospitalisation delay has mean 1. The matrix delay_hosp curve then
    # scales this to the desired mean delay in calendar days.
    onset_to_hospitalisation = make_gamma_sampler(
      mean = scalar_inputs$raw_onset_to_hospitalisation_mean,
      sd = scalar_inputs$raw_onset_to_hospitalisation_sd
    ),

    onset_to_death = make_gamma_sampler(
      mean = scalar_inputs$onset_to_death_mean,
      sd = scalar_inputs$onset_to_death_sd
    ),
    onset_to_recovery = make_gamma_sampler(
      mean = scalar_inputs$onset_to_recovery_mean,
      sd = scalar_inputs$onset_to_recovery_sd
    ),
    hospitalisation_to_death = make_gamma_sampler(
      mean = scalar_inputs$hospitalisation_to_death_mean,
      sd = scalar_inputs$hospitalisation_to_death_sd
    ),
    hospitalisation_to_recovery = make_gamma_sampler(
      mean = scalar_inputs$hospitalisation_to_recovery_mean,
      sd = scalar_inputs$hospitalisation_to_recovery_sd
    ),

    # Disease severity and healthcare seeking.
    prob_symptomatic = scalar_inputs$prob_symptomatic,
    prob_death_comm = scalar_inputs$prob_death_comm,
    prob_death_hosp = scalar_inputs$prob_death_hosp,

    # Conditional class probabilities.
    prob_hcw_cond_genPop_comm = scalar_inputs$prob_hcw_cond_genPop_comm,
    prob_hcw_cond_genPop_hospital = scalar_inputs$prob_hcw_cond_genPop_hospital,
    prob_hcw_cond_hcw_comm = scalar_inputs$prob_hcw_cond_hcw_comm,
    prob_hcw_cond_hcw_hospital = scalar_inputs$prob_hcw_cond_hcw_hospital,
    prob_hospital_cond_hcw_preAdm = scalar_inputs$prob_hospital_cond_hcw_preAdm,

    # Safe funerals are assumed to block funeral transmission fully for now.
    safe_funeral_efficacy = scalar_inputs$safe_funeral_efficacy,
    prob_hcw_cond_funeral_hcw = scalar_inputs$prob_hcw_cond_funeral_hcw,
    prob_hcw_cond_funeral_genPop = scalar_inputs$prob_hcw_cond_funeral_genPop,

    # Simulation controls.
    tf = scalar_inputs$tf,
    population = population,
    hcw_per_capita = hcw_per_capita,
    check_final_size = check_final_size,
    initial_immune = scalar_inputs$initial_immune,
    seeding_cases = seeding_cases,
    susceptible_deplete = scalar_inputs$susceptible_deplete
  )
}

# -------------------------------------------------------------------------
# Running and summarising simulations
# -------------------------------------------------------------------------

summarise_one_outbreak <- function(out, scenario_id, scenario_label, replicate, check_final_size) {
  tdf <- out$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]

  deaths <- !is.na(tdf$outcome) & tdf$outcome
  hcw <- !is.na(tdf$class) & tdf$class == "HCW"
  unsafe_funeral <- !is.na(tdf$funeral_safety) & tdf$funeral_safety == "unsafe"

  data.frame(
    scenario = scenario_id,
    scenario_label = scenario_label,
    replicate = replicate,
    n_cases = nrow(tdf),
    n_deaths = sum(deaths),
    n_hcw_cases = sum(hcw),
    n_hcw_deaths = sum(hcw & deaths),
    outbreak_duration_days = suppressWarnings(max(c(tdf$time_infection_absolute,
                                                    tdf$time_outcome_absolute), na.rm = TRUE)),
    max_generation = suppressWarnings(max(tdf$generation, na.rm = TRUE)),
    prop_hospitalised = mean(tdf$hospitalisation, na.rm = TRUE),
    prop_unsafe_funeral_among_deaths = ifelse(sum(deaths) > 0, sum(unsafe_funeral & deaths) / sum(deaths), NA_real_),
    hit_size_cap = nrow(tdf) >= check_final_size,
    stringsAsFactors = FALSE
  )
}

make_weekly_incidence <- function(out, scenario_id, scenario_label, replicate, bin_width = 7) {
  tdf <- out$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
  if (nrow(tdf) == 0L) {
    return(data.frame())
  }

  week <- floor(tdf$time_infection_absolute / bin_width) * bin_width
  tab <- as.data.frame(table(week), stringsAsFactors = FALSE)
  names(tab) <- c("week_start_day", "cases")
  tab$week_start_day <- as.numeric(tab$week_start_day)
  tab$cases <- as.integer(tab$cases)
  tab$scenario <- scenario_id
  tab$scenario_label <- scenario_label
  tab$replicate <- replicate
  tab[, c("scenario", "scenario_label", "replicate", "week_start_day", "cases")]
}

summarise_replicate_table <- function(summary_df) {
  q <- function(x, p) as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))

  data.frame(
    scenario = unique(summary_df$scenario)[1L],
    scenario_label = unique(summary_df$scenario_label)[1L],
    n_reps = nrow(summary_df),
    median_cases = q(summary_df$n_cases, 0.50),
    lo95_cases = q(summary_df$n_cases, 0.025),
    hi95_cases = q(summary_df$n_cases, 0.975),
    median_deaths = q(summary_df$n_deaths, 0.50),
    median_hcw_cases = q(summary_df$n_hcw_cases, 0.50),
    median_duration_days = q(summary_df$outbreak_duration_days, 0.50),
    prop_hit_size_cap = mean(summary_df$hit_size_cap),
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------------
# Weekly-incidence diagnostics
# -------------------------------------------------------------------------
# The raw weekly incidence table contains only replicate-weeks with at least one
# infection. If we summarise that table directly, extinct or zero-incidence
# replicates disappear from later weeks. That is useful for a conditional
# "among outbreaks with cases this week" diagnostic, but misleading as the main
# epidemic curve. The functions below therefore create:
#   1. a zero-filled unconditional curve, where every replicate contributes zero
#      after extinction or during weeks with no new infections;
#   2. a conditional active-only curve, which reproduces the old behaviour but is
#      explicitly labelled as conditional;
#   3. proportions with non-zero weekly incidence and not-yet-ended outbreaks.

summarise_weekly_values <- function(df, value_col = "cases") {
  if (nrow(df) == 0L) {
    return(data.frame())
  }

  by_week <- split(df[[value_col]], df$week_start_day)

  out <- do.call(
    rbind,
    lapply(names(by_week), function(w) {
      x <- by_week[[w]]
      data.frame(
        week_start_day = as.numeric(w),
        median = median(x, na.rm = TRUE),
        lo95 = unname(stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE)),
        hi95 = unname(stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE)),
        stringsAsFactors = FALSE
      )
    })
  )

  out <- out[order(out$week_start_day), , drop = FALSE]
  out$week_since_start <- out$week_start_day / 7
  out
}

make_weekly_diagnostics <- function(summary_df, weekly_df, bin_width = 7) {
  if (nrow(summary_df) == 0L) {
    return(list(
      zero_filled_summary = data.frame(),
      active_only_summary = data.frame(),
      activity_summary = data.frame(),
      zero_filled_full = data.frame()
    ))
  }

  scenario_id <- unique(summary_df$scenario)[1L]
  scenario_label <- unique(summary_df$scenario_label)[1L]

  duration_col <- if ("outbreak_duration_days" %in% names(summary_df)) {
    "outbreak_duration_days"
  } else if ("duration_days" %in% names(summary_df)) {
    "duration_days"
  } else {
    stop("Could not find outbreak-duration column in summary_df.", call. = FALSE)
  }

  finite_duration <- summary_df[[duration_col]][is.finite(summary_df[[duration_col]])]
  finite_week <- if (nrow(weekly_df) > 0L) weekly_df$week_start_day[is.finite(weekly_df$week_start_day)] else numeric(0)

  max_day <- max(c(finite_duration, finite_week, 0), na.rm = TRUE)
  max_week <- ceiling(max_day / bin_width) * bin_width

  full_grid <- expand.grid(
    replicate = unique(summary_df$replicate),
    week_start_day = seq(0, max_week, by = bin_width)
  )

  if (nrow(weekly_df) > 0L) {
    merge_cols <- c("replicate", "week_start_day", "cases")
    weekly_full <- merge(
      full_grid,
      weekly_df[, merge_cols, drop = FALSE],
      by = c("replicate", "week_start_day"),
      all.x = TRUE
    )
  } else {
    weekly_full <- full_grid
    weekly_full$cases <- NA_integer_
  }

  weekly_full$scenario <- scenario_id
  weekly_full$scenario_label <- scenario_label
  weekly_full$cases_zero_filled <- weekly_full$cases
  weekly_full$cases_zero_filled[is.na(weekly_full$cases_zero_filled)] <- 0L
  weekly_full$has_cases_this_week <- !is.na(weekly_full$cases) & weekly_full$cases > 0

  duration_lookup <- summary_df[, c("replicate", duration_col), drop = FALSE]
  names(duration_lookup)[2L] <- "outbreak_duration_days"
  weekly_full <- merge(weekly_full, duration_lookup, by = "replicate", all.x = TRUE)
  weekly_full$not_yet_ended <- weekly_full$week_start_day <= weekly_full$outbreak_duration_days

  zero_filled_summary <- summarise_weekly_values(
    weekly_full,
    value_col = "cases_zero_filled"
  )
  zero_filled_summary$scenario <- scenario_id
  zero_filled_summary$scenario_label <- scenario_label
  zero_filled_summary <- zero_filled_summary[, c(
    "scenario", "scenario_label", "week_start_day", "week_since_start", "median", "lo95", "hi95"
  )]

  active_only <- weekly_full[weekly_full$has_cases_this_week, , drop = FALSE]
  active_only_summary <- summarise_weekly_values(active_only, value_col = "cases")
  if (nrow(active_only_summary) > 0L) {
    active_only_summary$scenario <- scenario_id
    active_only_summary$scenario_label <- scenario_label
    active_only_summary <- active_only_summary[, c(
      "scenario", "scenario_label", "week_start_day", "week_since_start", "median", "lo95", "hi95"
    )]
  }

  activity_summary <- aggregate(
    cbind(has_cases_this_week, not_yet_ended) ~ week_start_day,
    data = weekly_full,
    FUN = mean
  )
  activity_summary$scenario <- scenario_id
  activity_summary$scenario_label <- scenario_label
  names(activity_summary) <- c(
    "week_start_day", "prop_with_cases_this_week", "prop_not_yet_ended",
    "scenario", "scenario_label"
  )
  activity_summary$week_since_start <- activity_summary$week_start_day / bin_width
  activity_summary <- activity_summary[, c(
    "scenario", "scenario_label", "week_start_day", "week_since_start",
    "prop_with_cases_this_week", "prop_not_yet_ended"
  )]

  list(
    zero_filled_summary = zero_filled_summary,
    active_only_summary = active_only_summary,
    activity_summary = activity_summary,
    zero_filled_full = weekly_full
  )
}

plot_weekly_summary <- function(summary_df, output_path, main, ylab = "Weekly infections") {
  if (nrow(summary_df) == 0L) {
    return(invisible(NULL))
  }

  x <- if ("week_since_start" %in% names(summary_df)) {
    summary_df$week_since_start
  } else {
    summary_df$week_start_day / 7
  }

  png(output_path, width = 1100, height = 700, res = 120)
  on.exit(dev.off(), add = TRUE)

  ylim_upper <- max(summary_df$hi95, summary_df$median, na.rm = TRUE)
  if (!is.finite(ylim_upper) || ylim_upper <= 0) {
    ylim_upper <- 1
  }

  plot(
    x,
    summary_df$median,
    type = "l",
    lwd = 2,
    ylim = c(0, ylim_upper),
    xlab = "Weeks since first infection",
    ylab = ylab,
    main = main
  )
  lines(x, summary_df$lo95, lty = 2)
  lines(x, summary_df$hi95, lty = 2)
  legend(
    "topright",
    legend = c("Median", "95% interval"),
    lty = c(1, 2),
    lwd = c(2, 1),
    bty = "n"
  )

  invisible(summary_df)
}

plot_activity_summary <- function(activity_df, output_path, main) {
  if (nrow(activity_df) == 0L) {
    return(invisible(NULL))
  }

  x <- if ("week_since_start" %in% names(activity_df)) {
    activity_df$week_since_start
  } else {
    activity_df$week_start_day / 7
  }

  png(output_path, width = 1100, height = 700, res = 120)
  on.exit(dev.off(), add = TRUE)

  plot(
    x,
    activity_df$prop_with_cases_this_week,
    type = "l",
    lwd = 2,
    ylim = c(0, 1),
    xlab = "Weeks since first infection",
    ylab = "Proportion of simulations",
    main = main
  )
  lines(x, activity_df$prop_not_yet_ended, lty = 2)
  legend(
    "topright",
    legend = c("With >=1 new case this week", "Outbreak not yet ended"),
    lty = c(1, 2),
    lwd = c(2, 1),
    bty = "n"
  )

  invisible(activity_df)
}

make_final_size_diagnostics <- function(summary_df, seeding_cases, check_final_size) {
  data.frame(
    scenario = unique(summary_df$scenario)[1L],
    scenario_label = unique(summary_df$scenario_label)[1L],
    n_reps = nrow(summary_df),
    prop_no_onward_transmission = mean(summary_df$n_cases <= seeding_cases, na.rm = TRUE),
    prop_10_cases_or_fewer = mean(summary_df$n_cases <= 10, na.rm = TRUE),
    prop_100_cases_or_more = mean(summary_df$n_cases >= 100, na.rm = TRUE),
    prop_1000_cases_or_more = mean(summary_df$n_cases >= 1000, na.rm = TRUE),
    prop_hit_size_cap = mean(summary_df$n_cases >= check_final_size, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

plot_final_size_histogram <- function(summary_df, output_path, main, check_final_size) {
  if (nrow(summary_df) == 0L) {
    return(invisible(NULL))
  }

  png(output_path, width = 1100, height = 700, res = 120)
  on.exit(dev.off(), add = TRUE)

  hist(
    summary_df$n_cases,
    breaks = "FD",
    xlab = "Total cases by end of simulation",
    ylab = "Number of replicate simulations",
    main = main
  )
  abline(v = check_final_size, lty = 2, lwd = 2)
  legend(
    "topright",
    legend = paste0("Simulation size cap = ", check_final_size),
    lty = 2,
    lwd = 2,
    bty = "n"
  )

  invisible(summary_df)
}

plot_final_size_histogram_log10 <- function(summary_df, output_path, main, check_final_size) {
  if (nrow(summary_df) == 0L) {
    return(invisible(NULL))
  }

  n_cases <- pmax(summary_df$n_cases, 1)

  png(output_path, width = 1100, height = 700, res = 120)
  on.exit(dev.off(), add = TRUE)

  hist(
    log10(n_cases),
    breaks = "FD",
    xlab = "Total cases by end of simulation, log10 scale",
    ylab = "Number of replicate simulations",
    main = main,
    xaxt = "n"
  )

  ticks <- c(1, 3, 10, 30, 100, 300, 1000, 3000, check_final_size)
  ticks <- sort(unique(ticks[ticks <= max(n_cases, check_final_size, na.rm = TRUE)]))
  axis(1, at = log10(ticks), labels = ticks)
  abline(v = log10(check_final_size), lty = 2, lwd = 2)
  legend(
    "topright",
    legend = paste0("Simulation size cap = ", check_final_size),
    lty = 2,
    lwd = 2,
    bty = "n"
  )

  invisible(summary_df)
}

make_cumulative_case_summary <- function(zero_filled_full) {
  if (nrow(zero_filled_full) == 0L) {
    return(data.frame())
  }

  x <- zero_filled_full[order(zero_filled_full$replicate, zero_filled_full$week_start_day), ]
  x$cumulative_cases <- ave(x$cases_zero_filled, x$replicate, FUN = cumsum)

  out <- summarise_weekly_values(x, value_col = "cumulative_cases")
  out$scenario <- unique(x$scenario)[1L]
  out$scenario_label <- unique(x$scenario_label)[1L]
  out <- out[, c(
    "scenario", "scenario_label", "week_start_day", "week_since_start", "median", "lo95", "hi95"
  )]
  out
}

plot_cumulative_case_summary <- function(cumulative_df, output_path, main) {
  plot_weekly_summary(
    summary_df = cumulative_df,
    output_path = output_path,
    main = main,
    ylab = "Cumulative cases"
  )
}

plot_zero_case_proportion <- function(activity_df, output_path, main) {
  if (nrow(activity_df) == 0L) {
    return(invisible(NULL))
  }

  x <- if ("week_since_start" %in% names(activity_df)) {
    activity_df$week_since_start
  } else {
    activity_df$week_start_day / 7
  }

  prop_zero_this_week <- 1 - activity_df$prop_with_cases_this_week
  prop_ended <- 1 - activity_df$prop_not_yet_ended

  png(output_path, width = 1100, height = 700, res = 120)
  on.exit(dev.off(), add = TRUE)

  plot(
    x,
    prop_zero_this_week,
    type = "l",
    lwd = 2,
    ylim = c(0, 1),
    xlab = "Weeks since first infection",
    ylab = "Proportion of simulations",
    main = main
  )
  lines(x, prop_ended, lty = 2)
  legend(
    "bottomright",
    legend = c("Zero new cases this week", "Outbreak ended"),
    lty = c(1, 2),
    lwd = c(2, 1),
    bty = "n"
  )

  invisible(activity_df)
}


# -------------------------------------------------------------------------
# Simpler epidemic-curve diagnostics for quick visual inspection
# -------------------------------------------------------------------------
# These are intended to be easier to read than median/quantile weekly-incidence
# curves. They show the distribution of cases over time using bar plots and
# overlay the proportion of simulations that have gone extinct/ended.

make_weekly_histogram_summary <- function(zero_filled_full, activity_summary) {
  if (nrow(zero_filled_full) == 0L) {
    return(data.frame())
  }

  by_week <- split(zero_filled_full$cases_zero_filled, zero_filled_full$week_start_day)
  out <- do.call(
    rbind,
    lapply(names(by_week), function(w) {
      x <- by_week[[w]]
      data.frame(
        week_start_day = as.numeric(w),
        week_since_start = as.numeric(w) / 7,
        mean_weekly_cases = mean(x, na.rm = TRUE),
        median_weekly_cases = median(x, na.rm = TRUE),
        lo95_weekly_cases = unname(stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE)),
        hi95_weekly_cases = unname(stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE)),
        prop_zero_cases_this_week = mean(x == 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )

  out <- out[order(out$week_start_day), , drop = FALSE]

  if (!is.null(activity_summary) && nrow(activity_summary) > 0L) {
    keep <- c("week_start_day", "prop_not_yet_ended", "prop_with_cases_this_week")
    activity_keep <- activity_summary[, intersect(keep, names(activity_summary)), drop = FALSE]
    out <- merge(out, activity_keep, by = "week_start_day", all.x = TRUE)
  }

  if (!"prop_not_yet_ended" %in% names(out)) {
    out$prop_not_yet_ended <- NA_real_
  }
  if (!"prop_with_cases_this_week" %in% names(out)) {
    out$prop_with_cases_this_week <- 1 - out$prop_zero_cases_this_week
  }

  out$prop_extinct_or_ended <- 1 - out$prop_not_yet_ended
  out$scenario <- unique(zero_filled_full$scenario)[1L]
  out$scenario_label <- unique(zero_filled_full$scenario_label)[1L]

  out[, c(
    "scenario", "scenario_label", "week_start_day", "week_since_start",
    "mean_weekly_cases", "median_weekly_cases", "lo95_weekly_cases", "hi95_weekly_cases",
    "prop_zero_cases_this_week", "prop_with_cases_this_week",
    "prop_not_yet_ended", "prop_extinct_or_ended"
  )]
}

plot_weekly_case_histogram_with_extinction <- function(hist_df, output_path, main) {
  if (nrow(hist_df) == 0L) {
    return(invisible(NULL))
  }

  png(output_path, width = 1200, height = 750, res = 120)
  on.exit(dev.off(), add = TRUE)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  y_cases_max <- max(hist_df$mean_weekly_cases, na.rm = TRUE)
  if (!is.finite(y_cases_max) || y_cases_max <= 0) {
    y_cases_max <- 1
  }

  bar_mid <- barplot(
    height = hist_df$mean_weekly_cases,
    names.arg = hist_df$week_since_start,
    space = 0,
    border = NA,
    xlab = "Weeks since first infection",
    ylab = "Mean weekly cases per simulation",
    main = main,
    ylim = c(0, y_cases_max * 1.15),
    axes = TRUE
  )

  # Keep the x-axis readable by labelling roughly every 4 weeks.
  axis_ticks <- seq(1, length(bar_mid), by = max(1, floor(length(bar_mid) / 12)))
  axis(1, at = bar_mid[axis_ticks], labels = hist_df$week_since_start[axis_ticks])

  par(new = TRUE)
  plot(
    bar_mid,
    hist_df$prop_extinct_or_ended,
    type = "l",
    lwd = 2,
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = c(0, 1)
  )
  axis(4)
  mtext("Proportion extinct/ended", side = 4, line = 3)
  legend(
    "topright",
    legend = c("Mean weekly cases", "Proportion extinct/ended"),
    lty = c(NA, 1),
    lwd = c(NA, 2),
    pch = c(15, NA),
    bty = "n"
  )

  invisible(hist_df)
}

plot_smooth_epidemic_curve <- function(hist_df, output_path, main, span = 0.35) {
  if (nrow(hist_df) == 0L) {
    return(invisible(NULL))
  }

  x <- hist_df$week_since_start
  y <- hist_df$mean_weekly_cases

  png(output_path, width = 1200, height = 750, res = 120)
  on.exit(dev.off(), add = TRUE)

  ylim_upper <- max(y, na.rm = TRUE)
  if (!is.finite(ylim_upper) || ylim_upper <= 0) {
    ylim_upper <- 1
  }

  plot(
    x,
    y,
    type = "n",
    ylim = c(0, ylim_upper * 1.15),
    xlab = "Weeks since first infection",
    ylab = "Mean weekly cases per simulation",
    main = main
  )

  points(x, y, pch = 16, cex = 0.55)

  enough_points <- sum(is.finite(x) & is.finite(y)) >= 5L && length(unique(x[is.finite(x)])) >= 5L
  if (enough_points) {
    fit <- try(stats::loess(y ~ x, span = span, degree = 1), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      pred <- predict(fit, newdata = data.frame(x = x))
      pred <- pmax(pred, 0)
      lines(x, pred, lwd = 3)
    } else {
      lines(x, y, lwd = 2)
    }
  } else {
    lines(x, y, lwd = 2)
  }

  legend(
    "topright",
    legend = c("Weekly mean", "Smoothed mean"),
    pch = c(16, NA),
    lty = c(NA, 1),
    lwd = c(NA, 3),
    bty = "n"
  )

  invisible(hist_df)
}

make_escape_probability_summary <- function(summary_df, thresholds = c(10, 50, 100, 500, 1000, 5000)) {
  if (nrow(summary_df) == 0L) {
    return(data.frame())
  }

  data.frame(
    scenario = unique(summary_df$scenario)[1L],
    scenario_label = unique(summary_df$scenario_label)[1L],
    threshold_cases = thresholds,
    prop_exceeding_or_equal = sapply(thresholds, function(th) mean(summary_df$n_cases >= th, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

make_key_scenario_summary <- function(summary_df, check_final_size) {
  q <- function(x, p) as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))

  data.frame(
    scenario = unique(summary_df$scenario)[1L],
    scenario_label = unique(summary_df$scenario_label)[1L],
    n_reps = nrow(summary_df),
    median_cases = q(summary_df$n_cases, 0.50),
    mean_cases = mean(summary_df$n_cases, na.rm = TRUE),
    q75_cases = q(summary_df$n_cases, 0.75),
    q90_cases = q(summary_df$n_cases, 0.90),
    q95_cases = q(summary_df$n_cases, 0.95),
    q975_cases = q(summary_df$n_cases, 0.975),
    prop_10_cases_or_fewer = mean(summary_df$n_cases <= 10, na.rm = TRUE),
    prop_100_cases_or_more = mean(summary_df$n_cases >= 100, na.rm = TRUE),
    prop_500_cases_or_more = mean(summary_df$n_cases >= 500, na.rm = TRUE),
    prop_1000_cases_or_more = mean(summary_df$n_cases >= 1000, na.rm = TRUE),
    prop_hit_size_cap = mean(summary_df$n_cases >= check_final_size, na.rm = TRUE),
    median_deaths = q(summary_df$n_deaths, 0.50),
    median_hcw_cases = q(summary_df$n_hcw_cases, 0.50),
    q95_hcw_cases = q(summary_df$n_hcw_cases, 0.95),
    median_duration_days = q(summary_df$outbreak_duration_days, 0.50),
    q95_duration_days = q(summary_df$outbreak_duration_days, 0.95),
    stringsAsFactors = FALSE
  )
}

make_time_varying_parameter_summary <- function(
    matrix = read_scenario_matrix(),
    etu_efficacy_baseline = DEFAULT_SCALAR_INPUTS$etu_efficacy_baseline
) {
  matrix$etu_efficacy <- clip01(
    etu_efficacy_baseline +
      (1 - etu_efficacy_baseline) * matrix$ipc_helper
  )

  matrix$hospital_quarantine_efficacy <- clip01(
    matrix$prop_etu * matrix$etu_efficacy +
      (1 - matrix$prop_etu) * matrix$ipc_helper
  )

  matrix$p_unsafe_funeral_hosp_mixed <- clip01(
    (1 - matrix$prop_etu) * matrix$prob_unsafe_funeral_hosp +
      matrix$prop_etu * matrix$prob_unsafe_funeral_etu
  )

  param_cols <- c(
    "prob_hosp", "delay_hosp", "prob_unsafe_funeral_comm",
    "prob_unsafe_funeral_hosp", "prob_unsafe_funeral_etu",
    "p_unsafe_funeral_hosp_mixed", "prop_etu", "ipc_helper",
    "etu_efficacy", "hospital_quarantine_efficacy"
  )

  pieces <- list()
  k <- 1L
  for (sc in unique(matrix$scenario)) {
    sc_df <- matrix[matrix$scenario == sc, , drop = FALSE]
    sc_df <- sc_df[order(sc_df$relative_day), , drop = FALSE]
    for (param in param_cols) {
      if (!param %in% names(sc_df)) next
      x <- sc_df[[param]]
      pieces[[k]] <- data.frame(
        scenario = sc,
        scenario_label = unique(sc_df$scenario_label)[1L],
        parameter = param,
        start_value = x[1L],
        end_value = x[length(x)],
        min_value = min(x, na.rm = TRUE),
        max_value = max(x, na.rm = TRUE),
        mean_value = mean(x, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  do.call(rbind, pieces)
}

run_scenario <- function(
    scenario_id,
    n_reps = as.integer(Sys.getenv("N_REPS", "100")),
    seed_base = 1000,
    output_dir = file.path(SCRIPT_DIR, "outputs", scenario_id),
    matrix_path = file.path(SCRIPT_DIR, "final_four_scenario_values.csv"),
    base_args = make_base_args(check_final_size = as.integer(Sys.getenv("CHECK_FINAL_SIZE", "30000"))),
    save_individual_trees = FALSE
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  matrix <- read_scenario_matrix(matrix_path)
  tv <- build_time_varying_args(scenario_id = scenario_id, matrix = matrix)
  scenario_label <- tv$scenario_label

  scenario_matrix <- tv$scenario_matrix
  tv$scenario_label <- NULL
  tv$scenario_matrix <- NULL

  args_without_seed <- c(base_args, tv)

  summaries <- vector("list", n_reps)
  weekly <- vector("list", n_reps)
  trees <- if (save_individual_trees) vector("list", n_reps) else NULL

  for (i in seq_len(n_reps)) {
    sim_seed <- seed_base + i
    message("[", scenario_id, "] replicate ", i, "/", n_reps, " seed=", sim_seed)

    out <- do.call(branching_process_main, c(args_without_seed, list(seed = sim_seed)))

    summaries[[i]] <- summarise_one_outbreak(
      out = out,
      scenario_id = scenario_id,
      scenario_label = scenario_label,
      replicate = i,
      check_final_size = base_args$check_final_size
    )

    weekly[[i]] <- make_weekly_incidence(
      out = out,
      scenario_id = scenario_id,
      scenario_label = scenario_label,
      replicate = i
    )

    if (save_individual_trees) {
      trees[[i]] <- out
    }
  }

  summary_df <- do.call(rbind, summaries)
  weekly_df <- do.call(rbind, weekly)
  quick_summary <- summarise_replicate_table(summary_df)

  write.csv(scenario_matrix,
            file = file.path(output_dir, paste0(scenario_id, "_time_varying_matrix_used.csv")),
            row.names = FALSE)

  scalar_defaults_used <- data.frame(
    parameter = names(DEFAULT_SCALAR_INPUTS),
    value = vapply(DEFAULT_SCALAR_INPUTS, function(x) paste(x, collapse = "; "), character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(scalar_defaults_used,
            file = file.path(output_dir, paste0(scenario_id, "_scalar_defaults_used.csv")),
            row.names = FALSE)
  write.csv(summary_df,
            file = file.path(output_dir, paste0(scenario_id, "_replicate_summary.csv")),
            row.names = FALSE)
  write.csv(weekly_df,
            file = file.path(output_dir, paste0(scenario_id, "_weekly_incidence_active_only_raw.csv")),
            row.names = FALSE)

  weekly_diagnostics <- make_weekly_diagnostics(
    summary_df = summary_df,
    weekly_df = weekly_df
  )

  weekly_histogram_summary <- make_weekly_histogram_summary(
    zero_filled_full = weekly_diagnostics$zero_filled_full,
    activity_summary = weekly_diagnostics$activity_summary
  )

  escape_probability_summary <- make_escape_probability_summary(summary_df)
  key_scenario_summary <- make_key_scenario_summary(
    summary_df = summary_df,
    check_final_size = base_args$check_final_size
  )

  write.csv(weekly_diagnostics$zero_filled_summary,
            file = file.path(output_dir, paste0(scenario_id, "_weekly_incidence_zero_filled_summary.csv")),
            row.names = FALSE)
  write.csv(weekly_diagnostics$active_only_summary,
            file = file.path(output_dir, paste0(scenario_id, "_weekly_incidence_active_only_summary.csv")),
            row.names = FALSE)
  write.csv(weekly_diagnostics$activity_summary,
            file = file.path(output_dir, paste0(scenario_id, "_weekly_activity_summary.csv")),
            row.names = FALSE)
  write.csv(weekly_histogram_summary,
            file = file.path(output_dir, paste0(scenario_id, "_weekly_case_histogram_summary.csv")),
            row.names = FALSE)
  write.csv(escape_probability_summary,
            file = file.path(output_dir, paste0(scenario_id, "_escape_probability_summary.csv")),
            row.names = FALSE)
  write.csv(key_scenario_summary,
            file = file.path(output_dir, paste0(scenario_id, "_key_scenario_summary.csv")),
            row.names = FALSE)

  cumulative_summary <- make_cumulative_case_summary(
    weekly_diagnostics$zero_filled_full
  )
  final_size_diagnostics <- make_final_size_diagnostics(
    summary_df = summary_df,
    seeding_cases = base_args$seeding_cases,
    check_final_size = base_args$check_final_size
  )

  write.csv(cumulative_summary,
            file = file.path(output_dir, paste0(scenario_id, "_cumulative_cases_summary.csv")),
            row.names = FALSE)
  write.csv(final_size_diagnostics,
            file = file.path(output_dir, paste0(scenario_id, "_final_size_diagnostics.csv")),
            row.names = FALSE)
  write.csv(quick_summary,
            file = file.path(output_dir, paste0(scenario_id, "_quick_summary.csv")),
            row.names = FALSE)

  # Main plot: unconditional weekly incidence. Replicates that have gone extinct
  # or have no infections in a given week contribute zeroes.
  plot_weekly_summary(
    summary_df = weekly_diagnostics$zero_filled_summary,
    output_path = file.path(output_dir, paste0(scenario_id, "_weekly_incidence.png")),
    main = paste0(scenario_label, "\nUnconditional incidence; zero-filled extinct/quiet weeks")
  )
  plot_weekly_summary(
    summary_df = weekly_diagnostics$zero_filled_summary,
    output_path = file.path(output_dir, paste0(scenario_id, "_weekly_incidence_zero_filled.png")),
    main = paste0(scenario_label, "\nUnconditional incidence; zero-filled extinct/quiet weeks")
  )
  plot_weekly_summary(
    summary_df = weekly_diagnostics$active_only_summary,
    output_path = file.path(output_dir, paste0(scenario_id, "_weekly_incidence_active_only.png")),
    main = paste0(scenario_label, "\nConditional incidence among replicate-weeks with cases")
  )
  plot_activity_summary(
    activity_df = weekly_diagnostics$activity_summary,
    output_path = file.path(output_dir, paste0(scenario_id, "_weekly_activity_summary.png")),
    main = paste0(scenario_label, "\nSimulation activity over time")
  )
  plot_zero_case_proportion(
    activity_df = weekly_diagnostics$activity_summary,
    output_path = file.path(output_dir, paste0(scenario_id, "_zero_case_proportion.png")),
    main = paste0(scenario_label, "\nProportion with zero weekly incidence / ended outbreaks")
  )
  plot_weekly_case_histogram_with_extinction(
    hist_df = weekly_histogram_summary,
    output_path = file.path(output_dir, paste0(scenario_id, "_weekly_case_histogram_with_extinction.png")),
    main = paste0(scenario_label, "\nWeekly case histogram with extinction line")
  )
  plot_smooth_epidemic_curve(
    hist_df = weekly_histogram_summary,
    output_path = file.path(output_dir, paste0(scenario_id, "_smooth_epidemic_curve.png")),
    main = paste0(scenario_label, "\nSmoothed mean weekly epidemic curve")
  )
  plot_cumulative_case_summary(
    cumulative_df = cumulative_summary,
    output_path = file.path(output_dir, paste0(scenario_id, "_cumulative_cases.png")),
    main = paste0(scenario_label, "\nCumulative cases over time; zero-filled extinct/quiet weeks")
  )
  plot_final_size_histogram(
    summary_df = summary_df,
    output_path = file.path(output_dir, paste0(scenario_id, "_final_size_histogram.png")),
    main = paste0(scenario_label, "\nDistribution of final outbreak sizes"),
    check_final_size = base_args$check_final_size
  )
  plot_final_size_histogram_log10(
    summary_df = summary_df,
    output_path = file.path(output_dir, paste0(scenario_id, "_final_size_histogram_log10.png")),
    main = paste0(scenario_label, "\nDistribution of final outbreak sizes, log10 scale"),
    check_final_size = base_args$check_final_size
  )

  if (save_individual_trees) {
    saveRDS(trees, file = file.path(output_dir, paste0(scenario_id, "_trees.rds")))
  }

  print(quick_summary)

  invisible(list(
    scenario_matrix = scenario_matrix,
    summary = summary_df,
    weekly_incidence_active_only_raw = weekly_df,
    weekly_diagnostics = weekly_diagnostics,
    weekly_histogram_summary = weekly_histogram_summary,
    escape_probability_summary = escape_probability_summary,
    key_scenario_summary = key_scenario_summary,
    quick_summary = quick_summary,
    output_dir = output_dir
  ))
}

abc_summarise <- function(out) {
  tdf <- out$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]

  if (nrow(tdf) == 0L) {
    return(c(n_cases = 0, n_deaths = 0, n_hcw_deaths = 0, duration = 0))   # <-- added n_cases = 0
  }

  deaths <- !is.na(tdf$outcome) & tdf$outcome
  hcw    <- !is.na(tdf$class) & tdf$class == "HCW"

  n_cases      <- nrow(tdf)
  n_deaths     <- sum(deaths)
  n_hcw_deaths <- sum(deaths & hcw)

  if (n_deaths == 0L) {
    duration <- 0
  } else {
    t_first_death <- min(tdf$time_outcome_absolute[deaths], na.rm = TRUE)
    t_last_event  <- max(tdf$time_outcome_absolute,         na.rm = TRUE)
    duration      <- max(t_last_event - t_first_death, 0)
  }

  c(n_cases = n_cases, n_deaths = n_deaths, n_hcw_deaths = n_hcw_deaths, duration = duration)
}
