# run_single_default_simulation.R
# -----------------------------------------------------------------------------
# Minimal single-run script for the FIBER Ebola branching-process model.
#
# Purpose:
#   Run ONE stochastic simulation using the default scalar parameters and the
#   selected time-varying scenario curves from 00_common_time_varying_scenario_setup.R.
#
# Expected folder structure:
#   run_single_default_simulation.R
#   00_common_time_varying_scenario_setup.R
#   final_four_scenario_values.csv
#   model_functions/
#     branching_process_main.R
#     complete_offspring_info.R
#     make_time_varying.R
#     offspring_function_funeral.R
#     offspring_function_genPop.R
#     offspring_function_hcw.R
#     resolve_time_varying.R
#     ...
# -----------------------------------------------------------------------------

rm(list = ls())

# -----------------------------------------------------------------------------
# User choices for this single run
# -----------------------------------------------------------------------------

SCENARIO_ID <- "worst_west_africa"
SEED <- 12345
OUTPUT_DIR <- file.path("outputs", "single_default_simulation")

# The run uses the default scalar values defined in 00_common_time_varying_scenario_setup.R.
# To inspect those values after sourcing 00_common, run:
#   DEFAULT_SCALAR_INPUTS

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

required_files <- c(
  "00_common_time_varying_scenario_setup.R",
  "final_four_scenario_values.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing required file(s): ", paste(missing_files, collapse = ", "),
    "\nSet the working directory to the run folder containing 00_common_time_varying_scenario_setup.R.",
    call. = FALSE
  )
}

if (!dir.exists("model_functions")) {
  stop("Missing required folder: model_functions/", call. = FALSE)
}

source("00_common_time_varying_scenario_setup.R")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Build one complete argument list
# -----------------------------------------------------------------------------

scenario_matrix <- read_scenario_matrix("final_four_scenario_values.csv")

if (!SCENARIO_ID %in% unique(scenario_matrix$scenario)) {
  stop(
    "SCENARIO_ID = '", SCENARIO_ID, "' not found in final_four_scenario_values.csv.\n",
    "Available scenarios: ", paste(unique(scenario_matrix$scenario), collapse = ", "),
    call. = FALSE
  )
}

base_args <- make_base_args(
  scalar_inputs = DEFAULT_SCALAR_INPUTS
)

scenario_args <- build_time_varying_args(
  scenario_id = SCENARIO_ID,
  matrix = scenario_matrix,
  etu_efficacy_baseline = DEFAULT_SCALAR_INPUTS$etu_efficacy_baseline
)

scenario_label <- scenario_args$scenario_label
scenario_args_for_model <- scenario_args[setdiff(names(scenario_args), c("scenario_label", "scenario_matrix"))]

model_args <- c(
  base_args,
  scenario_args_for_model,
  list(seed = SEED)
)

# Save the exact default scalar values used for this run.
scalar_defaults_df <- data.frame(
  parameter = names(DEFAULT_SCALAR_INPUTS),
  value = vapply(DEFAULT_SCALAR_INPUTS, function(x) paste(x, collapse = ";"), character(1)),
  stringsAsFactors = FALSE
)
write.csv(
  scalar_defaults_df,
  file.path(OUTPUT_DIR, "scalar_defaults_used.csv"),
  row.names = FALSE
)

# Save the selected scenario trajectory for traceability.
write.csv(
  scenario_args$scenario_matrix,
  file.path(OUTPUT_DIR, paste0(SCENARIO_ID, "_time_varying_matrix_used.csv")),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Run one stochastic simulation
# -----------------------------------------------------------------------------

message("Running one stochastic simulation")
message("Scenario: ", SCENARIO_ID, " (", scenario_label, ")")
message("Seed: ", SEED)
message("Seed cases: ", model_args$seeding_cases)
message("Final-size cap: ", model_args$check_final_size)

out <- do.call(branching_process_main, model_args)

# -----------------------------------------------------------------------------
# Summarise and save outputs
# -----------------------------------------------------------------------------

summary_df <- summarise_one_outbreak(
  out = out,
  scenario_id = SCENARIO_ID,
  scenario_label = scenario_label,
  replicate = 1L,
  check_final_size = model_args$check_final_size
)

weekly_df <- make_weekly_incidence(
  out = out,
  scenario_id = SCENARIO_ID,
  scenario_label = scenario_label,
  replicate = 1L
)

tdf <- out$tdf
tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]

saveRDS(out, file.path(OUTPUT_DIR, "single_default_simulation_output.rds"))
write.csv(tdf, file.path(OUTPUT_DIR, "single_default_simulation_transmission_tree.csv"), row.names = FALSE)
write.csv(summary_df, file.path(OUTPUT_DIR, "single_default_simulation_summary.csv"), row.names = FALSE)
write.csv(weekly_df, file.path(OUTPUT_DIR, "single_default_simulation_weekly_incidence.csv"), row.names = FALSE)

# Simple base-R incidence plot; avoids extra package dependencies.
if (nrow(weekly_df) > 0L) {
  png(
    filename = file.path(OUTPUT_DIR, "single_default_simulation_weekly_incidence.png"),
    width = 1100,
    height = 700,
    res = 120
  )
  barplot(
    height = weekly_df$cases,
    names.arg = weekly_df$week_start_day / 7,
    border = NA,
    main = paste0(scenario_label, ": single stochastic simulation"),
    xlab = "Weeks since first infection",
    ylab = "Weekly infections",
    las = 2
  )
  dev.off()
}

message("Done. Outputs written to: ", normalizePath(OUTPUT_DIR))
print(summary_df)
