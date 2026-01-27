# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`fiber` is an R package for simulating filovirus outbreaks (Ebola, Marburg) using a branching-process model. It models the impact of antivirals and other medical countermeasures (MCMs) on outbreak control.

## Common Commands

```r
# Load package for development
devtools::load_all()

# Run tests
devtools::test()

# Regenerate documentation from Roxygen2 comments
devtools::document()

# Check package for issues
devtools::check()

# Build package
devtools::build()

# Install from GitHub
devtools::install_github("petal-code/fiber")
```

## Architecture

### Core Simulation Flow

The main entry point is `branching_process_main()` in `R/branching_process_main.R`. It orchestrates the simulation by:

1. Creating seeding cases (initial infections)
2. Iteratively generating offspring (secondary infections) from active cases
3. Completing offspring information (symptoms, hospitalization, death/recovery)
4. Tracking transmission chains until outbreak ends or size limit reached

### Offspring Generation (Three Transmission Tiers)

Each infected person generates offspring based on their class:

- **`offspring_function_genPop()`** - General population parents; offspring can be infected in community or hospital
- **`offspring_function_hcw()`** - Healthcare worker parents; includes workplace transmission modeling
- **`offspring_function_funeral()`** - Unsafe funeral transmission from deceased cases

All use Negative Binomial distribution for offspring count and truncated Gamma for infection timing.

### Key Support Functions

- **`complete_offspring_info()`** (`R/complete_offspring_info.R`) - Fills in offspring details: symptomatic status, hospitalization, death/recovery outcomes, delay times
- **`helper_functions.R`** - Utilities including `rtrunc_gamma()`, probability calculations

### Medical Countermeasures (MCMs)

Three intervention mechanisms with efficacy parameters:
- Hospital quarantine (reduces post-admission transmission)
- PPE/IPC for HCWs (reduces pre-admission hospital transmission)
- Safe burial practices (reduces funeral transmission)

### Infection Locations

Cases are tracked by where infection occurred: `community`, `hospital`, or `funeral`

## Testing

Uses testthat (edition 3). Test files go in `tests/testthat/`. Currently no tests are implemented.

## Notes

- Function documentation uses Roxygen2 (`#'` comments generate `man/*.Rd` files)
- Analysis helper functions are in `inst/` (not exported)
- Project status: WIP - active development on `main_function_offspring_included` branch
