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

The NB mean (`mn_offspring_genPop`, `mn_offspring_hcw`, `mn_offspring_funeral`) may be a scalar or a
time-varying `function(t)` of absolute calendar time (e.g. via `make_time_varying()`). It is resolved
to a single value per parent immediately before the NB draw: genPop/HCW at the **parent's infection
time**, funeral at the **parent's death (outcome) time** (since the funeral occurs then). A constant
function reproduces the equivalent scalar run bit-for-bit under a fixed seed (no extra RNG draws).

### Key Support Functions

- **`complete_offspring_info()`** (`R/complete_offspring_info.R`) - Fills in offspring details: symptomatic status, hospitalization, death/recovery outcomes, delay times
- **`helper_functions.R`** - Utilities including `rtrunc_gamma()`, probability calculations

### Medical Countermeasures (MCMs)

Four intervention mechanisms with efficacy parameters:
- Hospital quarantine (reduces post-admission transmission)
- PPE/IPC for HCWs (reduces pre-admission hospital transmission, plus receiver-PPE for HCW recipients in hospital)
- Safe burial practices (reduces funeral transmission)
- **Obeldesivir (OBV) PEP** — post-exposure prophylaxis for exposed candidates in target class/locations (default: HCWs at hospital)

Hospital keep-probability composes the first two multiplicatively (Swiss-cheese):
`(1 − source_PPE) × (1 − receiver_PPE) × (1 − hospital_quarantine)`.
OBV PEP is then applied as an additional, independent thinning step on top of that.

### OBV PEP details

Implemented as a two-phase gate (`apply_obv_pep_gate()` in `R/helper_functions.R`):

1. **Pre-thinning phase**: among all candidate exposures in `obv_pep_target_class` × `obv_pep_target_locations`, draw treatment status (received, adherent) and resolve DPC. These populate the `pre_*` counters — interpretable as "Policy A: treat all contacts" doses.
2. **Post-thinning phase**: for candidates that *also* survived PPE/quarantine thinning, the `post_*` counters are populated as the subset of pre-thinning treated/adherent — interpretable as "Policy B: treat only PPE failures" doses. For kept-and-adherent candidates, `Bernoulli(efficacy(dpc))` decides whether infection is prevented.

A given candidate has *one* treatment status (consistency property): per-individual `post_* ⊆ pre_*` as sets.

Seven gate counters surface in `summarise_output()`:
- `n_obv_pep_pre_eligible`, `n_obv_pep_pre_treated`, `n_obv_pep_pre_adherent`
- `n_obv_pep_post_eligible`, `n_obv_pep_post_treated`, `n_obv_pep_post_adherent`
- `n_obv_pep_prevented`

Plus three tdf-based cohort counters of HCW cases who became cases despite being in the eligible / treated / adherent cohort (`n_obv_pep_eligible_cases`, `n_obv_pep_treated_cases`, `n_obv_pep_breakthroughs`). Only `n_obv_pep_breakthroughs` (received AND adhered AND still infected) is a clinical breakthrough.

Default efficacy curve `obv_pep_efficacy_from_dpc()` is NHP-derived (E0 = 0.82 at dpc=0, decays to 0 at dpc ≥ 10).

### Infection Locations

Cases are tracked by where infection occurred: `community`, `hospital`, or `funeral`

## Testing

Uses testthat (edition 3). Test files go in `tests/testthat/`. Current suites:
- `test-make_time_varying.R`
- `test-time_varying_hospitalisation.R`
- `test-time_varying_transmissibility.R`
- `test-ppe_thinning.R`
- `test-conditional_cfr.R`
- `test-obv_pep.R`
- `test-compute_reproduction_number.R`

## Notes

- Function documentation uses Roxygen2 (`#'` comments generate `man/*.Rd` files)
- Analysis helper functions are in `inst/` (not exported)
- Project status: WIP - active development on `main_function_offspring_included` branch
