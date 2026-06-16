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

Four intervention mechanisms. The three NPIs are each parameterised as a
time-varying **coverage** lever (probability the intervention reaches a given
exposure) times a fixed scalar **efficacy** (effect conditional on reaching it).
The intended fitting workflow pre-specifies the coverage curves from the
literature and varies the efficacies.

- **Hospital / ETU quarantine** (reduces post-admission transmission). Coverage =
  `prop_etu(t)`, the fraction of hospitalised cases managed in an ETU vs a general
  hospital. The post-admission quarantine efficacy is the mixture
  `prop_etu(t)·etu_efficacy + (1 − prop_etu(t))·general_hospital_quarantine_efficacy`,
  where `etu_efficacy` and `general_hospital_quarantine_efficacy` are fixed,
  independently-togglable scalars (no ordering enforced).
- **PPE/IPC for HCWs** (reduces pre-admission hospital transmission, plus
  receiver-PPE for HCW recipients in hospital). Each PPE layer thins by
  `ppe_coverage_hcw(t) · ppe_efficacy`, with `ppe_coverage_hcw` the time-varying
  coverage and `ppe_efficacy` a fixed scalar shared by the source and receiver
  layers.
- **Safe burial practices** (reduces funeral transmission). Coverage =
  `1 − p_unsafe_funeral_*(t)`; efficacy = `safe_funeral_efficacy` (fixed scalar).
- **Obeldesivir (OBV) PEP** — post-exposure prophylaxis for exposed candidates in
  target class/locations (default: HCWs at hospital). Pharmaceutical, not an NPI.

Hospital keep-probability composes the PPE and quarantine layers multiplicatively
(Swiss-cheese): `(1 − source_PPE) × (1 − receiver_PPE) × (1 − hospital_quarantine)`,
where each PPE factor is `ppe_coverage_hcw(t)·ppe_efficacy` and `hospital_quarantine`
is the ETU/general mixture above. OBV PEP is then applied as an additional,
independent thinning step on top of that.

### OBV PEP details

Implemented as a two-phase gate (`apply_obv_pep_gate()` in `R/helper_functions.R`):

1. **Pre-thinning phase**: among all candidate exposures in `obv_pep_target_class` × `obv_pep_target_locations`, draw treatment status (received, adherent) and resolve DPC. These populate the `pre_*` counters — interpretable as "Policy A: treat all contacts" doses.

   DPC (days post challenge/exposure to first dose, the argument to the fixed `efficacy(dpc)` curve) is anchored at each candidate's *own* infection time. `obv_pep_dpc` is resolved at the candidate's **absolute infection time** purely so the value can be *calendar*-time-varying (e.g. delays shrink as the response matures); it is the realised delay, not the calendar clock. By default DPC is deterministic (= the resolved `obv_pep_dpc`). Supplying `obv_pep_dpc_shape` (a fixed positive scalar) instead makes `obv_pep_dpc(t)` the **mean** of a per-recipient `Gamma(shape, scale = mean/shape)` draw (mean `obv_pep_dpc(t)`, variance `mean²/shape`, `CV = 1/√shape`), giving individual variation in time-to-treatment. `obv_pep_dpc_shape = NULL` (default) draws no Gamma and reproduces prior runs bit-for-bit; a zero mean degenerates to a point mass at DPC 0.
2. **Post-thinning phase**: for candidates that *also* survived PPE/quarantine thinning, the `post_*` counters are populated as the subset of pre-thinning treated/adherent — interpretable as "Policy B: treat only PPE failures" doses. For kept-and-adherent candidates, `Bernoulli(efficacy(dpc))` decides whether infection is prevented.

A given candidate has *one* treatment status (consistency property): per-individual `post_* ⊆ pre_*` as sets.

Seven gate counters surface in `summarise_output()`:
- `n_obv_pep_pre_eligible`, `n_obv_pep_pre_treated`, `n_obv_pep_pre_adherent`
- `n_obv_pep_post_eligible`, `n_obv_pep_post_treated`, `n_obv_pep_post_adherent`
- `n_obv_pep_prevented`

Plus three tdf-based cohort counters of HCW cases who became cases despite being in the eligible / treated / adherent cohort (`n_obv_pep_eligible_cases`, `n_obv_pep_treated_cases`, `n_obv_pep_breakthroughs`). Only `n_obv_pep_breakthroughs` (received AND adhered AND still infected) is a clinical breakthrough.

Plus a deferred counterfactual counter `n_obv_pep_prevented_deaths` (≤ `n_obv_pep_prevented`): the subset of prevented infections that *would have died* had they occurred. The infections OBV blocks are removed before `complete_offspring_info()` resolves any outcome, so their would-be death status is constructed after the simulation loop by replaying the stashed prevented infections through the **same** outcome model as realised cases (`complete_offspring_info()`: symptomatic → potential hospitalisation → community CFR → hospital second-chance) against a zero-time dummy parent. Like `n_obv_pep_prevented`, this is a *direct* count (the prevented index infections only, excluding their averted onward chains), and is therefore a **lower bound** on the total deaths the programme averts. The replay runs after the transmission tree is finalised and is skipped entirely when nothing was prevented, so it never perturbs the simulated trajectory's RNG stream — letting deaths averted by OBV be read off a single run instead of differencing a separate (stochastically divergent) no-OBV run.

The completed replay frame itself is returned as `out$prevented_completed` (alongside `out$tdf` and `out$sim_info`) — one row per averted index infection, carrying its counterfactual natural history (`time_infection_absolute` is the would-be infection time; `outcome` the would-be death/recovery status, whose sum is `prevented_deaths`). It is `NULL` when nothing was prevented. Because the replay uses a zero-time dummy parent, `parent`/`generation` are `NA` and `time_infection_relative` equals `time_infection_absolute`.

Built-in efficacy curve `obv_pep_efficacy_from_dpc()` is **flat by default** (`shape = "flat"`): efficacy is a constant `E0` at every DPC. This is a deliberately model-neutral default — the prior logistic constants were placeholders. Pass `shape = "logistic"` to opt into the NHP-derived scaled-logistic decay (E0 at dpc=0 → 0 by `dpc_zero`, cut to 0 past `max_dpc`); the normalised logistic **cannot** be made flat by any `d50/k/dpc_zero` (as `k→0` it tends to a *linear* decay, and `k=0` is degenerate), which is why flat is a distinct mode rather than a parameter setting. The curve is **separable**: `eff(dpc) = E0 · shape(dpc)` with `shape(0) = 1`, so `E0` is a pure vertical scale = the peak (dpc=0) efficacy. The `obv_pep_efficacy` argument selects the model: **NULL or a scalar** use the built-in curve (a scalar IS the constant/`E0`, so a bare scalar gives a flat efficacy); a **function(dpc)** is used as-is. The curve parameters (`shape`, `E0`, `d50`, `k`, `dpc_zero`, `max_dpc`) can be swept without writing a closure by passing `obv_pep_efficacy_args` (a named list, e.g. `list(shape = "logistic", d50 = 4)`) to `branching_process_main()`; it applies to the built-in curve only (NULL or scalar `obv_pep_efficacy`), errors if combined with a function `obv_pep_efficacy` / if it names `E0` while `obv_pep_efficacy` is a scalar / given an unknown name, and is validated up front (e.g. `E0 ∈ [0, 1]`) so a bad sweep point fails before the run.

### Infection Locations

Cases are tracked by where infection occurred: `community`, `hospital`, or `funeral`

## Testing

Uses testthat (edition 3). Test files go in `tests/testthat/`. Current suites:
- `test-make_time_varying.R`
- `test-time_varying_hospitalisation.R`
- `test-time_varying_transmissibility.R`
- `test-ppe_thinning.R`
- `test-hospital_quarantine_efficacy.R`
- `test-conditional_cfr.R`
- `test-obv_pep.R`
- `test-compute_reproduction_number.R`

## Notes

- Function documentation uses Roxygen2 (`#'` comments generate `man/*.Rd` files)
- Analysis helper functions are in `inst/` (not exported)
- Project status: WIP - active development on `main_function_offspring_included` branch
