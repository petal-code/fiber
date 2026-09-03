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

### Offspring Generation (Three Transmission Routes)

Transmission is **contact-first**: the Negative Binomial draw in each offspring function is the
number of **contacts** (exposure events), not infections. Each contact is assigned a risk tier and
then transmits or not. Each infected person generates contacts based on their class:

- **`offspring_function_genPop()`** - General population parents; contacts in community or hospital
- **`offspring_function_hcw()`** - Healthcare worker parents; includes workplace transmission modeling
- **`offspring_function_funeral()`** - Funeral contacts of deceased cases

All use Negative Binomial for the contact count and truncated Gamma for contact timing.

Each route is parameterised independently with `mn_contacts_*`, `overdisp_contacts_*` and
`baseline_risk_*`. The NB mean may be a scalar or a time-varying `function(t)` of absolute calendar
time (e.g. via `make_time_varying()`), resolved once per parent immediately before the NB draw:
genPop/HCW at the **parent's infection time**, funeral at the **parent's death (outcome) time**
(since the funeral occurs then).

A tier-`l` contact transmits with probability `baseline_risk(t) * relative_risk[l]`, resolved at
**that contact's own** calendar time. Contacts that transmit then pass through the intervention
layers below. Because tiers are drawn independently per contact, a contact's marginal transmission
probability is `baseline_risk * mean_relative_risk` whatever tier it lands in — so with a flat tier
structure and `baseline_risk = 1` the model reduces exactly to a direct offspring draw.

### Contact Risk Tiers (`R/contact_risk.R`)

`make_contact_risk()` builds the tier structure: `fractions` (tier shares, summing to 1),
`relative_risk`, per-tier `trace_prob`, and a `reference` tier. Relative risks are normalised by the
reference tier's own value, so `baseline_risk` is literally that tier's per-contact transmission
probability. Default is five equal, flat tiers with no tracing. `contact_risk_gradient()` builds a
log-spaced gradient from a single `ratio`, optionally with a `trace_prob_range`.

Structures are per route (`contact_risk_genPop` / `_hcw` / `_funeral`), each falling back to a shared
`contact_risk`. Typical use: set `contact_risk` once for genPop and HCW, override
`contact_risk_funeral`.

Two derived quantities matter. `mean_relative_risk` is the fraction-weighted mean and scales R0.
`case_weights` is the tier distribution among **realised infections** — proportional to
`fractions * relative_risk`, so cases over-represent high-risk tiers. Anything downstream that
depends on the tier a case was infected in (i.e. contact tracing) must average over `case_weights`,
not the raw fractions.

The top tier's probability must stay ≤ 1, so `baseline_risk * max_relative_risk <= 1`. This is
checked up front across the simulation horizon rather than clipped mid-run.

### R0 Calibration (`R/approx_r0.R`)

Single-type (genPop-dominant) approximation at t = 0, matching the ABC calibration scripts:

```
R0_direct  = mn_contacts_genPop  * baseline_risk_genPop  * rr_bar_genPop  * D
R0_funeral = mn_contacts_funeral * baseline_risk_funeral * rr_bar_funeral * F
D = 1 - hospital_quarantine_efficacy(0) * Q_g - isolation_efficacy * Q_iso
F = p_die_comm * (1 - safe_eff * (1 - p_unsafe_comm)) + p_die_hosp * (1 - safe_eff * (1 - p_unsafe_hosp))
```

`Q_g` is the expected fraction of generation-time mass falling after admission; `Q_iso` the fraction
between isolation onset and admission. The two windows are **disjoint**, so the terms are additive
with no double counting. Contact overdispersion never enters R0 (thinning an NB leaves its mean
unchanged), so superspreading can be dialled without disturbing the calibration.

`compute_r0_invariants()` runs the Monte Carlo and returns only efficacy-**independent** quantities,
so an ABC loop caches it once and recomputes the cheap closed-form `r0_direct_multiplier()` /
`r0_funeral_multiplier()` per particle. `solve_baseline_risk_for_r0()` inverts a target
`(R0, prop_funeral)` into baseline risks; `branching_process_main()` accepts `r0_target` /
`r0_prop_funeral` directly and reports the solved values in `sim_info`. Infeasible targets error with
the achievable ceiling named.

### Key Support Functions

- **`complete_offspring_info()`** (`R/complete_offspring_info.R`) - Fills in offspring details: symptomatic status, hospitalization, death/recovery outcomes, delay times
- **`helper_functions.R`** - Utilities including `rtrunc_gamma()`, probability calculations

### Contact Tracing and Isolation

Contact tracing is the channel by which the risk tiers drive the NPIs. Every contact is independently
marked traced with probability `trace_coverage(t) * trace_prob[tier]` — a time-varying programme
coverage lever times a fixed per-tier traceability, following the coverage × efficacy pattern used
throughout. Tracing is drawn for **every** contact, not just those that transmit, so the contact log
carries the true programme denominator.

Traced status travels with any contact that becomes a case. In `complete_offspring_info()` a traced,
symptomatic case then:

- enters **isolation** with probability `prob_isolate_given_traced(t)`, starting `onset_to_isolation()`
  after symptom onset (default: at onset). Isolation thins that case's **pre-admission** transmission
  by `isolation_efficacy` — for HCW parents this covers workplace as well as community contacts, since
  an isolated HCW is off work. Isolation runs up to admission, where hospital quarantine takes over.
- may be admitted more often (`prob_hospitalised_multiplier_traced`, capped at 1) and sooner
  (`hospitalisation_delay_factor_traced`, an extra multiplier on the admission delay).

All three default to no effect. Asymptomatic cases are not isolated: isolation triggers on symptom
onset in a monitored contact, not on the contact event — a quarantine-all-traced-contacts policy
would need a different trigger.

Note this is the model's **only** pre-admission transmission-reducing state; before it, everything
between infection and admission was unthinned.

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

   DPC (days post challenge/exposure to first dose, the argument to the fixed `efficacy(dpc)` curve) is anchored at each candidate's *own* infection time. `obv_pep_dpc` is resolved at the candidate's **absolute infection time** purely so the value can be *calendar*-time-varying (e.g. delays shrink as the response matures); it is the realised delay, not the calendar clock. By default DPC is deterministic (= the resolved `obv_pep_dpc`). Supplying `obv_pep_dpc_sd` (a fixed positive scalar) instead makes `obv_pep_dpc(t)` the **mean** of a per-recipient Gamma draw with that mean and a fixed standard deviation `obv_pep_dpc_sd`, reparameterised as `Gamma(shape = mean²/sd², scale = sd²/mean)` (mean `obv_pep_dpc(t)`, variance `sd²`, `CV = sd/mean`), giving individual variation in time-to-treatment. The spread is a **fixed absolute** sd even when the mean is time-varying (a fixed shape would instead hold the CV constant). `obv_pep_dpc_sd = NULL` (default) draws no Gamma and reproduces prior runs bit-for-bit; a zero mean degenerates to a point mass at DPC 0.
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

### Outputs

`branching_process_main()` returns `tdf`, `contact_log`, `prevented_completed` and `sim_info`.

`tdf` is the transmission tree as before, plus the risk tier of the contact that produced each case
(`contact_risk_level`, `contact_risk_category`), whether that contact was traced (`traced`), and the
case's isolation state (`isolated`, `time_isolation_relative`, `time_isolation_absolute`). Seed cases
have no tier and are never traced.

`contact_log` has one row per contact generated over the whole run — including contacts that never
became infections. Columns: `parent`, `case_id` (joins to `tdf$id` for contacts that became cases, NA
otherwise), `record_type` (`"contact"` / `"infection"`), `class`, `infection_location`,
`time_contact_relative`, `time_contact_absolute`, `contact_risk_level`, `contact_risk_category`,
`relative_risk`, `transmission_prob`, `traced`, and `blocked_by` — NA for infections, else
`"no_transmission"`, the route's intervention layer, or `"obv_pep"`. This is the denominator for
anything tracing- or prophylaxis-related.

`sim_info` additionally carries the resolved per-route risk structures, the baseline risks actually
used (solved from `r0_target` where applicable), and the R0 inversion diagnostics.

`summarise_output()` takes an optional `contact_log` argument and then reports contact counts by tier
and location, the realised per-tier attack rate, tracing/isolation counts, and the `blocked_by`
breakdown.

## Testing

Uses testthat (edition 3). Test files go in `tests/testthat/`. Current suites:
- `test-contact_risk.R`
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
