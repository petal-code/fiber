## Tests for the post-admission hospital-quarantine efficacy formulation.
##
## When `hospital_quarantine_efficacy` is not supplied directly, the efficacy is
## derived as a prop_etu(t)-weighted mixture of two fixed scalar efficacies:
##
##   hospital_quarantine_efficacy(t) = prop_etu(t) * etu_efficacy
##                                   + (1 - prop_etu(t)) * general_hospital_quarantine_efficacy
##
## The only time-variation enters through prop_etu(t); both efficacies are fixed
## scalars and are independently togglable (no ordering between them is enforced).
##
## The mixture is pinned by comparing the derived path against the equivalent
## direct `hospital_quarantine_efficacy` scalar. Neither path consumes RNG when
## resolving the efficacy, so equal effective efficacy under a fixed seed gives
## bit-identical draws.

hq_make_parent <- function(time_infection_absolute = 0,
                           time_to_hospitalisation = 0.001,
                           time_to_outcome = 30) {
  data.frame(
    id                             = 1L,
    class                          = "genPop",
    infection_location             = "community",
    parent                         = NA_integer_,
    generation                     = 1L,
    time_infection_relative        = 0,
    time_infection_absolute        = time_infection_absolute,
    incubation_period              = 5,
    symptomatic                    = TRUE,
    time_symptom_onset_relative    = 5,
    time_symptom_onset_absolute    = time_infection_absolute + 5,
    hospitalisation                = TRUE,
    time_hospitalisation_relative  = time_to_hospitalisation,
    time_hospitalisation_absolute  = time_infection_absolute + time_to_hospitalisation,
    outcome                        = TRUE,
    outcome_location               = "hospital",
    time_outcome_relative          = time_to_outcome,
    time_outcome_absolute          = time_infection_absolute + time_to_outcome,
    funeral_safety                 = "unsafe",
    n_offspring                    = NA_integer_,
    offspring_generated            = FALSE,
    stringsAsFactors               = FALSE
  )
}

## A hospitalised genPop source whose hospital recipients are all genPop
## (prob_hcw_cond_genPop_hospital = 0), so the PPE layer never fires and the only
## hospital thinning is the post-admission quarantine layer.
hq_run <- function(seed, ...) {
  set.seed(seed)
  offspring_function_genPop(
    parent_info                   = hq_make_parent(),
    mn_offspring_genPop           = 800,
    overdisp_offspring_genPop     = 400,
    Tg_shape_genPop               = 4,
    Tg_rate_genPop                = 1,
    ppe_coverage_hcw              = 0,
    ppe_efficacy                  = 1,
    prob_hcw_cond_genPop_comm     = 0,
    prob_hcw_cond_genPop_hospital = 0,
    ...
  )
}

test_that("derived ETU/general mixture equals the equivalent direct scalar (prop_etu = 0.5)", {
  ## 0.5 * 0.8 + 0.5 * 0.4 = 0.6
  res_derived <- hq_run(
    seed = 100,
    prop_etu = 0.5,
    etu_efficacy = 0.8,
    general_hospital_quarantine_efficacy = 0.4
  )
  res_direct <- hq_run(seed = 100, hospital_quarantine_efficacy = 0.6)

  expect_identical(res_derived$infection_location, res_direct$infection_location)
  expect_identical(res_derived$time_infection_relative, res_direct$time_infection_relative)
})

test_that("prop_etu = 1 reduces the mixture to the ETU efficacy", {
  res_derived <- hq_run(
    seed = 7,
    prop_etu = 1,
    etu_efficacy = 0.8,
    general_hospital_quarantine_efficacy = 0.4
  )
  res_direct <- hq_run(seed = 7, hospital_quarantine_efficacy = 0.8)
  expect_identical(res_derived$time_infection_relative, res_direct$time_infection_relative)
})

test_that("prop_etu = 0 reduces the mixture to the general-hospital efficacy", {
  res_derived <- hq_run(
    seed = 7,
    prop_etu = 0,
    etu_efficacy = 0.8,
    general_hospital_quarantine_efficacy = 0.4
  )
  res_direct <- hq_run(seed = 7, hospital_quarantine_efficacy = 0.4)
  expect_identical(res_derived$time_infection_relative, res_direct$time_infection_relative)
})

test_that("a constant function(t) prop_etu reproduces the scalar mixture", {
  res_scalar <- hq_run(
    seed = 42,
    prop_etu = 0.3,
    etu_efficacy = 0.9,
    general_hospital_quarantine_efficacy = 0.2
  )
  res_fn <- hq_run(
    seed = 42,
    prop_etu = function(t) 0.3,
    etu_efficacy = 0.9,
    general_hospital_quarantine_efficacy = 0.2
  )
  expect_identical(res_scalar$time_infection_relative, res_fn$time_infection_relative)
  expect_identical(res_scalar$class, res_fn$class)
})

test_that("higher ETU coverage yields stronger quarantine when etu_efficacy > general", {
  ## etu_efficacy (0.9) > general (0.1): more ETU coverage => fewer surviving
  ## hospital events. Compare expected survivors across many seeds.
  surv <- function(prop_etu_val) {
    vapply(1:40, function(s) {
      r <- hq_run(
        seed = s,
        prop_etu = prop_etu_val,
        etu_efficacy = 0.9,
        general_hospital_quarantine_efficacy = 0.1
      )
      sum(r$infection_location == "hospital")
    }, numeric(1))
  }
  expect_lt(mean(surv(1)), mean(surv(0)))
})

## --- Validation -----------------------------------------------------------

test_that("derived path rejects incomplete inputs (missing general efficacy)", {
  expect_error(
    hq_run(seed = 1, prop_etu = 0.5, etu_efficacy = 0.8),
    "general_hospital_quarantine_efficacy"
  )
})

test_that("derived path rejects an out-of-range etu_efficacy", {
  expect_error(
    hq_run(
      seed = 1,
      prop_etu = 0.5,
      etu_efficacy = 1.5,
      general_hospital_quarantine_efficacy = 0.4
    ),
    "etu_efficacy"
  )
})

test_that("branching_process_main errors when derived inputs are incomplete", {
  ## End-to-end: the top-level entry point requires either a direct efficacy or
  ## the full derived set (prop_etu + etu_efficacy + general_*).
  expect_error(
    branching_process_main(
      mn_offspring_genPop           = 1.5,
      overdisp_offspring_genPop     = 5,
      Tg_shape_genPop               = 4,
      Tg_rate_genPop                = 1,
      mn_offspring_hcw              = 1.5,
      overdisp_offspring_hcw        = 5,
      Tg_shape_hcw                  = 4,
      Tg_rate_hcw                   = 1,
      mn_offspring_funeral          = 1.5,
      overdisp_offspring_funeral    = 5,
      Tg_shape_funeral              = 10,
      Tg_rate_funeral               = 5,
      incubation_period             = function(n) rep(5, n),
      onset_to_hospitalisation      = function(n) rep(3, n),
      onset_to_death                = function(n) rep(7, n),
      onset_to_recovery             = function(n) rep(10, n),
      hospitalisation_to_death      = function(n) rep(4, n),
      hospitalisation_to_recovery   = function(n) rep(6, n),
      prob_symptomatic              = 1,
      prob_hospitalised_hcw         = 0.5,
      prob_hospitalised_genPop      = 0.5,
      prob_death_comm               = 0.5,
      prob_death_hosp               = 0.3,
      prob_hcw_cond_genPop_comm     = 0.01,
      prob_hcw_cond_genPop_hospital = 0.3,
      prob_hcw_cond_hcw_comm        = 0.01,
      prob_hcw_cond_hcw_hospital    = 0.3,
      prob_hospital_cond_hcw_preAdm = 0.5,
      ppe_coverage_hcw              = 0.5,
      ppe_efficacy                  = 0.5,
      prop_etu                      = 0.5,
      etu_efficacy                  = 0.8,
      ## general_hospital_quarantine_efficacy deliberately omitted
      p_unsafe_funeral_comm_hcw     = 0.5,
      p_unsafe_funeral_hosp_hcw     = 0.2,
      p_unsafe_funeral_comm_genPop  = 0.5,
      p_unsafe_funeral_hosp_genPop  = 0.2,
      safe_funeral_efficacy         = 1.0,
      prob_hcw_cond_funeral_hcw     = 0.3,
      prob_hcw_cond_funeral_genPop  = 0.1,
      population                    = 1000,
      check_final_size              = 50,
      seeding_cases                 = 1,
      seed                          = 1
    ),
    "general_hospital_quarantine_efficacy"
  )
})
