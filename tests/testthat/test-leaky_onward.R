## Tests for estimate_leaky_onward(): the post-hoc "leaky OBV" sensitivity
## estimator. Covers:
##   - r = 0 boundary (added = averted index cases, deaths = prevented_deaths,
##     zero onward) -- the clean nesting with the sterilizing model
##   - r > 0 produces onward transmission that increases with r (monotone)
##   - accounting identities (added = seeds + onward on every route)
##   - NULL / empty prevented_completed fast path
##   - reproducibility under a fixed seed and non-perturbation of the main run
##   - input validation
##   - sim_info$params round-trips the run's parameters

## --- helpers -----------------------------------------------------------

leaky_bpm_args <- function(...) {
  fixed <- function(v) function(n) rep(v, n)
  defaults <- list(
    mn_offspring_genPop           = 1.5,
    overdisp_offspring_genPop     = 0.5,
    Tg_shape_genPop               = 2,
    Tg_rate_genPop                = 0.15,
    mn_offspring_hcw              = 1.5,
    overdisp_offspring_hcw        = 0.5,
    Tg_shape_hcw                  = 2,
    Tg_rate_hcw                   = 0.15,
    mn_offspring_funeral          = 1.5,
    overdisp_offspring_funeral    = 0.5,
    Tg_shape_funeral              = 10,
    Tg_rate_funeral               = 5,
    incubation_period             = fixed(5),
    onset_to_hospitalisation      = fixed(3),
    onset_to_death                = fixed(7),
    onset_to_recovery             = fixed(10),
    hospitalisation_to_death      = fixed(7),
    hospitalisation_to_recovery   = fixed(10),
    prob_symptomatic              = 0.9,
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
    ppe_efficacy                  = 1,
    prop_etu                      = 1,
    etu_efficacy                  = 0.5,
    general_hospital_quarantine_efficacy = 0.5,
    p_unsafe_funeral_comm_hcw     = 0.5,
    p_unsafe_funeral_hosp_hcw     = 0.2,
    p_unsafe_funeral_comm_genPop  = 0.5,
    p_unsafe_funeral_hosp_genPop  = 0.2,
    safe_funeral_efficacy         = 1.0,
    prob_hcw_cond_funeral_hcw     = 0.05,
    prob_hcw_cond_funeral_genPop  = 0.05,
    population                    = 5000,
    hcw_per_capita                = 0.02,
    check_final_size              = 200,
    seeding_cases                 = 20,
    seed                          = 1L
  )
  utils::modifyList(defaults, list(...))
}

## Full-prevention, HCW-at-hospital-dominant scenario (mirrors test-obv_pep.R):
## engineered so the default HCW-at-hospital target reliably yields prevented
## infections to feed the estimator.
leaky_full_prevention_args <- function(...) {
  leaky_bpm_args(
    obv_pep_enabled               = TRUE,
    obv_pep_coverage              = 1,
    obv_pep_adherence             = 1,
    obv_pep_dpc                   = 0,
    obv_pep_efficacy              = 1,
    ppe_coverage_hcw              = 0,
    etu_efficacy                  = 0,
    general_hospital_quarantine_efficacy = 0,
    prob_hospitalised_hcw         = 0.9,
    prob_hospitalised_genPop      = 0.9,
    prob_hcw_cond_genPop_hospital = 0.9,
    prob_hcw_cond_hcw_hospital    = 0.9,
    seeding_cases                 = 20,
    ...
  )
}

## Return the first run (over a seed scan) that actually prevents something, so
## prevented_completed is non-trivial without pinning a fragile seed.
leaky_run_with_prevention <- function(seeds = 1:40) {
  for (sd in seeds) {
    res <- do.call(branching_process_main, leaky_full_prevention_args(seed = sd))
    if (summarise_output(res$tdf, sim_info = res$sim_info)$n_obv_pep_prevented > 0) {
      return(res)
    }
  }
  stop("no seed in 1:40 produced any prevented infections; check the scenario")
}

## Drop the estimator's mean-offspring parameters to a subcritical regime so
## forward-simulated trees stay small and the tests run fast.
subcritical <- function(params, mn = 0.6) {
  utils::modifyList(params, list(mn_offspring_genPop = mn,
                                 mn_offspring_hcw     = mn,
                                 mn_offspring_funeral = mn))
}

## --- sim_info$params round-trip ----------------------------------------

test_that("branching_process_main stashes a usable params bundle in sim_info", {
  res <- do.call(branching_process_main, leaky_bpm_args(seed = 3L))
  p <- res$sim_info$params
  expect_type(p, "list")
  ## A representative spread of the parameters the estimator relies on.
  for (nm in c("mn_offspring_genPop", "mn_offspring_hcw", "mn_offspring_funeral",
               "incubation_period", "onset_to_death", "prob_symptomatic",
               "prob_death_comm", "prob_death_hosp", "prop_etu", "etu_efficacy",
               "obv_pep_enabled", "obv_pep_target_locations", "safe_funeral_efficacy")) {
    expect_true(nm %in% names(p), info = nm)
  }
  expect_identical(p$mn_offspring_genPop, 1.5)
  ## This run leaves OBV at its default (disabled); the flag round-trips verbatim.
  expect_identical(p$obv_pep_enabled, FALSE)
})

test_that("adding sim_info$params does not perturb the trajectory", {
  args <- leaky_full_prevention_args(seed = 11L)
  r1 <- do.call(branching_process_main, args)
  r2 <- do.call(branching_process_main, args)
  expect_equal(r1$tdf, r2$tdf)
})

## --- r = 0 boundary -----------------------------------------------------

test_that("r = 0 adds back the index cases and their deaths but no onward chains", {
  res <- leaky_run_with_prevention()
  s   <- summarise_output(res$tdf, sim_info = res$sim_info)
  e   <- estimate_leaky_onward(res$prevented_completed, res$sim_info$params,
                               residual_transmissibility = 0,
                               n_replicates = 6, seed = 1L)
  ## Engaged individuals at r = 0 transmit nothing -> no descendants ...
  expect_true(all(e$per_replicate$onward_infections == 0))
  expect_true(all(e$per_replicate$onward_deaths == 0))
  ## ... but they still exist as cases and can die.
  expect_equal(e$n_seeds, s$n_obv_pep_prevented)
  expect_equal(e$seed_deaths, s$n_obv_pep_prevented_deaths)
  expect_true(all(e$per_replicate$added_infections == e$n_seeds))
  expect_true(all(e$per_replicate$added_deaths == e$seed_deaths))
})

## --- r > 0 onward transmission, monotone in r ---------------------------

test_that("onward transmission is positive and increases with residual transmissibility", {
  res <- leaky_run_with_prevention()
  sub <- subcritical(res$sim_info$params)
  e_lo <- estimate_leaky_onward(res$prevented_completed, sub, 0.4,
                                n_replicates = 40, seed = 2L)
  e_hi <- estimate_leaky_onward(res$prevented_completed, sub, 0.9,
                                n_replicates = 40, seed = 2L)
  expect_gt(e_lo$summary$onward_infections$mean, 0)
  expect_gt(e_hi$summary$onward_infections$mean, e_lo$summary$onward_infections$mean)
})

## --- accounting identities ---------------------------------------------

test_that("added counts decompose exactly into seeds + onward", {
  res <- leaky_run_with_prevention()
  sub <- subcritical(res$sim_info$params)
  e <- estimate_leaky_onward(res$prevented_completed, sub, 0.5,
                             n_replicates = 20, seed = 5L)
  pr <- e$per_replicate
  expect_true(all(pr$added_infections == e$n_seeds + pr$onward_infections))
  expect_true(all(pr$added_deaths == e$seed_deaths + pr$onward_deaths))
  ## Deaths never exceed infections, on either scope.
  expect_true(all(pr$onward_deaths <= pr$onward_infections))
  expect_true(all(pr$added_deaths <= pr$added_infections))
})

## --- reproducibility ----------------------------------------------------

test_that("estimate_leaky_onward is reproducible under a fixed seed", {
  res <- leaky_run_with_prevention()
  sub <- subcritical(res$sim_info$params)
  a <- estimate_leaky_onward(res$prevented_completed, sub, 0.6, n_replicates = 10, seed = 7L)
  b <- estimate_leaky_onward(res$prevented_completed, sub, 0.6, n_replicates = 10, seed = 7L)
  expect_equal(a$per_replicate, b$per_replicate)
  expect_equal(a$summary, b$summary)
})

## --- empty / NULL fast path --------------------------------------------

test_that("NULL or empty prevented_completed yields an all-zero result", {
  res <- do.call(branching_process_main, leaky_bpm_args(seed = 3L))  # OBV disabled
  params <- res$sim_info$params

  e_null <- estimate_leaky_onward(NULL, params, 0.5, n_replicates = 3)
  expect_equal(e_null$n_seeds, 0L)
  expect_equal(e_null$seed_deaths, 0L)
  expect_true(all(e_null$per_replicate$added_infections == 0))
  expect_true(all(e_null$per_replicate$onward_infections == 0))

  empty_pc <- res$prevented_completed  # NULL when nothing prevented
  e_empty <- estimate_leaky_onward(empty_pc, params, 0.5, n_replicates = 2)
  expect_equal(e_empty$n_seeds, 0L)
})

## --- return_trees -------------------------------------------------------

test_that("return_trees returns per-replicate frames flagged with engaged", {
  res <- leaky_run_with_prevention()
  sub <- subcritical(res$sim_info$params)
  e <- estimate_leaky_onward(res$prevented_completed, sub, 0.7,
                             n_replicates = 2, seed = 3L, return_trees = TRUE)
  expect_length(e$trees, 2L)
  for (i in seq_along(e$trees)) {
    tr <- e$trees[[i]]
    expect_true("engaged" %in% names(tr))
    ## Tree size matches the reported added infections; the seeds are all engaged.
    expect_equal(nrow(tr), e$per_replicate$added_infections[i])
    expect_gte(sum(tr$engaged), e$n_seeds)
  }
})

## --- input validation ---------------------------------------------------

test_that("estimate_leaky_onward validates its inputs", {
  res <- leaky_run_with_prevention()
  pc  <- res$prevented_completed
  params <- res$sim_info$params

  expect_error(estimate_leaky_onward(pc, params, residual_transmissibility = 1.5),
               "residual_transmissibility")
  expect_error(estimate_leaky_onward(pc, params, residual_transmissibility = -0.1),
               "residual_transmissibility")
  expect_error(estimate_leaky_onward(pc, params, 0.5, n_replicates = 0),
               "n_replicates")
  expect_error(estimate_leaky_onward(pc, params, 0.5, n_replicates = 2.5),
               "n_replicates")
  expect_error(estimate_leaky_onward(pc, "not a list", 0.5),
               "params")
  expect_error(estimate_leaky_onward(pc, params, 0.5, max_tree_size = 0),
               "max_tree_size")
})

test_that("hitting max_tree_size warns rather than silently truncating", {
  res <- leaky_run_with_prevention()
  ## Full (supercritical) means + a tiny cap guarantees the cap is hit.
  expect_warning(
    estimate_leaky_onward(res$prevented_completed, res$sim_info$params, 1.0,
                          n_replicates = 1, max_tree_size = 50, seed = 1L),
    "max_tree_size"
  )
})
