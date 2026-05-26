## Tests for the OBV PEP gate.
## Covers:
##   - The two-phase semantics of apply_obv_pep_gate()
##   - Disabled / zero-coverage no-op behaviour
##   - Full-prevention regression (coverage = adherence = efficacy = 1 zeroes
##     out HCW infections at hospital)
##   - Boundary behaviour of obv_pep_efficacy_from_dpc()
##   - Set-nesting invariants between the 7 counters

## --- helpers -----------------------------------------------------------

make_parent_info_genPop_hospitalised <- function(time_infection_absolute = 0,
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

obv_off_args <- function() {
  list(
    mn_offspring_genPop          = 500,
    overdisp_offspring_genPop    = 200,
    Tg_shape_genPop              = 4,
    Tg_rate_genPop               = 1,
    hospital_quarantine_efficacy = 0,
    ppe_efficacy_hcw             = 0,
    prob_hcw_cond_genPop_comm    = 0,
    prob_hcw_cond_genPop_hospital = 1
  )
}

## --- obv_pep_efficacy_from_dpc boundaries ------------------------------

test_that("obv_pep_efficacy_from_dpc returns E0 at dpc = 0", {
  expect_equal(obv_pep_efficacy_from_dpc(0), 0.82342697, tolerance = 1e-6)
})

test_that("obv_pep_efficacy_from_dpc is 0 at dpc_zero = 15", {
  expect_equal(obv_pep_efficacy_from_dpc(15), 0)
})

test_that("obv_pep_efficacy_from_dpc is 0 just past max_dpc = 10", {
  expect_equal(obv_pep_efficacy_from_dpc(10.01), 0)
})

test_that("obv_pep_efficacy_from_dpc is in [0, 1] across a sensible range", {
  vals <- obv_pep_efficacy_from_dpc(seq(0, 20, by = 0.5))
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("obv_pep_efficacy_from_dpc is monotone non-increasing on [0, 15]", {
  dpc <- seq(0, 15, by = 0.5)
  vals <- obv_pep_efficacy_from_dpc(dpc)
  expect_true(all(diff(vals) <= 1e-10))
})

## --- apply_obv_pep_gate no-op paths ------------------------------------

test_that("gate is a no-op when obv_pep_enabled = FALSE", {
  pre_thinning <- list(
    infection_location      = rep("hospital", 5),
    offspring_class         = rep("HCW", 5),
    infection_time_absolute = rep(10, 5)
  )
  set.seed(1)
  res <- apply_obv_pep_gate(
    pre_thinning = pre_thinning,
    kept_indices = 1:5,
    obv_pep_enabled = FALSE,
    obv_pep_coverage = 1,
    obv_pep_adherence = 1,
    obv_pep_efficacy = 1
  )
  expect_true(all(res$keep))
  expect_equal(res$num_treated$pre_eligible, 0L)
  expect_equal(res$num_treated$prevented, 0L)
  expect_true(all(!res$metadata$obv_pep_eligible))
  expect_true(all(is.na(res$metadata$obv_pep_dpc)))
})

test_that("gate is a no-op when obv_pep_coverage = 0", {
  pre_thinning <- list(
    infection_location      = rep("hospital", 5),
    offspring_class         = rep("HCW", 5),
    infection_time_absolute = rep(10, 5)
  )
  set.seed(1)
  res <- apply_obv_pep_gate(
    pre_thinning = pre_thinning,
    kept_indices = 1:5,
    obv_pep_enabled = TRUE,
    obv_pep_coverage = 0,
    obv_pep_adherence = 1,
    obv_pep_efficacy = 1
  )
  expect_true(all(res$keep))
  expect_equal(res$num_treated$pre_eligible, 5L)   ## eligibility is still recorded
  expect_equal(res$num_treated$pre_treated, 0L)
  expect_equal(res$num_treated$post_treated, 0L)
  expect_equal(res$num_treated$prevented, 0L)
})

## --- full-prevention regression ----------------------------------------

test_that("gate prevents all HCW-hospital candidates when cov = adh = eff = 1", {
  pre_thinning <- list(
    infection_location      = c("hospital", "hospital", "hospital", "community", "community"),
    offspring_class         = c("HCW", "HCW", "HCW", "HCW", "genPop"),
    infection_time_absolute = rep(10, 5)
  )
  set.seed(1)
  res <- apply_obv_pep_gate(
    pre_thinning = pre_thinning,
    kept_indices = 1:5,
    obv_pep_enabled = TRUE,
    obv_pep_coverage = 1,
    obv_pep_adherence = 1,
    obv_pep_dpc = 0,
    obv_pep_efficacy = 1,
    obv_pep_target_class = "HCW",
    obv_pep_target_locations = "hospital"
  )
  ## All three HCW-hospital rows must be flipped to FALSE (prevented).
  expect_equal(res$keep, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(res$num_treated$pre_eligible, 3L)
  expect_equal(res$num_treated$pre_treated,  3L)
  expect_equal(res$num_treated$pre_adherent, 3L)
  expect_equal(res$num_treated$post_eligible, 3L)
  expect_equal(res$num_treated$post_treated,  3L)
  expect_equal(res$num_treated$post_adherent, 3L)
  expect_equal(res$num_treated$prevented, 3L)
})

test_that("gate respects obv_pep_target_class (genPop community example)", {
  pre_thinning <- list(
    infection_location      = c("community", "community", "hospital", "hospital"),
    offspring_class         = c("genPop", "HCW", "genPop", "HCW"),
    infection_time_absolute = rep(10, 4)
  )
  set.seed(1)
  res <- apply_obv_pep_gate(
    pre_thinning = pre_thinning,
    kept_indices = 1:4,
    obv_pep_enabled = TRUE,
    obv_pep_coverage = 1,
    obv_pep_adherence = 1,
    obv_pep_dpc = 0,
    obv_pep_efficacy = 1,
    obv_pep_target_class = "genPop",
    obv_pep_target_locations = "community"
  )
  ## Only the genPop-community row (index 1) is eligible.
  expect_equal(res$num_treated$pre_eligible, 1L)
  expect_equal(res$keep, c(FALSE, TRUE, TRUE, TRUE))
})

## --- pre vs post nesting invariants ------------------------------------

test_that("pre/post counts are nested for any random draw", {
  set.seed(123)
  n_pre <- 200
  pre_thinning <- list(
    infection_location      = rep("hospital", n_pre),
    offspring_class         = sample(c("HCW", "genPop"), n_pre, TRUE),
    infection_time_absolute = rep(20, n_pre)
  )
  kept <- sort(sample.int(n_pre, n_pre %/% 2))
  res <- apply_obv_pep_gate(
    pre_thinning = pre_thinning,
    kept_indices = kept,
    obv_pep_enabled = TRUE,
    obv_pep_coverage = 0.6,
    obv_pep_adherence = 0.7,
    obv_pep_dpc = 2,
    obv_pep_efficacy = 0.5
  )
  nt <- res$num_treated
  ## Pre-thinning nesting
  expect_true(nt$pre_treated  <= nt$pre_eligible)
  expect_true(nt$pre_adherent <= nt$pre_treated)
  ## Post-thinning is a subset of pre-thinning
  expect_true(nt$post_eligible <= nt$pre_eligible)
  expect_true(nt$post_treated  <= nt$pre_treated)
  expect_true(nt$post_adherent <= nt$pre_adherent)
  ## Post-thinning is also nested with itself
  expect_true(nt$post_treated  <= nt$post_eligible)
  expect_true(nt$post_adherent <= nt$post_treated)
  ## prevented is bounded above by post_adherent (only adherent kept can be prevented)
  expect_true(nt$prevented <= nt$post_adherent)
  ## metadata length matches kept_indices length
  expect_equal(nrow(res$metadata), length(kept))
  expect_equal(length(res$keep), length(kept))
})

## --- integration through offspring_function_genPop ---------------------

test_that("offspring_function_genPop returns obv_pep_num_treated attribute (disabled = zero)", {
  parent <- make_parent_info_genPop_hospitalised()
  set.seed(7)
  out <- do.call(offspring_function_genPop,
                 c(list(parent_info = parent), obv_off_args()))
  nt <- attr(out, "obv_pep_num_treated", exact = TRUE)
  expect_false(is.null(nt))
  expect_true(all(unlist(nt) == 0))
})

test_that("offspring_function_genPop with full prevention zeroes out HCW hospital offspring", {
  parent <- make_parent_info_genPop_hospitalised()
  set.seed(7)
  args <- obv_off_args()
  args$obv_pep_enabled        <- TRUE
  args$obv_pep_coverage   <- 1
  args$obv_pep_adherence      <- 1
  args$obv_pep_dpc            <- 0
  args$obv_pep_efficacy       <- 1
  out <- do.call(offspring_function_genPop, c(list(parent_info = parent), args))
  ## No HCWs should remain at hospital in the realised offspring set.
  expect_equal(sum(out$class == "HCW" & out$infection_location == "hospital"), 0)
  ## genPop offspring at community / hospital are untouched.
  nt <- attr(out, "obv_pep_num_treated", exact = TRUE)
  expect_true(nt$pre_eligible >= nt$prevented)
  expect_true(nt$prevented == nt$post_adherent)
})

## --- summarise_output OBV field contract -------------------------------
## Regression guard: pin the set of OBV-related field names surfaced by
## summarise_output(). If anyone renames or drops these in the future,
## downstream analysis scripts will break -- this test catches it.

test_that("summarise_output() exposes the expected OBV PEP field names", {
  ## Build a tiny tdf that has just enough structure for summarise_output to run.
  ## We don't care about the simulation contents -- we're just checking output
  ## field names.
  tdf <- data.frame(
    id                            = 1L,
    class                         = "genPop",
    infection_location            = "community",
    parent                        = NA_integer_,
    generation                    = 1L,
    time_infection_relative       = 0,
    time_infection_absolute       = 0,
    incubation_period             = 5,
    symptomatic                   = TRUE,
    time_symptom_onset_relative   = 5,
    time_symptom_onset_absolute   = 5,
    hospitalisation               = FALSE,
    time_hospitalisation_relative = NA_real_,
    time_hospitalisation_absolute = NA_real_,
    outcome                       = TRUE,
    outcome_location              = "community",
    time_outcome_relative         = 15,
    time_outcome_absolute         = 15,
    funeral_safety                = "safe",
    obv_pep_eligible              = FALSE,
    obv_pep_received              = FALSE,
    obv_pep_adherent              = FALSE,
    obv_pep_dpc                   = NA_real_,
    n_offspring                   = 0L,
    offspring_generated           = TRUE,
    stringsAsFactors              = FALSE
  )
  attr(tdf, "obv_pep_num_treated") <- empty_obv_pep_num_treated()
  sim_info <- list(population = 100, hcw_per_capita = 0.05, hcw_total = 5,
                   obv_pep_enabled = FALSE,
                   obv_pep_num_treated = empty_obv_pep_num_treated())

  out <- summarise_output(tdf, sim_info = sim_info)

  expected_obv_fields <- c(
    ## seven gate counters
    "n_obv_pep_pre_eligible",
    "n_obv_pep_pre_treated",
    "n_obv_pep_pre_adherent",
    "n_obv_pep_post_eligible",
    "n_obv_pep_post_treated",
    "n_obv_pep_post_adherent",
    "n_obv_pep_prevented",
    ## three tdf-based cohort counters
    "n_obv_pep_eligible_cases",
    "n_obv_pep_treated_cases",
    "n_obv_pep_breakthroughs",
    ## derived ratio
    "prop_obv_pep_prevented_among_adherent"
  )
  missing_fields <- setdiff(expected_obv_fields, names(out))
  expect_equal(missing_fields, character(0),
               info = sprintf("summarise_output() missing fields: %s",
                              paste(missing_fields, collapse = ", ")))
})
