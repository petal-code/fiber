## Tests for the ring-vaccination gate (apply_ring_vax_gate) in isolation.
## Probabilities are pinned to 0/1 wherever possible so rbinom() outcomes are
## deterministic and the timing / eligibility / two-ring logic can be asserted
## exactly; stochastic cases use set.seed and check nesting invariants.

## --- helper ------------------------------------------------------------
make_brood <- function(class, infection_location, time_infection_relative) {
  data.frame(class = class,
             infection_location = infection_location,
             time_infection_relative = time_infection_relative,
             stringsAsFactors = FALSE)
}

## A genPop/community brood at a spread of infection times after the parent.
gp_brood <- function(rel) make_brood(rep("genPop", length(rel)),
                                     rep("community", length(rel)), rel)

## --- detection / empty paths ------------------------------------------

test_that("an asymptomatic parent starts no ring-1 campaign", {
  res <- apply_ring_vax_gate(
    gp_brood(c(1, 5, 10, 20)),
    parent_time_infection_absolute = 0,
    parent_symptomatic = FALSE,
    parent_time_symptom_onset_absolute = NA_real_,
    ring_vax_detection_prob = 1, ring_vax_trace_prob = 1, ring_vax_coverage = 1,
    ring_vax_efficacy_infection = 1, ring_vax_logistical_delay = 0,
    ring_vax_protection_delay = 0, ring_vax_n_rings = 1
  )
  expect_true(is.na(res$parent_detection_time))
  expect_true(all(res$keep))
  expect_true(all(!res$vaccinated))
})

test_that("empty brood returns zero-length per-child vectors but resolves detection", {
  res <- apply_ring_vax_gate(
    gp_brood(numeric(0)),
    parent_time_infection_absolute = 0,
    parent_symptomatic = TRUE,
    parent_time_symptom_onset_absolute = 5,
    ring_vax_detection_prob = 1, ring_vax_reporting_delay = 0
  )
  expect_equal(length(res$keep), 0L)
  expect_equal(res$parent_detection_time, 5)
  expect_equal(res$num_treated$prevented, 0L)
})

## --- ring-1 timing ----------------------------------------------------

test_that("ring 1 averts only contacts protected in time (full coverage/efficacy)", {
  ## onset 5, reporting 0 -> detection 5; logistical 0 -> v1 = 5; protection 0 -> p1 = 5.
  ## Children infected at rel >= 5 are reached and protected in time -> averted.
  res <- apply_ring_vax_gate(
    gp_brood(c(1, 3, 5, 7, 9, 11)),
    parent_time_infection_absolute = 0,
    parent_symptomatic = TRUE, parent_time_symptom_onset_absolute = 5,
    ring_vax_detection_prob = 1, ring_vax_reporting_delay = 0,
    ring_vax_trace_prob = 1, ring_vax_coverage = 1, ring_vax_efficacy_infection = 1,
    ring_vax_logistical_delay = 0, ring_vax_protection_delay = 0,
    ring_vax_start = 0, ring_vax_n_rings = 1
  )
  expect_equal(res$keep, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(res$num_treated$prevented, 4L)
  expect_equal(res$parent_detection_time, 5)
  expect_equal(nrow(res$prevented_info), 4L)
  expect_true(all(res$prevented_info$time_infection_absolute >= 5))
})

test_that("a vaccinated contact infected before protection is a non-breakthrough full transmitter", {
  ## v1 = 5, protection delay 10 -> p1 = 15. rel in [5, 15) is vaccinated-too-late;
  ## rel >= 15 is protected (and averted at efficacy 1).
  res <- apply_ring_vax_gate(
    gp_brood(c(1, 5, 10, 15, 20)),
    parent_time_infection_absolute = 0,
    parent_symptomatic = TRUE, parent_time_symptom_onset_absolute = 5,
    ring_vax_detection_prob = 1, ring_vax_trace_prob = 1, ring_vax_coverage = 1,
    ring_vax_efficacy_infection = 1, ring_vax_logistical_delay = 0,
    ring_vax_protection_delay = 10, ring_vax_n_rings = 1
  )
  expect_equal(res$vaccinated, c(FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_equal(res$keep,       c(TRUE,  TRUE, TRUE, FALSE, FALSE))
  ## The two vaccinated-too-late contacts are not breakthroughs (never protected).
  expect_equal(res$breakthrough, c(FALSE, FALSE, FALSE, FALSE, FALSE))
})

test_that("efficacy = 0 yields breakthroughs, not aversions", {
  res <- apply_ring_vax_gate(
    gp_brood(c(1, 3, 5, 7, 9, 11)),
    parent_time_infection_absolute = 0,
    parent_symptomatic = TRUE, parent_time_symptom_onset_absolute = 5,
    ring_vax_detection_prob = 1, ring_vax_trace_prob = 1, ring_vax_coverage = 1,
    ring_vax_efficacy_infection = 0, ring_vax_logistical_delay = 0,
    ring_vax_protection_delay = 0, ring_vax_n_rings = 1
  )
  expect_true(all(res$keep))
  expect_equal(res$breakthrough, c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_equal(res$vaccinated,   c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_equal(res$num_treated$prevented, 0L)
})

test_that("the campaign-start gate suppresses vaccination before ring_vax_start", {
  ## v1 = 5 but start = 100, so no child is reachable.
  res <- apply_ring_vax_gate(
    gp_brood(c(10, 50, 200)),
    parent_time_infection_absolute = 0,
    parent_symptomatic = TRUE, parent_time_symptom_onset_absolute = 5,
    ring_vax_detection_prob = 1, ring_vax_trace_prob = 1, ring_vax_coverage = 1,
    ring_vax_efficacy_infection = 1, ring_vax_logistical_delay = 0,
    ring_vax_protection_delay = 0, ring_vax_start = 100, ring_vax_n_rings = 1
  )
  expect_true(all(res$keep))
  expect_true(all(!res$vaccinated))
})

## --- eligibility ------------------------------------------------------

test_that("only target class/location contacts are eligible", {
  brood <- make_brood(c("genPop", "HCW", "genPop"),
                      c("community", "community", "hospital"),
                      c(10, 10, 10))
  res <- apply_ring_vax_gate(
    brood,
    parent_time_infection_absolute = 0,
    parent_symptomatic = TRUE, parent_time_symptom_onset_absolute = 1,
    ring_vax_detection_prob = 1, ring_vax_trace_prob = 1, ring_vax_coverage = 1,
    ring_vax_efficacy_infection = 1, ring_vax_logistical_delay = 0,
    ring_vax_protection_delay = 0, ring_vax_n_rings = 1,
    ring_vax_target_class = "genPop", ring_vax_target_locations = "community"
  )
  ## Only row 1 (genPop/community) is eligible and averted.
  expect_equal(res$keep, c(FALSE, TRUE, TRUE))
})

## --- ring 2 (grandparent) + require-intermediate-traced ---------------

test_that("ring 2 averts grandchildren when the intermediate was traced", {
  ## Parent asymptomatic (no ring 1) isolates ring 2. Grandparent detected at 5.
  res <- apply_ring_vax_gate(
    gp_brood(c(1, 8, 20)),
    parent_time_infection_absolute = 0,
    parent_symptomatic = FALSE, parent_time_symptom_onset_absolute = NA_real_,
    grandparent_detection_time = 5, grandparent_infection_time = -10,
    parent_traced_by_grandparent = TRUE,
    ring_vax_n_rings = 2, ring_vax_require_intermediate_traced = TRUE,
    ring_vax_trace_prob = 1, ring_vax_coverage = 1, ring_vax_efficacy_infection = 1,
    ring_vax_logistical_delay = 0, ring_vax_protection_delay = 0,
    ring_vax_ring2_delay_increment = 0, ring_vax_start = 0
  )
  ## v2 = 5; rel >= 5 averted.
  expect_equal(res$keep, c(TRUE, FALSE, FALSE))
})

test_that("ring 2 is blocked when the intermediate was not traced", {
  res <- apply_ring_vax_gate(
    gp_brood(c(1, 8, 20)),
    parent_time_infection_absolute = 0,
    parent_symptomatic = FALSE, parent_time_symptom_onset_absolute = NA_real_,
    grandparent_detection_time = 5, grandparent_infection_time = -10,
    parent_traced_by_grandparent = FALSE,
    ring_vax_n_rings = 2, ring_vax_require_intermediate_traced = TRUE,
    ring_vax_trace_prob = 1, ring_vax_coverage = 1, ring_vax_efficacy_infection = 1,
    ring_vax_logistical_delay = 0, ring_vax_protection_delay = 0, ring_vax_start = 0
  )
  expect_true(all(res$keep))
})

test_that("ring 2 is inactive when n_rings = 1", {
  res <- apply_ring_vax_gate(
    gp_brood(c(8, 20)),
    parent_time_infection_absolute = 0,
    parent_symptomatic = FALSE, parent_time_symptom_onset_absolute = NA_real_,
    grandparent_detection_time = 5, grandparent_infection_time = -10,
    parent_traced_by_grandparent = TRUE,
    ring_vax_n_rings = 1,
    ring_vax_trace_prob = 1, ring_vax_coverage = 1, ring_vax_efficacy_infection = 1,
    ring_vax_logistical_delay = 0, ring_vax_protection_delay = 0, ring_vax_start = 0
  )
  expect_true(all(res$keep))
})

test_that("a child failed by ring 2 can still be averted by the parent's ring 1", {
  ## Grandparent detected late (so ring-2 protection p2 lands after the child's
  ## infection), parent detected early (ring-1 protection p1 lands before it):
  ## the parent's own ring rescues the contact.
  ##   grandparent detection 30 -> v2 = p2 = 30; parent onset 2 -> v1 = p1 = 2.
  res <- apply_ring_vax_gate(
    gp_brood(10),
    parent_time_infection_absolute = 0,
    parent_symptomatic = TRUE, parent_time_symptom_onset_absolute = 2,
    grandparent_detection_time = 30, grandparent_infection_time = -20,
    parent_traced_by_grandparent = TRUE,
    ring_vax_n_rings = 2, ring_vax_require_intermediate_traced = TRUE,
    ring_vax_trace_prob = 1, ring_vax_coverage = 1, ring_vax_efficacy_infection = 1,
    ring_vax_logistical_delay = 0, ring_vax_protection_delay = 0, ring_vax_start = 0
  )
  ## Ring 2 too late (v2 = 30 > 10) but ring 1 in time (p1 = 2 <= 10) -> averted.
  expect_equal(res$keep, FALSE)
})

## --- coverage modes ---------------------------------------------------

test_that("independent vs shared coverage agree at coverage = 1", {
  args <- list(
    gp_brood(c(8, 20)),
    parent_time_infection_absolute = 0,
    parent_symptomatic = TRUE, parent_time_symptom_onset_absolute = 2,
    ring_vax_trace_prob = 1, ring_vax_coverage = 1, ring_vax_efficacy_infection = 1,
    ring_vax_logistical_delay = 0, ring_vax_protection_delay = 0, ring_vax_n_rings = 1
  )
  res_ind <- do.call(apply_ring_vax_gate, c(args, ring_vax_independent_coverage = TRUE))
  res_shr <- do.call(apply_ring_vax_gate, c(args, ring_vax_independent_coverage = FALSE))
  expect_equal(res_ind$keep, res_shr$keep)
  expect_true(all(!res_ind$keep))
})

## --- counter nesting (stochastic) -------------------------------------

test_that("gate counters are nested for arbitrary intermediate probabilities", {
  set.seed(1)
  rel <- runif(300, 5, 40)
  brood <- make_brood(sample(c("genPop", "HCW"), 300, TRUE),
                      rep("community", 300), rel)
  res <- apply_ring_vax_gate(
    brood,
    parent_time_infection_absolute = 0,
    parent_symptomatic = TRUE, parent_time_symptom_onset_absolute = 2,
    ring_vax_detection_prob = 1, ring_vax_trace_prob = 0.7, ring_vax_coverage = 0.6,
    ring_vax_efficacy_infection = 0.5, ring_vax_logistical_delay = 1,
    ring_vax_protection_delay = 1, ring_vax_n_rings = 1
  )
  nt <- res$num_treated
  expect_true(nt$prevented  <= nt$protected)
  expect_true(nt$protected  <= nt$vaccinated)
  expect_true(nt$vaccinated <= nt$traced)
  ## prevented_info rows match the prevented counter.
  expect_equal(nrow(res$prevented_info), nt$prevented)
  ## keep / metadata vectors are all brood-length.
  expect_equal(length(res$keep), nrow(brood))
  expect_equal(length(res$vaccinated), nrow(brood))
})

test_that("a constant function reproduces the scalar gate exactly", {
  args <- list(
    gp_brood(c(3, 7, 12, 25)),
    parent_time_infection_absolute = 0,
    parent_symptomatic = TRUE, parent_time_symptom_onset_absolute = 2,
    ring_vax_efficacy_infection = 1, ring_vax_logistical_delay = 2,
    ring_vax_protection_delay = 3, ring_vax_n_rings = 1
  )
  set.seed(99)
  a <- do.call(apply_ring_vax_gate,
               c(args, ring_vax_detection_prob = 1, ring_vax_trace_prob = 0.5,
                 ring_vax_coverage = 0.8))
  set.seed(99)
  b <- do.call(apply_ring_vax_gate,
               c(args, ring_vax_detection_prob = function(t) 1,
                 ring_vax_trace_prob = function(t) 0.5,
                 ring_vax_coverage = function(t) 0.8))
  expect_identical(a$keep, b$keep)
  expect_identical(a$vaccinated, b$vaccinated)
})

## --- end-to-end through branching_process_main ------------------------

rv_bpm_args <- function(...) {
  fixed <- function(v) function(n) rep(v, n)
  defaults <- list(
    mn_offspring_genPop           = 2.0,
    overdisp_offspring_genPop     = 0.5,
    Tg_shape_genPop               = 2,
    Tg_rate_genPop                = 0.15,
    mn_offspring_hcw              = 2.0,
    overdisp_offspring_hcw        = 0.5,
    Tg_shape_hcw                  = 2,
    Tg_rate_hcw                   = 0.15,
    mn_offspring_funeral          = 1.0,
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
    check_final_size              = 300,
    seeding_cases                 = 5,
    seed                          = 1L
  )
  utils::modifyList(defaults, list(...))
}

test_that("disabled ring vaccination is byte-identical to a default run (existing columns)", {
  base <- do.call(branching_process_main, rv_bpm_args(seed = 3L))
  off  <- do.call(branching_process_main, rv_bpm_args(seed = 3L, ring_vax_enabled = FALSE))
  shared <- setdiff(names(base$tdf), c("ring_vax_parent_detection_time",
                                       "ring_vax_parent_infection_time",
                                       "ring_vax_traced_by_parent", "ring_vaccinated",
                                       "ring_vax_breakthrough", "time_ring_vaccinated_absolute"))
  expect_equal(base$tdf[shared], off$tdf[shared])
  ## And the ring columns are all at their inert defaults when disabled.
  expect_true(all(!off$tdf$ring_vaccinated))
  expect_true(all(!off$tdf$ring_vax_breakthrough))
  expect_true(all(is.na(off$tdf$time_ring_vaccinated_absolute)))
  expect_null(off$ring_prevented_completed)
  s <- summarise_output(off$tdf, sim_info = off$sim_info)
  expect_equal(s$n_ring_vax_prevented, 0L)
  expect_equal(s$n_ring_vax_prevented_deaths, 0L)
})

test_that("ring vaccination runs end-to-end and reduces realised cases under full prevention", {
  ## Full, fast prevention: every symptomatic case detected, contacts always traced,
  ## vaccinated and protected instantly, vaccine fully efficacious.
  full <- rv_bpm_args(
    seed = 3L, ring_vax_enabled = TRUE, ring_vax_n_rings = 2,
    ring_vax_detection_prob = 1, ring_vax_reporting_delay = 0,
    ring_vax_trace_prob = 1, ring_vax_coverage = 1, ring_vax_efficacy_infection = 1,
    ring_vax_efficacy_transmission = 1, ring_vax_logistical_delay = 0,
    ring_vax_protection_delay = 0
  )
  res <- do.call(branching_process_main, full)
  s_on <- summarise_output(res$tdf, sim_info = res$sim_info)
  off  <- do.call(branching_process_main, rv_bpm_args(seed = 3L))
  s_off <- summarise_output(off$tdf, sim_info = off$sim_info)
  ## Ring vax should curtail the outbreak well below the no-intervention size.
  expect_lt(s_on$n_cases_total, s_off$n_cases_total)
  expect_gt(s_on$n_ring_vax_prevented, 0)
})

test_that("ring-vax counters nest and prevented_deaths is bounded and reproducible", {
  args <- rv_bpm_args(
    seed = 7L, ring_vax_enabled = TRUE, ring_vax_detection_prob = 0.8,
    ring_vax_trace_prob = 0.8, ring_vax_coverage = 0.7, ring_vax_efficacy_infection = 0.8,
    ring_vax_efficacy_transmission = 0.5, ring_vax_logistical_delay = 1,
    ring_vax_protection_delay = 2
  )
  r1 <- do.call(branching_process_main, args)
  r2 <- do.call(branching_process_main, args)
  s  <- summarise_output(r1$tdf, sim_info = r1$sim_info)
  expect_true(s$n_ring_vax_prevented  <= s$n_ring_vax_protected)
  expect_true(s$n_ring_vax_protected  <= s$n_ring_vax_vaccinated)
  expect_true(s$n_ring_vax_vaccinated <= s$n_ring_vax_traced)
  expect_true(s$n_ring_vax_prevented_deaths <= s$n_ring_vax_prevented)
  ## Deferred counterfactual is post-loop, so the tree is reproducible under a seed.
  expect_equal(r1$tdf, r2$tdf)
  expect_equal(s$n_ring_vax_prevented_deaths,
               summarise_output(r2$tdf, sim_info = r2$sim_info)$n_ring_vax_prevented_deaths)
  ## prevented_completed has one row per prevented infection; outcomes sum to deaths.
  if (!is.null(r1$ring_prevented_completed)) {
    expect_equal(nrow(r1$ring_prevented_completed), s$n_ring_vax_prevented)
    expect_equal(sum(r1$ring_prevented_completed$outcome, na.rm = TRUE),
                 s$n_ring_vax_prevented_deaths)
  }
})

test_that("summarise_output exposes the expected ring-vax field names", {
  off <- do.call(branching_process_main, rv_bpm_args(seed = 1L))
  out <- summarise_output(off$tdf, sim_info = off$sim_info)
  expected <- c("n_ring_vax_traced", "n_ring_vax_vaccinated", "n_ring_vax_protected",
                "n_ring_vax_prevented", "n_ring_vax_prevented_deaths",
                "n_ring_vaccinated_cases", "n_ring_vax_breakthroughs")
  expect_equal(setdiff(expected, names(out)), character(0))
})
