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
