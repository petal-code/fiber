## Tests for the contact-first transmission model: the risk-tier structure, the
## R0 approximation and its inversion, the contact log, and the contact-tracing /
## isolation pathway that lets the risk tiers drive the NPIs.

## --- helpers -----------------------------------------------------------

bpm_args <- function(...) {
  args <- list(
    mn_contacts_genPop        = 10,
    overdisp_contacts_genPop  = 0.5,
    baseline_risk_genPop      = 0.1,
    Tg_shape_genPop           = 4, Tg_rate_genPop = 0.5,
    mn_contacts_hcw           = 10,
    overdisp_contacts_hcw     = 0.5,
    Tg_shape_hcw              = 4, Tg_rate_hcw = 0.5,
    mn_contacts_funeral       = 15,
    overdisp_contacts_funeral = 1,
    baseline_risk_funeral     = 0.1,
    Tg_shape_funeral          = 20, Tg_rate_funeral = 10,
    incubation_period           = function(n) rgamma(n, shape = 4, rate = 0.5),
    onset_to_hospitalisation    = function(n) rgamma(n, shape = 2, rate = 0.5),
    onset_to_death              = function(n) rgamma(n, shape = 4, rate = 0.5),
    onset_to_recovery           = function(n) rgamma(n, shape = 6, rate = 0.5),
    hospitalisation_to_death    = function(n) rgamma(n, shape = 3, rate = 0.5),
    hospitalisation_to_recovery = function(n) rgamma(n, shape = 5, rate = 0.5),
    prob_symptomatic         = 0.9,
    prob_hospitalised_hcw    = 0.5,
    prob_hospitalised_genPop = 0.4,
    prob_death_comm = 0.7, prob_death_hosp = 0.5,
    prob_hcw_cond_genPop_comm = 0.02, prob_hcw_cond_genPop_hospital = 0.3,
    prob_hcw_cond_hcw_comm = 0.05,    prob_hcw_cond_hcw_hospital = 0.3,
    prob_hospital_cond_hcw_preAdm = 0.5,
    ppe_coverage_hcw = 0.5, ppe_efficacy = 0.7,
    prop_etu = 0.5, etu_efficacy = 0.9,
    general_hospital_quarantine_efficacy = 0.3,
    p_unsafe_funeral_comm_hcw = 0.5, p_unsafe_funeral_hosp_hcw = 0.2,
    p_unsafe_funeral_comm_genPop = 0.6, p_unsafe_funeral_hosp_genPop = 0.2,
    safe_funeral_efficacy = 0.8,
    prob_hcw_cond_funeral_hcw = 0.1, prob_hcw_cond_funeral_genPop = 0.05,
    population = 1e6, hcw_per_capita = 0.01,
    check_final_size = 300, seeding_cases = 5, seed = 1
  )
  utils::modifyList(args, list(...))
}

parent_genPop <- function(time_infection_absolute = 0,
                          time_to_hospitalisation = NA_real_,
                          time_to_outcome = 30,
                          isolated = FALSE,
                          time_isolation_relative = NA_real_) {
  data.frame(
    id                            = 1L,
    class                         = "genPop",
    infection_location            = "community",
    parent                        = NA_integer_,
    generation                    = 1L,
    time_infection_relative       = 0,
    time_infection_absolute       = time_infection_absolute,
    incubation_period             = 5,
    symptomatic                   = TRUE,
    time_symptom_onset_relative   = 5,
    time_symptom_onset_absolute   = time_infection_absolute + 5,
    hospitalisation               = !is.na(time_to_hospitalisation),
    time_hospitalisation_relative = time_to_hospitalisation,
    time_hospitalisation_absolute = time_infection_absolute + time_to_hospitalisation,
    outcome                       = TRUE,
    outcome_location              = "community",
    time_outcome_relative         = time_to_outcome,
    time_outcome_absolute         = time_infection_absolute + time_to_outcome,
    funeral_safety                = "unsafe",
    isolated                      = isolated,
    time_isolation_relative       = time_isolation_relative,
    n_offspring                   = NA_integer_,
    offspring_generated           = FALSE,
    stringsAsFactors              = FALSE
  )
}

## --- risk structure ----------------------------------------------------

test_that("make_contact_risk validates and derives its summary quantities", {
  r <- make_contact_risk()
  expect_s3_class(r, "fiber_contact_risk")
  expect_equal(r$n_levels, 5)
  expect_equal(r$mean_relative_risk, 1)
  expect_equal(r$max_relative_risk, 1)
  expect_equal(sum(r$case_weights), 1)
  expect_equal(r$trace_prob, rep(0, 5))

  ## Fractions must sum to 1, lengths must match, risks must be positive.
  expect_error(make_contact_risk(fractions = c(0.5, 0.4)), "sum to 1")
  expect_error(make_contact_risk(fractions = c(0.5, 0.5), relative_risk = c(1, 2, 3)),
               "same length")
  expect_error(make_contact_risk(fractions = c(0.5, 0.5), relative_risk = c(1, -2)),
               "strictly positive")
  expect_error(make_contact_risk(fractions = c(0.5, 0.5), relative_risk = c(1, 2),
                                 trace_prob = c(0.5, 1.5)), "\\[0, 1\\]")
  expect_error(make_contact_risk(reference = 9), "between 1 and 5")
})

test_that("relative risks are normalised so the reference tier is the baseline", {
  ## Writing the reference tier's risk as 4 rather than 1 must not change anything:
  ## the baseline risk parameter always describes the reference tier.
  a <- make_contact_risk(fractions = rep(0.25, 4), relative_risk = c(1, 2, 4, 8))
  b <- make_contact_risk(fractions = rep(0.25, 4), relative_risk = c(4, 8, 16, 32))
  expect_equal(a$relative_risk, b$relative_risk)
  expect_equal(a$mean_relative_risk, b$mean_relative_risk)

  ## Naming a different reference rescales relative to that tier.
  c3 <- make_contact_risk(fractions = rep(0.25, 4), relative_risk = c(1, 2, 4, 8),
                          reference = 3)
  expect_equal(c3$relative_risk[3], 1)
  expect_equal(c3$relative_risk, a$relative_risk / a$relative_risk[3])
})

test_that("case weights over-represent high-risk tiers", {
  r <- make_contact_risk(fractions = rep(0.25, 4), relative_risk = c(1, 2, 4, 8))
  ## Cases are drawn from contacts in proportion to fractions * relative risk.
  expect_equal(r$case_weights,
               r$fractions * r$relative_risk / r$mean_relative_risk)
  expect_true(all(diff(r$case_weights) > 0))
  ## The lowest tier is under-represented among cases relative to contacts.
  expect_lt(r$case_weights[1], r$fractions[1])
  expect_gt(r$case_weights[4], r$fractions[4])
})

test_that("contact_risk_gradient builds a log-spaced gradient", {
  g <- contact_risk_gradient(n_levels = 5, ratio = 16, trace_prob_range = c(0.1, 0.9))
  expect_equal(g$relative_risk[1], 1)
  expect_equal(g$relative_risk[5], 16)
  ## Log-spaced: successive ratios are constant.
  expect_equal(diff(log(g$relative_risk)), rep(log(16) / 4, 4))
  expect_equal(g$trace_prob, seq(0.1, 0.9, length.out = 5))
  ## ratio = 1 degenerates to a flat structure.
  expect_equal(contact_risk_gradient(ratio = 1)$relative_risk, rep(1, 5))
})

## --- R0 approximation and its inversion --------------------------------

test_that("solve_baseline_risk_for_r0 round-trips through approx_r0", {
  args <- bpm_args(contact_risk = contact_risk_gradient(5, ratio = 6),
                   contact_risk_funeral = contact_risk_gradient(5, ratio = 3))
  args$baseline_risk_genPop <- NULL
  args$baseline_risk_funeral <- NULL

  set.seed(11)
  fit <- solve_baseline_risk_for_r0(R0 = 2.2, args = args,
                                    proportion_transmission_from_funerals = 0.25,
                                    n = 20000, seed = 3)
  args$baseline_risk_genPop  <- fit$baseline_risk_genPop_required
  args$baseline_risk_funeral <- fit$baseline_risk_funeral_required

  fwd <- approx_r0(args, invariants = fit$invariants)
  expect_equal(fwd$R0, 2.2)
  ## The funeral share must land where it was asked to.
  expect_equal(fwd$R0_funeral / fwd$R0, 0.25)
})

test_that("R0 is linear in the mean contact number and the mean relative risk", {
  args <- bpm_args()
  inv <- compute_r0_invariants(args, n = 20000, seed = 5)

  r1 <- approx_r0(args, invariants = inv)
  r2 <- approx_r0(bpm_args(mn_contacts_genPop = 20), invariants = inv)
  ## Doubling contacts doubles the direct contribution.
  expect_equal(r2$R0_direct, 2 * r1$R0_direct)

  ## Doubling every relative risk is absorbed by the reference-tier normalisation,
  ## so it must NOT change R0; changing the SPREAD does.
  flat <- make_contact_risk(fractions = rep(0.2, 5), relative_risk = rep(2, 5))
  r3 <- approx_r0(bpm_args(contact_risk = flat), invariants = inv)
  expect_equal(r3$R0_direct, r1$R0_direct)

  graded <- contact_risk_gradient(5, ratio = 10)
  r4 <- approx_r0(bpm_args(contact_risk = graded), invariants = inv)
  expect_equal(r4$R0_direct, r1$R0_direct * graded$mean_relative_risk)
})

test_that("the contact overdispersion does not enter R0", {
  args <- bpm_args()
  inv <- compute_r0_invariants(args, n = 10000, seed = 7)
  expect_equal(approx_r0(bpm_args(overdisp_contacts_genPop = 0.1), invariants = inv)$R0,
               approx_r0(bpm_args(overdisp_contacts_genPop = 50),  invariants = inv)$R0)
})

test_that("an infeasible R0 target errors and names the achievable ceiling", {
  ## 5 contacts cannot deliver R0 = 20, whatever the baseline risk.
  args <- bpm_args(mn_contacts_genPop = 5)
  args$baseline_risk_genPop <- NULL
  args$baseline_risk_funeral <- NULL
  expect_error(
    solve_baseline_risk_for_r0(R0 = 20, args = args,
                               proportion_transmission_from_funerals = 0,
                               n = 5000, seed = 2),
    "not achievable"
  )
})

test_that("isolation reduces the direct multiplier D", {
  args_no_trace <- bpm_args()
  inv_no <- compute_r0_invariants(args_no_trace, n = 20000, seed = 9)
  ## With no tracing there is no isolated generation-time mass at all.
  expect_equal(inv_no$Q_iso, 0)

  args_trace <- bpm_args(
    contact_risk = contact_risk_gradient(5, ratio = 4, trace_prob_range = 0.8),
    trace_coverage = 1, prob_isolate_given_traced = 1
  )
  inv_tr <- compute_r0_invariants(args_trace, n = 20000, seed = 9)
  expect_gt(inv_tr$Q_iso, 0)

  D_off <- r0_direct_multiplier(inv_tr, 0.9, 0.3, isolation_efficacy = 0)
  D_on  <- r0_direct_multiplier(inv_tr, 0.9, 0.3, isolation_efficacy = 0.9)
  expect_lt(D_on, D_off)
})

## --- offspring-function behaviour --------------------------------------

test_that("baseline_risk scales the number of infections, not the contacts", {
  ## With a flat structure, halving the baseline risk should roughly halve infections
  ## while leaving the contact count untouched.
  run <- function(p, seed) {
    set.seed(seed)
    o <- offspring_function_genPop(
      parent_info = parent_genPop(),
      mn_contacts_genPop = 400, overdisp_contacts_genPop = 200,
      baseline_risk_genPop = p, Tg_shape_genPop = 4, Tg_rate_genPop = 1,
      prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
      ppe_coverage_hcw = 0, ppe_efficacy = 0,
      prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
    )
    c(contacts = nrow(attr(o, "contact_log")), infections = nrow(o))
  }
  a <- run(0.8, 100); b <- run(0.4, 100)
  expect_equal(a[["contacts"]], b[["contacts"]])   # same seed, same contact draw
  expect_gt(a[["infections"]], b[["infections"]])
  expect_lt(abs(b[["infections"]] / a[["infections"]] - 0.5), 0.15)
})

test_that("baseline_risk = 1 with flat tiers makes every contact an infection", {
  set.seed(21)
  o <- offspring_function_genPop(
    parent_info = parent_genPop(),
    mn_contacts_genPop = 200, overdisp_contacts_genPop = 100,
    baseline_risk_genPop = 1, Tg_shape_genPop = 4, Tg_rate_genPop = 1,
    prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
    ppe_coverage_hcw = 0, ppe_efficacy = 0,
    prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
  )
  log <- attr(o, "contact_log")
  expect_equal(nrow(o), nrow(log))
  expect_true(all(log$record_type == "infection"))
})

test_that("high-risk tiers transmit more often than low-risk tiers", {
  risk <- make_contact_risk(fractions = rep(0.5, 2), relative_risk = c(1, 5))
  set.seed(31)
  o <- offspring_function_genPop(
    parent_info = parent_genPop(),
    mn_contacts_genPop = 4000, overdisp_contacts_genPop = 2000,
    baseline_risk_genPop = 0.1, contact_risk_genPop = risk,
    Tg_shape_genPop = 4, Tg_rate_genPop = 1,
    prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
    ppe_coverage_hcw = 0, ppe_efficacy = 0,
    prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
  )
  log <- attr(o, "contact_log")
  rate <- tapply(log$record_type == "infection", log$contact_risk_level, mean)
  expect_lt(abs(rate[["1"]] - 0.1), 0.03)
  expect_lt(abs(rate[["2"]] - 0.5), 0.05)
})

test_that("a baseline risk that pushes the top tier above 1 is rejected", {
  risk <- make_contact_risk(fractions = rep(0.5, 2), relative_risk = c(1, 5))
  expect_error(
    offspring_function_genPop(
      parent_info = parent_genPop(),
      mn_contacts_genPop = 50, overdisp_contacts_genPop = 25,
      baseline_risk_genPop = 0.5, contact_risk_genPop = risk,
      Tg_shape_genPop = 4, Tg_rate_genPop = 1,
      prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
      ppe_coverage_hcw = 0, ppe_efficacy = 0,
      prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
    ),
    "exceeds 1"
  )
})

test_that("isolation thins a genPop parent's post-isolation community transmission", {
  run_iso <- function(eff, seed = 41) {
    set.seed(seed)
    nrow(offspring_function_genPop(
      parent_info = parent_genPop(isolated = TRUE, time_isolation_relative = 0.001),
      mn_contacts_genPop = 2000, overdisp_contacts_genPop = 1000,
      baseline_risk_genPop = 1, Tg_shape_genPop = 4, Tg_rate_genPop = 1,
      isolation_efficacy = eff,
      prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
      ppe_coverage_hcw = 0, ppe_efficacy = 0,
      prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
    ))
  }
  ## Isolating from (almost) time zero with full efficacy blocks essentially everything.
  expect_equal(run_iso(1), 0)
  expect_gt(run_iso(0), 1000)
  ## Half efficacy should keep roughly half.
  expect_lt(abs(run_iso(0.5) / run_iso(0) - 0.5), 0.1)
})

test_that("a parent who is not isolated is unaffected by isolation_efficacy", {
  run <- function(eff) {
    set.seed(51)
    nrow(offspring_function_genPop(
      parent_info = parent_genPop(isolated = FALSE),
      mn_contacts_genPop = 500, overdisp_contacts_genPop = 250,
      baseline_risk_genPop = 1, Tg_shape_genPop = 4, Tg_rate_genPop = 1,
      isolation_efficacy = eff,
      prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
      ppe_coverage_hcw = 0, ppe_efficacy = 0,
      prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
    ))
  }
  expect_equal(run(0), run(1))
})

test_that("tracing probability follows the tier and scales with coverage", {
  risk <- make_contact_risk(fractions = rep(0.5, 2), relative_risk = c(1, 1),
                            trace_prob = c(0.2, 0.8))
  traced_rates <- function(cov, seed = 61) {
    set.seed(seed)
    o <- offspring_function_genPop(
      parent_info = parent_genPop(),
      mn_contacts_genPop = 4000, overdisp_contacts_genPop = 2000,
      baseline_risk_genPop = 0.1, contact_risk_genPop = risk,
      trace_coverage = cov, Tg_shape_genPop = 4, Tg_rate_genPop = 1,
      prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
      ppe_coverage_hcw = 0, ppe_efficacy = 0,
      prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
    )
    log <- attr(o, "contact_log")
    tapply(log$traced, log$contact_risk_level, mean)
  }
  full <- traced_rates(1)
  expect_lt(abs(full[["1"]] - 0.2), 0.03)
  expect_lt(abs(full[["2"]] - 0.8), 0.03)

  half <- traced_rates(0.5)
  expect_lt(abs(half[["1"]] - 0.1), 0.03)
  expect_lt(abs(half[["2"]] - 0.4), 0.03)

  ## Tracing is drawn for every contact, not only those that transmit.
  expect_equal(length(full), 2L)
})

## --- contact log integrity ---------------------------------------------

test_that("the contact log accounts for every contact exactly once", {
  out <- do.call(branching_process_main,
                 bpm_args(contact_risk = contact_risk_gradient(5, ratio = 5,
                                                               trace_prob_range = c(0.1, 0.9)),
                          trace_coverage = 0.7))
  log <- out$contact_log
  expect_gt(nrow(log), 0)

  ## Every row is either an infection (with a case id and no block reason) or a
  ## non-infection (with a reason and no case id). No third state.
  inf <- log$record_type == "infection"
  expect_true(all(log$record_type %in% c("contact", "infection")))
  expect_false(anyNA(log$case_id[inf]))
  expect_true(all(is.na(log$blocked_by[inf])))
  expect_true(all(is.na(log$case_id[!inf])))
  expect_false(anyNA(log$blocked_by[!inf]))
  expect_true(all(log$blocked_by[!inf] %in%
                    c("no_transmission", "ppe_quarantine_isolation",
                      "safe_funeral", "obv_pep")))

  ## Case ids are unique and every one resolves to a row in the tree.
  expect_false(anyDuplicated(log$case_id[inf]) > 0)
  expect_true(all(log$case_id[inf] %in% out$tdf$id))

  ## The tier, class and traced status recorded on the tree must agree with the log.
  m <- match(log$case_id[inf], out$tdf$id)
  expect_identical(out$tdf$contact_risk_level[m], log$contact_risk_level[inf])
  expect_identical(out$tdf$contact_risk_category[m], log$contact_risk_category[inf])
  expect_identical(out$tdf$traced[m], log$traced[inf])
  expect_identical(out$tdf$class[m], log$class[inf])

  ## Contacts always outnumber infections when the baseline risk is below 1.
  expect_gt(nrow(log), sum(inf))
})

test_that("seed cases carry no risk tier and are never traced", {
  out <- do.call(branching_process_main, bpm_args(trace_coverage = 1))
  seeds <- out$tdf[out$tdf$generation == 1 & !is.na(out$tdf$time_infection_absolute), ]
  expect_true(all(is.na(seeds$contact_risk_level)))
  expect_true(all(!seeds$traced))
})

## --- tracing drives the NPIs -------------------------------------------

test_that("tracing plus isolation reduces onward transmission", {
  ## Measured as realised offspring per expanded case rather than final size: these
  ## parameters are supercritical, so both arms run into `check_final_size` and the
  ## final size says more about the cap than about transmission.
  traced_risk <- contact_risk_gradient(5, ratio = 5, trace_prob_range = c(0.5, 0.95))
  mean_offspring <- function(extra, seeds = 1:8) {
    vapply(seeds, function(s) {
      o <- do.call(branching_process_main,
                   utils::modifyList(bpm_args(contact_risk = traced_risk, seed = s), extra))
      done <- o$tdf[!is.na(o$tdf$time_infection_absolute) & o$tdf$offspring_generated, ]
      mean(done$n_offspring)
    }, numeric(1))
  }
  off <- mean_offspring(list(trace_coverage = 0))
  on  <- mean_offspring(list(trace_coverage = 1, prob_isolate_given_traced = 0.9,
                             isolation_efficacy = 0.9))
  expect_lt(mean(on), mean(off))
})

test_that("isolation is recorded only for traced, symptomatic cases", {
  out <- do.call(branching_process_main,
                 bpm_args(contact_risk = contact_risk_gradient(5, ratio = 4,
                                                               trace_prob_range = 0.9),
                          trace_coverage = 1, prob_isolate_given_traced = 1,
                          isolation_efficacy = 0.8))
  real <- out$tdf[!is.na(out$tdf$time_infection_absolute), ]
  expect_gt(sum(real$isolated), 0)
  ## Isolation implies traced and symptomatic, and carries a time.
  expect_true(all(real$traced[real$isolated]))
  expect_true(all(real$symptomatic[real$isolated]))
  expect_false(anyNA(real$time_isolation_relative[real$isolated]))
  expect_true(all(is.na(real$time_isolation_relative[!real$isolated])))
  ## With prob_isolate_given_traced = 1, every traced symptomatic case isolates.
  expect_true(all(real$isolated[real$traced & real$symptomatic]))
})

test_that("tracing can accelerate admission for traced cases", {
  args <- bpm_args(contact_risk = contact_risk_gradient(5, ratio = 4, trace_prob_range = 0.9),
                   trace_coverage = 1, hospitalisation_delay_factor_traced = 0.2,
                   check_final_size = 600)
  out <- do.call(branching_process_main, args)
  real <- out$tdf[!is.na(out$tdf$time_infection_absolute) & out$tdf$hospitalisation, ]
  ## Admission delay is measured from symptom onset.
  delay <- real$time_hospitalisation_relative - real$incubation_period
  expect_lt(mean(delay[real$traced]), mean(delay[!real$traced]))
})

test_that("tracing can raise the probability of admission for traced cases", {
  args <- bpm_args(contact_risk = contact_risk_gradient(5, ratio = 4, trace_prob_range = 0.5),
                   trace_coverage = 1, prob_hospitalised_multiplier_traced = 2,
                   check_final_size = 800)
  out <- do.call(branching_process_main, args)
  real <- out$tdf[!is.na(out$tdf$time_infection_absolute) & out$tdf$symptomatic, ]
  expect_gt(mean(real$hospitalisation[real$traced]),
            mean(real$hospitalisation[!real$traced]))
})

## --- main-function wiring ----------------------------------------------

test_that("branching_process_main solves baseline risks from r0_target", {
  out <- do.call(branching_process_main,
                 utils::modifyList(
                   bpm_args(r0_target = 1.6, r0_prop_funeral = 0.3,
                            r0_solve_n = 10000, r0_solve_seed = 4),
                   list(baseline_risk_genPop = NULL, baseline_risk_funeral = NULL)))
  expect_equal(out$sim_info$r0_target, 1.6)
  expect_true(is.finite(out$sim_info$baseline_risk_genPop))
  expect_true(is.finite(out$sim_info$baseline_risk_funeral))
  ## HCW parents inherit the genPop per-contact risk when not given their own.
  expect_equal(out$sim_info$baseline_risk_hcw, out$sim_info$baseline_risk_genPop)
})

test_that("supplying both r0_target and a baseline risk is an error", {
  expect_error(
    do.call(branching_process_main, bpm_args(r0_target = 2)),
    "not both"
  )
})

test_that("per-route risk structures are independent and reported back", {
  rg <- contact_risk_gradient(5, ratio = 3)
  rf <- contact_risk_gradient(4, ratio = 9)
  out <- do.call(branching_process_main,
                 bpm_args(contact_risk = rg, contact_risk_funeral = rf))
  ## genPop and HCW inherit the shared structure; funerals use their own.
  expect_equal(out$sim_info$contact_risk_genPop$relative_risk, rg$relative_risk)
  expect_equal(out$sim_info$contact_risk_hcw$relative_risk, rg$relative_risk)
  expect_equal(out$sim_info$contact_risk_funeral$relative_risk, rf$relative_risk)
  expect_equal(out$sim_info$contact_risk_funeral$n_levels, 4)

  ## Funeral contacts must be labelled from the funeral structure.
  log <- out$contact_log
  fun_cats <- unique(log$contact_risk_category[log$infection_location == "funeral"])
  expect_true(all(fun_cats %in% rf$labels))
})

test_that("an infeasible baseline risk fails before the simulation starts", {
  expect_error(
    do.call(branching_process_main,
            bpm_args(contact_risk = contact_risk_gradient(5, ratio = 10),
                     baseline_risk_genPop = 0.5)),
    "highest-risk tier"
  )
})

test_that("summarise_output reports contact and tracing counts", {
  out <- do.call(branching_process_main,
                 bpm_args(contact_risk = contact_risk_gradient(5, ratio = 5,
                                                               trace_prob_range = c(0.2, 0.9)),
                          trace_coverage = 0.8, prob_isolate_given_traced = 0.5,
                          isolation_efficacy = 0.7))
  s <- summarise_output(out$tdf, sim_info = out$sim_info, contact_log = out$contact_log)

  expect_equal(s$n_contacts_total, nrow(out$contact_log))
  expect_equal(s$n_contacts_infected, sum(out$contact_log$record_type == "infection"))
  expect_gt(s$contacts_per_case, 1)
  expect_gt(s$n_cases_traced, 0)
  expect_true(s$n_cases_isolated <= s$n_cases_traced)

  ## The per-tier attack rate must increase with the tier's relative risk.
  ar <- s$attack_rate_by_risk_tier
  expect_equal(length(ar), 5L)
  expect_gt(ar[[5]], ar[[1]])

  ## Without the log the contact fields are NA but the call still works.
  s2 <- summarise_output(out$tdf, sim_info = out$sim_info)
  expect_true(is.na(s2$n_contacts_total))
  expect_equal(s2$n_cases_traced, s$n_cases_traced)
})
