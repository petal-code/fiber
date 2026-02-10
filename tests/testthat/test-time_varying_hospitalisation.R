## Helper: build a minimal parent_info row matching the tdf structure
make_parent_info <- function(time_infection_absolute = 0) {
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
    hospitalisation                = FALSE,
    time_hospitalisation_relative  = NA_real_,
    time_hospitalisation_absolute  = NA_real_,
    outcome                        = TRUE,
    outcome_location               = "community",
    time_outcome_relative          = 12,
    time_outcome_absolute          = time_infection_absolute + 12,
    funeral_safety                 = "unsafe",
    n_offspring                    = NA_integer_,
    offspring_generated            = FALSE,
    stringsAsFactors               = FALSE
  )
}

## Helper: build a minimal offspring_dataframe with n genPop offspring
make_offspring_df <- function(n, classes = rep("genPop", n)) {
  data.frame(
    id                             = rep(NA_integer_, n),
    class                          = classes,
    infection_location             = rep("community", n),
    parent                         = rep(NA_integer_, n),
    generation                     = rep(NA_integer_, n),
    time_infection_relative        = rep(1, n),
    time_infection_absolute        = rep(NA_real_, n),
    incubation_period              = rep(NA_real_, n),
    symptomatic                    = rep(NA, n),
    time_symptom_onset_relative    = rep(NA_real_, n),
    time_symptom_onset_absolute    = rep(NA_real_, n),
    hospitalisation                = rep(FALSE, n),
    time_hospitalisation_relative  = rep(NA_real_, n),
    time_hospitalisation_absolute  = rep(NA_real_, n),
    outcome                        = rep(FALSE, n),
    outcome_location               = rep(NA_character_, n),
    time_outcome_relative          = rep(NA_real_, n),
    time_outcome_absolute          = rep(NA_real_, n),
    funeral_safety                 = rep(NA_character_, n),
    n_offspring                    = rep(NA_integer_, n),
    offspring_generated            = rep(FALSE, n),
    stringsAsFactors               = FALSE
  )
}

## Common delay distribution stubs (deterministic for testing)
fixed_incubation    <- function(n) rep(5, n)
fixed_hosp_delay    <- function(n) rep(2, n)
fixed_death_delay   <- function(n) rep(7, n)
fixed_recover_delay <- function(n) rep(10, n)

## ------------------------------------------------------------------
## resolve_time_varying unit tests
## ------------------------------------------------------------------
test_that("resolve_time_varying works with scalar", {
  expect_equal(resolve_time_varying(0.5, 10, "test_param"), 0.5)
})

test_that("resolve_time_varying works with function", {
  f <- function(t) 0.1 + 0.01 * t
  expect_equal(resolve_time_varying(f, 10, "test_param"), 0.2)
})

test_that("resolve_time_varying errors on character input", {
  expect_error(
    resolve_time_varying("bad", 10, "test_param"),
    "must be a single numeric value or a function"
  )
})

test_that("resolve_time_varying errors when function returns NA", {
  f <- function(t) NA_real_
  expect_error(
    resolve_time_varying(f, 10, "test_param"),
    "evaluated to an invalid value"
  )
})

## ------------------------------------------------------------------
## Scalar backward compatibility
## ------------------------------------------------------------------
test_that("complete_offspring_info works with scalar prob_hospitalised_genPop", {
  set.seed(42)
  parent <- make_parent_info(time_infection_absolute = 0)
  offspring <- make_offspring_df(20)

  result <- complete_offspring_info(
    parent_info                = parent,
    offspring_dataframe        = offspring,
    prob_symptomatic           = 0.8,
    prob_hospitalised_hcw      = 0.5,
    prob_hospitalised_genPop   = 0.4,
    prob_death_comm            = 0.5,
    prob_death_hosp            = 0.3,
    p_unsafe_funeral_comm_hcw    = 0.5,
    p_unsafe_funeral_hosp_hcw    = 0.2,
    p_unsafe_funeral_comm_genPop = 0.5,
    p_unsafe_funeral_hosp_genPop = 0.2,
    incubation_period          = fixed_incubation,
    onset_to_hospitalisation   = fixed_hosp_delay,
    hospitalisation_to_death   = fixed_death_delay,
    hospitalisation_to_recovery = fixed_recover_delay,
    onset_to_death             = fixed_death_delay,
    onset_to_recovery          = fixed_recover_delay
  )

  expect_equal(nrow(result), 20)
  expect_true(all(c("hospitalisation", "outcome", "symptomatic") %in% names(result)))
})

## ------------------------------------------------------------------
## Constant function yields same results as equivalent scalar
## ------------------------------------------------------------------
test_that("constant function(t) gives same results as scalar", {
  set.seed(123)
  parent <- make_parent_info(time_infection_absolute = 0)
  offspring1 <- make_offspring_df(30)

  result_scalar <- complete_offspring_info(
    parent_info                = parent,
    offspring_dataframe        = offspring1,
    prob_symptomatic           = 0.8,
    prob_hospitalised_hcw      = 0.5,
    prob_hospitalised_genPop   = 0.4,
    prob_death_comm            = 0.5,
    prob_death_hosp            = 0.3,
    p_unsafe_funeral_comm_hcw    = 0.5,
    p_unsafe_funeral_hosp_hcw    = 0.2,
    p_unsafe_funeral_comm_genPop = 0.5,
    p_unsafe_funeral_hosp_genPop = 0.2,
    incubation_period          = fixed_incubation,
    onset_to_hospitalisation   = fixed_hosp_delay,
    hospitalisation_to_death   = fixed_death_delay,
    hospitalisation_to_recovery = fixed_recover_delay,
    onset_to_death             = fixed_death_delay,
    onset_to_recovery          = fixed_recover_delay
  )

  set.seed(123)
  offspring2 <- make_offspring_df(30)

  result_fn <- complete_offspring_info(
    parent_info                = parent,
    offspring_dataframe        = offspring2,
    prob_symptomatic           = 0.8,
    prob_hospitalised_hcw      = 0.5,
    prob_hospitalised_genPop   = function(t) 0.4,
    prob_death_comm            = 0.5,
    prob_death_hosp            = 0.3,
    p_unsafe_funeral_comm_hcw    = 0.5,
    p_unsafe_funeral_hosp_hcw    = 0.2,
    p_unsafe_funeral_comm_genPop = 0.5,
    p_unsafe_funeral_hosp_genPop = 0.2,
    incubation_period          = fixed_incubation,
    onset_to_hospitalisation   = fixed_hosp_delay,
    hospitalisation_to_death   = fixed_death_delay,
    hospitalisation_to_recovery = fixed_recover_delay,
    onset_to_death             = fixed_death_delay,
    onset_to_recovery          = fixed_recover_delay
  )

  expect_identical(result_scalar$hospitalisation, result_fn$hospitalisation)
  expect_identical(result_scalar$symptomatic, result_fn$symptomatic)
  expect_identical(result_scalar$outcome, result_fn$outcome)
})

## ------------------------------------------------------------------
## Time-variation has effect on hospitalisation probability
## ------------------------------------------------------------------
test_that("time-varying prob_hospitalised_genPop changes outcomes", {
  ## All offspring are symptomatic; early ones get 0% hosp, late ones get high hosp
  parent_early <- make_parent_info(time_infection_absolute = 0)
  parent_late  <- make_parent_info(time_infection_absolute = 100)

  prob_hosp_fn <- function(t) ifelse(t < 10, 0.0, 0.8)

  set.seed(999)
  offspring_early <- make_offspring_df(50)
  result_early <- complete_offspring_info(
    parent_info                = parent_early,
    offspring_dataframe        = offspring_early,
    prob_symptomatic           = 1.0,
    prob_hospitalised_hcw      = 0.5,
    prob_hospitalised_genPop   = prob_hosp_fn,
    prob_death_comm            = 0.5,
    prob_death_hosp            = 0.3,
    p_unsafe_funeral_comm_hcw    = 0.5,
    p_unsafe_funeral_hosp_hcw    = 0.2,
    p_unsafe_funeral_comm_genPop = 0.5,
    p_unsafe_funeral_hosp_genPop = 0.2,
    incubation_period          = fixed_incubation,
    onset_to_hospitalisation   = fixed_hosp_delay,
    hospitalisation_to_death   = fixed_death_delay,
    hospitalisation_to_recovery = fixed_recover_delay,
    onset_to_death             = fixed_death_delay,
    onset_to_recovery          = fixed_recover_delay
  )

  set.seed(999)
  offspring_late <- make_offspring_df(50)
  result_late <- complete_offspring_info(
    parent_info                = parent_late,
    offspring_dataframe        = offspring_late,
    prob_symptomatic           = 1.0,
    prob_hospitalised_hcw      = 0.5,
    prob_hospitalised_genPop   = prob_hosp_fn,
    prob_death_comm            = 0.5,
    prob_death_hosp            = 0.3,
    p_unsafe_funeral_comm_hcw    = 0.5,
    p_unsafe_funeral_hosp_hcw    = 0.2,
    p_unsafe_funeral_comm_genPop = 0.5,
    p_unsafe_funeral_hosp_genPop = 0.2,
    incubation_period          = fixed_incubation,
    onset_to_hospitalisation   = fixed_hosp_delay,
    hospitalisation_to_death   = fixed_death_delay,
    hospitalisation_to_recovery = fixed_recover_delay,
    onset_to_death             = fixed_death_delay,
    onset_to_recovery          = fixed_recover_delay
  )

  ## Early parent (t=0): symptom onset at t=5 → prob_hosp = 0 → no hospitalisations
  expect_equal(sum(result_early$hospitalisation), 0)
  ## Late parent (t=100): symptom onset at t=105 → prob_hosp = 0.8 → some hospitalisations
  expect_true(sum(result_late$hospitalisation) > 0)
})

## ------------------------------------------------------------------
## Delay factor (scalar) works
## ------------------------------------------------------------------
test_that("hospitalisation_delay_factor scales hospitalisation times", {
  parent <- make_parent_info(time_infection_absolute = 0)

  set.seed(555)
  offspring1 <- make_offspring_df(50)
  result_base <- complete_offspring_info(
    parent_info                = parent,
    offspring_dataframe        = offspring1,
    prob_symptomatic           = 1.0,
    prob_hospitalised_hcw      = 0.5,
    prob_hospitalised_genPop   = 0.8,
    prob_death_comm            = 0.5,
    prob_death_hosp            = 0.3,
    p_unsafe_funeral_comm_hcw    = 0.5,
    p_unsafe_funeral_hosp_hcw    = 0.2,
    p_unsafe_funeral_comm_genPop = 0.5,
    p_unsafe_funeral_hosp_genPop = 0.2,
    incubation_period          = fixed_incubation,
    onset_to_hospitalisation   = fixed_hosp_delay,
    hospitalisation_delay_factor = 1.0,
    hospitalisation_to_death   = fixed_death_delay,
    hospitalisation_to_recovery = fixed_recover_delay,
    onset_to_death             = fixed_death_delay,
    onset_to_recovery          = fixed_recover_delay
  )

  set.seed(555)
  offspring2 <- make_offspring_df(50)
  result_scaled <- complete_offspring_info(
    parent_info                = parent,
    offspring_dataframe        = offspring2,
    prob_symptomatic           = 1.0,
    prob_hospitalised_hcw      = 0.5,
    prob_hospitalised_genPop   = 0.8,
    prob_death_comm            = 0.5,
    prob_death_hosp            = 0.3,
    p_unsafe_funeral_comm_hcw    = 0.5,
    p_unsafe_funeral_hosp_hcw    = 0.2,
    p_unsafe_funeral_comm_genPop = 0.5,
    p_unsafe_funeral_hosp_genPop = 0.2,
    incubation_period          = fixed_incubation,
    onset_to_hospitalisation   = fixed_hosp_delay,
    hospitalisation_delay_factor = 2.0,
    hospitalisation_to_death   = fixed_death_delay,
    hospitalisation_to_recovery = fixed_recover_delay,
    onset_to_death             = fixed_death_delay,
    onset_to_recovery          = fixed_recover_delay
  )

  ## With deterministic delays of 2, factor=2 gives delay of 4.
  ## Hospitalisation times (relative) should differ accordingly.
  ## Some offspring that are hospitalised with factor=1 may not be hospitalised with factor=2
  ## because the longer delay means community outcome might precede hospitalisation.
  ## For those hospitalised in both, the relative hosp time should be larger in the scaled case.
  hosp_base <- which(result_base$hospitalisation)
  hosp_scaled <- which(result_scaled$hospitalisation)

  ## The same random draws are used, so hospitalisation decisions are identical
  ## but the timing differs — potentially hospitalised offspring may lose their
  ## hospitalisation if the scaled delay exceeds community outcome time.
  ## With fixed delays: onset_to_hosp = 2 (base) vs 4 (scaled),
  ## onset_to_death/recovery = 7/10. Both 2 and 4 < 7, so hospitalisation
  ## should be the same set.
  expect_equal(hosp_base, hosp_scaled)

  ## For hospitalised offspring, check that the onset-to-hospitalisation
  ## component is 2× larger in the scaled run.
  ## time_hospitalisation_relative = incubation_period + hosp_delay
  ## incubation = 5 in both cases, so hosp_delay = time_hosp_rel - 5
  if (length(hosp_base) > 0) {
    delay_base   <- result_base$time_hospitalisation_relative[hosp_base] - 5
    delay_scaled <- result_scaled$time_hospitalisation_relative[hosp_scaled] - 5
    expect_equal(delay_scaled, delay_base * 2)
  }
})

## ------------------------------------------------------------------
## Delay factor as function of time
## ------------------------------------------------------------------
test_that("hospitalisation_delay_factor as function varies with time", {
  ## Early parent: factor = 3; late parent: factor = 1
  delay_factor_fn <- function(t) ifelse(t < 10, 3.0, 1.0)

  parent_early <- make_parent_info(time_infection_absolute = 0)
  parent_late  <- make_parent_info(time_infection_absolute = 100)

  set.seed(777)
  offspring1 <- make_offspring_df(30)
  result_early <- complete_offspring_info(
    parent_info                = parent_early,
    offspring_dataframe        = offspring1,
    prob_symptomatic           = 1.0,
    prob_hospitalised_hcw      = 0.5,
    prob_hospitalised_genPop   = 0.8,
    prob_death_comm            = 0.5,
    prob_death_hosp            = 0.3,
    p_unsafe_funeral_comm_hcw    = 0.5,
    p_unsafe_funeral_hosp_hcw    = 0.2,
    p_unsafe_funeral_comm_genPop = 0.5,
    p_unsafe_funeral_hosp_genPop = 0.2,
    incubation_period          = fixed_incubation,
    onset_to_hospitalisation   = fixed_hosp_delay,
    hospitalisation_delay_factor = delay_factor_fn,
    hospitalisation_to_death   = fixed_death_delay,
    hospitalisation_to_recovery = fixed_recover_delay,
    onset_to_death             = fixed_death_delay,
    onset_to_recovery          = fixed_recover_delay
  )

  set.seed(777)
  offspring2 <- make_offspring_df(30)
  result_late <- complete_offspring_info(
    parent_info                = parent_late,
    offspring_dataframe        = offspring2,
    prob_symptomatic           = 1.0,
    prob_hospitalised_hcw      = 0.5,
    prob_hospitalised_genPop   = 0.8,
    prob_death_comm            = 0.5,
    prob_death_hosp            = 0.3,
    p_unsafe_funeral_comm_hcw    = 0.5,
    p_unsafe_funeral_hosp_hcw    = 0.2,
    p_unsafe_funeral_comm_genPop = 0.5,
    p_unsafe_funeral_hosp_genPop = 0.2,
    incubation_period          = fixed_incubation,
    onset_to_hospitalisation   = fixed_hosp_delay,
    hospitalisation_delay_factor = delay_factor_fn,
    hospitalisation_to_death   = fixed_death_delay,
    hospitalisation_to_recovery = fixed_recover_delay,
    onset_to_death             = fixed_death_delay,
    onset_to_recovery          = fixed_recover_delay
  )

  ## Early: symptom onset at t=5 → factor=3 → hosp delay = 6
  ## Late:  symptom onset at t=105 → factor=1 → hosp delay = 2
  hosp_early <- which(result_early$hospitalisation)
  hosp_late  <- which(result_late$hospitalisation)

  if (length(hosp_early) > 0) {
    delay_early <- result_early$time_hospitalisation_relative[hosp_early] - 5
    expect_true(all(delay_early == 6))  # 2 * 3
  }
  if (length(hosp_late) > 0) {
    delay_late <- result_late$time_hospitalisation_relative[hosp_late] - 5
    expect_true(all(delay_late == 2))  # 2 * 1
  }
})

## ------------------------------------------------------------------
## Validation errors in complete_offspring_info
## ------------------------------------------------------------------
test_that("complete_offspring_info errors on invalid prob_hospitalised_genPop", {
  parent <- make_parent_info()
  offspring <- make_offspring_df(5)

  expect_error(
    complete_offspring_info(
      parent_info                = parent,
      offspring_dataframe        = offspring,
      prob_symptomatic           = 0.8,
      prob_hospitalised_hcw      = 0.5,
      prob_hospitalised_genPop   = "bad",
      prob_death_comm            = 0.5,
      prob_death_hosp            = 0.3,
      p_unsafe_funeral_comm_hcw    = 0.5,
      p_unsafe_funeral_hosp_hcw    = 0.2,
      p_unsafe_funeral_comm_genPop = 0.5,
      p_unsafe_funeral_hosp_genPop = 0.2,
      incubation_period          = fixed_incubation,
      onset_to_hospitalisation   = fixed_hosp_delay,
      hospitalisation_to_death   = fixed_death_delay,
      hospitalisation_to_recovery = fixed_recover_delay,
      onset_to_death             = fixed_death_delay,
      onset_to_recovery          = fixed_recover_delay
    ),
    "prob_hospitalised_genPop"
  )
})

test_that("complete_offspring_info errors on invalid hospitalisation_delay_factor", {
  parent <- make_parent_info()
  offspring <- make_offspring_df(5)

  expect_error(
    complete_offspring_info(
      parent_info                = parent,
      offspring_dataframe        = offspring,
      prob_symptomatic           = 0.8,
      prob_hospitalised_hcw      = 0.5,
      prob_hospitalised_genPop   = 0.4,
      prob_death_comm            = 0.5,
      prob_death_hosp            = 0.3,
      p_unsafe_funeral_comm_hcw    = 0.5,
      p_unsafe_funeral_hosp_hcw    = 0.2,
      p_unsafe_funeral_comm_genPop = 0.5,
      p_unsafe_funeral_hosp_genPop = 0.2,
      incubation_period          = fixed_incubation,
      onset_to_hospitalisation   = fixed_hosp_delay,
      hospitalisation_delay_factor = -1,
      hospitalisation_to_death   = fixed_death_delay,
      hospitalisation_to_recovery = fixed_recover_delay,
      onset_to_death             = fixed_death_delay,
      onset_to_recovery          = fixed_recover_delay
    ),
    "hospitalisation_delay_factor"
  )
})

test_that("complete_offspring_info errors when function returns NA", {
  parent <- make_parent_info()
  offspring <- make_offspring_df(5)

  expect_error(
    complete_offspring_info(
      parent_info                = parent,
      offspring_dataframe        = offspring,
      prob_symptomatic           = 1.0,
      prob_hospitalised_hcw      = 0.5,
      prob_hospitalised_genPop   = function(t) NA_real_,
      prob_death_comm            = 0.5,
      prob_death_hosp            = 0.3,
      p_unsafe_funeral_comm_hcw    = 0.5,
      p_unsafe_funeral_hosp_hcw    = 0.2,
      p_unsafe_funeral_comm_genPop = 0.5,
      p_unsafe_funeral_hosp_genPop = 0.2,
      incubation_period          = fixed_incubation,
      onset_to_hospitalisation   = fixed_hosp_delay,
      hospitalisation_to_death   = fixed_death_delay,
      hospitalisation_to_recovery = fixed_recover_delay,
      onset_to_death             = fixed_death_delay,
      onset_to_recovery          = fixed_recover_delay
    ),
    "evaluated to an invalid value"
  )
})
