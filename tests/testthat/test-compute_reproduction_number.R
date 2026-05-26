## Tests for compute_reproduction_number().

## --- helpers -----------------------------------------------------------

## Build a minimal tdf for testing. cases are at times t = seq(...).
## n_offspring is per case; parent links them in a simple chain (each case
## is parent of the next, except the first which has parent = NA).
make_test_tdf <- function(times, n_offspring,
                          offspring_generated = rep(TRUE, length(times)),
                          parents = NULL) {
  n <- length(times)
  if (is.null(parents)) parents <- c(NA_integer_, seq_len(n - 1L))
  data.frame(
    id                      = seq_len(n),
    parent                  = parents,
    time_infection_absolute = times,
    n_offspring             = n_offspring,
    offspring_generated     = offspring_generated,
    stringsAsFactors        = FALSE
  )
}

## --- input validation --------------------------------------------------

test_that("compute_reproduction_number errors on missing columns", {
  bad_tdf <- data.frame(id = 1, time_infection_absolute = 0)
  expect_error(compute_reproduction_number(bad_tdf),
               "missing required column")
})

test_that("compute_reproduction_number errors on bad bin_width", {
  tdf <- make_test_tdf(times = 0:9, n_offspring = rep(1, 10))
  expect_error(compute_reproduction_number(tdf, bin_width = 0),
               "must be a single positive numeric")
  expect_error(compute_reproduction_number(tdf, bin_width = c(1, 2)),
               "must be a single positive numeric")
})

test_that("compute_reproduction_number errors when no cases remain", {
  tdf <- make_test_tdf(times = 0:4, n_offspring = rep(1, 5),
                       offspring_generated = rep(FALSE, 5))
  ## Two conditions are raised: a warning (about the exclusions) AND an error
  ## (because nothing is left). Check them separately to avoid testthat 3.0
  ## nesting surprises.
  expect_error(
    suppressWarnings(compute_reproduction_number(tdf)),
    "No cases remain"
  )
  expect_warning(
    tryCatch(compute_reproduction_number(tdf),
             error = function(e) invisible(NULL)),
    "5 individual"
  )
})

## --- exclusion warning --------------------------------------------------

test_that("compute_reproduction_number warns about excluded cases", {
  tdf <- make_test_tdf(times = c(0, 1, 2, 3),
                       n_offspring = c(2, 3, 1, NA),
                       offspring_generated = c(TRUE, TRUE, TRUE, FALSE))
  expect_warning(
    out <- compute_reproduction_number(tdf, bin_width = 1, min_cases_per_bin = 1,
                                       type = "case"),
    "1 individual\\(s\\) had offspring_generated == FALSE"
  )
  ## The excluded case (id 4) should not contribute to any bin
  expect_equal(sum(out$n_cases), 3)
})

## --- case R correctness -----------------------------------------------

test_that("case R equals mean(n_offspring) per bin", {
  ## 10 cases at t = 0..9, n_offspring = 2 for first 5 cases, 3 for next 5.
  ## With bin_width = 5 we expect two bins: [0,5) with mean 2 and [5,10) with mean 3.
  tdf <- make_test_tdf(times = 0:9, n_offspring = c(rep(2, 5), rep(3, 5)))
  out <- compute_reproduction_number(tdf, bin_width = 5,
                                     min_cases_per_bin = 1, type = "case")
  expect_equal(nrow(out), 2)
  expect_equal(out$R_case, c(2, 3))
  expect_equal(out$n_cases, c(5L, 5L))
})

test_that("bins below min_cases_per_bin get NA for R_case", {
  ## 5 cases in bin [0,5), 1 case in bin [5,10).
  tdf <- make_test_tdf(times = c(0, 1, 2, 3, 4, 7),
                       n_offspring = c(1, 1, 2, 2, 1, 5))
  out <- compute_reproduction_number(tdf, bin_width = 5,
                                     min_cases_per_bin = 5, type = "case")
  ## First bin has 5 cases (meets threshold), second has 1 (below).
  expect_false(is.na(out$R_case[1]))
  expect_true(is.na(out$R_case[2]))
})

test_that("case R can be requested alone (type = 'case')", {
  tdf <- make_test_tdf(times = 0:9, n_offspring = rep(2, 10))
  out <- compute_reproduction_number(tdf, bin_width = 5,
                                     min_cases_per_bin = 1, type = "case")
  expect_true("R_case" %in% names(out))
  expect_false("R_instantaneous" %in% names(out))
})

## --- instantaneous R: structure and provided GI -----------------------

test_that("instantaneous R can be requested alone (type = 'instantaneous')", {
  tdf <- make_test_tdf(times = 0:9, n_offspring = rep(1, 10))
  out <- compute_reproduction_number(tdf, bin_width = 1,
                                     min_cases_per_bin = 1,
                                     type = "instantaneous",
                                     generation_interval = c(1))
  expect_true("R_instantaneous" %in% names(out))
  expect_false("R_case" %in% names(out))
})

test_that("instantaneous R with delta-function GI at lag 1 equals I(t)/I(t-1)", {
  ## Build incidence 1, 2, 4, 8, 16 over bins 0..4 with bin_width = 1.
  ## n_offspring is irrelevant for instantaneous R.
  inc <- c(1, 2, 4, 8, 16)
  times <- unlist(lapply(seq_along(inc), function(i) rep(i - 1L, inc[i])))
  tdf <- make_test_tdf(times = times,
                       n_offspring = rep(0, length(times)),
                       parents = c(NA, rep(1, length(times) - 1L)))
  ## GI pmf at lag 1 only (i.e. all transmission happens with exactly 1 bin lag).
  out <- compute_reproduction_number(tdf, bin_width = 1,
                                     min_cases_per_bin = 1,
                                     type = "instantaneous",
                                     generation_interval = c(1))
  ## R(t) = I(t) / (I(t-1) * 1) for t >= 1. First bin has no predecessor -> NA.
  expected_R <- c(NA, 2, 2, 2, 2)
  expect_equal(out$R_instantaneous, expected_R, tolerance = 1e-9)
})

test_that("instantaneous R with provided GI normalises non-pmf inputs", {
  ## GI = c(2, 4) is not a pmf (sums to 6); should normalise to c(1/3, 2/3).
  inc <- c(3, 6, 12)
  times <- unlist(lapply(seq_along(inc), function(i) rep(i - 1L, inc[i])))
  tdf <- make_test_tdf(times = times,
                       n_offspring = rep(0, length(times)),
                       parents = c(NA, rep(1, length(times) - 1L)))
  out <- compute_reproduction_number(tdf, bin_width = 1,
                                     min_cases_per_bin = 1,
                                     type = "instantaneous",
                                     generation_interval = c(2, 4))
  ## At t=2 (bin index 3 in grid), denom = I(1)*1/3 + I(0)*2/3 = 6/3 + 3*2/3 = 4
  ## R(2) = I(2)/denom = 12 / 4 = 3
  expect_equal(out$R_instantaneous[3], 3, tolerance = 1e-9)
})

test_that("instantaneous R rejects malformed generation_interval", {
  tdf <- make_test_tdf(times = 0:9, n_offspring = rep(1, 10))
  expect_error(compute_reproduction_number(tdf, type = "instantaneous",
                                           generation_interval = c(-1, 1)),
               "non-negative")
  expect_error(compute_reproduction_number(tdf, type = "instantaneous",
                                           generation_interval = c(0, 0)),
               "positive sum")
})

## --- instantaneous R: empirical GI estimation -------------------------

test_that("empirical GI is estimated from parent-child gaps when not supplied", {
  ## 10-case chain with each gap = 2 bin_widths (so empirical GI pmf = c(0, 1)).
  tdf <- make_test_tdf(times = seq(0, 18, by = 2), n_offspring = rep(1, 10))
  ## Should NOT error and should compute some R_inst values.
  out <- compute_reproduction_number(tdf, bin_width = 1,
                                     min_cases_per_bin = 1,
                                     type = "instantaneous")
  ## R_inst[1] is NA (no predecessor); R_inst[2] is NA (only 1 incidence,
  ## but more importantly the GI requires lag 2). R_inst[3] = I(2)/(I(0)*1)
  ## with empirical GI pmf c(0, 1) -> R(2) = 1/1 = 1.
  expect_true(any(!is.na(out$R_instantaneous)))
})

test_that("compute_reproduction_number warns when fewer than 5 parent-child gaps", {
  ## Only 3 parent-child links (4 cases in a chain).
  tdf <- make_test_tdf(times = c(0, 5, 10, 15), n_offspring = c(1, 1, 1, 0))
  expect_warning(
    compute_reproduction_number(tdf, bin_width = 5, min_cases_per_bin = 1,
                                type = "instantaneous"),
    "Only 3 parent-child gap"
  )
})

## --- output structure -------------------------------------------------

test_that("output covers a regular grid from first to last bin", {
  ## 3 cases at t = 0 and 2 at t = 20 (bin_width = 5).
  ## Expected grid: 0, 5, 10, 15, 20 (5 bins).
  tdf <- make_test_tdf(times = c(0, 0, 0, 20, 20),
                       n_offspring = c(1, 1, 1, 0, 0))
  out <- compute_reproduction_number(tdf, bin_width = 5,
                                     min_cases_per_bin = 1, type = "case")
  expect_equal(out$time, c(0, 5, 10, 15, 20))
  expect_equal(out$n_cases, c(3L, 0L, 0L, 0L, 2L))
})

test_that("type='both' returns both R columns", {
  tdf <- make_test_tdf(times = 0:9, n_offspring = rep(1, 10))
  out <- compute_reproduction_number(tdf, bin_width = 1,
                                     min_cases_per_bin = 1,
                                     type = "both",
                                     generation_interval = c(1))
  expect_true(all(c("R_case", "R_instantaneous") %in% names(out)))
})
