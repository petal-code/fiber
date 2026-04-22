test_that("make_time_varying returns a function", {
  f <- make_time_varying(c(0, 10), c(0.1, 0.5))
  expect_true(is.function(f))
  expect_s3_class(f, "time_varying_fn")
})

test_that("linear interpolation works correctly", {
  f <- make_time_varying(c(0, 10), c(0.1, 0.5))
  expect_equal(f(0), 0.1)
  expect_equal(f(5), 0.3)
  expect_equal(f(10), 0.5)
})

test_that("linear interpolation clamps outside range", {
  f <- make_time_varying(c(0, 10), c(0.1, 0.5))
  expect_equal(f(-5), 0.1)
  expect_equal(f(100), 0.5)
})

test_that("constant/step method works correctly", {
  f <- make_time_varying(c(0, 10, 20), c(0.1, 0.3, 0.5), method = "constant")
  expect_equal(f(5), 0.1)
  expect_equal(f(10), 0.3)
  expect_equal(f(15), 0.3)
})

test_that("constant method clamps outside range", {
  f <- make_time_varying(c(0, 10, 20), c(0.1, 0.3, 0.5), method = "constant")
  expect_equal(f(-5), 0.1)
  expect_equal(f(100), 0.5)
})

test_that("validation: non-numeric times errors", {
  expect_error(
    make_time_varying(c("a", "b"), c(0.1, 0.5)),
    "`times` must be a numeric vector"
  )
})

test_that("validation: non-numeric values errors", {
  expect_error(
    make_time_varying(c(0, 10), c("a", "b")),
    "`values` must be a numeric vector"
  )
})

test_that("validation: mismatched lengths errors", {
  expect_error(
    make_time_varying(c(0, 10, 20), c(0.1, 0.5)),
    "same length"
  )
})

test_that("validation: non-increasing times errors", {
  expect_error(
    make_time_varying(c(10, 5), c(0.1, 0.5)),
    "non-decreasing"
  )
})

test_that("validation: length-1 inputs error", {
  expect_error(
    make_time_varying(c(0), c(0.1)),
    "length >= 2"
  )
})

test_that("validation: invalid method errors", {
  expect_error(
    make_time_varying(c(0, 10), c(0.1, 0.5), method = "cubic"),
    "'arg' should be one of"
  )
})
