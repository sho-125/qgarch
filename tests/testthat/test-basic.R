test_that("moment diagnostics return named statistics", {
  out <- qgarch_moments(rnorm(50))
  expect_true(is.list(out))
  expect_true(all(c("mean", "variance", "skewness", "excess_kurtosis") %in% names(out$estimate)))
})

test_that("default starts return a non-empty list", {
  starts <- qgarch_default_starts(rnorm(100), model = "free")
  expect_true(is.list(starts))
  expect_true(length(starts) >= 1)
})
