test_that("moment diagnostics return named statistics", {
  data("us_monthly", package = "qgarch", envir = environment())
  x <- us_monthly$ER
  
  out <- qgarch_moments(x)
  
  expect_s3_class(out, "qgarch_moments")
  expect_true(is.list(out))
  expect_true(all(c("mean", "variance", "skewness", "excess_kurtosis") %in% names(out$estimate)))
  expect_true(all(c("mean", "variance", "skewness", "excess_kurtosis") %in% names(out$std_error)))
  expect_equal(out$n, length(x))
})

test_that("default starts return a non-empty list with named parameters", {
  data("us_monthly", package = "qgarch", envir = environment())
  x <- us_monthly$ER
  
  starts <- qgarch_default_starts(x, model = "free")
  
  expect_true(is.list(starts))
  expect_true(length(starts) >= 1)
  expect_true(is.numeric(starts[[1]]))
  expect_true(all(c("omega", "alpha1", "b", "beta1", "mu", "gamma", "lambda") %in% names(starts[[1]])))
})

test_that("qgarch_fit returns a usable qgarch object for zero model", {
  data("us_monthly", package = "qgarch", envir = environment())
  x <- us_monthly$ER
  
  fit <- qgarch_fit(x, model = "zero")
  
  expect_s3_class(fit, "qgarch")
  expect_equal(fit$model, "zero")
  expect_true(is.numeric(fit$loglik))
  expect_true(is.numeric(fit$nll))
  expect_equal(length(fit$x), length(x))
  expect_equal(length(fit$fitted), length(x))
  expect_equal(length(fit$residuals), length(x))
  expect_equal(length(fit$standardized_residuals), length(x))
  expect_equal(length(fit$sigma2), length(x))
  expect_true(all(c("omega", "alpha1", "b", "beta1", "mu", "gamma") %in% names(fit$coefficients)))
})

test_that("basic S3 methods work for qgarch objects", {
  data("us_monthly", package = "qgarch", envir = environment())
  x <- us_monthly$ER
  
  fit <- qgarch_fit(x, model = "zero")
  
  expect_true(is.numeric(coef(fit)))
  expect_true(is.matrix(vcov(fit)) || all(is.na(vcov(fit))))
  expect_true(is.numeric(fitted(fit)))
  expect_true(is.numeric(residuals(fit)))
  expect_true(is.numeric(residuals(fit, type = "standardized")))
  expect_s3_class(logLik(fit), "logLik")
  expect_s3_class(summary(fit), "summary.qgarch")
  
  fcst <- predict(fit, n.ahead = 3)
  expect_true(is.data.frame(fcst))
  expect_equal(nrow(fcst), 3)
  expect_true(all(c("step", "mean", "sigma2", "sigma") %in% names(fcst)))
})

test_that("likelihood ratio test returns a data frame with expected columns", {
  data("us_monthly", package = "qgarch", envir = environment())
  x <- us_monthly$ER
  
  fit_restricted <- qgarch_fit(x, model = "zero")
  fit_unrestricted <- qgarch_fit(x, model = "free")
  
  lr <- suppressWarnings(qgarch_lr_test(fit_restricted, fit_unrestricted))
  
  expect_s3_class(lr, "qgarch_lr_test")
  expect_true(is.data.frame(lr))
  expect_true(all(c("statistic", "df", "p.value", "restricted_model", "unrestricted_model") %in% names(lr)))
  expect_equal(nrow(lr), 1)
})