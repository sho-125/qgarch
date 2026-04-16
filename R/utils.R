qgarch_moments <- function(x, scale_mean_variance = 1) {
  x <- .validate_input_series(x)
  nn <- length(x)

  out <- list(
    estimate = c(
      mean = mean(x) * scale_mean_variance,
      variance = stats::var(x) * scale_mean_variance,
      skewness = .skewness(x),
      excess_kurtosis = .excess_kurtosis(x)
    ),
    std_error = c(
      mean = stats::sd(x) / sqrt(nn) * scale_mean_variance,
      variance = sqrt(2 * stats::var(x)^2 / (nn - 1)) * scale_mean_variance,
      skewness = sqrt(6 * nn * (nn - 1) / ((nn - 2) * (nn + 1) * (nn + 3))),
      excess_kurtosis = 2 * sqrt(6 * nn * (nn - 1) / ((nn - 2) * (nn + 1) * (nn + 3))) * sqrt((nn^2 - 1) / ((nn - 3) * (nn + 5)))
    ),
    n = nn
  )
  class(out) <- "qgarch_moments"
  out
}

qgarch_lr_test <- function(restricted, unrestricted, df = NULL) {
  if (!inherits(restricted, "qgarch") || !inherits(unrestricted, "qgarch")) {
    stop("Both inputs must be 'qgarch' objects.", call. = FALSE)
  }
  if (is.null(df)) {
    df <- length(unrestricted$coefficients) - length(restricted$coefficients)
  }
  stat <- 2 * (restricted$nll - unrestricted$nll)
  pval <- stats::pchisq(stat, df = df, lower.tail = FALSE)
  out <- data.frame(
    statistic = stat,
    df = df,
    p.value = pval,
    restricted_model = restricted$model,
    unrestricted_model = unrestricted$model,
    row.names = NULL
  )
  class(out) <- c("qgarch_lr_test", class(out))
  out
}
