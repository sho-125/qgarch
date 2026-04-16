#' Compute sample moments for a series
#'
#' Computes the sample mean, variance, skewness, and excess kurtosis
#' for a numeric series, along with simple standard errors.
#'
#' @param x A numeric vector or a one-column data frame containing the series.
#' @param scale_mean_variance A finite scalar used to rescale the reported mean
#'   and variance and their standard errors.
#'
#' @return A list with elements:
#' \describe{
#'   \item{estimate}{Named numeric vector of moment estimates.}
#'   \item{std_error}{Named numeric vector of standard errors.}
#'   \item{n}{Number of usable observations.}
#' }
#' The returned object has class `"qgarch_moments"`.
#'
#' @export
qgarch_moments <- function(x, scale_mean_variance = 1) {
  if (!is.numeric(scale_mean_variance) ||
      length(scale_mean_variance) != 1L ||
      !is.finite(scale_mean_variance)) {
    stop("'scale_mean_variance' must be a single finite numeric value.", call. = FALSE)
  }
  
  x <- .validate_input_series(x)
  nn <- length(x)
  vx <- stats::var(x)
  sx <- stats::sd(x)
  
  out <- list(
    estimate = c(
      mean = mean(x) * scale_mean_variance,
      variance = vx * scale_mean_variance,
      skewness = .skewness(x),
      excess_kurtosis = .excess_kurtosis(x)
    ),
    std_error = c(
      mean = sx / sqrt(nn) * scale_mean_variance,
      variance = sqrt(2 * vx^2 / (nn - 1)) * scale_mean_variance,
      skewness = sqrt(6 * nn * (nn - 1) / ((nn - 2) * (nn + 1) * (nn + 3))),
      excess_kurtosis =
        2 * sqrt(6 * nn * (nn - 1) / ((nn - 2) * (nn + 1) * (nn + 3))) *
        sqrt((nn^2 - 1) / ((nn - 3) * (nn + 5)))
    ),
    n = nn
  )
  
  class(out) <- "qgarch_moments"
  out
}

#' Likelihood ratio test for nested qgarch models
#'
#' Compares two fitted `qgarch` models using the likelihood ratio test.
#'
#' @param restricted A fitted `qgarch` object representing the restricted model.
#' @param unrestricted A fitted `qgarch` object representing the unrestricted model.
#' @param df Degrees of freedom for the test. If `NULL`, it is computed as the
#'   difference in the number of estimated coefficients.
#'
#' @return A data frame with the likelihood ratio statistic, degrees of freedom,
#'   p-value, and the model names. The returned object has class
#'   `c("qgarch_lr_test", "data.frame")`.
#'
#' @export
qgarch_lr_test <- function(restricted, unrestricted, df = NULL) {
  if (!inherits(restricted, "qgarch") || !inherits(unrestricted, "qgarch")) {
    stop("Both inputs must be 'qgarch' objects.", call. = FALSE)
  }
  
  if (!is.numeric(restricted$nll) || length(restricted$nll) != 1L || !is.finite(restricted$nll)) {
    stop("The restricted model does not contain a finite negative log-likelihood.", call. = FALSE)
  }
  if (!is.numeric(unrestricted$nll) || length(unrestricted$nll) != 1L || !is.finite(unrestricted$nll)) {
    stop("The unrestricted model does not contain a finite negative log-likelihood.", call. = FALSE)
  }
  
  if (is.null(df)) {
    df <- length(unrestricted$coefficients) - length(restricted$coefficients)
  }
  
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 0) {
    stop("'df' must be a single positive number.", call. = FALSE)
  }
  
  stat <- 2 * (restricted$nll - unrestricted$nll)
  
  if (!is.finite(stat)) {
    stop("Likelihood-ratio statistic is not finite.", call. = FALSE)
  }
  
  if (stat < 0) {
    warning(
      "Negative LR statistic encountered; setting statistic to 0. ",
      "This usually indicates numerical optimization error or that the ",
      "unrestricted fit did not improve on the restricted fit.",
      call. = FALSE
    )
    stat <- 0
  }
  
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