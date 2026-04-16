`%||%` <- function(x, y) if (is.null(x)) y else x

#' Extract coefficients from a qgarch model
#'
#' @param object A fitted `qgarch` object.
#' @param ... Additional arguments. Supports `type = "estimated"` or
#'   `type = "full"`.
#'
#' @return A named numeric vector of coefficients.
#' @export
coef.qgarch <- function(object, ...) {
  dots <- list(...)
  type <- dots$type %||% "estimated"
  type <- match.arg(type, c("estimated", "full"))
  if (type == "estimated") object$coefficients else object$coefficients_full
}

#' Variance-covariance matrix for a qgarch model
#'
#' @param object A fitted `qgarch` object.
#' @param ... Unused.
#'
#' @return A variance-covariance matrix.
#' @export
vcov.qgarch <- function(object, ...) {
  object$vcov
}

#' Fitted values from a qgarch model
#'
#' @param object A fitted `qgarch` object.
#' @param ... Unused.
#'
#' @return A numeric vector of fitted values.
#' @export
fitted.qgarch <- function(object, ...) {
  object$fitted
}

#' Residuals from a qgarch model
#'
#' @param object A fitted `qgarch` object.
#' @param type Type of residuals to return: `"raw"`, `"standardized"`,
#'   or `"eta"`.
#' @param ... Unused.
#'
#' @return A numeric vector of residuals.
#' @export
residuals.qgarch <- function(object, type = c("raw", "standardized", "eta"), ...) {
  type <- match.arg(type)
  switch(
    type,
    raw = object$residuals,
    standardized = object$standardized_residuals,
    eta = object$eta
  )
}

#' Log-likelihood for a qgarch model
#'
#' @param object A fitted `qgarch` object.
#' @param ... Unused.
#'
#' @return An object of class `"logLik"`.
#' @export
logLik.qgarch <- function(object, ...) {
  out <- object$loglik
  attr(out, "df") <- length(object$coefficients)
  attr(out, "nobs") <- length(object$x)
  class(out) <- "logLik"
  out
}

#' Forecast from a qgarch model
#'
#' @param object A fitted `qgarch` object.
#' @param n.ahead Number of periods ahead to forecast.
#' @param ... Unused.
#'
#' @return A data frame with forecast horizon, conditional mean,
#'   conditional variance, and conditional standard deviation.
#' @export
predict.qgarch <- function(object, n.ahead = 1L, ...) {
  n.ahead <- as.integer(n.ahead)
  if (n.ahead < 1L) {
    stop("'n.ahead' must be at least 1.", call. = FALSE)
  }
  
  co <- object$coefficients
  omega <- co["omega"]
  alpha <- co["alpha"]
  b <- co["b"]
  beta <- co["beta"]
  mu <- co["mu"]
  gamma <- co["gamma"]
  
  sigma_fcst <- numeric(n.ahead)
  sigma_fcst[1] <- omega +
    alpha * (utils::tail(object$eta, 1) - b)^2 +
    beta * utils::tail(object$sigma2, 1)
  
  if (n.ahead > 1L) {
    for (h in 2:n.ahead) {
      sigma_fcst[h] <- omega + alpha * b^2 + (alpha + beta) * sigma_fcst[h - 1L]
    }
  }
  
  data.frame(
    step = seq_len(n.ahead),
    mean = mu + gamma * sigma_fcst,
    sigma2 = sigma_fcst,
    sigma = sqrt(sigma_fcst),
    row.names = NULL
  )
}

#' Summarize a qgarch model
#'
#' @param object A fitted `qgarch` object.
#' @param ... Unused.
#'
#' @return An object of class `"summary.qgarch"`.
#' @export
summary.qgarch <- function(object, ...) {
  est <- object$coefficients
  se <- object$standard_errors
  z <- est / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  
  coef_table <- cbind(
    Estimate = est,
    Std.Error = se,
    z.value = z,
    p.value = p
  )
  
  if (!is.null(object$implied_lambda)) {
    coef_table <- rbind(
      coef_table,
      lambda = c(
        Estimate = object$implied_lambda,
        Std.Error = object$implied_lambda_se,
        z.value = object$implied_lambda / object$implied_lambda_se,
        p.value = 2 * stats::pnorm(
          abs(object$implied_lambda / object$implied_lambda_se),
          lower.tail = FALSE
        )
      )
    )
  }
  
  structure(
    list(
      model = object$model,
      coefficients = coef_table,
      loglik = object$loglik,
      nobs = length(object$x),
      convergence = object$convergence,
      unconditional_variance = object$unconditional_variance,
      rho = object$rho
    ),
    class = "summary.qgarch"
  )
}

#' Print a qgarch model
#'
#' @param x A fitted `qgarch` object.
#' @param digits Number of digits to print.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.qgarch <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("QGARCH fit\n")
  cat("  model:", x$model, "\n")
  cat("  observations:", length(x$x), "\n")
  cat("  log-likelihood:", format(round(x$loglik, digits), nsmall = digits), "\n")
  cat("  convergence code:", x$convergence, "\n\n")
  print(round(x$coefficients_full, digits))
  invisible(x)
}

#' Print a qgarch summary
#'
#' @param x An object of class `"summary.qgarch"`.
#' @param digits Number of digits to print.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.summary.qgarch <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("QGARCH summary\n")
  cat("  model:", x$model, "\n")
  cat("  observations:", x$nobs, "\n")
  cat("  log-likelihood:", format(round(x$loglik, digits), nsmall = digits), "\n")
  cat("  convergence code:", x$convergence, "\n")
  cat(
    "  unconditional variance:",
    format(round(x$unconditional_variance, digits), nsmall = digits),
    "\n"
  )
  if (!identical(x$model, "zero")) {
    cat("  rho used in restricted mapping:", x$rho, "\n")
  }
  cat("\nCoefficients:\n")
  print(round(x$coefficients, digits))
  invisible(x)
}

#' Plot a qgarch model
#'
#' @param x A fitted `qgarch` object.
#' @param which Which plot to show: `"sigma2"` or `"standardized"`.
#' @param ... Additional graphical arguments passed to [graphics::plot()].
#'
#' @return The input object, invisibly.
#' @export
plot.qgarch <- function(x, which = c("sigma2", "standardized"), ...) {
  which <- match.arg(which)
  
  if (which == "sigma2") {
    graphics::plot(
      x$sigma2,
      type = "l",
      xlab = "Time",
      ylab = expression(sigma[t]^2),
      ...
    )
  } else {
    graphics::plot(
      x$standardized_residuals,
      type = "h",
      xlab = "Time",
      ylab = "Standardized residual",
      ...
    )
    graphics::abline(h = 0, lty = 2)
  }
  
  invisible(x)
}