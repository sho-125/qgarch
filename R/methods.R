`%||%` <- function(x, y) if (is.null(x)) y else x

.get_alpha_beta <- function(object) {
  co <- object$coefficients
  arch_order <- object$arch_order %||% 1L
  garch_order <- object$garch_order %||% 1L
  
  alpha_names <- paste0("alpha", seq_len(arch_order))
  beta_names <- paste0("beta", seq_len(garch_order))
  
  if (all(alpha_names %in% names(co))) {
    alpha <- unname(co[alpha_names])
    names(alpha) <- alpha_names
  } else if ("alpha" %in% names(co) && arch_order == 1L) {
    alpha <- unname(co["alpha"])
    names(alpha) <- "alpha1"
  } else {
    stop("Could not find ARCH coefficients in object.", call. = FALSE)
  }
  
  if (all(beta_names %in% names(co))) {
    beta <- unname(co[beta_names])
    names(beta) <- beta_names
  } else if ("beta" %in% names(co) && garch_order == 1L) {
    beta <- unname(co["beta"])
    names(beta) <- "beta1"
  } else {
    stop("Could not find GARCH coefficients in object.", call. = FALSE)
  }
  
  list(
    alpha = alpha,
    beta = beta,
    arch_order = arch_order,
    garch_order = garch_order
  )
}

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

#' Forecast from a generalized qgarch(m, n) model
#'
#' @param object A fitted `qgarch` object.
#' @param n.ahead Number of periods ahead to forecast.
#' @param ... Unused.
#'
#' @return A data frame with forecast horizon, conditional mean,
#'   conditional variance, and conditional standard deviation.
#' @export
predict.qgarch <- function(object, n.ahead = 1L, ...) {
  if (!is.numeric(n.ahead) || length(n.ahead) != 1L || !is.finite(n.ahead)) {
    stop("'n.ahead' must be a single finite number.", call. = FALSE)
  }
  
  n.ahead <- as.integer(n.ahead)
  if (n.ahead < 1L) {
    stop("'n.ahead' must be at least 1.", call. = FALSE)
  }
  
  blocks <- .get_alpha_beta(object)
  alpha <- blocks$alpha
  beta <- blocks$beta
  arch_order <- blocks$arch_order
  garch_order <- blocks$garch_order
  
  co <- object$coefficients
  omega <- unname(co["omega"])
  b <- unname(co["b"])
  mu <- unname(co["mu"])
  gamma <- unname(co["gamma"])
  
  eta_hist <- object$eta
  sigma_hist <- object$sigma2
  uncond <- object$unconditional_variance %||% stats::var(object$x, na.rm = TRUE)
  
  sigma_fcst <- numeric(n.ahead)
  
  for (h in seq_len(n.ahead)) {
    arch_term <- 0
    for (j in seq_len(arch_order)) {
      idx <- h - j
      if (idx <= 0L) {
        hist_pos <- length(eta_hist) + idx
        eta_lag <- if (hist_pos >= 1L) eta_hist[hist_pos] else 0
        arch_piece <- (eta_lag - b)^2
      } else {
        arch_piece <- sigma_fcst[idx] + b^2
      }
      arch_term <- arch_term + alpha[j] * arch_piece
    }
    
    garch_term <- 0
    for (k in seq_len(garch_order)) {
      idx <- h - k
      if (idx <= 0L) {
        hist_pos <- length(sigma_hist) + idx
        sigma_lag <- if (hist_pos >= 1L) sigma_hist[hist_pos] else uncond
      } else {
        sigma_lag <- sigma_fcst[idx]
      }
      garch_term <- garch_term + beta[k] * sigma_lag
    }
    
    sigma_fcst[h] <- omega + arch_term + garch_term
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
  
  if (is.null(se) || length(se) != length(est)) {
    se <- rep(NA_real_, length(est))
    names(se) <- names(est)
  }
  
  z <- rep(NA_real_, length(est))
  ok <- is.finite(est) & is.finite(se) & se > 0
  z[ok] <- est[ok] / se[ok]
  
  p <- rep(NA_real_, length(est))
  p[ok] <- 2 * stats::pnorm(abs(z[ok]), lower.tail = FALSE)
  
  coef_table <- cbind(
    Estimate = est,
    Std.Error = se,
    z.value = z,
    p.value = p
  )
  
  if (!is.null(object$implied_lambda)) {
    lambda_z <- NA_real_
    lambda_p <- NA_real_
    
    if (is.finite(object$implied_lambda_se) && object$implied_lambda_se > 0) {
      lambda_z <- object$implied_lambda / object$implied_lambda_se
      lambda_p <- 2 * stats::pnorm(abs(lambda_z), lower.tail = FALSE)
    }
    
    coef_table <- rbind(
      coef_table,
      lambda = c(
        Estimate = object$implied_lambda,
        Std.Error = object$implied_lambda_se,
        z.value = lambda_z,
        p.value = lambda_p
      )
    )
  }
  
  structure(
    list(
      model = object$model,
      arch_order = object$arch_order %||% 1L,
      garch_order = object$garch_order %||% 1L,
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
  arch_order <- x$arch_order %||% 1L
  garch_order <- x$garch_order %||% 1L
  
  cat("QGARCH fit\n")
  cat("  model:", x$model, "\n")
  cat("  order:", sprintf("(%d, %d)", arch_order, garch_order), "\n")
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
  cat("  order:", sprintf("(%d, %d)", x$arch_order, x$garch_order), "\n")
  cat("  observations:", x$nobs, "\n")
  cat("  log-likelihood:", format(round(x$loglik, digits), nsmall = digits), "\n")
  cat("  convergence code:", x$convergence, "\n")
  cat(
    "  unconditional variance:",
    format(round(x$unconditional_variance, digits), nsmall = digits),
    "\n"
  )
  if (identical(x$model, "restricted") && !is.null(x$rho)) {
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