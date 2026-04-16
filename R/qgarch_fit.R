# Internal penalty used by likelihood functions.
.qgarch_big_penalty <- 1e200

.validate_input_series <- function(x) {
  if (is.data.frame(x) && ncol(x) == 1L) {
    x <- x[[1L]]
  }
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 10L) {
    stop("'x' must contain at least 10 finite observations.", call. = FALSE)
  }
  x
}

.unconditional_variance <- function(omega, alpha, b, beta) {
  den <- 1 - alpha - beta
  if (!is.finite(den) || den <= 0 || omega <= 0 || alpha < 0 || beta < 0) {
    return(NA_real_)
  }
  (omega + alpha * b^2) / den
}

.safe_invert <- function(mat, tol = 1e-7) {
  out <- tryCatch(solve(mat), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(out)) {
    out <- MASS::ginv(mat, tol = tol)
  }
  out
}

.diag_se <- function(mat) {
  inv <- .safe_invert(mat)
  se <- sqrt(pmax(diag(inv), 0))
  list(se = se, vcov = inv)
}

.skewness <- function(x) {
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  mean((x - m)^3) / s^3
}

.excess_kurtosis <- function(x) {
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  mean((x - m)^4) / s^4 - 3
}

.make_named_start <- function(model, values) {
  model <- match.arg(model, c("zero", "restricted", "free", "threshold"))
  nm <- switch(
    model,
    zero = c("omega", "alpha", "b", "beta", "mu", "gamma"),
    restricted = c("omega", "alpha", "b", "beta", "mu", "gamma"),
    free = c("omega", "alpha", "b", "beta", "mu", "gamma", "lambda"),
    threshold = c("omega", "alpha", "b", "beta", "mu", "gamma", "lambda1", "lambda2")
  )
  stats::setNames(as.numeric(values), nm)
}

qgarch_default_starts <- function(x, model = c("zero", "restricted", "free", "threshold")) {
  model <- match.arg(model)
  x <- .validate_input_series(x)
  vx <- stats::var(x)
  sx <- stats::sd(x)
  mx <- mean(x)
  omega_base <- max(vx * 0.05, 1e-6)
  b_small <- sx * 0.10
  base <- switch(
    model,
    zero = list(
      c(omega_base, 0.05, b_small, 0.90, mx, 0.00),
      c(omega_base * 0.5, 0.10, 0.00, 0.85, mx, 0.05),
      c(omega_base * 1.5, 0.08, b_small * 0.5, 0.80, mx, -0.05)
    ),
    restricted = list(
      c(omega_base, 0.05, b_small, 0.90, mx, 0.05),
      c(omega_base * 0.5, 0.10, 0.00, 0.85, mx, 0.10),
      c(omega_base * 1.5, 0.08, b_small * 0.5, 0.80, mx, -0.02)
    ),
    free = list(
      c(omega_base, 0.05, b_small, 0.90, mx, 0.00, 0.005),
      c(omega_base * 0.5, 0.10, 0.00, 0.85, mx, 0.05, 0.010),
      c(omega_base * 1.5, 0.08, b_small * 0.5, 0.80, mx, -0.05, 0.001)
    ),
    threshold = list(
      c(omega_base, 0.05, b_small, 0.90, mx, 0.00, 0.005, 0.000),
      c(omega_base * 0.5, 0.10, 0.00, 0.85, mx, 0.05, 0.010, 0.005),
      c(omega_base * 1.5, 0.08, b_small * 0.5, 0.80, mx, -0.05, 0.001, -0.001)
    )
  )
  lapply(base, function(z) .make_named_start(model, z))
}

.solve_eta_quadratic <- function(obs, sigma2, mu, gamma, lambda, kappa, tol = 1e-10) {
  if (abs(lambda) <= tol) {
    return((obs - mu - gamma * sigma2) / kappa)
  }
  rhs <- obs - mu - (gamma + lambda) * sigma2
  disc <- kappa^2 - 4 * lambda * rhs
  if (!is.finite(disc) || disc < 0) {
    return(NA_real_)
  }
  (kappa - sqrt(disc)) / (2 * lambda)
}

.qgarch_nll <- function(theta, x, model, rho = 1, indicator = NULL, return_details = FALSE) {
  n <- length(x)
  eta <- numeric(n)
  sigma2 <- numeric(n)
  logl <- 0

  if (model == "zero") {
    omega <- theta[1]
    alpha <- theta[2]
    b <- theta[3]
    beta <- theta[4]
    mu <- theta[5]
    gamma <- theta[6]
    lambda_vec <- rep(0, n)
    kappa <- 1
  } else if (model == "restricted") {
    omega <- theta[1]
    alpha <- theta[2]
    b <- theta[3]
    beta <- theta[4]
    mu <- theta[5]
    gamma <- theta[6]
    lambda <- rho * gamma * alpha / (1 - rho * (alpha + beta))
    lambda_vec <- rep(lambda, n)
    kappa <- 1 + 2 * lambda * b
  } else if (model == "free") {
    omega <- theta[1]
    alpha <- theta[2]
    b <- theta[3]
    beta <- theta[4]
    mu <- theta[5]
    gamma <- theta[6]
    lambda <- theta[7]
    lambda_vec <- rep(lambda, n)
    kappa <- 1 + 2 * lambda * b
  } else if (model == "threshold") {
    omega <- theta[1]
    alpha <- theta[2]
    b <- theta[3]
    beta <- theta[4]
    mu <- theta[5]
    gamma <- theta[6]
    lambda1 <- theta[7]
    lambda2 <- theta[8]
    lambda_vec <- lambda1 + lambda2 * indicator
    kappa <- 1 + 2 * (lambda1 + mean(indicator) * lambda2) * b
  } else {
    stop("Unknown model.", call. = FALSE)
  }

  uncond <- .unconditional_variance(omega, alpha, b, beta)
  if (!is.finite(uncond)) {
    return(.qgarch_big_penalty)
  }

  for (i in seq_len(n)) {
    if (i == 1L) {
      sigma2[i] <- uncond
      if (model == "zero") {
        eta[i] <- x[i] - mu - gamma * sigma2[i]
      } else if (model == "restricted") {
        eta[i] <- .solve_eta_quadratic(x[i], sigma2[i], mu, gamma, lambda_vec[i], kappa)
      } else {
        eta[i] <- (x[i] - mu - gamma * sigma2[i]) / kappa
      }
    } else {
      sigma2[i] <- omega + alpha * (eta[i - 1L] - b)^2 + beta * sigma2[i - 1L]
      if (model == "zero") {
        eta[i] <- x[i] - mu - gamma * sigma2[i]
      } else {
        eta[i] <- .solve_eta_quadratic(x[i], sigma2[i], mu, gamma, lambda_vec[i], kappa)
      }
    }

    if (!is.finite(sigma2[i]) || sigma2[i] <= 0 || !is.finite(eta[i])) {
      return(.qgarch_big_penalty)
    }

    if (model == "zero") {
      logl <- logl + log((1 / sqrt(2 * pi * sigma2[i])) * exp(-(eta[i]^2) / (2 * sigma2[i])))
    } else {
      jac <- kappa - 2 * lambda_vec[i] * eta[i]
      if (!is.finite(jac) || jac <= 0) {
        return(.qgarch_big_penalty)
      }
      logl <- logl - log(jac) - 0.5 * log(2 * pi * sigma2[i]) - 0.5 * (eta[i]^2 / sigma2[i])
    }
  }

  nll <- -logl
  if (return_details) {
    return(list(nll = nll, eta = eta, sigma2 = sigma2, uncond_var = uncond, lambda = lambda_vec, kappa = kappa))
  }
  nll
}

.restricted_lambda_delta_se <- function(alpha, beta, gamma, vcov_sub, rho = 1) {
  den <- 1 - rho * (alpha + beta)
  grad <- c(
    rho * alpha / den,
    gamma * rho * (1 - rho * beta) / (den^2),
    gamma * rho^2 * alpha / (den^2)
  )
  se2 <- as.numeric(t(grad) %*% vcov_sub %*% grad)
  sqrt(max(se2, 0))
}

.prepare_start_list <- function(start, x, model) {
  if (is.null(start)) {
    return(qgarch_default_starts(x, model = model))
  }
  if (is.list(start) && !is.null(names(start))) {
    return(list(.make_named_start(model, unname(unlist(start)))))
  }
  if (is.list(start) && length(start) > 0 && is.numeric(start[[1L]])) {
    return(lapply(start, function(z) .make_named_start(model, z)))
  }
  if (is.numeric(start)) {
    return(list(.make_named_start(model, start)))
  }
  stop("'start' must be NULL, a numeric vector, or a list of numeric vectors.", call. = FALSE)
}

.build_fit_object <- function(best_fit, x, model, rho, indicator, starts, call) {
  details <- .qgarch_nll(best_fit$estimate, x = x, model = model, rho = rho, indicator = indicator, return_details = TRUE)
  hessian <- best_fit$hessian
  inv <- .diag_se(hessian)

  est_names <- names(best_fit$estimate)
  names(inv$se) <- est_names
  rownames(inv$vcov) <- est_names
  colnames(inv$vcov) <- est_names

  coef_full <- best_fit$estimate
  full_se <- inv$se

  if (model == "zero") {
    coef_full <- c(coef_full, lambda = 0)
    full_se <- c(full_se, lambda = NA_real_)
  }

  implied_lambda <- NULL
  implied_lambda_se <- NULL

  if (model == "restricted") {
    alpha <- best_fit$estimate["alpha"]
    beta <- best_fit$estimate["beta"]
    gamma <- best_fit$estimate["gamma"]
    implied_lambda <- rho * gamma * alpha / (1 - rho * (alpha + beta))
    idx <- c("gamma", "alpha", "beta")
    implied_lambda_se <- .restricted_lambda_delta_se(
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      vcov_sub = inv$vcov[idx, idx, drop = FALSE],
      rho = rho
    )
    coef_full <- c(coef_full, lambda = implied_lambda)
    full_se <- c(full_se, lambda = implied_lambda_se)
  }

  structure(
    list(
      call = call,
      model = model,
      coefficients = best_fit$estimate,
      coefficients_full = coef_full,
      standard_errors = inv$se,
      standard_errors_full = full_se,
      vcov = inv$vcov,
      loglik = -best_fit$minimum,
      nll = best_fit$minimum,
      eta = details$eta,
      sigma2 = details$sigma2,
      fitted = best_fit$estimate["mu"] + best_fit$estimate["gamma"] * details$sigma2,
      residuals = x - (best_fit$estimate["mu"] + best_fit$estimate["gamma"] * details$sigma2),
      standardized_residuals = details$eta / sqrt(details$sigma2),
      lambda_t = details$lambda,
      kappa = details$kappa,
      unconditional_variance = details$uncond_var,
      hessian = hessian,
      convergence = best_fit$code,
      iterations = best_fit$iterations,
      minimum = best_fit$minimum,
      x = x,
      indicator = indicator,
      rho = rho,
      starts_tried = starts,
      implied_lambda = implied_lambda,
      implied_lambda_se = implied_lambda_se
    ),
    class = "qgarch"
  )
}

qgarch_fit <- function(x,
                       model = c("zero", "restricted", "free", "threshold"),
                       threshold_indicator = NULL,
                       start = NULL,
                       rho = 1,
                       steptol = 1e-10,
                       typsize = 0.1,
                       print.level = 0,
                       hessian = TRUE) {
  model <- match.arg(model)
  x <- .validate_input_series(x)

  if (model == "threshold") {
    if (is.null(threshold_indicator)) {
      stop("'threshold_indicator' must be supplied when model = 'threshold'.", call. = FALSE)
    }
    indicator <- as.numeric(threshold_indicator)
    if (length(indicator) != length(x)) {
      stop("'threshold_indicator' must have the same length as 'x'.", call. = FALSE)
    }
    indicator <- ifelse(indicator != 0, 1, 0)
  } else {
    indicator <- NULL
  }

  starts <- .prepare_start_list(start = start, x = x, model = model)
  fits <- list()
  values <- rep(NA_real_, length(starts))

  for (i in seq_along(starts)) {
    fit_i <- try(
      stats::nlm(
        f = .qgarch_nll,
        p = unname(starts[[i]]),
        x = x,
        model = model,
        rho = rho,
        indicator = indicator,
        hessian = hessian,
        steptol = steptol,
        typsize = rep(typsize, length(starts[[i]])),
        print.level = print.level
      ),
      silent = TRUE
    )
    if (!inherits(fit_i, "try-error") && is.finite(fit_i$minimum) && fit_i$minimum < .qgarch_big_penalty) {
      names(fit_i$estimate) <- names(starts[[i]])
      fits[[i]] <- fit_i
      values[i] <- fit_i$minimum
    }
  }

  if (!any(is.finite(values))) {
    stop("No converged solution was found. Try supplying custom starting values.", call. = FALSE)
  }

  best_fit <- fits[[which.min(values)]]
  .build_fit_object(best_fit = best_fit, x = x, model = model, rho = rho, indicator = indicator, starts = starts, call = match.call())
}
