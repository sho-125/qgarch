# Internal penalty used by likelihood functions.
.qgarch_big_penalty <- 1e200

.validate_input_series <- function(x) {
  if (is.data.frame(x) && ncol(x) == 1L) {
    x <- x[[1L]]
  }
  
  x <- as.numeric(x)
  
  if (length(x) < 10L) {
    stop("'x' must contain at least 10 observations.", call. = FALSE)
  }
  
  if (anyNA(x) || any(!is.finite(x))) {
    stop("'x' must contain only finite, non-missing values.", call. = FALSE)
  }
  
  x
}

.validate_orders <- function(arch_order, garch_order) {
  if (!is.numeric(arch_order) || length(arch_order) != 1L || !is.finite(arch_order)) {
    stop("'arch_order' must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(garch_order) || length(garch_order) != 1L || !is.finite(garch_order)) {
    stop("'garch_order' must be a single positive integer.", call. = FALSE)
  }
  
  arch_order <- as.integer(arch_order)
  garch_order <- as.integer(garch_order)
  
  if (arch_order < 1L || garch_order < 1L) {
    stop("'arch_order' and 'garch_order' must both be at least 1.", call. = FALSE)
  }
  
  list(arch_order = arch_order, garch_order = garch_order)
}

.unconditional_variance <- function(omega, alpha, b, beta) {
  alpha_sum <- sum(alpha)
  beta_sum <- sum(beta)
  den <- 1 - alpha_sum - beta_sum
  
  if (!is.finite(den) || den <= 0 || !is.finite(omega) || omega <= 0) {
    return(NA_real_)
  }
  if (any(!is.finite(alpha)) || any(alpha < 0)) {
    return(NA_real_)
  }
  if (any(!is.finite(beta)) || any(beta < 0)) {
    return(NA_real_)
  }
  if (!is.finite(b)) {
    return(NA_real_)
  }
  
  (omega + alpha_sum * b^2) / den
}

.clean_sym_matrix <- function(mat) {
  if (is.null(mat) || !is.matrix(mat)) {
    return(NULL)
  }
  0.5 * (mat + t(mat))
}

.safe_invert <- function(mat, tol = 1e-7) {
  out <- tryCatch(
    solve(mat),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  
  if (is.null(out)) {
    out <- MASS::ginv(mat, tol = tol)
  }
  
  out
}

.regularize_pd <- function(mat, rel_floor = 1e-8, abs_floor = 1e-10) {
  mat <- .clean_sym_matrix(mat)
  
  if (is.null(mat) || any(!is.finite(mat))) {
    return(NULL)
  }
  
  ee <- tryCatch(
    eigen(mat, symmetric = TRUE),
    error = function(e) NULL
  )
  if (is.null(ee)) {
    return(NULL)
  }
  
  max_eval <- max(abs(ee$values), na.rm = TRUE)
  floor_val <- max(max_eval * rel_floor, abs_floor)
  
  vals <- ee$values
  vals[!is.finite(vals)] <- floor_val
  vals[vals < floor_val] <- floor_val
  
  out <- ee$vectors %*% (vals * t(ee$vectors))
  .clean_sym_matrix(out)
}

.make_se_from_vcov <- function(vcov_mat, tol = 1e-12) {
  vcov_mat <- .clean_sym_matrix(vcov_mat)
  
  if (is.null(vcov_mat) || any(!is.finite(vcov_mat))) {
    return(list(se = NULL, vcov = NULL))
  }
  
  dv <- diag(vcov_mat)
  se <- rep(NA_real_, length(dv))
  ok <- is.finite(dv) & dv > tol
  se[ok] <- sqrt(dv[ok])
  
  list(se = se, vcov = vcov_mat)
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

.lag_weights <- function(k, pattern = c("geom", "equal", "front"), decay = 0.65) {
  pattern <- match.arg(pattern)
  
  w <- switch(
    pattern,
    geom = decay^(seq_len(k) - 1L),
    equal = rep(1, k),
    front = rev(seq_len(k))
  )
  
  w / sum(w)
}

.parameter_names <- function(model, arch_order, garch_order) {
  model <- match.arg(model, c("zero", "restricted", "free", "threshold"))
  
  base <- c(
    "omega",
    paste0("alpha", seq_len(arch_order)),
    "b",
    paste0("beta", seq_len(garch_order)),
    "mu",
    "gamma"
  )
  
  extra <- switch(
    model,
    zero = character(0),
    restricted = character(0),
    free = "lambda",
    threshold = c("lambda1", "lambda2")
  )
  
  c(base, extra)
}

.make_named_start <- function(model, arch_order, garch_order, values) {
  nm <- .parameter_names(model, arch_order = arch_order, garch_order = garch_order)
  values <- as.numeric(values)
  
  if (length(values) != length(nm)) {
    stop(
      "Starting values have length ", length(values),
      " but model/order combination requires ", length(nm), ".",
      call. = FALSE
    )
  }
  
  stats::setNames(values, nm)
}

.target_omega <- function(target_var, alpha_total, b, beta_total) {
  out <- target_var * (1 - alpha_total - beta_total) - alpha_total * b^2
  max(out, 1e-8)
}

.compose_start_values <- function(model,
                                  arch_order,
                                  garch_order,
                                  target_var,
                                  alpha_total,
                                  b,
                                  beta_total,
                                  mu,
                                  gamma,
                                  lambda = NULL,
                                  lambda1 = NULL,
                                  lambda2 = NULL,
                                  alpha_pattern = c("geom", "equal", "front"),
                                  beta_pattern = c("geom", "equal", "front")) {
  alpha_pattern <- match.arg(alpha_pattern)
  beta_pattern <- match.arg(beta_pattern)
  
  alpha <- alpha_total * .lag_weights(arch_order, pattern = alpha_pattern)
  beta <- beta_total * .lag_weights(garch_order, pattern = beta_pattern)
  omega <- .target_omega(target_var, alpha_total = sum(alpha), b = b, beta_total = sum(beta))
  
  out <- c(
    omega,
    alpha,
    b,
    beta,
    mu,
    gamma
  )
  
  if (model == "free") {
    out <- c(out, lambda)
  } else if (model == "threshold") {
    out <- c(out, lambda1, lambda2)
  }
  
  .make_named_start(model, arch_order, garch_order, out)
}

.restricted_lambda_value <- function(alpha, beta, gamma, rho = 1) {
  alpha_sum <- sum(alpha)
  beta_sum <- sum(beta)
  den <- 1 - rho * (alpha_sum + beta_sum)
  
  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  
  rho * gamma * alpha_sum / den
}

#' Default starting values for qgarch estimation
#'
#' Creates a compact but order-adaptive set of starting values for the selected
#' qgarch model.
#'
#' @param x A numeric vector or one-column data frame containing the observed
#'   series.
#' @param model Character string specifying the model variant. Must be one of
#'   `"zero"`, `"restricted"`, `"free"`, or `"threshold"`.
#' @param arch_order Positive integer giving the ARCH lag order.
#' @param garch_order Positive integer giving the GARCH lag order.
#'
#' @return A list of named numeric vectors containing candidate starting values.
#' @export
qgarch_default_starts <- function(x,
                                  model = c("zero", "restricted", "free", "threshold"),
                                  arch_order = 1L,
                                  garch_order = 1L) {
  model <- match.arg(model)
  orders <- .validate_orders(arch_order, garch_order)
  arch_order <- orders$arch_order
  garch_order <- orders$garch_order
  
  x <- .validate_input_series(x)
  
  vx <- stats::var(x)
  sx <- stats::sd(x)
  mx <- mean(x)
  b_small <- max(sx * 0.10, 1e-4)
  total_order <- arch_order + garch_order
  
  profiles <- list(
    list(alpha_total = 0.08, beta_total = 0.84, b = b_small,        mu = mx, gamma = 0.00),
    list(alpha_total = 0.10, beta_total = 0.82, b = 0.00,           mu = mx, gamma = 0.05),
    list(alpha_total = 0.06, beta_total = 0.88, b = 0.75 * b_small, mu = mx, gamma = 0.15),
    list(alpha_total = 0.12, beta_total = 0.75, b = 0.50 * b_small, mu = mx, gamma = -0.02)
  )
  
  if (total_order <= 4L) {
    pattern_pairs <- list(
      c("geom", "geom"),
      c("equal", "equal")
    )
  } else {
    pattern_pairs <- list(
      c("geom", "geom")
    )
  }
  
  starts <- list()
  
  for (prof in profiles) {
    for (pp in pattern_pairs) {
      alpha_pattern <- pp[1]
      beta_pattern <- pp[2]
      
      base <- .compose_start_values(
        model = "restricted",
        arch_order = arch_order,
        garch_order = garch_order,
        target_var = vx,
        alpha_total = prof$alpha_total,
        b = prof$b,
        beta_total = prof$beta_total,
        mu = prof$mu,
        gamma = prof$gamma,
        alpha_pattern = alpha_pattern,
        beta_pattern = beta_pattern
      )
      
      if (model %in% c("zero", "restricted")) {
        starts[[length(starts) + 1L]] <- .make_named_start(model, arch_order, garch_order, base)
        
      } else if (model == "free") {
        alpha_idx <- paste0("alpha", seq_len(arch_order))
        beta_idx <- paste0("beta", seq_len(garch_order))
        
        lam_restricted <- .restricted_lambda_value(
          alpha = base[alpha_idx],
          beta = base[beta_idx],
          gamma = base["gamma"],
          rho = 1
        )
        
        lambda_grid <- unique(c(0, lam_restricted, 0.25))
        lambda_grid <- lambda_grid[is.finite(lambda_grid)]
        
        for (lam in lambda_grid) {
          tmp <- c(base, lambda = lam)
          starts[[length(starts) + 1L]] <- .make_named_start(model, arch_order, garch_order, tmp)
        }
        
      } else if (model == "threshold") {
        lam1_grid <- c(0, 0.25)
        lam2_grid <- c(0, 0.10)
        
        for (lam1 in lam1_grid) {
          for (lam2 in lam2_grid) {
            tmp <- c(base, lambda1 = lam1, lambda2 = lam2)
            starts[[length(starts) + 1L]] <- .make_named_start(model, arch_order, garch_order, tmp)
          }
        }
      }
    }
  }
  
  dedup_key <- vapply(
    starts,
    function(z) paste(signif(unname(z), 8), collapse = "|"),
    character(1)
  )
  
  starts[!duplicated(dedup_key)]
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

.qgarch_nll <- function(theta,
                        x,
                        model,
                        arch_order,
                        garch_order,
                        rho = 1,
                        indicator = NULL,
                        return_details = FALSE) {
  n <- length(x)
  eta <- numeric(n)
  sigma2 <- numeric(n)
  loglik_i <- numeric(n)
  
  idx <- 1L
  omega <- theta[idx]
  idx <- idx + 1L
  
  alpha <- theta[idx:(idx + arch_order - 1L)]
  idx <- idx + arch_order
  
  b <- theta[idx]
  idx <- idx + 1L
  
  beta <- theta[idx:(idx + garch_order - 1L)]
  idx <- idx + garch_order
  
  mu <- theta[idx]
  idx <- idx + 1L
  
  gamma <- theta[idx]
  idx <- idx + 1L
  
  if (!is.finite(omega) || omega <= 0) {
    return(if (return_details) NULL else .qgarch_big_penalty)
  }
  if (any(!is.finite(alpha)) || any(alpha < 0)) {
    return(if (return_details) NULL else .qgarch_big_penalty)
  }
  if (any(!is.finite(beta)) || any(beta < 0)) {
    return(if (return_details) NULL else .qgarch_big_penalty)
  }
  if (sum(alpha) + sum(beta) >= 0.9995) {
    return(if (return_details) NULL else .qgarch_big_penalty)
  }
  
  if (model == "zero") {
    lambda_vec <- rep(0, n)
  } else if (model == "restricted") {
    lam <- .restricted_lambda_value(alpha, beta, gamma, rho = rho)
    if (!is.finite(lam)) {
      return(if (return_details) NULL else .qgarch_big_penalty)
    }
    lambda_vec <- rep(lam, n)
  } else if (model == "free") {
    lam <- theta[idx]
    if (!is.finite(lam) || abs(lam) > 50) {
      return(if (return_details) NULL else .qgarch_big_penalty)
    }
    lambda_vec <- rep(lam, n)
  } else if (model == "threshold") {
    lambda1 <- theta[idx]
    lambda2 <- theta[idx + 1L]
    if (!is.finite(lambda1) || !is.finite(lambda2) || max(abs(c(lambda1, lambda2))) > 50) {
      return(if (return_details) NULL else .qgarch_big_penalty)
    }
    lambda_vec <- lambda1 + lambda2 * indicator
  } else {
    stop("Unknown model.", call. = FALSE)
  }
  
  kappa_vec <- 1 + 2 * lambda_vec * b
  
  uncond <- .unconditional_variance(omega, alpha, b, beta)
  if (!is.finite(uncond)) {
    return(if (return_details) NULL else .qgarch_big_penalty)
  }
  
  for (i in seq_len(n)) {
    arch_term <- 0
    for (j in seq_len(arch_order)) {
      eta_lag <- if (i - j > 0L) eta[i - j] else 0
      arch_term <- arch_term + alpha[j] * (eta_lag - b)^2
    }
    
    garch_term <- 0
    for (k in seq_len(garch_order)) {
      sigma_lag <- if (i - k > 0L) sigma2[i - k] else uncond
      garch_term <- garch_term + beta[k] * sigma_lag
    }
    
    sigma2[i] <- omega + arch_term + garch_term
    
    if (model == "zero") {
      eta[i] <- x[i] - mu - gamma * sigma2[i]
    } else {
      eta[i] <- .solve_eta_quadratic(
        obs = x[i],
        sigma2 = sigma2[i],
        mu = mu,
        gamma = gamma,
        lambda = lambda_vec[i],
        kappa = kappa_vec[i]
      )
    }
    
    if (!is.finite(sigma2[i]) || sigma2[i] <= 0 || !is.finite(eta[i])) {
      return(if (return_details) NULL else .qgarch_big_penalty)
    }
    
    if (model == "zero") {
      loglik_i[i] <- -0.5 * log(2 * pi * sigma2[i]) - 0.5 * (eta[i]^2 / sigma2[i])
    } else {
      jac <- kappa_vec[i] - 2 * lambda_vec[i] * eta[i]
      if (!is.finite(jac) || jac <= 0) {
        return(if (return_details) NULL else .qgarch_big_penalty)
      }
      loglik_i[i] <- -log(jac) - 0.5 * log(2 * pi * sigma2[i]) - 0.5 * (eta[i]^2 / sigma2[i])
    }
  }
  
  nll <- -sum(loglik_i)
  
  if (return_details) {
    return(list(
      nll = nll,
      eta = eta,
      sigma2 = sigma2,
      uncond_var = uncond,
      lambda = lambda_vec,
      kappa = kappa_vec,
      loglik_contrib = loglik_i
    ))
  }
  
  nll
}

.numeric_score_matrix <- function(theta,
                                  x,
                                  model,
                                  arch_order,
                                  garch_order,
                                  rho,
                                  indicator,
                                  eps = 1e-5) {
  p <- length(theta)
  base <- .qgarch_nll(
    theta,
    x = x,
    model = model,
    arch_order = arch_order,
    garch_order = garch_order,
    rho = rho,
    indicator = indicator,
    return_details = TRUE
  )
  
  if (is.null(base)) {
    return(NULL)
  }
  
  n <- length(base$loglik_contrib)
  score <- matrix(NA_real_, nrow = n, ncol = p)
  
  for (j in seq_len(p)) {
    step <- eps * max(1, abs(theta[j]))
    
    theta_up <- theta
    theta_dn <- theta
    theta_up[j] <- theta_up[j] + step
    theta_dn[j] <- theta_dn[j] - step
    
    up <- .qgarch_nll(
      theta_up,
      x = x,
      model = model,
      arch_order = arch_order,
      garch_order = garch_order,
      rho = rho,
      indicator = indicator,
      return_details = TRUE
    )
    dn <- .qgarch_nll(
      theta_dn,
      x = x,
      model = model,
      arch_order = arch_order,
      garch_order = garch_order,
      rho = rho,
      indicator = indicator,
      return_details = TRUE
    )
    
    if (is.null(up) || is.null(dn)) {
      return(NULL)
    }
    
    score[, j] <- (up$loglik_contrib - dn$loglik_contrib) / (2 * step)
  }
  
  colnames(score) <- names(theta)
  score
}

.compute_vcov <- function(theta,
                          raw_hessian,
                          x,
                          model,
                          arch_order,
                          garch_order,
                          rho,
                          indicator,
                          vcov_type = c("auto", "sandwich", "hessian", "opg", "none")) {
  vcov_type <- match.arg(vcov_type)
  est_names <- names(theta)
  
  empty <- list(
    vcov = matrix(
      NA_real_,
      nrow = length(theta),
      ncol = length(theta),
      dimnames = list(est_names, est_names)
    ),
    se = stats::setNames(rep(NA_real_, length(theta)), est_names),
    method = "none"
  )
  
  if (vcov_type == "none") {
    return(empty)
  }
  
  make_result <- function(vcov_mat, method) {
    se_info <- .make_se_from_vcov(vcov_mat)
    if (is.null(se_info$se) || !any(is.finite(se_info$se))) {
      return(NULL)
    }
    
    names(se_info$se) <- est_names
    rownames(se_info$vcov) <- est_names
    colnames(se_info$vcov) <- est_names
    
    list(vcov = se_info$vcov, se = se_info$se, method = method)
  }
  
  H_inv <- NULL
  
  if (vcov_type %in% c("auto", "sandwich", "hessian")) {
    H_raw <- .clean_sym_matrix(raw_hessian)
    
    H_use <- NULL
    if (!is.null(H_raw) && all(is.finite(H_raw))) {
      H_use <- H_raw
    } else {
      H_num <- tryCatch(
        stats::optimHess(
          par = unname(theta),
          fn = .qgarch_nll,
          x = x,
          model = model,
          arch_order = arch_order,
          garch_order = garch_order,
          rho = rho,
          indicator = indicator
        ),
        error = function(e) NULL,
        warning = function(w) NULL
      )
      H_num <- .clean_sym_matrix(H_num)
      
      if (!is.null(H_num) && all(is.finite(H_num))) {
        H_use <- H_num
      }
    }
    
    if (!is.null(H_use)) {
      H_reg <- .regularize_pd(H_use)
      H_inv <- .clean_sym_matrix(.safe_invert(H_reg))
    }
    
    if (!is.null(H_inv) && all(is.finite(H_inv))) {
      if (vcov_type == "hessian") {
        out <- make_result(H_inv, "hessian")
        return(if (is.null(out)) empty else out)
      }
      
      if (vcov_type == "auto") {
        out <- make_result(H_inv, "hessian")
        if (!is.null(out)) {
          return(out)
        }
      }
    }
  }
  
  if (vcov_type %in% c("auto", "sandwich", "opg")) {
    score_mat <- .numeric_score_matrix(
      theta = unname(theta),
      x = x,
      model = model,
      arch_order = arch_order,
      garch_order = garch_order,
      rho = rho,
      indicator = indicator
    )
    
    S <- NULL
    if (!is.null(score_mat) && all(is.finite(score_mat))) {
      S <- .regularize_pd(crossprod(score_mat))
    }
    
    if (vcov_type == "opg") {
      if (!is.null(S) && all(is.finite(S))) {
        S_inv <- .clean_sym_matrix(.safe_invert(S))
        out <- make_result(S_inv, "opg")
        return(if (is.null(out)) empty else out)
      }
      return(empty)
    }
    
    if (vcov_type == "sandwich") {
      if (!is.null(H_inv) && !is.null(S) && all(is.finite(H_inv)) && all(is.finite(S))) {
        out <- make_result(.clean_sym_matrix(H_inv %*% S %*% H_inv), "sandwich")
        return(if (is.null(out)) empty else out)
      }
      return(empty)
    }
    
    if (vcov_type == "auto") {
      if (!is.null(H_inv) && !is.null(S) && all(is.finite(H_inv)) && all(is.finite(S))) {
        out <- make_result(.clean_sym_matrix(H_inv %*% S %*% H_inv), "sandwich")
        if (!is.null(out)) {
          return(out)
        }
      }
      
      if (!is.null(S) && all(is.finite(S))) {
        S_inv <- .clean_sym_matrix(.safe_invert(S))
        out <- make_result(S_inv, "opg")
        if (!is.null(out)) {
          return(out)
        }
      }
    }
  }
  
  empty
}

.restricted_lambda_delta_se <- function(alpha, beta, gamma, vcov_sub, rho = 1) {
  alpha_sum <- sum(alpha)
  beta_sum <- sum(beta)
  den <- 1 - rho * (alpha_sum + beta_sum)
  
  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  
  grad <- c(
    gamma = rho * alpha_sum / den,
    stats::setNames(
      rep(rho * gamma * (1 - rho * beta_sum) / (den^2), length(alpha)),
      names(alpha)
    ),
    stats::setNames(
      rep(rho^2 * gamma * alpha_sum / (den^2), length(beta)),
      names(beta)
    )
  )
  
  grad <- grad[rownames(vcov_sub)]
  se2 <- as.numeric(t(grad) %*% vcov_sub %*% grad)
  
  if (!is.finite(se2) || se2 <= 0) {
    return(NA_real_)
  }
  
  sqrt(se2)
}

.deduplicate_start_list <- function(starts) {
  if (length(starts) <= 1L) {
    return(starts)
  }
  
  key <- vapply(
    starts,
    function(z) paste(signif(unname(z), 8), collapse = "|"),
    character(1)
  )
  
  starts[!duplicated(key)]
}

.prepare_start_list <- function(start, x, model, arch_order, garch_order) {
  defaults <- qgarch_default_starts(
    x,
    model = model,
    arch_order = arch_order,
    garch_order = garch_order
  )
  
  if (is.null(start)) {
    return(defaults)
  }
  
  user_starts <- NULL
  
  if (is.numeric(start)) {
    user_starts <- list(.make_named_start(model, arch_order, garch_order, start))
  } else if (is.list(start) && length(start) > 0L && all(vapply(start, is.numeric, logical(1)))) {
    user_starts <- lapply(start, function(z) .make_named_start(model, arch_order, garch_order, z))
  } else {
    stop("'start' must be NULL, a numeric vector, or a list of numeric vectors.", call. = FALSE)
  }
  
  .deduplicate_start_list(c(user_starts, defaults))
}

.is_better_fit <- function(candidate, incumbent) {
  if (is.null(incumbent)) {
    return(TRUE)
  }
  
  cand_min <- candidate$minimum
  inc_min <- incumbent$minimum
  
  if (cand_min < inc_min - 1e-8) {
    return(TRUE)
  }
  if (abs(cand_min - inc_min) <= 1e-8 && candidate$code < incumbent$code) {
    return(TRUE)
  }
  
  FALSE
}

.refine_fit <- function(fit_obj,
                        x,
                        model,
                        arch_order,
                        garch_order,
                        rho,
                        indicator,
                        steptol,
                        iterlim,
                        print.level,
                        hessian) {
  p0 <- unname(fit_obj$estimate)
  typsize2 <- pmax(abs(p0), 0.05)
  
  refit <- try(
    stats::nlm(
      f = .qgarch_nll,
      p = p0,
      x = x,
      model = model,
      arch_order = arch_order,
      garch_order = garch_order,
      rho = rho,
      indicator = indicator,
      hessian = hessian,
      steptol = steptol / 10,
      typsize = typsize2,
      iterlim = max(iterlim, 500L),
      print.level = print.level
    ),
    silent = TRUE
  )
  
  if (inherits(refit, "try-error") || !is.finite(refit$minimum) || refit$minimum >= .qgarch_big_penalty) {
    return(fit_obj)
  }
  
  names(refit$estimate) <- names(fit_obj$estimate)
  
  if (.is_better_fit(refit, fit_obj)) {
    return(refit)
  }
  
  fit_obj
}

.build_fit_object <- function(best_fit,
                              x,
                              model,
                              arch_order,
                              garch_order,
                              rho,
                              indicator,
                              starts,
                              call,
                              vcov_type) {
  details <- .qgarch_nll(
    best_fit$estimate,
    x = x,
    model = model,
    arch_order = arch_order,
    garch_order = garch_order,
    rho = rho,
    indicator = indicator,
    return_details = TRUE
  )
  
  est_names <- names(best_fit$estimate)
  
  vcov_info <- .compute_vcov(
    theta = best_fit$estimate,
    raw_hessian = best_fit$hessian,
    x = x,
    model = model,
    arch_order = arch_order,
    garch_order = garch_order,
    rho = rho,
    indicator = indicator,
    vcov_type = vcov_type
  )
  
  vcov_mat <- vcov_info$vcov
  se_vec <- vcov_info$se
  
  names(se_vec) <- est_names
  rownames(vcov_mat) <- est_names
  colnames(vcov_mat) <- est_names
  
  coef_full <- best_fit$estimate
  full_se <- se_vec
  
  if (model == "zero") {
    coef_full <- c(coef_full, lambda = 0)
    full_se <- c(full_se, lambda = NA_real_)
  }
  
  implied_lambda <- NULL
  implied_lambda_se <- NULL
  
  if (model == "restricted") {
    alpha_idx <- paste0("alpha", seq_len(arch_order))
    beta_idx <- paste0("beta", seq_len(garch_order))
    
    alpha <- best_fit$estimate[alpha_idx]
    beta <- best_fit$estimate[beta_idx]
    gamma <- best_fit$estimate["gamma"]
    
    implied_lambda <- .restricted_lambda_value(alpha, beta, gamma, rho = rho)
    
    needed <- c("gamma", alpha_idx, beta_idx)
    if (all(needed %in% rownames(vcov_mat))) {
      implied_lambda_se <- .restricted_lambda_delta_se(
        alpha = alpha,
        beta = beta,
        gamma = gamma,
        vcov_sub = vcov_mat[needed, needed, drop = FALSE],
        rho = rho
      )
    } else {
      implied_lambda_se <- NA_real_
    }
    
    coef_full <- c(coef_full, lambda = implied_lambda)
    full_se <- c(full_se, lambda = implied_lambda_se)
  }
  
  structure(
    list(
      call = call,
      model = model,
      arch_order = arch_order,
      garch_order = garch_order,
      coefficients = best_fit$estimate,
      coefficients_full = coef_full,
      standard_errors = se_vec,
      standard_errors_full = full_se,
      vcov = vcov_mat,
      vcov_type = vcov_info$method,
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
      hessian = best_fit$hessian,
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

#' Fit generalized qgarch(m, n) models
#'
#' Fits QGARCH-in-mean models using nonlinear minimization of the negative
#' log-likelihood. Four variants are supported: a zero-lambda model, a
#' restricted-lambda model, a free-lambda model, and a threshold model with
#' state-dependent lambda.
#'
#' @param x A numeric vector or one-column data frame containing the observed
#'   series. The series must contain only finite, non-missing values.
#' @param model Character string specifying the model variant. Must be one of
#'   `"zero"`, `"restricted"`, `"free"`, or `"threshold"`.
#' @param arch_order Positive integer giving the ARCH lag order `m`.
#' @param garch_order Positive integer giving the GARCH lag order `n`.
#' @param threshold_indicator Optional threshold indicator used only when
#'   `model = "threshold"`. Must have the same length as `x`. Nonzero values
#'   are converted to 1 and zero values to 0.
#' @param start Optional starting values. May be `NULL`, a numeric vector, or a
#'   list of numeric vectors.
#' @param rho Scalar used in the restricted-lambda mapping.
#' @param steptol Step tolerance passed to [stats::nlm()].
#' @param typsize Typical size passed to [stats::nlm()]. A single value is
#'   repeated to the appropriate length.
#' @param iterlim Maximum number of iterations passed to [stats::nlm()].
#' @param print.level Print level passed to [stats::nlm()].
#' @param hessian Logical; should the Hessian be returned by [stats::nlm()]?
#' @param vcov_type Character string controlling standard-error estimation.
#'   One of `"auto"`, `"sandwich"`, `"hessian"`, `"opg"`, or `"none"`.
#'   In `"auto"` mode, Hessian-based standard errors are tried first.
#'
#' @return An object of class `"qgarch"`.
#' @export
qgarch_fit <- function(x,
                       model = c("zero", "restricted", "free", "threshold"),
                       arch_order = 1L,
                       garch_order = 1L,
                       threshold_indicator = NULL,
                       start = NULL,
                       rho = 1,
                       steptol = 1e-10,
                       typsize = 0.1,
                       iterlim = 300L,
                       print.level = 0,
                       hessian = TRUE,
                       vcov_type = c("auto", "sandwich", "hessian", "opg", "none")) {
  model <- match.arg(model)
  vcov_type <- match.arg(vcov_type)
  
  orders <- .validate_orders(arch_order, garch_order)
  arch_order <- orders$arch_order
  garch_order <- orders$garch_order
  
  x <- .validate_input_series(x)
  
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho)) {
    stop("'rho' must be a single finite numeric value.", call. = FALSE)
  }
  
  if (!is.numeric(typsize) || length(typsize) != 1L || !is.finite(typsize) || typsize <= 0) {
    stop("'typsize' must be a single positive finite numeric value.", call. = FALSE)
  }
  
  if (!is.numeric(steptol) || length(steptol) != 1L || !is.finite(steptol) || steptol <= 0) {
    stop("'steptol' must be a single positive finite numeric value.", call. = FALSE)
  }
  
  if (!is.numeric(iterlim) || length(iterlim) != 1L || !is.finite(iterlim) || iterlim < 1) {
    stop("'iterlim' must be a single positive integer.", call. = FALSE)
  }
  iterlim <- as.integer(iterlim)
  
  if (!is.numeric(print.level) || length(print.level) != 1L || !is.finite(print.level)) {
    stop("'print.level' must be a single finite numeric value.", call. = FALSE)
  }
  
  if (!is.logical(hessian) || length(hessian) != 1L || is.na(hessian)) {
    stop("'hessian' must be either TRUE or FALSE.", call. = FALSE)
  }
  
  if (model == "threshold") {
    if (is.null(threshold_indicator)) {
      stop("'threshold_indicator' must be supplied when model = 'threshold'.", call. = FALSE)
    }
    
    indicator <- as.numeric(threshold_indicator)
    
    if (length(indicator) != length(x)) {
      stop("'threshold_indicator' must have the same length as 'x'.", call. = FALSE)
    }
    
    if (anyNA(indicator) || any(!is.finite(indicator))) {
      stop("'threshold_indicator' must contain only finite, non-missing values.", call. = FALSE)
    }
    
    indicator <- ifelse(indicator != 0, 1, 0)
  } else {
    indicator <- NULL
  }
  
  starts <- .prepare_start_list(
    start = start,
    x = x,
    model = model,
    arch_order = arch_order,
    garch_order = garch_order
  )
  
  best_fit <- NULL
  
  for (i in seq_along(starts)) {
    fit_i <- try(
      stats::nlm(
        f = .qgarch_nll,
        p = unname(starts[[i]]),
        x = x,
        model = model,
        arch_order = arch_order,
        garch_order = garch_order,
        rho = rho,
        indicator = indicator,
        hessian = hessian,
        steptol = steptol,
        typsize = rep(typsize, length(starts[[i]])),
        iterlim = iterlim,
        print.level = print.level
      ),
      silent = TRUE
    )
    
    if (!inherits(fit_i, "try-error") &&
        is.finite(fit_i$minimum) &&
        fit_i$minimum < .qgarch_big_penalty) {
      names(fit_i$estimate) <- names(starts[[i]])
      if (.is_better_fit(fit_i, best_fit)) {
        best_fit <- fit_i
      }
    }
  }
  
  if (is.null(best_fit)) {
    stop("No converged solution was found. Try supplying custom starting values.", call. = FALSE)
  }
  
  best_fit <- .refine_fit(
    fit_obj = best_fit,
    x = x,
    model = model,
    arch_order = arch_order,
    garch_order = garch_order,
    rho = rho,
    indicator = indicator,
    steptol = steptol,
    iterlim = iterlim,
    print.level = print.level,
    hessian = hessian
  )
  
  .build_fit_object(
    best_fit = best_fit,
    x = x,
    model = model,
    arch_order = arch_order,
    garch_order = garch_order,
    rho = rho,
    indicator = indicator,
    starts = starts,
    call = match.call(),
    vcov_type = vcov_type
  )
}