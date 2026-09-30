.arm_revision_loglik <- function(theta, sensitivity, specificity,
                                 history, arm) {
  idx <- history$A == arm & history$Q == 1L
  if (!any(idx)) return(0)
  z <- history$QZ[idx]
  reviewed <- history$R[idx] == 1L
  y <- history$RY[idx]
  q <- .clip_probability(
    sensitivity * theta + (1 - specificity) * (1 - theta)
  )
  value <- sum(ifelse(z[!reviewed] == 1L, log(q), log1p(-q)))
  if (any(reviewed)) {
    yr <- y[reviewed]
    zr <- z[reviewed]
    value <- value + sum(ifelse(yr == 1L, log(theta), log1p(-theta)))
    value <- value + sum(ifelse(
      yr == 1L,
      ifelse(zr == 1L, log(sensitivity), log1p(-sensitivity)),
      ifelse(zr == 0L, log(specificity), log1p(-specificity))
    ))
  }
  value
}

.fit_revision_arm <- function(history, arm, validation_counts) {
  hc <- validation_counts
  start_se <- (hc$se_success + 0.5) /
    (hc$se_success + hc$se_failure + 1)
  start_sp <- (hc$sp_success + 0.5) /
    (hc$sp_success + hc$sp_failure + 1)
  returned <- history$A == arm & history$R == 1L
  start_theta <- (sum(history$RY[returned], na.rm = TRUE) + 0.5) /
    (sum(returned) + 1)
  objective <- function(eta) {
    parameter <- .clip_probability(stats::plogis(eta), 1e-8, 1 - 1e-8)
    theta <- parameter[1L]
    sensitivity <- parameter[2L]
    specificity <- parameter[3L]
    loglik <- .arm_revision_loglik(
      theta, sensitivity, specificity, history, arm
    ) +
      hc$se_success * log(sensitivity) +
      hc$se_failure * log1p(-sensitivity) +
      hc$sp_success * log(specificity) +
      hc$sp_failure * log1p(-specificity)
    -loglik
  }
  fit <- stats::optim(
    stats::qlogis(c(start_theta, start_se, start_sp)),
    objective, method = "BFGS", hessian = TRUE,
    control = list(reltol = 1e-10, maxit = 300)
  )
  parameter <- stats::plogis(fit$par)
  covariance_eta <- tryCatch(
    solve(fit$hessian),
    error = function(e) matrix(Inf, 3L, 3L)
  )
  gradient_theta <- c(parameter[1L] * (1 - parameter[1L]), 0, 0)
  variance_theta <- as.numeric(
    t(gradient_theta) %*% covariance_eta %*% gradient_theta
  )
  if (!is.finite(variance_theta) || variance_theta <= 0) {
    variance_theta <- Inf
  }
  c(
    theta = parameter[1L],
    variance = variance_theta,
    sensitivity = parameter[2L],
    specificity = parameter[3L],
    converged = as.integer(fit$convergence == 0L)
  )
}

.revision_effect_summary <- function(history, validation) {
  if (!nrow(history) || !all(0:1 %in% history$A)) {
    return(list(theta = c(`0` = NA_real_, `1` = NA_real_),
                variance = c(`0` = NA_real_, `1` = NA_real_),
                z = NA_real_, fits = NULL))
  }
  has_provisional <- vapply(0:1, function(arm) {
    any(history$A == arm & history$Q == 1L)
  }, logical(1))
  if (!all(has_provisional)) {
    return(list(theta = c(`0` = NA_real_, `1` = NA_real_),
                variance = c(`0` = NA_real_, `1` = NA_real_),
                z = NA_real_, fits = NULL))
  }
  fits <- lapply(0:1, function(arm) {
    .fit_revision_arm(
      history, arm, validation[validation$arm == arm, , drop = FALSE]
    )
  })
  theta <- c(`0` = unname(fits[[1L]]["theta"]),
             `1` = unname(fits[[2L]]["theta"]))
  variance <- c(`0` = unname(fits[[1L]]["variance"]),
                `1` = unname(fits[[2L]]["variance"]))
  total_variance <- sum(variance)
  z <- if (is.finite(total_variance) && total_variance > 0) {
    unname((theta["1"] - theta["0"]) / sqrt(total_variance))
  } else {
    NA_real_
  }
  list(theta = theta, variance = variance, z = z, fits = fits)
}

