source("simulation/allocation_constraints.R")

HISTORICAL_RAR_METHODS <- c(
  "CR", "Full-oracle-target", "Safe-Y-only", "Safe-Naive-Z",
  "H-RAR", "Known-RM-RAR"
)

historical_revision_counts <- function(config = NULL) {
  if (!is.null(config) && !is.null(config$historical_revision_counts)) {
    return(config$historical_revision_counts)
  }

  data.frame(
    arm = 0:1,
    se_success = c(71, 144), se_failure = c(19, 6),
    sp_success = c(80, 25), sp_failure = c(9, 5)
  )
}

clamp_unit <- function(x, tolerance = 1e-10) {
  pmin(1 - tolerance, pmax(tolerance, x))
}

arm_observed_loglik <- function(theta, sensitivity, specificity, history,
                                arm) {
  idx <- history$A == arm & history$Q == 1
  if (!any(idx)) return(0)
  z <- history$QZ[idx]
  reviewed <- history$R[idx] == 1
  y <- history$RY[idx]

  q <- clamp_unit(
    sensitivity * theta + (1 - specificity) * (1 - theta)
  )
  loglik <- sum(ifelse(z[!reviewed] == 1, log(q), log1p(-q)))

  if (any(reviewed)) {
    yr <- y[reviewed]
    zr <- z[reviewed]
    loglik <- loglik + sum(ifelse(yr == 1, log(theta), log1p(-theta)))
    loglik <- loglik + sum(ifelse(
      yr == 1,
      ifelse(zr == 1, log(sensitivity), log1p(-sensitivity)),
      ifelse(zr == 0, log(specificity), log1p(-specificity))
    ))
  }
  loglik
}

fit_historical_revision_arm <- function(history, arm,
                                        known_revision = FALSE,
                                        config = canonical_config()) {
  counts <- historical_revision_counts(config)
  hc <- counts[counts$arm == arm, , drop = FALSE]
  start_se <- (hc$se_success + 0.5) /
    (hc$se_success + hc$se_failure + 1)
  start_sp <- (hc$sp_success + 0.5) /
    (hc$sp_success + hc$sp_failure + 1)
  idx_y <- history$A == arm & history$R == 1
  start_theta <- smoothed_mean(history$RY[idx_y])

  if (known_revision) {
    sensitivity <- config$sensitivity[as.character(arm)]
    specificity <- config$specificity[as.character(arm)]
    objective <- function(eta) {
      theta <- clamp_unit(stats::plogis(eta), tolerance = 1e-8)
      -arm_observed_loglik(theta, sensitivity, specificity, history, arm)
    }
    fit <- stats::optim(
      stats::qlogis(start_theta), objective,
      method = "BFGS", hessian = TRUE,
      control = list(reltol = 1e-10, maxit = 200)
    )
    theta <- stats::plogis(fit$par)
    var_eta <- if (is.finite(fit$hessian[1, 1]) && fit$hessian[1, 1] > 0) {
      1 / fit$hessian[1, 1]
    } else Inf
    return(c(
      theta = theta,
      variance = theta^2 * (1 - theta)^2 * var_eta,
      sensitivity = unname(sensitivity),
      specificity = unname(specificity),
      converged = as.integer(fit$convergence == 0)
    ))
  }

  objective <- function(eta) {
    parameter <- clamp_unit(stats::plogis(eta), tolerance = 1e-8)
    theta <- parameter[1]
    sensitivity <- parameter[2]
    specificity <- parameter[3]
    loglik <- arm_observed_loglik(
      theta, sensitivity, specificity, history, arm
    )
    loglik <- loglik +
      hc$se_success * log(sensitivity) +
      hc$se_failure * log1p(-sensitivity) +
      hc$sp_success * log(specificity) +
      hc$sp_failure * log1p(-specificity)
    -loglik
  }
  fit <- stats::optim(
    stats::qlogis(c(start_theta, start_se, start_sp)), objective,
    method = "BFGS", hessian = TRUE,
    control = list(reltol = 1e-10, maxit = 300)
  )
  parameter <- stats::plogis(fit$par)
  covariance_eta <- tryCatch(
    solve(fit$hessian),
    error = function(e) matrix(Inf, 3, 3)
  )
  gradient_theta <- c(parameter[1] * (1 - parameter[1]), 0, 0)
  variance <- as.numeric(t(gradient_theta) %*% covariance_eta %*%
                           gradient_theta)
  if (!is.finite(variance) || variance <= 0) variance <- Inf
  c(
    theta = parameter[1], variance = variance,
    sensitivity = parameter[2], specificity = parameter[3],
    converged = as.integer(fit$convergence == 0)
  )
}

historical_revision_effect_z <- function(method, history, config) {
  known <- method == "Known-RM-RAR"
  fits <- lapply(0:1, function(a) {
    fit_historical_revision_arm(history, a, known_revision = known,
                                config = config)
  })
  variance <- fits[[1]]["variance"] + fits[[2]]["variance"]
  if (!is.finite(variance) || variance <= 0) return(NA_real_)
  unname((fits[[2]]["theta"] - fits[[1]]["theta"]) / sqrt(variance))
}

historical_revision_means <- function(history, config) {
  fits <- lapply(0:1, function(a) {
    fit_historical_revision_arm(
      history, a, known_revision = FALSE, config = config
    )
  })
  c(`0` = unname(fits[[1]]["theta"]),
    `1` = unname(fits[[2]]["theta"]))
}

wilson_interval <- function(x, n, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  phat <- x / n
  denominator <- 1 + z^2 / n
  centre <- (phat + z^2 / (2 * n)) / denominator
  half_width <- z / denominator * sqrt(
    phat * (1 - phat) / n + z^2 / (4 * n^2)
  )
  c(lower = centre - half_width, upper = centre + half_width)
}

analyse_binomial_final <- function(dat, config) {
  idx1 <- dat$A == 1
  idx0 <- dat$A == 0
  n1 <- sum(idx1)
  n0 <- sum(idx0)
  x1 <- sum(dat$Y[idx1])
  x0 <- sum(dat$Y[idx0])
  p1 <- x1 / n1
  p0 <- x0 / n0
  estimate <- p1 - p0
  se <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
  pooled <- (x1 + x0) / (n1 + n0)
  score_se <- sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n0))
  interval1 <- wilson_interval(x1, n1)
  interval0 <- wilson_interval(x0, n0)
  newcombe <- c(
    lower = unname(interval1["lower"] - interval0["upper"]),
    upper = unname(interval1["upper"] - interval0["lower"])
  )
  truth <- unname(diff(config$p_y))
  c(
    binomial_estimate = estimate,
    binomial_se = se,
    binomial_reject = as.integer(estimate / score_se >
                                   config$final_critical_value),
    binomial_cover = as.integer(abs(estimate - truth) <=
                                  stats::qnorm(0.975) * se),
    newcombe_cover = as.integer(
      truth >= newcombe["lower"] && truth <= newcombe["upper"]
    )
  )
}

next_historical_rar_probability <- function(method, history, month, config,
                                            gate, gamma = 2) {
  if (method == "CR") {
    return(list(pi = 0.5, certified_direction = 0L))
  }
  if (method == "Full-oracle-target") {
    return(list(
      pi = unname(oracle_welfare_target(config)["allocation_1"]),
      certified_direction = 0L
    ))
  }
  if (month <= gate$burn_months || nrow(history) == 0) {
    return(list(pi = 0.5, certified_direction = 0L))
  }

  if (method == "H-DBCD") {
    means <- historical_revision_means(history, config)
    target <- sqrt(clamp_unit(means["1"])) /
      (sqrt(clamp_unit(means["1"])) + sqrt(clamp_unit(means["0"])))
    return(list(
      pi = dbcd_update(
        mean(history$A), target, gamma,
        config$power_envelope_lower, config$power_envelope_upper
      ),
      certified_direction = 0L
    ))
  }

  z <- if (method %in% c("Safe-Y-only", "Safe-Naive-Z")) {
    safe_effect_z(method, history, config)
  } else {
    historical_revision_effect_z(method, history, config)
  }
  direction <- if (!is.finite(z)) {
    0L
  } else if (z > gate$upper_z_boundary) {
    1L
  } else if (z < -gate$lower_z_boundary) {
    -1L
  } else {
    0L
  }
  target <- if (direction == 1L) {
    config$power_envelope_upper
  } else if (direction == -1L) {
    config$power_envelope_lower
  } else {
    0.5
  }
  list(
    pi = dbcd_update(
      mean(history$A), target, gamma,
      config$power_envelope_lower, config$power_envelope_upper
    ),
    certified_direction = direction
  )
}

simulate_historical_rar_method <- function(method, stream, config, gate,
                                           assignment_seed) {
  set.seed(assignment_seed)
  n <- config$n
  dat <- data.frame(
    id = seq_len(n), entry_month = stream$entry_month,
    A = NA_integer_, pi = NA_real_, Q = NA_integer_, QZ = NA_integer_,
    R = NA_integer_, RY = NA_integer_, Y = NA_integer_
  )
  monthly_direction <- integer(config$n_months)

  for (month in seq_len(config$n_months)) {
    current <- which(dat$entry_month == month)
    previous <- which(dat$entry_month < month & is.finite(dat$A))
    if (length(previous)) {
      info <- canonical_information_at(stream, month, dat$A)
      dat$Q[previous] <- info$Q[previous]
      dat$QZ[previous] <- info$QZ[previous]
      dat$R[previous] <- info$R[previous]
      dat$RY[previous] <- info$RY[previous]
    }
    history <- dat[previous, , drop = FALSE]
    decision <- next_historical_rar_probability(
      method, history, month, config, gate
    )
    monthly_direction[month] <- decision$certified_direction
    dat$pi[current] <- decision$pi
    dat$A[current] <- stats::rbinom(length(current), 1, decision$pi)
    ids <- dat$id[current]
    dat$Y[current] <- ifelse(
      dat$A[current] == 1, stream$y1[ids], stream$y0[ids]
    )
  }

  analysed <- analyse_oracle_trial(dat, config)
  binomial <- analyse_binomial_final(dat, config)
  for (label in names(binomial)) analysed[[label]] <- unname(binomial[label])
  analysed$certified_upper <- as.integer(any(monthly_direction == 1L))
  analysed$certified_lower <- as.integer(any(monthly_direction == -1L))
  analysed$first_certified_upper_month <- if (any(monthly_direction == 1L)) {
    which(monthly_direction == 1L)[1]
  } else NA_integer_
  analysed$first_certified_lower_month <- if (any(monthly_direction == -1L)) {
    which(monthly_direction == -1L)[1]
  } else NA_integer_
  cbind(method = method, analysed, stringsAsFactors = FALSE)
}

run_historical_revision_rar_check <- function(n_rep = 500L,
                                              seed = 20261114L,
                                              config = canonical_config(),
                                              eta = 1,
                                              methods = HISTORICAL_RAR_METHODS) {
  gate <- safe_gate_config(config, eta)
  rows <- vector("list", n_rep * length(methods))
  k <- 1L
  for (replicate in seq_len(n_rep)) {
    stream <- generate_canonical_stream(config, seed + 10000L * replicate)
    for (j in seq_along(methods)) {
      rows[[k]] <- simulate_historical_rar_method(
        methods[j], stream, config, gate,
        assignment_seed = seed + 10000L * replicate + j
      )
      rows[[k]]$replicate <- replicate
      k <- k + 1L
    }
  }
  list(gate = gate, results = do.call(rbind, rows))
}

summarise_historical_revision_rar_check <- function(results,
                                                    config = canonical_config()) {
  method_order <- HISTORICAL_RAR_METHODS[
    HISTORICAL_RAR_METHODS %in% unique(results$method)
  ]
  ans <- summarise_oracle_check(results, config, method_order = method_order)
  certification <- do.call(rbind, lapply(split(results, results$method), function(x) {
    data.frame(
      method = x$method[1],
      probability_certified_upper = mean(x$certified_upper),
      probability_certified_lower = mean(x$certified_lower),
      median_first_certified_upper_month =
        stats::median(x$first_certified_upper_month, na.rm = TRUE),
      median_first_certified_lower_month =
        stats::median(x$first_certified_lower_month, na.rm = TRUE)
    )
  }))
  rownames(certification) <- NULL
  ans <- merge(ans, certification, by = "method", all.x = TRUE, sort = FALSE)
  binomial <- do.call(rbind, lapply(split(results, results$method), function(x) {
    data.frame(
      method = x$method[1],
      binomial_bias = mean(x$binomial_estimate) - diff(config$p_y),
      binomial_empirical_se = stats::sd(x$binomial_estimate),
      binomial_mean_se = mean(x$binomial_se),
      binomial_coverage = mean(x$binomial_cover),
      newcombe_coverage = mean(x$newcombe_cover),
      binomial_power = mean(x$binomial_reject)
    )
  }))
  rownames(binomial) <- NULL
  ans <- merge(ans, binomial, by = "method", all.x = TRUE, sort = FALSE)
  ans[match(method_order, ans$method), ]
}
