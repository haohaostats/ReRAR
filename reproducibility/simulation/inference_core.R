source("simulation/canonical_dgm.R")

ORACLE_METHODS <- c(
  "CR", "Full-oracle-target", "Y-only", "Naive-Z", "Revision-oracle"
)

clip_probability <- function(x, lower, upper) pmin(upper, pmax(lower, x))

smoothed_mean <- function(x, prior = 0.5) {
  x <- x[is.finite(x)]
  (sum(x) + prior) / (length(x) + 2 * prior)
}

dbcd_update <- function(current, target, gamma, lower, upper) {
  current <- clip_probability(current, 0.02, 0.98)
  target <- clip_probability(target, lower, upper)
  numerator <- target * (target / current)^gamma
  denominator <- numerator + (1 - target) * ((1 - target) / (1 - current))^gamma
  clip_probability(numerator / denominator, lower, upper)
}

constrained_welfare_target <- function(p0, p1, config) {
  if (p1 > p0) return(config$power_envelope_upper)
  if (p1 < p0) return(config$power_envelope_lower)
  0.5
}

true_revision_regression <- function(config) {
  ans <- matrix(NA_real_, 2, 2,
                dimnames = list(A = c("0", "1"), Z = c("0", "1")))
  for (a in 0:1) {
    p <- config$p_y[as.character(a)]
    se <- config$sensitivity[as.character(a)]
    sp <- config$specificity[as.character(a)]
    q1 <- se * p + (1 - sp) * (1 - p)
    ans[as.character(a), "1"] <- se * p / q1
    ans[as.character(a), "0"] <- (1 - se) * p / (1 - q1)
  }
  ans
}

estimate_online_means <- function(method, history, config) {
  p <- c(0.5, 0.5)
  m_true <- true_revision_regression(config)
  for (a in 0:1) {
    idx_a <- history$A == a
    if (method == "Y-only") {
      signal <- history$RY[idx_a]
    } else if (method == "Naive-Z") {
      signal <- history$QZ[idx_a]
    } else if (method == "Revision-oracle") {
      q <- history$Q[idx_a] == 1
      r <- history$R[idx_a] == 1
      z <- history$QZ[idx_a]
      signal <- rep(NA_real_, sum(idx_a))
      signal[q] <- m_true[cbind(as.character(rep(a, sum(q))), as.character(z[q]))]
      signal[r] <- history$RY[idx_a][r]
    } else {
      stop("Unknown online estimator: ", method)
    }
    p[a + 1] <- smoothed_mean(signal)
  }
  p
}

next_oracle_probability <- function(method, history, month, config,
                                    burn_months = 5L, gamma = 2) {
  if (method == "CR") return(0.5)
  if (method == "Full-oracle-target") {
    return(unname(oracle_welfare_target(config)["allocation_1"]))
  }
  if (month <= burn_months || nrow(history) == 0) return(0.5)
  p <- estimate_online_means(method, history, config)
  target <- constrained_welfare_target(p[1], p[2], config)
  dbcd_update(mean(history$A), target, gamma,
              config$power_envelope_lower, config$power_envelope_upper)
}

analyse_oracle_trial <- function(dat, config) {
  w1 <- dat$A / dat$pi
  w0 <- (1 - dat$A) / (1 - dat$pi)
  theta1 <- sum(w1 * dat$Y) / sum(w1)
  theta0 <- sum(w0 * dat$Y) / sum(w0)
  estimate <- theta1 - theta0




  influence <- w1 * (dat$Y - theta1) / mean(w1) -
    w0 * (dat$Y - theta0) / mean(w0)
  se <- stats::sd(influence) / sqrt(nrow(dat))
  critical <- config$final_critical_value
  pi_by_month <- split(dat$pi, dat$entry_month)
  monthly_pi <- vapply(pi_by_month, mean, numeric(1))
  upper_hit <- vapply(
    pi_by_month,
    function(x) any(x >= config$power_envelope_upper - 1e-10),
    logical(1)
  )
  lower_hit <- vapply(
    pi_by_month,
    function(x) any(x <= config$power_envelope_lower + 1e-10),
    logical(1)
  )
  first_upper_month <- if (any(upper_hit)) which(upper_hit)[1] else NA_integer_
  first_lower_month <- if (any(lower_hit)) which(lower_hit)[1] else NA_integer_





  working0 <- unname(config$design_p_y["0"])
  working1 <- unname(config$design_p_y["1"])
  aipw_score <- (working1 - working0) +
    w1 * (dat$Y - working1) - w0 * (dat$Y - working0)
  aipw_estimate <- mean(aipw_score)
  aipw_se <- stats::sd(aipw_score) / sqrt(nrow(dat))
  truth <- unname(diff(config$p_y))





  stabilized_arm <- function(arm, working) {
    indicator <- as.integer(dat$A == arm)
    propensity <- if (arm == 1L) dat$pi else 1 - dat$pi
    h <- sqrt(propensity)
    correction <- numeric(length(indicator))
    observed <- indicator == 1L
    correction[observed] <-
      (dat$Y[observed] - working) / propensity[observed]
    score <- working + correction
    estimate_arm <- sum(h * score) / sum(h)
    variance_arm <- sum(h^2 * (score - estimate_arm)^2) / sum(h)^2
    c(estimate = estimate_arm, variance = variance_arm)
  }
  stabilized0 <- stabilized_arm(0L, working0)
  stabilized1 <- stabilized_arm(1L, working1)
  stabilized_estimate <-
    unname(stabilized1["estimate"] - stabilized0["estimate"])
  stabilized_se <- sqrt(
    unname(stabilized1["variance"] + stabilized0["variance"])
  )





  contrast_h <- sqrt(dat$pi * (1 - dat$pi))
  contrast_stabilized_estimate <-
    sum(contrast_h * aipw_score) / sum(contrast_h)
  contrast_stabilized_variance <-
    sum(contrast_h^2 *
          (aipw_score - contrast_stabilized_estimate)^2) /
    sum(contrast_h)^2
  contrast_stabilized_se <- sqrt(contrast_stabilized_variance)





  predictable_working <- matrix(
    NA_real_, nrow = nrow(dat), ncol = 2L,
    dimnames = list(NULL, c("0", "1"))
  )
  has_availability <- "Y_available_month" %in% names(dat)
  for (i in seq_len(nrow(dat))) {
    prior <- if (i > 1L) seq_len(i - 1L) else integer(0)
    available <- if (has_availability && length(prior)) {
      prior[dat$Y_available_month[prior] <= dat$entry_month[i]]
    } else {
      integer(0)
    }
    for (arm in 0:1) {
      idx <- available[dat$A[available] == arm]
      predictable_working[i, as.character(arm)] <- if (length(idx)) {
        mean(dat$Y[idx])
      } else if (arm == 1L) {
        working1
      } else {
        working0
      }
    }
  }
  predictable_score <-
    predictable_working[, "1"] - predictable_working[, "0"] +
    w1 * (dat$Y - predictable_working[, "1"]) -
    w0 * (dat$Y - predictable_working[, "0"])
  predictable_h <- sqrt(dat$pi * (1 - dat$pi))
  predictable_stabilized_estimate <-
    sum(predictable_h * predictable_score) / sum(predictable_h)
  predictable_stabilized_variance <-
    sum(predictable_h^2 *
          (predictable_score - predictable_stabilized_estimate)^2) /
    sum(predictable_h)^2
  predictable_stabilized_se <- sqrt(predictable_stabilized_variance)





  idx1 <- dat$A == 1L
  idx0 <- dat$A == 0L
  n1 <- sum(idx1)
  n0 <- sum(idx0)
  x1 <- sum(dat$Y[idx1])
  x0 <- sum(dat$Y[idx0])
  likelihood_p1 <- x1 / n1
  likelihood_p0 <- x0 / n0
  likelihood_estimate <- likelihood_p1 - likelihood_p0
  likelihood_se <- sqrt(
    likelihood_p1 * (1 - likelihood_p1) / n1 +
      likelihood_p0 * (1 - likelihood_p0) / n0
  )
  pooled <- (x1 + x0) / (n1 + n0)
  likelihood_score_se <- sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n0))
  wilson <- function(x, n, level = 0.95) {
    z <- stats::qnorm(1 - (1 - level) / 2)
    phat <- x / n
    denominator <- 1 + z^2 / n
    centre <- (phat + z^2 / (2 * n)) / denominator
    half_width <- z / denominator * sqrt(
      phat * (1 - phat) / n + z^2 / (4 * n^2)
    )
    c(lower = centre - half_width, upper = centre + half_width)
  }
  likelihood_interval1 <- wilson(x1, n1)
  likelihood_interval0 <- wilson(x0, n0)
  likelihood_newcombe <- c(
    lower = unname(
      likelihood_interval1["lower"] - likelihood_interval0["upper"]
    ),
    upper = unname(
      likelihood_interval1["upper"] - likelihood_interval0["lower"]
    )
  )
  data.frame(
    estimate = estimate,
    se = se,
    reject = as.integer(estimate / se > critical),
    cover = as.integer(abs(estimate - diff(config$p_y)) <= stats::qnorm(0.975) * se),
    allocation_1 = mean(dat$A),
    mean_randomization_probability = mean(dat$pi),
    final_randomization_probability = tail(dat$pi, 1),
    fraction_months_at_upper = mean(upper_hit),
    first_upper_month = first_upper_month,
    fraction_months_at_lower = mean(lower_hit),
    first_lower_month = first_lower_month,
    expected_failures = sum(ifelse(dat$A == 1, 1 - config$p_y["1"],
                                   1 - config$p_y["0"])),
    realized_failures = sum(1 - dat$Y),
    prespecified_aipw_estimate = aipw_estimate,
    prespecified_aipw_se = aipw_se,
    prespecified_aipw_reject = as.integer(
      aipw_estimate / aipw_se > critical
    ),
    prespecified_aipw_cover = as.integer(
      abs(aipw_estimate - truth) <= stats::qnorm(0.975) * aipw_se
    ),
    stabilized_aipw_estimate = stabilized_estimate,
    stabilized_aipw_se = stabilized_se,
    stabilized_aipw_reject = as.integer(
      stabilized_estimate / stabilized_se > critical
    ),
    stabilized_aipw_cover = as.integer(
      abs(stabilized_estimate - truth) <=
        stats::qnorm(0.975) * stabilized_se
    ),
    contrast_stabilized_estimate = contrast_stabilized_estimate,
    contrast_stabilized_se = contrast_stabilized_se,
    contrast_stabilized_reject = as.integer(
      contrast_stabilized_estimate / contrast_stabilized_se > critical
    ),
    contrast_stabilized_cover = as.integer(
      abs(contrast_stabilized_estimate - truth) <=
        stats::qnorm(0.975) * contrast_stabilized_se
    ),
    predictable_stabilized_estimate = predictable_stabilized_estimate,
    predictable_stabilized_se = predictable_stabilized_se,
    predictable_stabilized_reject = as.integer(
      predictable_stabilized_estimate / predictable_stabilized_se > critical
    ),
    predictable_stabilized_cover = as.integer(
      abs(predictable_stabilized_estimate - truth) <=
        stats::qnorm(0.975) * predictable_stabilized_se
    ),
    likelihood_estimate = likelihood_estimate,
    likelihood_se = likelihood_se,
    likelihood_reject = as.integer(
      likelihood_estimate / likelihood_score_se > critical
    ),
    likelihood_cover = as.integer(
      abs(likelihood_estimate - truth) <=
        stats::qnorm(0.975) * likelihood_se
    ),
    likelihood_newcombe_cover = as.integer(
      truth >= likelihood_newcombe["lower"] &&
        truth <= likelihood_newcombe["upper"]
    )
  )
}

simulate_oracle_method <- function(method, stream, config, assignment_seed) {
  set.seed(assignment_seed)
  n <- config$n
  dat <- data.frame(
    id = seq_len(n), entry_month = stream$entry_month,
    A = NA_integer_, pi = NA_real_, Q = NA_integer_, QZ = NA_integer_,
    R = NA_integer_, RY = NA_integer_, Y = NA_integer_
  )

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
    pi_month <- next_oracle_probability(method, history, month, config)
    dat$pi[current] <- pi_month
    dat$A[current] <- stats::rbinom(length(current), 1, pi_month)
    ids <- dat$id[current]
    dat$Y[current] <- ifelse(dat$A[current] == 1, stream$y1[ids], stream$y0[ids])
  }
  cbind(method = method, analyse_oracle_trial(dat, config), stringsAsFactors = FALSE)
}

run_oracle_mechanism_check <- function(n_rep = 500L, seed = 20260928L,
                                       config = canonical_config()) {
  rows <- vector("list", n_rep * length(ORACLE_METHODS))
  k <- 1L
  for (r in seq_len(n_rep)) {
    stream <- generate_canonical_stream(config, seed + 10000L * r)
    for (j in seq_along(ORACLE_METHODS)) {
      rows[[k]] <- simulate_oracle_method(
        ORACLE_METHODS[j], stream, config,
        assignment_seed = seed + 10000L * r + j
      )
      rows[[k]]$replicate <- r
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

summarise_oracle_check <- function(results, config = canonical_config(),
                                   method_order = ORACLE_METHODS) {
  safe_median <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) stats::median(x) else NA_real_
  }
  pieces <- split(results, results$method)
  ans <- do.call(rbind, lapply(pieces, function(x) data.frame(
    method = x$method[1], n_rep = nrow(x),
    allocation_1 = mean(x$allocation_1),
    mean_randomization_probability = mean(x$mean_randomization_probability),
    final_randomization_probability = mean(x$final_randomization_probability),
    fraction_months_at_upper = mean(x$fraction_months_at_upper),
    median_first_upper_month = safe_median(x$first_upper_month),
    probability_ever_at_upper = mean(is.finite(x$first_upper_month)),
    fraction_months_at_lower = mean(x$fraction_months_at_lower),
    median_first_lower_month = safe_median(x$first_lower_month),
    probability_ever_at_lower = mean(is.finite(x$first_lower_month)),
    expected_failures = mean(x$expected_failures),
    realized_failures = mean(x$realized_failures),
    bias = mean(x$estimate) - diff(config$p_y),
    empirical_se = stats::sd(x$estimate),
    mean_se = mean(x$se),
    coverage = mean(x$cover),
    power = mean(x$reject)
  )))
  rownames(ans) <- NULL
  if (all(c(
    "prespecified_aipw_estimate", "prespecified_aipw_se",
    "prespecified_aipw_reject", "prespecified_aipw_cover"
  ) %in% names(results))) {
    aipw <- do.call(rbind, lapply(pieces, function(x) data.frame(
      method = x$method[1],
      prespecified_aipw_bias =
        mean(x$prespecified_aipw_estimate) - diff(config$p_y),
      prespecified_aipw_empirical_se =
        stats::sd(x$prespecified_aipw_estimate),
      prespecified_aipw_mean_se = mean(x$prespecified_aipw_se),
      prespecified_aipw_coverage = mean(x$prespecified_aipw_cover),
      prespecified_aipw_power = mean(x$prespecified_aipw_reject)
    )))
    rownames(aipw) <- NULL
    ans <- merge(ans, aipw, by = "method", all.x = TRUE, sort = FALSE)
  }
  if (all(c(
    "stabilized_aipw_estimate", "stabilized_aipw_se",
    "stabilized_aipw_reject", "stabilized_aipw_cover"
  ) %in% names(results))) {
    stabilized <- do.call(rbind, lapply(pieces, function(x) data.frame(
      method = x$method[1],
      stabilized_aipw_bias =
        mean(x$stabilized_aipw_estimate) - diff(config$p_y),
      stabilized_aipw_empirical_se =
        stats::sd(x$stabilized_aipw_estimate),
      stabilized_aipw_mean_se = mean(x$stabilized_aipw_se),
      stabilized_aipw_coverage = mean(x$stabilized_aipw_cover),
      stabilized_aipw_power = mean(x$stabilized_aipw_reject)
    )))
    rownames(stabilized) <- NULL
    ans <- merge(ans, stabilized, by = "method", all.x = TRUE, sort = FALSE)
  }
  if (all(c(
    "contrast_stabilized_estimate", "contrast_stabilized_se",
    "contrast_stabilized_reject", "contrast_stabilized_cover"
  ) %in% names(results))) {
    contrast_stabilized <- do.call(rbind, lapply(pieces, function(x) data.frame(
      method = x$method[1],
      contrast_stabilized_bias =
        mean(x$contrast_stabilized_estimate) - diff(config$p_y),
      contrast_stabilized_empirical_se =
        stats::sd(x$contrast_stabilized_estimate),
      contrast_stabilized_mean_se = mean(x$contrast_stabilized_se),
      contrast_stabilized_coverage = mean(x$contrast_stabilized_cover),
      contrast_stabilized_power = mean(x$contrast_stabilized_reject)
    )))
    rownames(contrast_stabilized) <- NULL
    ans <- merge(
      ans, contrast_stabilized, by = "method", all.x = TRUE, sort = FALSE
    )
  }
  if (all(c(
    "predictable_stabilized_estimate", "predictable_stabilized_se",
    "predictable_stabilized_reject", "predictable_stabilized_cover"
  ) %in% names(results))) {
    predictable_stabilized <- do.call(rbind, lapply(pieces, function(x) {
      data.frame(
        method = x$method[1],
        predictable_stabilized_bias =
          mean(x$predictable_stabilized_estimate) - diff(config$p_y),
        predictable_stabilized_empirical_se =
          stats::sd(x$predictable_stabilized_estimate),
        predictable_stabilized_mean_se =
          mean(x$predictable_stabilized_se),
        predictable_stabilized_coverage =
          mean(x$predictable_stabilized_cover),
        predictable_stabilized_power =
          mean(x$predictable_stabilized_reject)
      )
    }))
    rownames(predictable_stabilized) <- NULL
    ans <- merge(
      ans, predictable_stabilized, by = "method", all.x = TRUE, sort = FALSE
    )
  }
  if (all(c(
    "likelihood_estimate", "likelihood_se", "likelihood_reject",
    "likelihood_cover", "likelihood_newcombe_cover"
  ) %in% names(results))) {
    likelihood <- do.call(rbind, lapply(pieces, function(x) data.frame(
      method = x$method[1],
      likelihood_bias = mean(x$likelihood_estimate) - diff(config$p_y),
      likelihood_empirical_se = stats::sd(x$likelihood_estimate),
      likelihood_mean_se = mean(x$likelihood_se),
      likelihood_coverage = mean(x$likelihood_cover),
      likelihood_newcombe_coverage = mean(x$likelihood_newcombe_cover),
      likelihood_power = mean(x$likelihood_reject)
    )))
    rownames(likelihood) <- NULL
    ans <- merge(ans, likelihood, by = "method", all.x = TRUE, sort = FALSE)
  }
  ans[match(method_order, ans$method), ]
}
