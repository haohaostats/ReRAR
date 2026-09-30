








canonical_config <- function() {


  p_y <- c(`0` = 0.35, `1` = 0.50)






  sensitivity <- c(`0` = (90 - 19) / 90, `1` = (150 - 6) / 150)
  specificity <- c(`0` = (89 - 9) / 89, `1` = (30 - 5) / 30)

  list(



    n = 380L,
    n_months = 20L,
    patients_per_month = 19L,
    p_y = p_y,
    design_p_y = p_y,
    sensitivity = sensitivity,
    specificity = specificity,
    investigator_delay_months = 2,
    central_delay_median_months = 3.05,
    central_delay_q90_months = 7.5,
    alpha_one_sided = 0.025,
    final_critical_value = 2.00,
    target_power = 0.80,
    positivity_lower = 0.15,
    positivity_upper = 0.85,



    power_envelope_lower = 0.39,
    power_envelope_upper = 0.62
  )
}

central_delay_sdlog <- function(config) {
  log(config$central_delay_q90_months / config$central_delay_median_months) /
    stats::qnorm(0.90)
}

revision_quantities <- function(config = canonical_config()) {
  rows <- lapply(0:1, function(a) {
    p <- config$p_y[as.character(a)]
    se <- config$sensitivity[as.character(a)]
    sp <- config$specificity[as.character(a)]
    p_z <- se * p + (1 - sp) * (1 - p)
    discordance <- (1 - se) * p + (1 - sp) * (1 - p)
    data.frame(
      arm = a, final_response = p, sensitivity = se, specificity = sp,
      provisional_response = p_z, discordance = discordance
    )
  })
  do.call(rbind, rows)
}

fixed_allocation_power <- function(rho, config = canonical_config()) {
  p0 <- config$p_y["0"]
  p1 <- config$p_y["1"]
  variance <- p1 * (1 - p1) / (config$n * rho) +
    p0 * (1 - p0) / (config$n * (1 - rho))
  unname(stats::pnorm((p1 - p0) / sqrt(variance) -
                        config$final_critical_value))
}

exact_fixed_wald_power <- function(n1, config = canonical_config(),
                                   p_y = config$p_y) {
  n1 <- as.integer(n1)
  n0 <- config$n - n1
  if (n1 <= 0L || n0 <= 0L) return(NA_real_)

  x1 <- 0:n1
  x0 <- 0:n0
  phat1 <- x1 / n1
  phat0 <- x0 / n0
  difference <- outer(phat1, phat0, "-")
  se <- sqrt(
    outer(phat1 * (1 - phat1) / (n1 - 1), rep(1, length(phat0))) +
      outer(rep(1, length(phat1)), phat0 * (1 - phat0) / (n0 - 1))
  )
  reject <- difference / se > config$final_critical_value
  reject[is.na(reject)] <- FALSE
  probability <- outer(
    stats::dbinom(x1, n1, p_y["1"]),
    stats::dbinom(x0, n0, p_y["0"])
  )
  unname(sum(probability * reject))
}

power_envelope_audit <- function(config = canonical_config()) {
  n1 <- seq_len(config$n - 1L)
  power <- vapply(n1, exact_fixed_wald_power, numeric(1), config = config)
  locked_lower_n1 <- round(config$n * config$power_envelope_lower)
  locked_upper_n1 <- round(config$n * config$power_envelope_upper)
  data.frame(
    n1 = n1,
    allocation_1 = n1 / config$n,
    exact_power = power,
    feasible = power >= config$target_power,
    inside_locked_envelope = n1 >= locked_lower_n1 & n1 <= locked_upper_n1
  )
}

oracle_welfare_target <- function(config = canonical_config()) {
  target <- if (config$p_y["1"] > config$p_y["0"]) {
    config$power_envelope_upper
  } else if (config$p_y["1"] < config$p_y["0"]) {
    config$power_envelope_lower
  } else {
    0.5
  }
  n1 <- round(config$n * target)
  expected_failures <- config$n * (
    target * (1 - config$p_y["1"]) +
      (1 - target) * (1 - config$p_y["0"])
  )
  c(allocation_1 = unname(target),
    power = exact_fixed_wald_power(n1, config),
    expected_failures = unname(expected_failures))
}

revision_induced_reversal_interval <- function(config = canonical_config()) {
  p0 <- config$p_y["0"]
  q0 <- config$sensitivity["0"] * p0 +
    (1 - config$specificity["0"]) * (1 - p0)
  intercept1 <- 1 - config$specificity["1"]
  slope1 <- config$sensitivity["1"] + config$specificity["1"] - 1
  threshold <- (q0 - intercept1) / slope1
  c(lower_exclusive = unname(threshold), upper_exclusive = unname(p0))
}

generate_canonical_stream <- function(config = canonical_config(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(config$n == config$n_months * config$patients_per_month)
  n <- config$n
  entry <- rep(seq_len(config$n_months), each = config$patients_per_month)

  generate_arm <- function(a) {
    p <- config$p_y[as.character(a)]
    se <- config$sensitivity[as.character(a)]
    sp <- config$specificity[as.character(a)]
    y <- stats::rbinom(n, 1, p)
    z_prob <- ifelse(y == 1, se, 1 - sp)
    z <- stats::rbinom(n, 1, z_prob)
    incremental_delay <- stats::rlnorm(
      n,
      meanlog = log(config$central_delay_median_months),
      sdlog = central_delay_sdlog(config)
    )
    data.frame(y = y, z = z, incremental_delay = incremental_delay)
  }

  arm0 <- generate_arm(0)
  arm1 <- generate_arm(1)
  list(
    entry_month = entry,
    z_available_month = entry + config$investigator_delay_months,
    y0 = arm0$y,
    z0 = arm0$z,
    y0_available_month = entry + config$investigator_delay_months + arm0$incremental_delay,
    y1 = arm1$y,
    z1 = arm1$z,
    y1_available_month = entry + config$investigator_delay_months + arm1$incremental_delay
  )
}

canonical_information_at <- function(stream, month, assigned_arm) {
  ids <- seq_along(assigned_arm)
  y <- ifelse(assigned_arm == 1, stream$y1, stream$y0)
  z <- ifelse(assigned_arm == 1, stream$z1, stream$z0)
  y_time <- ifelse(
    assigned_arm == 1, stream$y1_available_month, stream$y0_available_month
  )
  q <- stream$z_available_month <= month
  r <- y_time <= month
  data.frame(
    id = ids, entry_month = stream$entry_month, A = assigned_arm,
    Q = as.integer(q), QZ = ifelse(q, z, NA_integer_),
    R = as.integer(r), RY = ifelse(r, y, NA_integer_)
  )
}
