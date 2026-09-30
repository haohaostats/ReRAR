source("simulation/inference_core.R")





fl_delayed_assign <- function(c1, c0, a0 = 0, b0 = 0, a, b, low, up,
                              tie_tolerance = 1e-10) {
  stopifnot(length(a) == length(b), low < up)
  n_stage <- length(a)
  if (!n_stage) return(list(allocation = numeric(), objective = Inf))




  active <- a + b > tie_tolerance
  answer <- rep(0.5, n_stage)
  if (!any(active)) {
    return(list(allocation = answer, objective = Inf))
  }

  original_active <- which(active)
  aa <- a[active]
  bb <- b[active]
  ratio <- ifelse(bb > tie_tolerance, aa / bb, Inf)
  stage_index <- order(ratio, seq_along(ratio))
  aa <- aa[stage_index]
  bb <- bb[stage_index]
  j_total <- length(aa)





  max_information1 <- a0 + sum(up * aa)
  max_information0 <- b0 + sum((1 - low) * bb)
  if (max_information1 <= tie_tolerance ||
      max_information0 <= tie_tolerance) {
    return(list(allocation = answer, objective = Inf))
  }

  candidates <- matrix(NA_real_, nrow = j_total, ncol = j_total)
  for (i in seq_len(j_total)) {
    for (j in seq_len(j_total)) {
      candidates[i, j] <- up * (j > i) + low * (j < i)
    }
  }
  objective <- rep(Inf, j_total)

  for (i in seq_len(j_total)) {
    objective_function <- function(x) {
      candidate <- candidates[i, ]
      candidate[i] <- x
      denominator1 <- a0 + sum(candidate * aa)
      denominator0 <- b0 + sum((1 - candidate) * bb)
      if (denominator1 <= 0 || denominator0 <= 0) return(Inf)
      c1 / denominator1 + c0 / denominator0
    }
    optimum <- stats::optimize(
      objective_function, interval = c(low, up), maximum = FALSE
    )
    candidates[i, i] <- optimum$minimum
    objective[i] <- optimum$objective
  }

  chosen <- which(objective <= min(objective) + tie_tolerance)[1]
  sorted_allocation <- candidates[chosen, ]
  active_allocation <- numeric(j_total)
  active_allocation[stage_index] <- sorted_allocation
  answer[original_active] <- active_allocation
  list(allocation = answer, objective = objective[chosen])
}

observed_delay_cdf <- function(dat, month, arm, max_delay,
                               conservative = TRUE) {
  mass <- numeric(max_delay + 1L)
  estimable_max <- max(0L, month - 1L)
  for (delay in 0:min(max_delay, estimable_max)) {
    risk <- dat$A == arm & dat$entry_month <= month - delay
    denominator <- sum(risk, na.rm = TRUE)
    mass[delay + 1L] <- if (denominator) {
      sum(risk & dat$delay_stage == delay, na.rm = TRUE) / denominator
    } else 0
  }
  cdf <- pmin(1, cummax(cumsum(mass)))
  if (max_delay > estimable_max) {
    tail_value <- if (conservative) cdf[estimable_max + 1L] else 1
    cdf[(estimable_max + 2L):(max_delay + 1L)] <- tail_value
  }
  cdf
}

fl_cara_variance_estimate <- function(dat, arm, prior = 0.5) {
  y <- dat$RY[dat$A == arm & dat$R == 1]
  y <- y[is.finite(y)]
  p <- (sum(y) + prior) / (length(y) + 2 * prior)
  max(0.01, p * (1 - p))
}

fl_cara_next_probability <- function(history, month, config,
                                     burn_months = 5L) {
  if (month <= burn_months || !nrow(history)) return(0.5)
  final_month <- config$n_months
  future_months <- (month + 1L):final_month
  if (!length(future_months) || future_months[1] > final_month) return(0.5)

  max_delay <- final_month - 1L
  rho1 <- observed_delay_cdf(history, month, 1, max_delay)
  rho0 <- observed_delay_cdf(history, month, 0, max_delay)
  stage_weight <- rep(1 / final_month, final_month)
  horizon <- final_month - seq_len(final_month)
  coefficient1 <- stage_weight * rho1[horizon + 1L]
  coefficient0 <- stage_weight * rho0[horizon + 1L]

  realised_allocation <- vapply(seq_len(month), function(stage) {
    values <- history$A[history$entry_month == stage]
    if (length(values)) mean(values) else 0.5
  }, numeric(1))
  a0 <- sum(coefficient1[seq_len(month)] * realised_allocation)
  b0 <- sum(coefficient0[seq_len(month)] * (1 - realised_allocation))

  result <- fl_delayed_assign(
    c1 = fl_cara_variance_estimate(history, 1),
    c0 = fl_cara_variance_estimate(history, 0),
    a0 = a0,
    b0 = b0,
    a = coefficient1[future_months],
    b = coefficient0[future_months],
    low = config$power_envelope_lower,
    up = config$power_envelope_upper
  )
  clip_probability(
    result$allocation[1],
    config$power_envelope_lower,
    config$power_envelope_upper
  )
}

simulate_fl_cara_method <- function(stream, config, assignment_seed,
                                    burn_months = 5L,
                                    return_trace = FALSE) {
  set.seed(assignment_seed)
  n <- config$n
  dat <- data.frame(
    id = seq_len(n), entry_month = stream$entry_month,
    A = NA_integer_, pi = NA_real_, Q = NA_integer_, QZ = NA_integer_,
    R = NA_integer_, RY = NA_integer_, Y = NA_integer_,
    y_available_month = NA_real_, delay_stage = NA_integer_
  )
  traces <- if (return_trace) vector("list", config$n_months) else NULL

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
    pi_month <- fl_cara_next_probability(
      history, month, config, burn_months = burn_months
    )
    dat$pi[current] <- pi_month
    dat$A[current] <- stats::rbinom(length(current), 1, pi_month)
    ids <- dat$id[current]
    arm1 <- dat$A[current] == 1
    dat$Y[current] <- ifelse(arm1, stream$y1[ids], stream$y0[ids])
    dat$y_available_month[current] <- ifelse(
      arm1, stream$y1_available_month[ids], stream$y0_available_month[ids]
    )
    dat$delay_stage[current] <- ceiling(
      dat$y_available_month[current] - dat$entry_month[current]
    )
    if (return_trace) {
      traces[[month]] <- monthly_trace_row(
        "FL-CARA", month, dat, current, pi_month
      )
    }
  }
  summary <- cbind(
    method = "FL-CARA",
    analyse_oracle_trial(dat, config),
    stringsAsFactors = FALSE
  )
  if (!return_trace) return(summary)
  list(summary = summary, trace = do.call(rbind, traces))
}

validate_fl_delayed_assign <- function(seed = 20261120L, grid_step = 0.01) {
  set.seed(seed)
  a <- stats::runif(3, 0.1, 1)
  b <- stats::runif(3, 0.1, 1)
  c1 <- 0.22
  c0 <- 0.18
  a0 <- 0.2
  b0 <- 0.2
  fit <- fl_delayed_assign(c1, c0, a0, b0, a, b, 0.39, 0.62)
  grid <- seq(0.39, 0.62, by = grid_step)
  brute <- expand.grid(e1 = grid, e2 = grid, e3 = grid)
  value <- apply(brute, 1, function(e) {
    c1 / (a0 + sum(a * e)) + c0 / (b0 + sum(b * (1 - e)))
  })
  data.frame(
    algorithm_objective = fit$objective,
    brute_grid_objective = min(value),
    excess_over_grid = fit$objective - min(value)
  )
}
