source("simulation/inference_core.R")

SAFE_ORACLE_METHODS <- c(
  "CR", "Full-oracle-target", "Safe-Y-only", "Safe-Naive-Z",
  "Safe-Revision-oracle"
)

safe_gate_config <- function(config = canonical_config(), eta = 1,
                             burn_months = 5L) {
  n_looks <- config$n_months - burn_months
  reversal <- revision_induced_reversal_interval(config)




  max_harm_upper <- config$n * (config$power_envelope_upper - 0.5) *
    (config$p_y["0"] - reversal["lower_exclusive"])
  max_harm_lower <- config$n * (0.5 - config$power_envelope_lower) *
    abs(diff(config$design_p_y))

  directional_error_upper <- min(0.5, eta / max_harm_upper)
  directional_error_lower <- min(0.5, eta / max_harm_lower)

  list(
    eta = eta,
    burn_months = burn_months,
    n_looks = n_looks,
    max_harm_upper = unname(max_harm_upper),
    max_harm_lower = unname(max_harm_lower),
    directional_error_upper = unname(directional_error_upper),
    directional_error_lower = unname(directional_error_lower),
    upper_z_boundary = stats::qnorm(
      1 - directional_error_upper / n_looks
    ),
    lower_z_boundary = stats::qnorm(
      1 - directional_error_lower / n_looks
    )
  )
}

safe_signal <- function(method, history, arm, config) {
  idx <- history$A == arm
  if (method == "Safe-Y-only") return(history$RY[idx])
  if (method == "Safe-Naive-Z") return(history$QZ[idx])
  if (method == "Safe-Revision-oracle") {
    m_true <- true_revision_regression(config)
    q <- history$Q[idx] == 1
    r <- history$R[idx] == 1
    z <- history$QZ[idx]
    signal <- rep(NA_real_, sum(idx))
    signal[q] <- m_true[cbind(
      as.character(rep(arm, sum(q))), as.character(z[q])
    )]
    signal[r] <- history$RY[idx][r]
    return(signal)
  }
  stop("Unknown safe online estimator: ", method)
}

safe_effect_z <- function(method, history, config, min_per_arm = 8L) {
  signals <- lapply(0:1, function(a) {
    x <- safe_signal(method, history, a, config)
    x[is.finite(x)]
  })
  n <- lengths(signals)
  if (any(n < min_per_arm)) return(NA_real_)
  variance <- vapply(signals, stats::var, numeric(1))
  se <- sqrt(sum(variance / n))
  if (!is.finite(se) || se <= 0) return(NA_real_)
  (mean(signals[[2]]) - mean(signals[[1]])) / se
}

next_safe_probability <- function(method, history, month, config, gate,
                                  gamma = 2) {
  if (method == "CR") return(0.5)
  if (method == "Full-oracle-target") {
    return(unname(oracle_welfare_target(config)["allocation_1"]))
  }
  if (month <= gate$burn_months || nrow(history) == 0) return(0.5)

  z <- safe_effect_z(method, history, config)
  target <- if (!is.finite(z)) {
    0.5
  } else if (z > gate$upper_z_boundary) {
    config$power_envelope_upper
  } else if (z < -gate$lower_z_boundary) {
    config$power_envelope_lower
  } else {
    0.5
  }

  dbcd_update(
    mean(history$A), target, gamma,
    config$power_envelope_lower, config$power_envelope_upper
  )
}

simulate_safe_oracle_method <- function(method, stream, config, gate,
                                        assignment_seed) {
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
    pi_month <- next_safe_probability(method, history, month, config, gate)
    dat$pi[current] <- pi_month
    dat$A[current] <- stats::rbinom(length(current), 1, pi_month)
    ids <- dat$id[current]
    dat$Y[current] <- ifelse(
      dat$A[current] == 1, stream$y1[ids], stream$y0[ids]
    )
  }
  cbind(method = method, analyse_oracle_trial(dat, config), stringsAsFactors = FALSE)
}

run_safe_oracle_check <- function(n_rep = 1000L, seed = 20261003L,
                                  config = canonical_config(), eta = 1,
                                  methods = SAFE_ORACLE_METHODS) {
  gate <- safe_gate_config(config, eta)
  rows <- vector("list", n_rep * length(methods))
  k <- 1L
  for (r in seq_len(n_rep)) {
    stream <- generate_canonical_stream(config, seed + 10000L * r)
    for (j in seq_along(methods)) {
      rows[[k]] <- simulate_safe_oracle_method(
        methods[j], stream, config, gate,
        assignment_seed = seed + 10000L * r + j
      )
      rows[[k]]$replicate <- r
      k <- k + 1L
    }
  }
  list(
    gate = gate,
    results = do.call(rbind, rows)
  )
}

summarise_safe_oracle_check <- function(results, config = canonical_config()) {
  method_order <- SAFE_ORACLE_METHODS[SAFE_ORACLE_METHODS %in% unique(results$method)]
  summarise_oracle_check(results, config, method_order = method_order)
}
