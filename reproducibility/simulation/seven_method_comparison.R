source("simulation/revision_likelihood.R")
source("simulation/fl_cara_canonical.R")
source("simulation/kim_idbcd_canonical.R")
source("simulation/welfare_rerar_core.R")

SEVEN_METHODS <- c(
  "CR", "Y-DBCD", "FL-CARA", "I-DBCD", "SP-RAR", "P-RAR", "ReRAR"
)

monthly_trace_row <- function(method, month, dat, current, pi_month,
                              target_1 = NA_real_, evidence_z = NA_real_) {
  assigned <- is.finite(dat$A)
  data.frame(
    method = method,
    month = month,
    randomization_probability = mean(pi_month),
    target_1 = target_1,
    evidence_z = evidence_z,
    cumulative_allocation_1 = mean(dat$A[assigned]),
    available_z = sum(dat$Q == 1, na.rm = TRUE),
    available_y = sum(dat$R == 1, na.rm = TRUE),
    enrolled = sum(assigned),
    stringsAsFactors = FALSE
  )
}

with_optional_trace <- function(summary, traces, return_trace) {
  if (!return_trace) return(summary)
  list(summary = summary, trace = do.call(rbind, traces))
}

ethical_response_target <- function(p0, p1, config) {
  target <- sqrt(clamp_unit(p1)) /
    (sqrt(clamp_unit(p1)) + sqrt(clamp_unit(p0)))
  clip_probability(
    target, config$power_envelope_lower, config$power_envelope_upper
  )
}

competitor_working_means <- function(method, history) {
  means <- c(`0` = 0.5, `1` = 0.5)

  for (a in 0:1) {
    idx <- history$A == a
    if (method == "Y-DBCD") {
      signal <- history$RY[idx]
    } else if (method == "P-RAR") {
      signal <- history$QZ[idx]
    } else if (method == "SP-RAR") {
      signal <- history$QZ[idx]
      reviewed <- history$R[idx] == 1
      signal[reviewed] <- history$RY[idx][reviewed]
    } else {
      stop("Unknown comparator: ", method)
    }
    means[as.character(a)] <- smoothed_mean(signal)
  }
  means
}

next_competitor_probability <- function(method, history, month, config,
                                        burn_months = 5L, gamma = 2) {
  if (method == "CR" || month <= burn_months || !nrow(history)) return(0.5)
  means <- competitor_working_means(method, history)
  target <- ethical_response_target(means["0"], means["1"], config)
  dbcd_update(
    mean(history$A), target, gamma,
    config$power_envelope_lower, config$power_envelope_upper
  )
}

simulate_standard_competitor <- function(method, stream, config,
                                         assignment_seed,
                                         burn_months = 5L,
                                         return_trace = FALSE) {
  set.seed(assignment_seed)
  n <- config$n
  dat <- data.frame(
    id = seq_len(n), entry_month = stream$entry_month,
    A = NA_integer_, pi = NA_real_, Q = NA_integer_, QZ = NA_integer_,
    R = NA_integer_, RY = NA_integer_, Y = NA_integer_
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
    pi_month <- next_competitor_probability(
      method, history, month, config, burn_months = burn_months
    )
    dat$pi[current] <- pi_month
    dat$A[current] <- stats::rbinom(length(current), 1, pi_month)
    ids <- dat$id[current]
    dat$Y[current] <- ifelse(
      dat$A[current] == 1, stream$y1[ids], stream$y0[ids]
    )
    if (return_trace) {
      traces[[month]] <- monthly_trace_row(
        method, month, dat, current, pi_month
      )
    }
  }
  summary <- cbind(
    method = method, analyse_oracle_trial(dat, config), stringsAsFactors = FALSE
  )
  with_optional_trace(summary, traces, return_trace)
}

simulate_seven_method <- function(method, stream, config, gate,
                                  assignment_seed, return_trace = FALSE) {
  if (method == "FL-CARA") {
    return(simulate_fl_cara_method(
      stream, config, assignment_seed, return_trace = return_trace
    ))
  }
  if (method == "I-DBCD") {
    return(simulate_kim_idbcd_method(
      stream, config, assignment_seed, return_trace = return_trace
    ))
  }
  if (method == "ReRAR") {
    result <- simulate_welfare_rerar_method(
      stream, config, assignment_seed, return_trace = return_trace
    )
    if (return_trace) {
      result$summary$method <- "ReRAR"
      result$trace$method <- "ReRAR"
    } else {
      result$method <- "ReRAR"
    }
    return(result)
  }
  if (method == "ReRAR-H-DBCD") {

    internal_method <- "H-DBCD"
    result <- simulate_historical_rar_method(
      internal_method, stream, config, gate, assignment_seed
    )
    result$method <- method
    return(result)
  }
  simulate_standard_competitor(
    method, stream, config, assignment_seed,
    return_trace = return_trace
  )
}

run_seven_method_comparison <- function(n_rep = 1000L, seed = 20261121L,
                                        config = canonical_config(), eta = 1,
                                        methods = SEVEN_METHODS) {
  gate <- safe_gate_config(config, eta)
  rows <- vector("list", n_rep * length(methods))
  k <- 1L
  for (replicate in seq_len(n_rep)) {
    stream <- generate_canonical_stream(config, seed + 10000L * replicate)
    for (j in seq_along(methods)) {
      rows[[k]] <- simulate_seven_method(
        methods[j], stream, config, gate,
        assignment_seed = seed + 10000L * replicate + j
      )
      rows[[k]]$replicate <- replicate
      k <- k + 1L
    }
  }
  columns <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(columns, names(row))
    for (label in missing) row[[label]] <- NA
    row[columns]
  })
  list(gate = gate, results = do.call(rbind, rows))
}

summarise_seven_method_comparison <- function(results,
                                              config = canonical_config()) {
  method_order <- c(SEVEN_METHODS, "ReRAR-H-DBCD")
  method_order <- method_order[method_order %in% unique(results$method)]
  summarise_oracle_check(results, config, method_order = method_order)
}
