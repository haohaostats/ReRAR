



historical_revision_summary <- function(history, config) {
  if (!nrow(history) || !all(0:1 %in% history$A)) {
    return(list(theta = c(`0` = NA_real_, `1` = NA_real_), z = NA_real_))
  }
  has_q <- vapply(0:1, function(a) {
    any(history$A == a & history$Q == 1, na.rm = TRUE)
  }, logical(1))
  if (!all(has_q)) {
    return(list(theta = c(`0` = NA_real_, `1` = NA_real_), z = NA_real_))
  }
  fits <- lapply(0:1, function(a) {
    fit_historical_revision_arm(
      history, a, known_revision = FALSE, config = config
    )
  })
  theta <- c(`0` = unname(fits[[1]]["theta"]),
             `1` = unname(fits[[2]]["theta"]))
  variance <- unname(fits[[1]]["variance"] + fits[[2]]["variance"])
  z <- if (is.finite(variance) && variance > 0) {
    unname((theta["1"] - theta["0"]) / sqrt(variance))
  } else {
    NA_real_
  }
  list(theta = theta, z = z)
}

welfare_rerar_target <- function(z, config) {
  if (length(z) != 1L || is.na(z)) return(0.5)
  evidence <- 2 * stats::pnorm(abs(z)) - 1
  if (z >= 0) {
    0.5 + evidence * (config$power_envelope_upper - 0.5)
  } else {
    0.5 - evidence * (0.5 - config$power_envelope_lower)
  }
}

simulate_welfare_rerar_method <- function(stream, config, assignment_seed,
                                          return_trace = FALSE) {
  set.seed(assignment_seed)
  n <- config$n
  dat <- data.frame(
    id = seq_len(n), entry_month = stream$entry_month,
    A = NA_integer_, pi = NA_real_, Q = NA_integer_, QZ = NA_integer_,
    R = NA_integer_, RY = NA_integer_, Y = NA_integer_,
    Y_available_month = NA_real_
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
    information_history <- dat[previous, , drop = FALSE]
    evidence <- historical_revision_summary(information_history, config)
    target1 <- welfare_rerar_target(evidence$z, config)


    pi_month <- if (!length(previous)) {
      0.5
    } else {
      dbcd_update(
        mean(dat$A[previous]), target1, gamma = 2,
        config$power_envelope_lower, config$power_envelope_upper
      )
    }
    dat$pi[current] <- pi_month
    dat$A[current] <- stats::rbinom(length(current), 1L, pi_month)
    ids <- dat$id[current]
    dat$Y[current] <- ifelse(
      dat$A[current] == 1L, stream$y1[ids], stream$y0[ids]
    )
    dat$Y_available_month[current] <- ifelse(
      dat$A[current] == 1L,
      stream$y1_available_month[ids], stream$y0_available_month[ids]
    )
    if (return_trace) {
      traces[[month]] <- monthly_trace_row(
        "ReRAR", month, dat, current, pi_month,
        target_1 = target1, evidence_z = evidence$z
      )
    }
  }
  summary <- cbind(
    method = "Welfare-ReRAR",
    analyse_oracle_trial(dat, config),
    stringsAsFactors = FALSE
  )
  if (!return_trace) return(summary)
  list(summary = summary, trace = do.call(rbind, traces))
}

validate_welfare_rerar <- function(config = canonical_config()) {
  stopifnot(abs(welfare_rerar_target(0, config) - 0.5) < 1e-12)
  stopifnot(welfare_rerar_target(Inf, config) ==
              config$power_envelope_upper)
  stopifnot(welfare_rerar_target(-Inf, config) ==
              config$power_envelope_lower)
  stopifnot(welfare_rerar_target(8, config) <= config$power_envelope_upper)
  stopifnot(welfare_rerar_target(-8, config) >= config$power_envelope_lower)
  stopifnot(welfare_rerar_target(1, config) > 0.5)
  stopifnot(welfare_rerar_target(-1, config) < 0.5)
  invisible(TRUE)
}
