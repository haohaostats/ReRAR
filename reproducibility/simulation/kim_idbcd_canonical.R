source("simulation/inference_core.R")







kim_urn_target <- function(theta1, theta2) {

  denominator <- (1 - theta1) + (1 - theta2)
  if (!is.finite(denominator) || denominator <= 0) return(0.5)
  (1 - theta2) / denominator
}

kim_dbcd_probability <- function(current1, target1) {

  if (current1 <= 0) return(1)
  if (current1 >= 1) return(0)
  numerator1 <- target1 * (target1 / current1)^2
  target2 <- 1 - target1
  current2 <- 1 - current1
  numerator2 <- target2 * (target2 / current2)^2
  numerator1 / (numerator1 + numerator2)
}

kim_conditional_cell_means <- function(history) {

  means <- matrix(
    NA_real_, nrow = 2L, ncol = 2L,
    dimnames = list(A = c("0", "1"), Z = c("0", "1"))
  )
  for (a in 0:1) {
    for (z in 0:1) {
      idx <- which(
        history$A == a & history$Q == 1 & is.finite(history$QZ) &
          history$QZ == z & history$R == 1
      )
      if (length(idx)) means[as.character(a), as.character(z)] <-
        mean(history$RY[idx])
    }
  }
  means
}

kim_delay_adaptive_means <- function(history, theta0 = c(`0` = 0.5, `1` = 0.5)) {



  cell_means <- kim_conditional_cell_means(history)
  answer <- c(`0` = NA_real_, `1` = NA_real_)
  for (a in 0:1) {
    idx <- which(history$A == a & history$Q == 1)
    if (!length(idx)) next
    yhat <- numeric(length(idx))
    for (j in seq_along(idx)) {
      i <- idx[j]
      if (history$R[i] == 1) {
        yhat[j] <- history$RY[i]
      } else {
        z <- as.character(history$QZ[i])
        yhat[j] <- cell_means[as.character(a), z]
      }
    }
    if (all(is.finite(yhat))) {
      answer[as.character(a)] <-
        (sum(yhat) + theta0[as.character(a)]) / (length(yhat) + 1)
    }
  }
  answer
}

kim_idbcd_next_probability <- function(history,
                                       theta0 = c(`0` = 0.5, `1` = 0.5)) {
  if (!nrow(history) || !all(0:1 %in% history$A)) return(0.5)
  theta <- kim_delay_adaptive_means(history, theta0 = theta0)



  if (!all(is.finite(theta))) return(0.5)
  target1 <- kim_urn_target(theta["1"], theta["0"])
  kim_dbcd_probability(mean(history$A), target1)
}

simulate_kim_idbcd_method <- function(stream, config, assignment_seed,
                                      theta0 = c(`0` = 0.5, `1` = 0.5),
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




    for (i in current) {
      assigned <- which(dat$id < i & is.finite(dat$A))
      history <- dat[assigned, , drop = FALSE]
      pi_i <- kim_idbcd_next_probability(history, theta0 = theta0)
      dat$pi[i] <- pi_i
      dat$A[i] <- stats::rbinom(1L, 1L, pi_i)
      dat$Y[i] <- if (dat$A[i] == 1L) stream$y1[i] else stream$y0[i]
      dat$Y_available_month[i] <- if (dat$A[i] == 1L) {
        stream$y1_available_month[i]
      } else {
        stream$y0_available_month[i]
      }
    }
    if (return_trace) {
      traces[[month]] <- monthly_trace_row(
        "I-DBCD", month, dat, current, dat$pi[current]
      )
    }
  }
  summary <- cbind(
    method = "I-DBCD",
    analyse_oracle_trial(dat, config),
    stringsAsFactors = FALSE
  )
  if (!return_trace) return(summary)
  list(summary = summary, trace = do.call(rbind, traces))
}

validate_kim_idbcd_formulas <- function() {
  stopifnot(abs(kim_urn_target(0.5, 0.2) - 0.8 / 1.3) < 1e-12)
  stopifnot(abs(kim_dbcd_probability(0.5, 0.6) -
                  (0.6 * (0.6 / 0.5)^2) /
                  (0.6 * (0.6 / 0.5)^2 + 0.4 * (0.4 / 0.5)^2)) < 1e-12)
  example <- data.frame(
    A = c(1L, 1L, 0L, 0L), Q = 1L, QZ = c(0L, 0L, 1L, 1L),
    R = c(1L, 0L, 1L, 0L), RY = c(1L, NA, 0L, NA)
  )
  theta <- kim_delay_adaptive_means(example)
  stopifnot(abs(theta["1"] - 2.5 / 3) < 1e-12)
  stopifnot(abs(theta["0"] - 0.5 / 3) < 1e-12)
  invisible(TRUE)
}
