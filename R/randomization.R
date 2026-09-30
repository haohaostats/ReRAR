.welfare_target <- function(z, lower, upper) {
  if (length(z) != 1L || !is.finite(z)) {
    if (is.infinite(z) && z > 0) return(upper)
    if (is.infinite(z) && z < 0) return(lower)
    return(0.5)
  }
  evidence <- 2 * stats::pnorm(abs(z)) - 1
  if (z >= 0) {
    0.5 + evidence * (upper - 0.5)
  } else {
    0.5 - evidence * (0.5 - lower)
  }
}

.dbcd_probability <- function(current_allocation, target, gamma,
                              lower, upper) {
  current_allocation <- .clip_probability(current_allocation, 0.02, 0.98)
  target <- .clip_probability(target, lower, upper)
  numerator <- target * (target / current_allocation)^gamma
  denominator <- numerator +
    (1 - target) * ((1 - target) / (1 - current_allocation))^gamma
  .clip_probability(numerator / denominator, lower, upper)
}

rerar_update <- function(history, validation, gamma = 2,
                         lower = 0.39, upper = 0.62) {
  if (!inherits(history, "rerar_history")) {
    stop("history must be created by rerar_history().", call. = FALSE)
  }
  if (!inherits(validation, "rerar_validation")) {
    stop("validation must be created by rerar_validation().", call. = FALSE)
  }
  .validate_bounds(lower, upper)
  if (length(gamma) != 1L || !is.finite(gamma) || gamma < 0) {
    stop("gamma must be one non-negative finite number.", call. = FALSE)
  }
  evidence <- .revision_effect_summary(history, validation)
  target <- .welfare_target(evidence$z, lower, upper)
  probability <- if (!nrow(history)) {
    0.5
  } else {
    .dbcd_probability(mean(history$A), target, gamma, lower, upper)
  }
  out <- list(
    probability = probability,
    target = target,
    z = evidence$z,
    theta = evidence$theta,
    variance = evidence$variance,
    gamma = gamma,
    bounds = c(lower = lower, upper = upper)
  )
  class(out) <- "rerar_update"
  out
}

print.rerar_update <- function(x, ...) {
  cat("ReRAR allocation update\n")
  cat(sprintf("  next-batch probability: %.4f\n", x$probability))
  cat(sprintf("  target allocation:      %.4f\n", x$target))
  cat(sprintf("  evidence z:             %s\n",
              if (is.finite(x$z)) sprintf("%.3f", x$z) else "not available"))
  invisible(x)
}

rerar_randomize <- function(n, update, seed = NULL) {
  if (length(n) != 1L || !is.finite(n) || n < 1 || n != floor(n)) {
    stop("n must be one positive integer.", call. = FALSE)
  }
  if (!inherits(update, "rerar_update")) {
    stop("update must be returned by rerar_update().", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  data.frame(
    A = stats::rbinom(n, 1L, update$probability),
    pi = rep(update$probability, n)
  )
}

