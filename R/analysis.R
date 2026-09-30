rerar_analyze <- function(A, Y, pi, working0 = 0.5, working1 = 0.5,
                          alternative = c("greater", "less", "two.sided"),
                          conf.level = 0.95) {
  alternative <- match.arg(alternative)
  n <- length(A)
  if (length(Y) != n || length(pi) != n || n < 2L) {
    stop("A, Y, and pi must have the same length of at least two.",
         call. = FALSE)
  }
  .validate_binary_arm(A, "A", require_both = TRUE)
  if (anyNA(Y)) stop("All final outcomes Y must be observed.", call. = FALSE)
  .validate_binary_outcome(Y, "Y")
  if (anyNA(pi) || any(!is.finite(pi)) || any(pi <= 0 | pi >= 1)) {
    stop("pi must contain probabilities strictly between 0 and 1.",
         call. = FALSE)
  }
  if (length(working0) == 1L) working0 <- rep(working0, n)
  if (length(working1) == 1L) working1 <- rep(working1, n)
  if (length(working0) != n || length(working1) != n ||
      anyNA(working0) || anyNA(working1)) {
    stop("Working means must be complete scalars or vectors of length n.",
         call. = FALSE)
  }
  score <- working1 - working0 +
    A / pi * (Y - working1) -
    (1 - A) / (1 - pi) * (Y - working0)
  h <- sqrt(pi * (1 - pi))
  estimate <- sum(h * score) / sum(h)
  variance <- sum(h^2 * (score - estimate)^2) / sum(h)^2
  standard_error <- sqrt(variance)
  z_value <- estimate / standard_error
  p_value <- switch(
    alternative,
    greater = stats::pnorm(z_value, lower.tail = FALSE),
    less = stats::pnorm(z_value),
    two.sided = 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)
  )
  critical <- stats::qnorm(1 - (1 - conf.level) / 2)
  out <- list(
    estimate = estimate,
    standard_error = standard_error,
    z_value = z_value,
    p_value = p_value,
    conf.int = estimate + c(-1, 1) * critical * standard_error,
    alternative = alternative,
    conf.level = conf.level,
    n = n
  )
  class(out) <- "rerar_analysis"
  out
}

print.rerar_analysis <- function(x, ...) {
  cat("ReRAR final analysis\n")
  cat(sprintf("  risk difference: %.4f\n", x$estimate))
  cat(sprintf("  standard error:  %.4f\n", x$standard_error))
  cat(sprintf("  p-value (%s): %s\n", x$alternative,
              format.pval(x$p_value, digits = 4)))
  cat(sprintf("  %.1f%% CI: [%.4f, %.4f]\n",
              100 * x$conf.level, x$conf.int[1L], x$conf.int[2L]))
  invisible(x)
}

