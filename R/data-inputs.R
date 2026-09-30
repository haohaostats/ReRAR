rerar_validation <- function(arm, se_success, se_failure,
                             sp_success, sp_failure) {
  x <- data.frame(
    arm = arm,
    se_success = se_success,
    se_failure = se_failure,
    sp_success = sp_success,
    sp_failure = sp_failure
  )
  .validate_binary_arm(x$arm, "arm", require_both = TRUE)
  if (nrow(x) != 2L || anyDuplicated(x$arm)) {
    stop("Supply exactly one validation row for each of arms 0 and 1.",
         call. = FALSE)
  }
  count_names <- names(x)[-1L]
  for (name in count_names) {
    .validate_counts(x[[name]], name)
  }
  if (any(x$se_success + x$se_failure == 0) ||
      any(x$sp_success + x$sp_failure == 0)) {
    stop("Sensitivity and specificity totals must be positive in each arm.",
         call. = FALSE)
  }
  x <- x[order(x$arm), , drop = FALSE]
  rownames(x) <- NULL
  class(x) <- c("rerar_validation", "data.frame")
  x
}

rerar_history <- function(A, Z = rep(NA_integer_, length(A)),
                          Y = rep(NA_integer_, length(A)),
                          id = seq_along(A)) {
  n <- length(A)
  if (length(Z) != n || length(Y) != n || length(id) != n) {
    stop("A, Z, Y, and id must have the same length.", call. = FALSE)
  }
  .validate_binary_arm(A, "A")
  .validate_binary_outcome(Z, "Z")
  .validate_binary_outcome(Y, "Y")
  if (any(!is.na(Y) & is.na(Z))) {
    stop("A returned final outcome Y requires a visible provisional outcome Z.",
         call. = FALSE)
  }
  x <- data.frame(
    id = id,
    A = as.integer(A),
    Q = as.integer(!is.na(Z)),
    QZ = as.integer(Z),
    R = as.integer(!is.na(Y)),
    RY = as.integer(Y)
  )
  class(x) <- c("rerar_history", "data.frame")
  x
}

