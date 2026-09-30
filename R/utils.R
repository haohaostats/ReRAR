.clip_probability <- function(x, lower = 1e-10, upper = 1 - 1e-10) {
  pmin(upper, pmax(lower, x))
}

.validate_binary_arm <- function(x, name, require_both = FALSE) {
  if (anyNA(x) || !all(x %in% 0:1)) {
    stop(name, " must contain only 0 and 1.", call. = FALSE)
  }
  if (require_both && !all(0:1 %in% x)) {
    stop(name, " must contain both 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}

.validate_binary_outcome <- function(x, name) {
  observed <- x[!is.na(x)]
  if (!all(observed %in% 0:1)) {
    stop(name, " must contain only 0, 1, and NA.", call. = FALSE)
  }
  invisible(TRUE)
}

.validate_counts <- function(x, name) {
  if (anyNA(x) || any(!is.finite(x)) || any(x < 0) || any(x != floor(x))) {
    stop(name, " must contain non-negative integer counts.", call. = FALSE)
  }
  invisible(TRUE)
}

.validate_bounds <- function(lower, upper) {
  if (length(lower) != 1L || length(upper) != 1L ||
      !is.finite(lower) || !is.finite(upper) ||
      lower <= 0 || lower >= 0.5 || upper <= 0.5 || upper >= 1) {
    stop("Require 0 < lower < 0.5 < upper < 1.", call. = FALSE)
  }
  invisible(TRUE)
}

