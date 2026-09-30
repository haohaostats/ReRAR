test_that("the example separates published counts from synthetic records", {
  e <- new.env()
  utils::data("exam_example", package = "ReRAR", envir = e)
  x <- e$exam_example
  expect_equal(x$validation$se_success, c(36, 131))
  expect_equal(x$validation$sp_failure, c(8, 25))
  expect_equal(nrow(x$trial), 154L)
  expect_equal(nrow(x$final), 330L)
  expect_equal(length(unique(x$final$id)), 330L)
  expect_match(x$metadata$participant_origin, "Synthetic")
  expect_true(all(x$final$working0 == 0.5 & x$final$working1 == 0.5))
})

test_that("every stored assignment probability uses only visible prior data", {
  e <- new.env()
  utils::data("exam_example", package = "ReRAR", envir = e)
  x <- e$exam_example
  for (month in seq_len(15L)) {
    previous <- x$final[x$final$entry_month < month, , drop = FALSE]
    z <- ifelse(previous$Z_available <= month, previous$Z, NA_integer_)
    y <- ifelse(previous$Y_available <= month, previous$Y, NA_integer_)
    h <- rerar_history(previous$A, z, y, previous$id)
    update <- rerar_update(h, x$validation, lower = 0.20, upper = 0.78)
    stored <- x$final$pi[x$final$entry_month == month]
    expect_equal(stored, rep(update$probability, length(stored)),
                 tolerance = 1e-7)
    if (month == 8L) expect_equal(h, x$trial)
  }
  fit <- with(x$final, rerar_analyze(A, Y, pi, working0, working1))
  expect_true(is.finite(fit$estimate))
  expect_true(is.finite(fit$standard_error))
})
