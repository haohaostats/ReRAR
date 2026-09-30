validation_fixture <- function() {
  rerar_validation(
    arm = 0:1,
    se_success = c(36, 131), se_failure = c(19, 19),
    sp_success = c(48, 44), sp_failure = c(8, 25)
  )
}

test_that("uninformative updates remain at equal allocation", {
  history <- rerar_history(integer(), integer(), integer())
  update <- rerar_update(history, validation_fixture())
  expect_equal(update$probability, 0.5)
  expect_equal(update$target, 0.5)
})

test_that("updates are bounded and randomization stores probabilities", {
  history <- rerar_history(
    A = rep(0:1, each = 12),
    Z = c(rep(c(0, 1), 6), rep(1, 12)),
    Y = c(rep(c(0, 1), 6), rep(c(1, 1, 0), 4))
  )
  update <- rerar_update(history, validation_fixture())
  expect_true(is.finite(update$probability))
  expect_gte(update$probability, 0.39)
  expect_lte(update$probability, 0.62)
  next_batch <- rerar_randomize(10, update, seed = 2026)
  expect_equal(nrow(next_batch), 10)
  expect_equal(next_batch$pi, rep(update$probability, 10))
})

