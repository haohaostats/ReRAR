test_that("balanced AIPW reproduces the observed risk difference", {
  A <- rep(0:1, each = 4)
  Y <- c(0, 0, 1, 1, 0, 1, 1, 1)
  fit <- rerar_analyze(A, Y, rep(0.5, 8))
  expect_equal(fit$estimate, mean(Y[A == 1]) - mean(Y[A == 0]))
  expect_true(is.finite(fit$standard_error))
  expect_length(fit$conf.int, 2)
})

