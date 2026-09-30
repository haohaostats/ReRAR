test_that("input constructors validate and standardize data", {
  validation <- rerar_validation(
    arm = 0:1,
    se_success = c(36, 131), se_failure = c(19, 19),
    sp_success = c(48, 44), sp_failure = c(8, 25)
  )
  expect_s3_class(validation, "rerar_validation")
  history <- rerar_history(c(0, 1), c(1, 0), c(1, NA))
  expect_s3_class(history, "rerar_history")
  expect_equal(history$R, c(1L, 0L))
  expect_error(rerar_history(0, NA, 1), "requires")
})

