context("Check strata")
library(CERMBlidar)

test_that("check_strata requires a data frame with columns: name, lower, upper", {
  # input is not a data frame
  expect_error(check_strata(matrix(0, nrow = 3, ncol = 3)))

  # input does not have expected columns
  expect_error(check_strata(data.frame(foo = "foo", lower = 1, other = 2)))
})


test_that("check_strata does not care about colname case or column order", {
  strata <- data.frame(LOWER = 1, name = "foo", UpPeR = 2)
  res <- check_strata(strata)
  expect_identical(colnames(res), c("name", "lower", "upper"))
})


test_that("check_strata drops extra columns", {
  strata <- data.frame(name = c("a", "b"),
                       lower = c(0, 5),
                       upper = c(5, 10),
                       extra = c("foo", "bar"))

  res <- check_strata(strata)
  expect_identical(colnames(res), c("name", "lower", "upper"))
})


test_that("check_strata detects overlapping layers", {
  # layers mid and upper overlap
  strata <- data.frame(name = c("lower", "mid", "upper"),
                       lower = c(1, 5, 8),
                       upper = c(5, 10, Inf))

  expect_error(check_strata(strata))
})


test_that("check_strata accepts gaps between layers", {
  # gap between lower and mid layers
  strata <- data.frame(name = c("lower", "mid", "upper"),
                       lower = c(1, 5, 10),
                       upper = c(2, 10, Inf))

  expect_identical(check_strata(strata), strata)
})


test_that("check_strata orders layers by upper height", {
  strata <- data.frame(name = c("c", "d", "a", "b"),
                       lower = c(5, 10, 0, 1),
                       upper = c(10, Inf, 1, 5))

  expected <- dplyr::arrange(strata, upper)

  expect_equivalent(check_strata(strata), expected)
})
