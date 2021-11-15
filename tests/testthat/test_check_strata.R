context("Check strata")
library(CERMBlidar)

test_that("check_strata requires a data frame with columns: label, lowerht, upperht", {
  # input is not a data frame
  expect_error(check_strata(matrix(0, nrow = 3, ncol = 3)))

  # input does not have expected columns
  expect_error(check_strata(data.frame(foo = "foo", lowerht = 1, other = 2)))
})


test_that("check_strata does not care about colname case or column order", {
  strata <- data.frame(LOWERht = 1, Label = "foo", UpPeRhT = 2)
  res <- check_strata(strata)
  expect_true( all(c("lowerht", "upperht", "label") %in% colnames(res)) )
})


test_that("check_strata retains any extra columns", {
  strata <- data.frame(label = c("a", "b"),
                       lowerht = c(0, 5),
                       upperht = c(5, 10),
                       extra1 = c("foo", "bar"),
                       extra2 = c(100, 200))

  res <- check_strata(strata)
  expect_true( all(colnames(strata) %in% colnames(res)) )
})


test_that("check_strata detects overlapping layers", {
  # layers mid and upper overlap
  strata <- data.frame(label = c("lower", "mid", "upper"),
                       lowerht = c(1, 5, 8),
                       upperht = c(5, 10, Inf))

  expect_error(check_strata(strata))
})


test_that("check_strata accepts gaps between layers", {
  # gap between lower and mid layers
  strata <- data.frame(label = c("lower", "mid", "upper"),
                       lowerht = c(1, 5, 10),
                       upperht = c(2, 10, Inf))

  expect_identical(check_strata(strata), strata)
})


test_that("check_strata orders layers by upper height", {
  strata <- data.frame(label = c("c", "d", "a", "b"),
                       lowerht = c(5, 10, 0, 1),
                       upperht = c(10, Inf, 1, 5))

  expected <- dplyr::arrange(strata, upperht)

  expect_equivalent(check_strata(strata), expected)
})
