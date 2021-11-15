context("Transform strata")
library(CERMBlidar)

test_that("get_strata_transform with valid strata", {
  s_from <- data.frame(label = paste0("layer", 1:10),
                       lowerht = c(-999, 2:10 / 2),
                       upperht = c(2:10 / 2, 999))

  s_to <- data.frame(label = c("lower", "upper"),
                     lowerht = c(-999, 2.5),
                     upperht = c(2.5, 999))

  xt <- get_strata_transform(s_from, s_to)

  expected_index <- c(1,1,1,1,2,2,2,2,2,2)
  expect_equal(xt$to_index, expected_index)
  expect_equal(xt$to_label, s_to$label[expected_index])
})


test_that("get_strata_transform with mis-aligned strata", {

  # Input strata have cut-points on whole metres
  s_from <- data.frame(label = paste0("input_layer", 1:10),
                       lowerht = c(-999.0, 1:9),
                       upperht = c(1:9, 999.0))

  # Destination strata have a cut-point on 2.5m
  s_to <- data.frame(label = paste0("dest_layer", 1:4),
                     lowerht = c(-999, 1, 2.5, 5),
                     upperht = c(1, 2.5, 5, 999))

  testthat::expect_error( get_strata_transform(s_from, s_to), regexp = "do not align" )
})
