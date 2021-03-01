context("Strata cover")
library(CERMBlidar)
library(raster)

test_that("Input object must be a RasterStack with more than one layer", {
  r <- raster(matrix(0, ncol = 10, nrow = 10))

  # Passing a RasterLayer instead of a RasterStack
  expect_error(get_stratum_cover(r))

  # Passing a RasterStack with a single layer
  expect_error( get_stratum_cover( stack(list(r)) ) )
})

test_that("Input RasterStack with layer names must have '*ground*'", {
  r <- raster(matrix(0, ncol = 10, nrow = 10))

  # Stack with no named ground layer
  rr <- stack(list(stratum1 = r, stratum2 = r, stratum3 = r))

  expect_error(get_stratum_cover(rr))
})
