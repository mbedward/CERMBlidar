context("Strata cover")
library(CERMBlidar)
library(raster)

test_that("Input RasterStack must have a ground layer", {
  r <- raster(matrix(0, ncol = 10, nrow = 10))

  # Stack with no ground layer
  rr <- stack(list(stratum1 = r, stratum2 = r, stratum3 = r))

  expect_error(get_stratum_cover(rr))
})
