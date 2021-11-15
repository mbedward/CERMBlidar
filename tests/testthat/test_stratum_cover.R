context("Stratum cover")
library(CERMBlidar)
library(raster)
library(terra)

test_that("Input object must be a raster with more than one layer", {
  m <- matrix(0, ncol = 10, nrow = 10)

  r <- raster(m)

  # Passing a RasterLayer instead of a RasterStack
  expect_error(get_stratum_cover(r))

  # Passing a RasterStack with a single layer
  expect_error( get_stratum_cover( stack(list(r)) ) )

  # Passing a SpatRast with a single layer
  r <- terra::rast(m)
  expect_error( get_stratum_cover())
})

test_that("Input RasterStack with layer names must have '*ground*'", {
  m <- matrix(0, ncol = 10, nrow = 10)

  # Stack with no named ground layer
  r <- raster(m)
  rr <- stack(list(stratum1 = r, stratum2 = r, stratum3 = r))
  expect_error(get_stratum_cover(rr))

  # SpatRast with no named ground layer
  m <- array(0, dim = c(10, 10, 3))
  rr <- terra::rast(m, nlyrs = 3)
  names(rr) <- c("foo", "bar", "baz")
  expect_error(get_stratum_cover(rr))
})

test_that("Cover is calculated correctly", {
  Ntests <- 100
  nrows <- 10
  ncols <- 10

  fn_counts <- function(max_count) {
    matrix(sample(0:max_count, size = nrows * ncols, replace = TRUE),
           nrow = nrows, ncol = ncols)
  }

  fn_make_rast <- function(ground, shrub, tree) {
    r <- lapply(list(ground, shrub, tree), terra::rast)
    r <- terra::rast(r)
    names(r) <- c("ground", "shrub", "tree")

    r
  }

  shrub_ok <- logical(Ntests)
  tree_ok <- logical(Ntests)

  for (i in 1:Ntests) {
    ground_counts <- fn_counts(20)
    shrub_counts <- fn_counts(40)
    tree_counts <- fn_counts(80)

    expected_shrub_cover <- shrub_counts / (shrub_counts + ground_counts)
    expected_tree_cover <- tree_counts / (tree_counts + shrub_counts + ground_counts)

    rcounts <- fn_make_rast(ground_counts, shrub_counts, tree_counts)
    rcover <- get_stratum_cover(rcounts)

    shrub_ok[i] <- all.equal(terra::values(rcover[[1]], mat=FALSE), as.vector(t(expected_shrub_cover)))
    tree_ok[i] <- all.equal(terra::values(rcover[[2]], mat=FALSE), as.vector(t(expected_tree_cover)))
  }

  expect_equal(sum(shrub_ok), Ntests)
  expect_equal(sum(tree_ok), Ntests)

})

