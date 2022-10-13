# Check that an argument is a scalar. If it is a vector, issue
# a warning and return just the first value
.as_scalar <- function(x) {
  if (length(x) == 1) x
  else {
    nm <- deparse(substitute(x))
    if (length(x) == 0) {
      stop("Expected a value for ", nm)
    }
    else {
      warning("Expected a single value for ", nm, ". Only using first value.")
      x[1]
    }
  }
}


# Check that a boolean scalar argument is valid. Convert NULL to FALSE.
.as_boolean <- function(x) {
  nm <- deparse(substitute(x))
  if (is.null(x)) {
    warning("Treating NULL value for ", nm, " as FALSE")
    FALSE
  }
  else {
    if (!is.logical(x)) {
      stop("Expected a boolean value for ", nm)
    }
    .as_scalar(x)
  }
}


# Check if a file name has a .zip extension
.is_zipped <- function(path) {
  stringr::str_detect(path, "\\.zip\\s*$")
}


# Ensure a raster-type object is a terra::SpatRaster, doing conversion
# from raster package objects if necessary. If the input is something other than
# a raster::Raster* or terra::SpatRaster object, issue an error message.
#
.as_spat_raster <- function(r) {
  if (.is_terra_object(r)) {
    # nothing to do - just return input
    r
  } else if (.is_raster_object(r)) {
    r <- terra::rast(r)
  } else {
    objname <- deparse(substitute(r))
    msg <- glue::glue("{objname} should be a SpatRaster (terra package) object or a Raster* (raster package) object")
    stop(msg)
  }
}


# Check if x is a raster package Raster* object
.is_raster_object <- function(x) {
  grepl("^Raster(Layer|Stack|Brick)", class(x))
}


# Check if x is a terra package SpatRaster object
.is_terra_object <- function(x) {
  inherits(x, "SpatRaster")
}


